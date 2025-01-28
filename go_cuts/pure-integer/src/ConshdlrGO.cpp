//
// Created by Lukas Schuermann
//

#include <iostream>
#include "ConshdlrGO.h"

#define DEFAULT_MAXROUNDS             1 /**< maximal number of GMI separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        25 /**< maximal number of GMI separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          10 /**< maximal number of GMI cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT      50 /**< maximal number of GMI cuts separated per separation round in root node */

#define DEFAULT_AWAY              0.01 /**< minimal fractionality of a basic variable in order to try GMI cut - default */
#define DEFAULT_MAX_DYN          1.0e+6 /**< maximal valid range max(|weights|)/min(|weights|) of cut coefficients - default */
#define DEFAULT_MAX_SUPP_ABS        100 /**< maximum cut support - absolute value in the formula - default */
#define DEFAULT_MAX_SUPP_REL       0.05 /**< maximum cut support - relative value in the formula - default */
#define CONSHDLR_NAME   "GO"


/**
 * UPWARDS = What if x = 1
 * DOWNWARDS = What if x = 0
 **/
enum CutDir
{
    UPWARDS = -1,
    DOWNWARDS = 1
};

static
bool getCPCFromRowAggr(
        SCIP*                 scip,               /**< pointer to the SCIP environment */
        ConshdlrGO*          objcpc,             /**< pointer to the conshdlr object */
        int                   ncols,              /**< number of columns in the LP */
        int                   nrows,              /**< number of rows in the LP */
        SCIP_COL**            cols,               /**< columns of the LP */
        SCIP_ROW**            rows,               /**< rows of the LP */
        SCIP_AggrRow*         aggrRow,
        int                   basisind,           /**< index of basic variable */
        const SCIP_Real*      binvrow,            /**< row of the basis inverse */
        SCIP_Real*            cutcoefs,           /**< array for cut elements in sparse format - must be of size ncols */
        int*                  cutind,             /**< array for indices of nonzero cut coefficients - must be of size ncols */
        int*                  cutnz,              /**< pointer to store number of nonzero elements in the cut */
        SCIP_Real*            cutrhs,             /**< pointer to store cut rhs */
        SCIP_Real*            workcoefs,          /**< working array of size ncols, allocated by caller for efficiency */
        const CutDir          cutdir,             /**< direction of the cut */
        const SCIP_Real       factor,             /**< factor for the reduced cost comparison */
        bool isLocal
){
    SCIP_COL* col;
    SCIP_ROW* row;
    SCIP_Var* var;
    double rowElem;

    assert(scip != nullptr);
    assert(cols != nullptr);
    assert(rows != nullptr);
    assert(binvrow != nullptr);
    assert(cutcoefs != nullptr);
    assert(cutind != nullptr);
    assert(cutnz != nullptr);
    assert(cutrhs != nullptr);
    assert(workcoefs != nullptr);

    /* clear cutcoefs and initialize cutcoefs
     * Downwards: (x_i - lb) >= 1 --> x_i >= 1 + lb
     * Upwards: (ub - x_i) >= 1 --> - x_i >= 1 - ub
     * */
    BMSclearMemoryArray(workcoefs, ncols);
    var = SCIPcolGetVar(cols[basisind]);
    *cutrhs = cutdir == DOWNWARDS ? 1 + SCIPvarGetLbLocal(var) : 1 - SCIPvarGetUbLocal(var);
    workcoefs[SCIPcolGetLPPos(cols[basisind])] = cutdir;

    std::vector<double> max_factors_vars = std::vector<double>();
    std::vector<double> max_factors_slack = std::vector<double>();

    double red_costs;
    *cutnz = 0;
    /* check which original variables join the cut */
    for(int i = 0; i < aggrRow->nnz; i++)
    {
        /* get variable and check if not fixed */
        var = SCIPgetVars(scip)[aggrRow->inds[i]];
        if(SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)))
            continue;

        SCIP_Real QUAD(rowelem_quad);
        col = SCIPvarGetCol(var);
        QUAD_ARRAY_LOAD(rowelem_quad, aggrRow->vals, aggrRow->inds[i]);

        switch (SCIPcolGetBasisStatus(col)){
            case SCIP_BASESTAT_LOWER:
                /* Check if element in simplex row has correct sign */
                rowElem = QUAD_TO_DBL(rowelem_quad);
                if(SCIPisPositive(scip, rowElem * cutdir)){
                    /* check if column can be excluded due to optimality condition */
                    red_costs = SCIPgetVarRedcost(scip, var);
                    if(SCIPisLE(scip, factor * rowElem, red_costs))
                        continue;

                    /* Should not happen if checked earlier */
                    if(!SCIPvarIsIntegral(var)){
                        return false;
                    }

                    /* set coefficient */
                    workcoefs[SCIPcolGetLPPos(col)] = 1;
                    if(!SCIPvarIsBinary(var)){
                        *cutrhs += SCIPvarGetLbLocal(var);
                    }
                    (*cutnz)++;
                }
                continue;
            case SCIP_BASESTAT_UPPER:
                /* Check if element in simplex row has correct sign */
                rowElem = QUAD_TO_DBL(rowelem_quad);
                if(SCIPisPositive(scip, -cutdir * rowElem)){
                    /* check if column can be excluded due to optimality condition */
                    red_costs = SCIPgetVarRedcost(scip, var);
                    max_factors_vars.push_back(red_costs / rowElem);
                    if(SCIPisGE(scip, factor * rowElem, red_costs))
                        continue;

                    /* Should not happen if checked earlier */
                    if(!SCIPvarIsIntegral(var)){
                        return false;
                    }

                    /* set coefficient and rhs */
                    workcoefs[SCIPcolGetLPPos(col)] = -1;
                    if(!SCIPvarIsBinary(var)){
                        *cutrhs -= SCIPvarGetUbLocal(var);
                    }else{
                        *cutrhs -= 1;
                    }
                    (*cutnz)++;
                }
                continue;
            default:
                /* In any other case: skip */
                continue;
        }
    }

    if(isLocal && objcpc->maxsupp_ < *cutnz)
        return false;

    if(!isLocal && ncols * objcpc->root_supp_rel_ < *cutnz)
        return false;

    /* check which slack variables join the cut */
    for(int c = 0; c < nrows; ++c)
    {
        row = rows[c];
        assert(row != nullptr);
        /* If the slack variable is fixed, we can ignore this cut coefficient. */
        if( SCIPisFeasZero(scip, SCIProwGetRhs(row) - SCIProwGetLhs(row)) )
            continue;
        /* Get simplex tableau element. */
        switch ( SCIProwGetBasisStatus(row) )
        {
            case SCIP_BASESTAT_LOWER:
                /* check if coefficient has correct sign */
                if(SCIPisPositive(scip, -cutdir * binvrow[SCIProwGetLPPos(row)]))
                {
                    /* Take element if nonbasic at lower bound.
                    * Then Slack Variable is at upper bound (\leq 0) -> flip it */

                    /* we might exclude the variable due to the gap */
                    /* would not have negative reduced cost with the changed objective function */
                    if(SCIPisGE(scip, factor * binvrow[SCIProwGetLPPos(row)], -SCIProwGetDualsol(row)))
                    {
                        continue;
                    }
                    /* If there is a non-integral slack variable, the cut derivation does not work -> abort */
                    if(!SCIProwIsIntegral(row)){
                        return false;
                    }
                    /* we do not actually add slack variable to cut but its row */
                    *cutrhs += (SCIProwGetLhs(row) - SCIProwGetConstant(row));
                    SCIP_COL** rowcols = SCIProwGetCols(row);
                    SCIP_Real* rowvals = SCIProwGetVals(row);
		            for(int i = 0; i < SCIProwGetNLPNonz(row); ++i )
                    {
                        if(SCIPisLT(scip, SCIPcolGetLb(rowcols[i]), SCIPcolGetUb(rowcols[i])))
                        {
                            workcoefs[SCIPcolGetLPPos(rowcols[i])] += rowvals[i];
                        }else// if(SCIPisEQ(scip, SCIPcolGetLb(rowcols[i]), 1))
                        {/* if fixed on 1, we have to adjust the rhs */
                            *cutrhs -= rowvals[i] * SCIPcolGetPrimsol(rowcols[i]);
                        }
                    }
                }
                continue;
            case SCIP_BASESTAT_UPPER:
                if(SCIPisPositive(scip, cutdir * binvrow[SCIProwGetLPPos(row)]))
                {
                    /* we might exclude the variable due to the gap */
                    /* would not have negative reduced cost with the changed objective function */
                    if(SCIPisLE(scip, factor * binvrow[SCIProwGetLPPos(row)], -SCIProwGetDualsol(row)))
                    {
                        continue;
                    }
                    /* No non-integral rows allowed */
                    if(!SCIProwIsIntegral(row)){
                        return false;
                    }
                    /* we do not actually add slack variable to cut but its row */
                    *cutrhs -= (SCIProwGetRhs(row) - SCIProwGetConstant(row));
                    SCIP_COL** rowcols = SCIProwGetCols(row);
                    SCIP_Real* rowvals = SCIProwGetVals(row);
                    for(int i = 0; i < SCIProwGetNLPNonz(row); ++i )
                    {
                        if(SCIPisLT(scip, SCIPcolGetLb(rowcols[i]), SCIPcolGetUb(rowcols[i])))
                        {
                            workcoefs[SCIPcolGetLPPos(rowcols[i])] -= rowvals[i];
                        }else
                        {/* if fixed on 1, we have to adjust the rhs */
                            *cutrhs += rowvals[i] * SCIPcolGetPrimsol(rowcols[i]);
                        }
                    }
                }
                continue;
            default:
                /* Basic variable: skip */
                continue;
        }
    }

    /* set cut coefficients */
    *cutnz = 0;
    for(int c = 0; c < ncols; ++c )
    {
        col = cols[c];
        assert(col != nullptr);
        auto i = SCIPcolGetLPPos(col);
        assert(0 <= i);

        if(!SCIPisZero(scip, workcoefs[i]))
        {
            cutcoefs[*cutnz] = workcoefs[i];
            cutind[*cutnz] = c;
            (*cutnz)++;
        }
    }

    if(isLocal && objcpc->maxsupp_ < *cutnz)
        return false;

    if(!isLocal && ncols * objcpc->root_supp_rel_ < *cutnz)
        return false;

    return true;
}


static
SCIP_RETCODE addCPCToModel(
        SCIP*           scip,
        SCIP_Conshdlr*  conshdlr,
        SCIP_Result*    result,
        SCIP_Col**      cols,
        SCIP_Real       cutlhs,
        CutDir          cutdir,
        bool            cutislocal,
        SCIP_Real*      cutcoefs,
        const int*      cutind,
        int             cutnz,
        int*            ncuts,
        int             c
){
    char cutname[SCIP_MAXSTRLEN];
    SCIP_ROW* cut;
    int j;

    /* construct cut name */
    if( c >= 0)
        (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "cpc%d_x%d_d%d", SCIPgetNLPs(scip), c, cutdir);
    else
        (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "cpc%d_s%d_d%d", SCIPgetNLPs(scip), -c-1, cutdir);

    /* create empty cut */
    SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &cut, conshdlr, cutname, cutlhs, SCIPinfinity(scip), cutislocal, FALSE, FALSE) );

    /* cache the row extension and only flush them if the cut gets added */
    SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

    /* collect all non-zero coefficients */
    for( j = 0; j < cutnz; ++j )
    {
        SCIP_CALL( SCIPaddVarToRow(scip, cut, SCIPcolGetVar(cols[cutind[j]]), cutcoefs[j]) );
    }

    if( SCIProwGetNNonz(cut) == 0 )
    {
        assert(SCIPisFeasNegative(scip, cutlhs));
        SCIPdebugMsg(scip, " -> GMI cut detected infeasibility with cut 0 <= %f.\n", cutlhs);
        *result = SCIP_CUTOFF;
//                    break;
    }

    /* Only take efficacious cuts, except for cuts with one non-zero coefficient (= bound
     * changes); the latter cuts will be handeled internally in sepastore. */
    if( SCIProwGetNNonz(cut) == 1 || SCIPisCutEfficacious(scip, nullptr, cut) )
    {
        SCIP_Bool infeasible;

        SCIPdebugMsg(scip, " -> found GMI cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f).\n",
                     cutname, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                     SCIPgetCutEfficacy(scip, nullptr, cut),
                     SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                     SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));

        /* flush all changes before adding the cut */
        SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

        SCIP_CALL( SCIPaddRow(scip, cut, FALSE, &infeasible) );

        /* add global cuts that are not implicit bound changes to the cut pool */
        if( ! cutislocal && SCIProwGetNNonz(cut) > 1 )
        {
            SCIP_CALL( SCIPaddPoolCut(scip, cut) );
        }

        if ( infeasible )
            *result = SCIP_CUTOFF;
        else
            *result = SCIP_SEPARATED;
        (*ncuts)++;
    }

    /* release the row */
    SCIP_CALL( SCIPreleaseRow(scip, &cut) );

    return SCIP_OKAY;
}

/**@name Callback methods
 *
 * @{
 */
/** frees specific constraint data */

/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(ConshdlrGO::scip_trans)
{   /*lint --e{715}*/
    SCIP_CONSDATA* sourcedata;
    SCIP_CONSDATA* targetdata;

    assert(conshdlr != nullptr);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), "CPC") == 0);
    assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
    assert(sourcecons != nullptr);
    assert(targetcons != nullptr);

    sourcedata = SCIPconsGetData(sourcecons);
    assert(sourcedata != nullptr);

    /* create target constraint */
    SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                              SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
                              SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
                              SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
                              SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

    return SCIP_OKAY;
}


SCIP_DECL_CONSSEPALP(ConshdlrGO::scip_sepalp){
    SCIP_VAR** vars;
    SCIP_COL** cols;
    SCIP_ROW** rows;
    SCIP_Real* binvrow;
    SCIP_Real* cutcoefs;
    SCIP_Real* workcoefs;
    SCIP_Real cutlhs;
    SCIP_Real absgap;
    SCIP_Real relgap;
    int* cutind;
    int* basisind;
    int nvars;
    int ncols;
    int nrows;
    int cutnz;
    int i, c;
    double away;
    bool cutislocal;

    *result = SCIP_DIDNOTRUN;

    /* max depth exceeded */
    if(maxDepth_ != -1 && maxDepth_ < SCIPgetDepth(scip))
    {
        return SCIP_OKAY;
    }
    /* Only call separator, if an optimal LP solution is at hand. */
    if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
        return SCIP_OKAY;

    if(SCIPisInfinity(scip, SCIPgetLocalDualbound(scip)))
        return SCIP_OKAY;

    /* Only call separator, if the LP solution is basic. */
    if( ! SCIPisLPSolBasic(scip) )
        return SCIP_OKAY;

    /* Only call separator, if there are fractional variables. */
    if( SCIPgetNLPBranchCands(scip) == 0 )
        return SCIP_OKAY;

    /* call frequency */
    if(SCIPgetDepth(scip) > 10){
        if(SCIPgetDepth(scip) % (2 * freq_) != 0){
            return SCIP_OKAY;
        }
    }

    /* Check for integrality of the model */
    if(!checkedIfInt_){
        if(!isModelInteger(scip)){
            return SCIP_OKAY;
        }
    }

    if(descending_supp_ && SCIPgetDepth(scip) > 1){
        maxsupp_ = max_supp_abs_ + int(max_supp_rel_ * SCIPgetNLPCols(scip));
        maxsupp_ = std::min(maxsupp_, int(0.1*SCIPgetNLPCols(scip))) / SCIPgetDepth(scip);
        maxsupp_ = std::max(maxsupp_, 1);
    }

    if(SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == lastnode_)
    {
        ncallsnode_++;
        if((SCIPgetDepth(scip) == 0 && maxroundsroot_ >= 0 && ncallsnode_ > maxroundsroot_) ||
           (SCIPgetDepth(scip) > 0 && maxrounds_ >= 0 && ncallsnode_ > maxrounds_))
            return SCIP_OKAY;
    }else{
        lastnode_ = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
        ncallsnode_ = 1;
    }

    /* get variables data */
    SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, nullptr, nullptr, nullptr, nullptr) );

    /* get LP data */
    SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
    SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

    /* exit if LP is trivial */
    if( ncols == 0 || nrows == 0 )
        return SCIP_OKAY;

    assert(SCIPisLE(scip, SCIPgetLocalDualbound(scip), SCIPgetPrimalbound(scip)));
    /* check local absolute gap */
    if(SCIPgetPrimalbound(scip) - SCIPgetLocalDualbound(scip) < SCIPinfinity(scip))
    {
        absgap = SCIPgetPrimalbound(scip) - SCIPgetLocalDualbound(scip);
        if(SCIPisZero(scip, SCIPgetLocalDualbound(scip)))
            relgap = -1;
        else
            relgap = absgap / abs(SCIPgetLocalDualbound(scip));

        assert(SCIPgetPrimalbound(scip) < SCIPinfinity(scip));
        assert(0 < absgap && absgap < SCIPinfinity(scip));
    }else{
        return SCIP_OKAY;
    }

    /* we do not expect good cuts for high gaps */
    int max_nofind = max_noFind_high_;
    int max_ntries = max_nTries_high_;

    if(relgap > rel_gap_low_)
    {
        max_nofind = max_noFind_low_;
        max_ntries = max_nTries_low_;
    }

    *result = SCIP_DIDNOTFIND;

    /* allocate temporary memory */
    SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, ncols) );
    SCIP_CALL( SCIPallocBufferArray(scip, &workcoefs, ncols) );
    SCIP_CALL( SCIPallocBufferArray(scip, &cutind, ncols) );
    SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );
    SCIP_CALL( SCIPallocBufferArray(scip, &binvrow, nrows) );

    /* get basis indices */
    SCIP_CALL( SCIPgetLPBasisInd(scip, basisind) );

    if(SCIPgetDepth(scip) == 0 ) {
        away = away_ / 3;
        cutislocal = false;
    }
    else {
        away = away_;
        cutislocal = true;
    }

    int ntries = 0;
    nRuns_++;

    /* For all basic columns belonging to integer variables, try to generate a CPC cut. */
    for( i = 0; i < nrows && ! SCIPisStopped(scip) && *result != SCIP_CUTOFF; ++i )
    {
        if(!cutislocal && ntries > SCIPgetNLPCols(scip) / 3){
            break;
        }
        if(ntries > max_nofind && cutislocal){
            break;
        }
        if(ntries > max_ntries && cutislocal){
            break;
        }

        SCIP_Bool tryrow;
        SCIP_Real primsol;

        tryrow = FALSE;
        c = basisind[i];
        primsol = SCIP_INVALID;
        SCIP_Var* var = nullptr;
        if( c >= 0 )
        {
            var = SCIPcolGetVar(cols[c]);
            if(!SCIPvarIsIntegral(var)){
                continue;
            }
            assert(c < ncols);
            assert(cols[c] != nullptr);

            primsol = SCIPcolGetPrimsol(cols[c]);

            /* shall be not too close to integer */
            if( (SCIPfeasFrac(scip, primsol) >= away) && (SCIPfeasFrac(scip, primsol) <= 1.0 - away) )
            {
                /* can happen if variable got fixed in earlier iteration */
                if(SCIPisLT(scip, primsol, SCIPvarGetLbLocal(var)) || SCIPisGT(scip, primsol, SCIPvarGetUbLocal(var))){
                    continue;
                }
                /* must either be < lb + 1 or > ub - 1 */
                if(SCIPisLT(scip, primsol, SCIPvarGetLbLocal(var) + 1) ||
                   SCIPisGT(scip, primsol, SCIPvarGetUbLocal(var) - 1)){
                    tryrow = TRUE;
                }
            }
        }
        /* found feasible basic variable */
        if (tryrow)
        {
            assert(var != nullptr);
            ntries++;
            SCIP_Real factor;
            bool success = false;

            /* get the row of B^-1 for this basic integer variable with fractional solution value */
            SCIP_CALL( SCIPgetLPBInvRow(scip, i, binvrow, nullptr, nullptr) );

            /* aggregate the simplex tableau row */
            SCIP_AggrRow* aggrrow;
            SCIP_Bool allowlocal = false;
            SCIP_Bool success2;
            SCIP_CALL( SCIPaggrRowCreate(scip, &aggrrow) );
            SCIP_CALL( SCIPaggrRowSumRows(scip, aggrrow, binvrow, nullptr, -1,
                                          true, allowlocal, 2, 1000000, &success2) );

            /* create a CPC out of the simplex tableau row */
            /* downwards */
            /* for constraint type DOWNWARDS the LP-value needs to be lower than lb + 1 */
            if(SCIPisLT(scip, primsol, SCIPvarGetLbLocal(var) + 1))
            {
                factor = absgap / SCIPfeasFrac(scip, primsol);
                success = getCPCFromRowAggr(scip, this, ncols, nrows, cols, rows, aggrrow, c, binvrow, cutcoefs, cutind,
                                        &cutnz, &cutlhs, workcoefs, DOWNWARDS, factor, cutislocal);

                if(success)
                {
                    /* potentially just change variable bound */
                    if(cutnz == 1){
                        cutlhs = std::ceil(cutlhs / std::fabs(cutcoefs[0]));
                        cutcoefs[0] = SCIPisPositive(scip, cutcoefs[0]) ? 1.0 : -1.0;
                    }
                    SCIP_CALL(addCPCToModel(scip, conshdlr, result, cols, cutlhs, DOWNWARDS, cutislocal, cutcoefs,
                                            cutind, cutnz, &nCuts_, c));
                    if(!cutislocal){
                        nCutsRoot_++;
                    }
                }

            }

            /* IF no success at DOWNWARDS, try UPWARDS:
             * for constraint type UPWARDS the LP-value needs to be higher than ub - 1 */
            if(!success && SCIPisGT(scip, primsol, SCIPvarGetUbLocal(SCIPcolGetVar(cols[c])) - 1))
            {
                factor = - absgap / (1 - SCIPfeasFrac(scip, primsol));
                success = getCPCFromRowAggr(scip, this, ncols, nrows, cols, rows, aggrrow, c, binvrow, cutcoefs, cutind,
                                            &cutnz, &cutlhs, workcoefs, UPWARDS, factor, cutislocal);
                if(success)
                {
                    /* potentially just change variable bound */
                    if(cutnz == 1){
                        cutlhs = std::ceil(cutlhs / std::fabs(cutcoefs[0]));
                        cutcoefs[0] = SCIPisPositive(scip, cutcoefs[0]) ? 1.0 : -1.0;
                    }
                    SCIP_CALL(addCPCToModel(scip, conshdlr, result, cols, cutlhs, UPWARDS, cutislocal, cutcoefs,
                                            cutind, cutnz, &nCuts_, c));
                    if(!cutislocal){
                        nCutsRoot_++;
                    }
                }
            }
            SCIPaggrRowFree(scip, &aggrrow);
        }
    }

    /* free temporary memory */
    SCIPfreeBufferArray(scip, &binvrow);
    SCIPfreeBufferArray(scip, &basisind);
    SCIPfreeBufferArray(scip, &workcoefs);
    SCIPfreeBufferArray(scip, &cutcoefs);
    SCIPfreeBufferArray(scip, &cutind);

    return SCIP_OKAY;

}

/** checks if model contains only variables with integer values in feasible IP solutions */
bool ConshdlrGO::isModelInteger(SCIP* scip) {
    maxsupp_ = max_supp_abs_ + int(max_supp_rel_ * SCIPgetNLPCols(scip));
    maxsupp_ = std::min(maxsupp_, int(0.1*SCIPgetNLPCols(scip)));

    checkedIfInt_ = true;
    start_gap_ = SCIPgetGap(scip);

    /* Check if all variable are integral */
    for(int i = 0; i < SCIPgetNVars(scip); i++)
    {
        if(!SCIPvarIsIntegral(SCIPgetVars(scip)[i]))
        {
            return false;
        }
    }
    return true;
}

ConshdlrGO::ConshdlrGO(SCIP *scip)
: ObjConshdlr(scip, "CPC", "counterpart cuts conshdlr", 50000, 400, 0,
                                                      1, 0, 0, 0, true, FALSE, FALSE,
                                                      SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_ALWAYS)
{
    checkedIfInt_ = false;
    lastnode_ = -1;
    ncallsnode_ = 0;
    nRuns_ = 0;
    nCuts_ = 0;
    start_gap_ = 0.0;

    maxrounds_ = DEFAULT_MAXROUNDS;
    maxroundsroot_ = DEFAULT_MAXROUNDSROOT;
    maxsepacuts_ = DEFAULT_MAXSEPACUTS;
    maxsepacutsroot_ = DEFAULT_MAXSEPACUTSROOT;
    away_ = DEFAULT_AWAY;
    maxdynamism_ = DEFAULT_MAX_DYN;
    max_supp_abs_ = DEFAULT_MAX_SUPP_ABS;
    max_supp_rel_ = DEFAULT_MAX_SUPP_REL;
}

