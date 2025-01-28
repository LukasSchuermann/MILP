//
// Created by lukasschuermann on 29.09.23.
//

#include <iostream>

#include "scip/cons_linear.h"
#include "ConshdlrGO.h"
#include "objscip/objscip.h"

#define DEFAULT_MAXROUNDS             1 /**< maximal number of GMI separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        25 /**< maximal number of GMI separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          10 /**< maximal number of GMI cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT      50 /**< maximal number of GMI cuts separated per separation round in root node */

#define DEFAULT_AWAY              0.01 /**< minimal fractionality of a basic variable in order to try GMI cut - default */
#define DEFAULT_MAX_DYN          1.0e+6 /**< maximal valid range max(|weights|)/min(|weights|) of cut coefficients - default */
#define DEFAULT_MAX_SUPP_ABS        100 /**< maximum cut support - absolute value in the formula - default */
#define DEFAULT_MAX_SUPP_REL       0.05 /**< maximum cut support - relative value in the formula - default */

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
bool checkNumerics(
        SCIP*               scip,
        SCIP_Real*          cutcoefs,
        int                 cutnz
){
    for(int i = 0; i < cutnz; i++){
        if(SCIPisZero(scip, cutcoefs[i])){
            cutcoefs[i] = 0;
        }else if(!SCIPisZero(scip, SCIPfeasFrac(scip, cutcoefs[i])) && SCIPisLT(scip, SCIPfeasFrac(scip, cutcoefs[i]), 1e-3)){
            return false;
        }
    }
    return true;
}

/** Calculate MIR cut based on GO-Cut idea */
static
bool getGOMIRFromRow(
        SCIP*                 scip,               /**< pointer to the SCIP environment */
        ConshdlrGO*          objcpc,             /**< pointer to the conshdlr object */
        int                   ncols,              /**< number of columns in the LP */
        int                   nrows,              /**< number of rows in the LP */
        SCIP_COL**            cols,               /**< columns of the LP */
        SCIP_ROW**            rows,               /**< rows of the LP */
        int                   basisind,           /**< index of basic variable */
        const SCIP_Real*      binvrow,            /**< row of the basis inverse */
        SCIP_Real*            cutcoefs,           /**< array for cut elements in sparse format - must be of size ncols */
        int*                  cutind,             /**< array for indices of nonzero cut coefficients - must be of size ncols */
        int*                  cutnz,              /**< pointer to store number of nonzero elements in the cut */
        SCIP_Real*            cutrhs,             /**< pointer to store cut rhs */
        SCIP_Real*            workcoefs,          /**< working array of size ncols, allocated by caller for efficiency */
        const SCIP_Real       bigR,               /**< factor for the reduced cost comparison */
        const SCIP_Real       smallG,
        SCIP_AggrRow*         aggrRow,
        double*               mir_effica,
        double                scalar,
        bool                  zeroOne,
        bool                  isLocal
){
    SCIP_VAR *var;
    SCIP_COL *col;
    SCIP_ROW *row;
    SCIP_Real factor;
    SCIP_Real rowelem;
    int i;
    int c;
    double red_cost;
    std::vector<int> lower_vars;
    std::vector<int> upper_vars;
    std::vector<int> slack_vars;
    double primsol = SCIPcolGetPrimsol(cols[basisind]);

    assert(scip != nullptr);
    assert(cols != nullptr);
    assert(rows != nullptr);
    assert(binvrow != nullptr);
    assert(cutcoefs != nullptr);
    assert(cutind != nullptr);
    assert(cutnz != nullptr);
    assert(cutrhs != nullptr);
    assert(workcoefs != nullptr);



    /* clear cutcoefs and initialize cutcoefs */
    BMSclearMemoryArray(workcoefs, ncols);
    BMSclearMemoryArray(cutcoefs, ncols); //TODO: only need it now
    *cutrhs = smallG / scalar + (std::floor(primsol) * bigR) / scalar;
    workcoefs[SCIPcolGetLPPos(cols[basisind])] = bigR / scalar;


    /* check which original variables join the cut */
    for (i = 0; i < aggrRow->nnz; i++) {
        SCIP_Real QUAD(val);
        var = SCIPgetVars(scip)[aggrRow->inds[i]];
        col = SCIPvarGetCol(var);
        QUAD_ARRAY_LOAD(val, aggrRow->vals, aggrRow->inds[i]);

        /* ignore fixed variables */
        if (SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)))
            continue;
        /* check for simplex basic status */
        switch (SCIPcolGetBasisStatus(col)) {
            /* just take the coefficient if variable at lower bound */
            case SCIP_BASESTAT_LOWER:
                /* Take element if nonbasic at lower bound.
                  * reduced costs >= 0
                  * --> need a_ij > 0
                  * */
                /* Calculate \bar{r}_j */
                red_cost = SCIPgetVarRedcost(scip, var) - QUAD_TO_DBL(val) * bigR;
                /* new reduced costs have to be < 0 to obtain a positive coefficient */
                if (SCIPisNegative(scip, red_cost)) {
                    factor = -(red_cost / scalar);
                    workcoefs[SCIPcolGetLPPos(col)] = factor;
                    *cutrhs += factor * SCIPvarGetLbLocal(var);
                }
                break;
                /* variables need to be flipped, either here or in SCIPcalcMIR call */
            case SCIP_BASESTAT_UPPER:
                /* Flip element if nonbasic at upper bound.
                 * EDIT: We do not flip the element, but let SCIP do the job
                 * NOTICE: We need to include the rhs for ignored variables, since they must be flipped
                 *
                 * reduced costs <= 0
                 * --> need a_ij < 0
                 * */
                red_cost = SCIPgetVarRedcost(scip, var) - QUAD_TO_DBL(val) * bigR;
                /* new reduced costs have to be > 0 since these variables are going to get complemented by SCIP */
                if (SCIPisPositive(scip, red_cost)) {
                    factor = -(red_cost / scalar);
                    workcoefs[SCIPcolGetLPPos(col)] = factor;
                    *cutrhs += factor * SCIPvarGetUbLocal(var);
                }
                break;
            default:
                /* Basic variable: skip */
                continue;
        }
    }

    *cutnz = 0;
    for (c = 0; c < ncols; ++c) {
        col = cols[c];
        assert(col != nullptr);
        i = SCIPcolGetLPPos(col);
        int var_ind = SCIPcolGetVarProbindex(col);
        assert(0 <= i);

        if (!SCIPisZero(scip, workcoefs[i])) {
            cutcoefs[*cutnz] = -workcoefs[i];
            cutind[*cutnz] = var_ind;
            (*cutnz)++;
        }
    }

    if(isLocal && objcpc->maxsupp_ < *cutnz)
        return false;

    if(!isLocal && ncols * objcpc->root_supp_rel_ < *cutnz)
        return false;

    SCIP_AggrRow *mir_aggr;
    SCIP_CALL(SCIPaggrRowCreate(scip, &mir_aggr));
    SCIP_CALL(SCIPaggrRowAddCustomCons(scip, mir_aggr, cutind, cutcoefs, *cutnz, -*cutrhs, 1.0, 1, FALSE));


    std::vector<int> slack_inds = std::vector<int>();
    std::vector<double> slack_weights = std::vector<double>();
    std::vector<bool> slack_dir = std::vector<bool>();
    /* check which slack variables join the cut */
    for (c = 0; c < nrows; ++c) {
        row = rows[c];
        assert(row != nullptr);
        bool add_slack = false;
        factor = 0;
        /* If the slack variable is fixed, we can ignore this cut coefficient. */
        if (SCIPisFeasZero(scip, SCIProwGetRhs(row) - SCIProwGetLhs(row)))
            continue;
        double bVal = binvrow[SCIProwGetLPPos(row)];

        /* reduced costs of slack variable = -dual_var of row */
        switch (SCIProwGetBasisStatus(row)) {
            case SCIP_BASESTAT_LOWER:
                if (SCIPisGE(scip, bVal, 0))
                    continue;

                assert(!SCIPisSumNegative(scip, SCIProwGetDualsol(row)));

                assert(!SCIPisInfinity(scip, -SCIProwGetLhs(row)));

                /* calculate new dual value \bar{y}_j */
                factor = SCIProwGetDualsol(row) + bigR * bVal;
                /* if dual variable still feasible (>= 0) continue */
                if (SCIPisGE(scip, factor, 0))
                    continue;
                /* scale variable for cut coefficient
                 * and invert it since aggr_row is a <= cut */
                factor /= -scalar;

                slack_weights.push_back(factor);
                slack_inds.push_back(SCIProwGetLPPos(row));
                slack_dir.push_back(false);
                break;
            case SCIP_BASESTAT_UPPER:
                if (SCIPisLE(scip, bVal, 0))
                    continue;

                assert(!SCIPisSumPositive(scip, SCIProwGetDualsol(row)));
                /* calculate new dual value \bar{y}_j */
                factor = SCIProwGetDualsol(row) + bigR * bVal;
                /* if dual variable still feasible (<= 0) continue */
                if (SCIPisLE(scip, factor, 0))
                    continue;
                /* scale variable for cut coefficient
                 * and invert it since aggr_row is a <= cut */
                factor /= -scalar;

                slack_weights.push_back(factor);
                slack_inds.push_back(SCIProwGetLPPos(row));
                slack_dir.push_back(true);

                break;
            case SCIP_BASESTAT_ZERO:
                /* Nonbasic free variable at zero: cut coefficient is zero, skip */
                SCIPdebugMsg(scip, "Free nonbasic slack variable, this should not happen!\n");
                continue;
            case SCIP_BASESTAT_BASIC:
                continue;
            default:
                assert(false);
        }
    }

    /* add slack variables to aggr-row */
    int nSlack = slack_inds.size();
    SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &mir_aggr->rowsinds, mir_aggr->rowssize, nSlack));
    SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &mir_aggr->slacksign, mir_aggr->rowssize, nSlack));
    SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &mir_aggr->rowweights, mir_aggr->rowssize, nSlack));
    mir_aggr->nrows += nSlack;
    mir_aggr->rowssize += nSlack;
    for (i = 0; i < nSlack; i++) {
        if (!slack_dir[i]) {
            mir_aggr->rowsinds[i] = slack_inds[i];
            mir_aggr->slacksign[i] = -1;
        } else {
            mir_aggr->rowsinds[i] = slack_inds[i];
            mir_aggr->slacksign[i] = +1;
        }
        mir_aggr->rowweights[i] = slack_weights[i];
    }

    double scale;
    int cutrank = -1;
    SCIP_Bool cutislocal = false, cutsuccess = false;

    /* calculate MIR cut */
    /* set scale of MIR */
    scale = zeroOne ? 1 / (objcpc->scale_zo_ * bigR) : 1 / (objcpc->scale_gen_ * bigR);

    SCIP_CALL(
            SCIPcalcMIR(scip, nullptr, TRUE, 0.999, FALSE, FALSE, FALSE, nullptr, nullptr, 0.001, 0.999, scale, mir_aggr,
                        cutcoefs, cutrhs, cutind, cutnz, mir_effica, &cutrank, &cutislocal, &cutsuccess));

    SCIPaggrRowFree(scip, &mir_aggr);

    if(isLocal && objcpc->maxsupp_ < *cutnz)
        return false;

    if(!isLocal && ncols * objcpc->root_supp_rel_ < *cutnz)
        return false;

    if(cutsuccess)
        return checkNumerics(scip, cutcoefs, *cutnz);
    return false;
}

static
SCIP_RETCODE addCPCToModel(
        SCIP*           scip,
        SCIP_Conshdlr*  conshdlr,
        SCIP_Result*    result,
        SCIP_Col**      cols,
        SCIP_Real       cutlhs,
        SCIP_Real       cutrhs,
        CutDir          cutdir,
        bool            cutislocal,
        SCIP_Real*      cutcoefs,
        const int*      cutind,
        int             cutnz,
        int*            ncuts,
        int             c,
        double*         effica
){
    char cutname[SCIP_MAXSTRLEN];
    SCIP_ConsData*  consdata;
    SCIP_ROW* cut;
    int j;

    /* construct cut name */
    if( c >= 0)
        (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "go%d_x%d_d%d", SCIPgetNLPs(scip), c, cutdir);
    else
        (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "go%d_s%d_d%d", SCIPgetNLPs(scip), -c-1, cutdir);

    /* create empty cut */
    SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &cut, conshdlr, cutname, cutlhs, cutrhs, cutislocal, FALSE, FALSE) );

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
        SCIPdebugMsg(scip, " -> GO cut detected infeasibility with cut 0 <= %f.\n", cutlhs);
        *result = SCIP_CUTOFF;
    }
    *effica = SCIPgetCutEfficacy(scip, nullptr, cut);
    /* Only take efficacious cuts, except for cuts with one non-zero coefficient (= bound
     * changes); the latter cuts will be handeled internally in sepastore. */
    if( SCIProwGetNNonz(cut) == 1 || SCIPisCutEfficacious(scip, nullptr, cut) )
    {
        SCIP_Bool infeasible;

        SCIPdebugMsg(scip, " -> found GO cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f).\n",
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

//    /* create constraint data for target constraint */
//    SCIP_CALL( consdataCreate(scip, &targetdata,
//                              sourcedata->customer, sourcedata->day, sourcedata->type, sourcedata->node) );

    /* create target constraint */
    SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                              SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
                              SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
                              SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
                              SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

    return SCIP_OKAY;
}

void ConshdlrGO::initConshdlr(SCIP* scip) {
    isInit_ = true;
    maxsupp_ = max_supp_abs_ + int(max_supp_rel_ * SCIPgetNLPCols(scip));
    maxsupp_ = std::min(maxsupp_, int(0.1*SCIPgetNLPCols(scip)));
}


SCIP_DECL_CONSSEPALP(ConshdlrGO::scip_sepalp){
    SCIP_VAR** vars;
    SCIP_COL** cols;
    SCIP_ROW** rows;
    SCIP_Real* binvrow;
    SCIP_Real* cutcoefs;
    SCIP_Real* workcoefs;
    SCIP_Real cutrhs;
    SCIP_Real absgap;
    SCIP_Real relgap;
    int* cutind;
    int* basisind;
    int nvars;
    int ncols;
    int nrows;
    int ncalls;
    int maxsepacuts;
    int ncuts;
    int cutnz;
    int i, j, c;
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
    {
        return SCIP_OKAY;
    }

    if(SCIPisInfinity(scip, SCIPgetLocalDualbound(scip)))
    {
        return SCIP_OKAY;
    }

    /* Only call separator, if the LP solution is basic. */
    if( ! SCIPisLPSolBasic(scip) )
    {
        return SCIP_OKAY;
    }

    /* Only call separator, if there are fractional variables. */
    if( SCIPgetNLPBranchCands(scip) == 0 )
    {
        return SCIP_OKAY;
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

    if(SCIPgetDepth(scip) > 10){
        if(SCIPgetDepth(scip) % (2 * freq_) != 0){
            return SCIP_OKAY;
        }
    }

    if(!isInit_){
        initConshdlr(scip);
    }

    /* smaller support the deeper into branching tree */
    if(descending_supp_ && SCIPgetDepth(scip) > 0){
        maxsupp_ = max_supp_abs_ + int(max_supp_rel_ * SCIPgetNLPCols(scip));
        maxsupp_ = std::min(maxsupp_, int(0.1*SCIPgetNLPCols(scip))) / SCIPgetDepth(scip);
        maxsupp_ = std::max(maxsupp_, 2);
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
    }else
    {
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

    /* get the maximal number of cuts allowed in a separation round */
    if(SCIPgetDepth(scip) == 0 ) {
        maxsepacuts = std::min(int(SCIPgetNLPCols(scip) / 4), maxsepacutsroot_);
        away = away_ / 3;
        cutislocal = false;
    }
    else {
        maxsepacuts = std::min(int(SCIPgetNLPCols(scip) / 20), maxsepacuts_);
        away = away_;
        cutislocal = true;
    }

    if( maxsepacuts == -1 )
        maxsepacuts = INT_MAX;

    ncuts = 0;
    int ntries = 0;
    double effica = 0;
    nRuns_++;

    /* For all basic columns belonging to integer variables, try to generate a CPC cut. */
    for( i = 0; i < nrows && ncuts < maxsepacuts && ! SCIPisStopped(scip) && *result != SCIP_CUTOFF; ++i )
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
        bool tryrow = FALSE;
        SCIP_Real primsol;

        tryrow = false;
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
                tryrow = TRUE;
            }
        }
        if (tryrow)
        {
            ntries++;
            bool success = false;
            double intDist = primsol - std::floor(primsol);
            double smallG = absgap;
            double bigR = (absgap + smallG) / intDist;
            double mir_effica = 0;

            /* get the row of B^-1 for this basic integer variable with fractional solution value */
            SCIP_CALL( SCIPgetLPBInvRow(scip, i, binvrow, nullptr, nullptr) );

            /* aggregate the simplex tableau row */
            SCIP_AggrRow* aggrrow;
            SCIP_Bool allowlocal = false;
            SCIP_Bool success2;
            SCIP_CALL( SCIPaggrRowCreate(scip, &aggrrow) );
            SCIP_CALL( SCIPaggrRowSumRows(scip, aggrrow, binvrow, nullptr, -1,
                                          true, allowlocal, 2, 1000000, &success2) );

            success = getGOMIRFromRow(scip, this, ncols, nrows, cols, rows, c, binvrow, cutcoefs, cutind,
                                      &cutnz, &cutrhs, workcoefs, bigR, smallG, aggrrow, &mir_effica, 1,
                                      (primsol < 1 && primsol > 0), cutislocal);
            if(success)
            {
                /* potentially just change variable bound */
                SCIP_CALL(addCPCToModel(scip, conshdlr, result, cols, -SCIPinfinity(scip), cutrhs, DOWNWARDS, cutislocal, cutcoefs,
                                        cutind, cutnz, &nCuts_, c, &effica));
                assert(SCIPisGT(scip, effica, 0));

                if(!cutislocal){
                    nCutsRoot_++;
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

SCIP_DECL_CONSPRINT(ConshdlrGO::scip_print)
{
    return SCIP_OKAY;
}

ConshdlrGO::ConshdlrGO(SCIP *scip)
: ObjConshdlr(scip, "CPC", "counterpart cuts conshdlr", 50000, 400, 0,
                                                      1, 0, 0, 0, true, FALSE, FALSE,
                                                      SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_ALWAYS)
{
    checkedIfInt_ = false;
    modelIsInt_ = false;
    lastnode_ = -1;
    ncallsnode_ = 0;
    skipped_slack_ = 0;
    checked_slack_ = 0;
    nRuns_ = 0;
    nCuts_ = 0;
    start_gap_ = 0.0;
    n_fixed_ = 0;


    maxrounds_ = DEFAULT_MAXROUNDS;
    maxroundsroot_ = DEFAULT_MAXROUNDSROOT;
    maxsepacuts_ = DEFAULT_MAXSEPACUTS;
    maxsepacutsroot_ = DEFAULT_MAXSEPACUTSROOT;
    away_ = DEFAULT_AWAY;
    maxdynamism_ = DEFAULT_MAX_DYN;
    max_supp_abs_ = DEFAULT_MAX_SUPP_ABS;
    max_supp_rel_ = DEFAULT_MAX_SUPP_REL;
}

