//
// Created by Lukas Schuermann
//

#include <iostream>

#include <algorithm>
#include "scip/cons_linear.h"
#include "ConshdlrRCFC.h"
#include "objscip/objscip.h"
#include "scip/scip_probing.h"

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
SCIP_RETCODE probingNoMove(
        SCIP*                   scip,
        ConshdlrRCFC*           obj,
        SCIP_Var*               var,
        int                     index,
        double                  gap_abs,
        double                  varLPVal,
        bool*                   success,
        double                  bound,
        double                  maxSteep
){
    double steepness = gap_abs/(varLPVal - bound);
    if(abs(steepness) > maxSteep * obj->factor_){
        return SCIP_OKAY;
    }

    SCIP_Bool cutoff, lperror;

    int k;
    double local_lpVal = varLPVal;

    SCIP_LPI* lpi;
    SCIP_CALL( SCIPgetLPI(scip, &lpi) );

    int nIters = 0;
    SCIPstartProbing(scip);
    SCIPchgVarObjProbing(scip, var, SCIPvarGetObj(var) + steepness);
    k = obj->iter_init_;

    /* as long as we did not leave old opt */
    while (SCIPisEQ(scip, varLPVal, local_lpVal))
    {
        if(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL)
            break;

        SCIPsolveProbingLP(scip, k, &lperror, &cutoff);
        nIters += k;
        k = MAX((int) (obj->iter_perc_ * nIters), 1);
        local_lpVal = SCIPlpiColGetNewLPvalRCFC(lpi, index);
        if(nIters > 500){
            break;
        }
    }

    /* if x^* is still optimal, fix variable */
    if(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL && SCIPisEQ(scip, varLPVal, local_lpVal)){
        *success = true;
    }

    SCIPendProbing(scip);

    /* update lowest fail */
    if(!*success)
    {
        assert(SCIPisPositive(scip, steepness));
        obj->lowest_fail_down_ = std::min(obj->lowest_fail_down_, abs(steepness));
    }

    return SCIP_OKAY;

}

static
SCIP_RETCODE exactBranching(
        SCIP*                   scip,
        SCIP_Var*               var,
        bool*                   success,
        double                  bound,
        CutDir                  cutDir
){
    SCIP_Bool cutoff, lperror;

    SCIPstartProbing(scip);
    if(cutDir == DOWNWARDS){
        SCIPchgVarUbProbing(scip, var, bound);
    }else{
        SCIPchgVarLbProbing(scip, var, bound);
    }
    SCIPsolveProbingLP(scip, -1, &lperror, &cutoff);


    std::cout << SCIPgetLPObjval(scip) << " vs " << SCIPgetCutoffbound(scip) << " ";
    std::cout << SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) << std::endl;
    if(SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)))
    {
        assert(cutoff);
        assert(success != nullptr);
        *success = true;
    }

    SCIPendProbing(scip);
    return SCIP_OKAY;
}
static
SCIP_RETCODE RCFC_Extended(
        SCIP*                   scip,
        ConshdlrRCFC*           obj,
        SCIP_Col*               col,
        int                     index,
        double                  gap_abs,
        double                  varLPVal,
        bool*                   success,
        double                  bound,
        CutDir                  cutDir,
        double                  maxSteep
){
    SCIP_Bool cutoff, lperror;
    SCIP_Var* var = SCIPcolGetVar(col);
    SCIP_LPI* lpi;
    bool tooSteep = false;
    int k;
    int nLPIterations;

    double local_gap = gap_abs;
    double local_lpVal = varLPVal;
    double origObj = SCIPvarGetObj(var);

    double orig_gap = local_gap;

    double oldValue;

    double steepness = gap_abs / (varLPVal - bound);
    double init_steepness = steepness;

    SCIP_CALL( SCIPgetLPI(scip, &lpi) );

    /* enter probing state */
    SCIPstartProbing(scip);
    /* move from one optimum to the optimum for the next objective function */
    while (SCIPisSumGT(scip, cutDir * (local_lpVal - bound), 0.0)) {
        /* if we stay at the current x^* we can fix the variable */
        oldValue = local_lpVal;
        /* if objective change is too steep --> abort (no success expected) */
        if(abs(steepness) > maxSteep * obj->factor_){
            *success = false;
            break;
        }

        /* do a few simplex iterations at a time */
        k = obj->iter_init_;
        nLPIterations = 0;

        /* rotate objective function */
        SCIPchgVarObjProbing(scip, var, origObj + steepness);

        /* First iteration */
        SCIPsolveProbingLP(scip, k, &lperror, &cutoff);
        local_lpVal = SCIPlpiColGetNewLPvalRCFC(lpi, index);
        /* Abort Probing once x_i = 0 */
        while (SCIPisSumGT(scip, cutDir * (local_lpVal - bound), 0.0)) {
            /* calculate the gap corresponding to the original objective function */
            orig_gap = SCIPgetCutoffbound(scip) - (SCIPgetLPObjval(scip) - steepness * local_lpVal);

            /* if optimal was found, abort loop */
            if(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL){
                break;
            }

            /* abort early if objective function is too steep relative to lp value of x_i */
            if(abs(orig_gap/(local_lpVal - bound)) > maxSteep * obj->factor_){
                tooSteep = true;
                break;
            }

            /* no optimum was found -> increase iteration number and continue solving */
            nLPIterations += k;
            k = MAX((int) (obj->iter_perc_ * nLPIterations), obj->iter_init_);
            SCIPsolveProbingLP(scip, k, &lperror, &cutoff);

            local_lpVal = SCIPlpiColGetNewLPvalRCFC(lpi, index);
            /* safety measure */
            if(nLPIterations > 100000){
                tooSteep = true;
                break;
            }
        }
        if(SCIPisSumLE(scip, cutDir * (local_lpVal - bound), 0) || tooSteep)
            break;
        assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);
        assert(SCIPvarGetLPSol(var) == SCIPlpiColGetNewLPvalRCFC(lpi, index));

        /* did not move after objective rotation -> can fix variable */
        if (SCIPisEQ(scip, local_lpVal, oldValue)){
            *success = true;
            break;
        }

        /* calculate new objective shift/steepness */
        local_gap = SCIPgetCutoffbound(scip) - (SCIPgetLPObjval(scip) - steepness * local_lpVal);
        assert(SCIPisPositive(scip, local_gap));
        steepness = (local_gap / (local_lpVal - bound));
    }

    SCIPendProbing(scip);

    /* update lowest fail */
    if(!*success)
    {
        if(cutDir == DOWNWARDS){
            assert(SCIPisPositive(scip, init_steepness));
            obj->lowest_fail_down_ = std::min(obj->lowest_fail_down_, abs(init_steepness));

        }else{
            assert(SCIPisNegative(scip, init_steepness));
            obj->lowest_fail_up_ = std::min(obj->lowest_fail_up_, abs(init_steepness));
        }
    }

    return SCIP_OKAY;
}

static
SCIP_RETCODE checkFixing(
        SCIP*                   scip,
        ConshdlrRCFC*           obj,
        std::vector<int>        &checkLPIcols,
        SCIP_COL**              lpicols,
        SCIP_RESULT*            result,
        bool                    isLocal
){
    std::vector<std::pair<SCIP_Var*, double>> fixableVars_ub;
    std::vector<std::pair<SCIP_Var*, double>> fixableVars_lb;
    int nFixable = 0;

    double varLPVal;
    bool success;

    /* relax steepness condition after last call */
    obj->lowest_fail_down_ = obj->lowest_fail_down_ * obj->fail_reset_;
    obj->lowest_fail_up_ = obj->lowest_fail_up_ * obj->fail_reset_;

    obj->nCalls_++;

    int nChecked = 0;
    for(auto index : checkLPIcols)
    {
        double bound_down;
        double bound_up;
        SCIP_Col* col = lpicols[index];
        SCIP_Var* var = SCIPcolGetVar(col);

        /* shifted LP value */
        varLPVal = SCIPcolGetPrimsol(col);

        success = false;
        double gap_abs = SCIPgetCutoffbound(scip) - SCIPgetLPObjval(scip);

        /** early abort conditions */
        /* found maximum number of fixable variables */
        if(nFixable > obj->max_nFixVar_rel_ * SCIPgetNVars(scip) || (isLocal && nFixable > obj->max_nFix_) ||
                (!isLocal && nFixable > obj->max_nFix_root_)){
            break;
        }
        /* maximum number of tries without any success */
        if(nFixable == 0 && isLocal && obj->max_noFind_ < nChecked){
            break;
        }
        /* maximum number of tries */
        if(isLocal && nChecked > obj->max_nTries_){
            break;
        }
        /* in case the LP value changed during probing */
        if((SCIPfeasFrac(scip, SCIPcolGetPrimsol(col)) < obj->away_) ||
        (SCIPfeasFrac(scip, SCIPcolGetPrimsol(col)) > 1.0 - obj->away_)){
            continue;
        }
        if(!isLocal && nChecked > (int) checkLPIcols.size() / 4 && nFixable == 0){
            break;
        }

        nChecked++;
        obj->nCheckedVars_++;
        bound_down = std::ceil(varLPVal) - 1;

        if(SCIPisLE(scip, SCIPvarGetLbLocal(var), bound_down)){
            double maxSteep = obj->lowest_fail_down_;
            /** call probing method */
            if(obj->probingMode_ == PROBING_HEURISTIC){
                SCIP_CALL(probingNoMove(scip, obj, var, index, gap_abs, varLPVal, &success, bound_down, maxSteep) );
            }else if(obj->probingMode_ == PROBING_EXTENDED){
                SCIP_CALL( RCFC_Extended(scip, obj, col, index, gap_abs, varLPVal, &success, bound_down,
                                         DOWNWARDS, maxSteep) );
            }else{
                assert(obj->probingMode_ == PROBING_STRONGBRANCHING);
                SCIP_CALL( exactBranching(scip, var, &success, bound_down, DOWNWARDS) );
            }
        }


        if(success)
        {
            fixableVars_ub.emplace_back(var, bound_down + 1);
            nFixable++;

            if(!SCIPisZero(scip, SCIPfeasFrac(scip, varLPVal))){
                obj->nFixed_frac_++;
            }else{
                continue;
            }

            if(!obj->withCutoff_){
                continue;
            }


            bound_up = std::floor(varLPVal) + 1;

            if(SCIPisLT(scip, SCIPvarGetUbLocal(var), bound_up)){
                continue;
            }

            double maxSteep = obj->lowest_fail_up_;
            success = false;
            SCIP_CALL( RCFC_Extended(scip, obj, col, index, gap_abs, varLPVal, &success, bound_up,
                                     UPWARDS, maxSteep) );
            if(success)
            {
                *result = SCIP_CUTOFF;
                obj->nCutoffs_++;
                SCIPchgVarLb(scip, var, bound_down + 1);
                SCIPchgVarUb(scip, var, bound_up - 1);
                return SCIP_OKAY;
            }
        }else
        {
            if(!obj->bothDir_){
                continue;
            }
            if(SCIPisZero(scip, SCIPfeasFrac(scip, varLPVal))){
                continue;
            }

            bound_up = std::floor(varLPVal) + 1;
            if(SCIPisLT(scip, SCIPvarGetUbLocal(var), bound_up)){
                continue;
            }
            double maxSteep = obj->lowest_fail_up_;
            SCIP_CALL( RCFC_Extended(scip, obj, col, index, gap_abs, varLPVal, &success, bound_up,
                                     UPWARDS, maxSteep) );
            if(success)
            {
                if(!SCIPisZero(scip, SCIPfeasFrac(scip, varLPVal))){
                    obj->nFixed_frac_++;
                }
                nFixable++;
                fixableVars_lb.emplace_back(var, bound_up - 1);
            }
        }
        if(SCIPgetSolvingTime(scip) > obj->time_limit_ - 10)
            break;
    }

    if(!fixableVars_ub.empty() || !fixableVars_lb.empty())
    {
        obj->nSuccess_++;
        obj->nFixed_ += nFixable;
        if(SCIPgetDepth(scip) == 0){
            obj->nFixed_Root_ += nFixable;
        }

        if(*result != SCIP_CUTOFF)
            *result = SCIP_REDUCEDDOM;
        for(auto& pair : fixableVars_ub)
        {
            SCIPchgVarLb(scip, pair.first, pair.second);
        }
        for(auto& pair : fixableVars_lb)
        {
            SCIPchgVarUb(scip, pair.first, pair.second);
        }
    }

    return SCIP_OKAY;
}

/**@name Callback methods
 *
 * @{
 */

/** constraint enforcing method of constraint handler for LP solutions */
SCIP_DECL_CONSENFOLP(ConshdlrRCFC::scip_enfolp)
{
    if(SCIPinProbing(scip))
        return SCIP_OKAY;
    if(SCIPinDive(scip))
        return SCIP_OKAY;
    if(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_NOTSOLVED || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE)
        return SCIP_OKAY;
    /* max depth reached */
    if(maxDepth_ >= 0 && SCIPgetDepth(scip) > maxDepth_)
        return SCIP_OKAY;
    if(SCIPgetDepth(scip) > 0 && SCIPgetDepth(scip) % freq_ != 0)
        return SCIP_OKAY;
    if(SCIPisInfinity(scip, SCIPgetCutoffbound(scip)))
        return SCIP_OKAY;

    if(noRoot_ && SCIPgetDepth(scip) == 0){
        return SCIP_OKAY;
    }

    if(useTime_){
        if(SCIPconshdlrGetEnfoLPTime(conshdlr) / SCIPgetSolvingTime(scip) > timePerc_){
            return SCIP_OKAY;
        }
    }else{
        if(double(nLpIter_) / double(SCIPgetNLPIterations(scip)) > timePerc_){
            return SCIP_OKAY;
        }
    }


    /* do not call this method too often at the same branching node */
    if(lastNode_ == SCIPnodeGetNumber(SCIPgetFocusNode(scip)))
    {
        if((lastNode_ == 1 && nCalls_Node_ >= max_nCallsRoot_) ||
            (lastNode_ != 1 && nCalls_Node_ >= max_nCalls_))
        {
            nCalls_Node_ = 0;
            return SCIP_OKAY;
        }
        nCalls_Node_++;
    }else{
        nCalls_Node_ = 1;
    }

    /* SCIP can not handle it if we are in probing mode after time limit */
    if(SCIPgetSolvingTime(scip) > time_limit_ - 10)
        return SCIP_OKAY;

    lastNode_ = SCIPnodeGetNumber(SCIPgetFocusNode(scip));

    double gap_rel = abs((SCIPgetCutoffbound(scip) - SCIPgetLPObjval(scip)) / SCIPgetLPObjval(scip));

    if(gap_rel > validGap_)
        return SCIP_OKAY;
    assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

    /** Get all active LPI-columns */
    SCIP_COL** lpicols = SCIPlukGetlpiCols(scip);
    std::vector<int> checkLPIcols;

    /* collect candidates */
    for(int i = 0; i < SCIPlukGetNlpiCols(scip); i++)
    {
        SCIP_Col* col = lpicols[i];
        if(!SCIPvarIsIntegral(SCIPcolGetVar(col)))
            continue;
        if(SCIPisEQ(scip, SCIPcolGetLb(col), SCIPcolGetUb(col)))
            continue;
        if((SCIPfeasFrac(scip, SCIPcolGetPrimsol(col)) >= away_) && (SCIPfeasFrac(scip, SCIPcolGetPrimsol(col)) <= 1.0 - away_))
        {
            checkLPIcols.push_back(i);
        }
    }
    std::sort(checkLPIcols.begin(), checkLPIcols.end(), [lpicols](auto &left, auto &right){
        return SCIPcolGetPrimsol(lpicols[left]) - SCIPcolGetUb(lpicols[left]) > SCIPcolGetPrimsol(lpicols[right]) - SCIPcolGetUb(lpicols[left]);
    });

    *result = SCIP_FEASIBLE;

    auto oldIter = SCIPgetNLPIterations(scip);


    auto start_time_ = std::chrono::high_resolution_clock::now();
    if(!checkLPIcols.empty())
        SCIP_CALL(checkFixing(scip, this, checkLPIcols, lpicols, result, !(SCIPgetDepth(scip) == 0)));

    auto end_time_ = std::chrono::high_resolution_clock::now();
    propTime_ += std::chrono::duration_cast<std::chrono::microseconds>(end_time_ - start_time_).count();
    nLpIter_ += (SCIPgetNLPIterations(scip) - oldIter);

    return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(ConshdlrRCFC::scip_trans)
{   /*lint --e{715}*/
    SCIP_CONSDATA* sourcedata;
    SCIP_CONSDATA* targetdata;

    assert(conshdlr != nullptr);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), "RCFC") == 0);
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



ConshdlrRCFC::ConshdlrRCFC(SCIP *scip, int mode, bool withCutoff, bool bothDir, double timePerc, bool noRoot, bool useTime)
: ObjConshdlr(scip, "RCFC", "counterpart cuts conshdlr", 5000, 400, 0,
                                                      4, 0, 0, 0, TRUE, FALSE, FALSE,
                                                      SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_ALWAYS),
                                                      withCutoff_(withCutoff),
                                                      bothDir_(bothDir),
                                                      timePerc_(timePerc),
                                                      noRoot_(noRoot),
                                                      useTime_(useTime)
{
    /* statistics initialization */
    nLpIter_ = 0;
    nCutoffs_ = 0;
    nSuccess_ = 0;
    nFixed_ = 0;
    nFixed_frac_ = 0;
    nCalls_ = 0;
    nCheckedVars_ = 0;
    lastNode_ = -1;

    lowest_fail_down_ = SCIPinfinity(scip);
    lowest_fail_up_ = SCIPinfinity(scip);
    probingMode_ = mode;

    factor_ = 0.85;
    validGap_ = 1.0;
    max_nCalls_ = 1;
    max_nCallsRoot_ = 5;
    nCalls_Node_ = 0;
    max_nFixVar_rel_ = 0.1;
    max_nFix_ = 5;
    max_nFix_root_ = 15;
    iter_init_ = 5;
    iter_perc_ = 0.5;

    max_noFind_ = 5;
    max_nTries_ = 20;

    SCIPgetRealParam(scip, "limits/time", &time_limit_);

    away_ = 0.2;
}

