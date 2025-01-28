//
// Created by Lukas Schuermann
//

#ifndef VRP_ConshdlrCPC_H
#define VRP_ConshdlrCPC_H

#include "objscip/objscip.h"
#include "scip/cons_linear.h"
#include "scip/scip.h"
#include "vector"

class ConshdlrGO : public scip::ObjConshdlr
{
public:
    std::vector<double> dualvals_;
    int nRuns_;
    double start_gap_;
    std::vector<double> fixed_gaps_;
    int nCuts_ = 0;
    int nCutsRoot_ = 0;
    int             max_supp_abs_;         /**< maximum cut support - absolute value in the formula */
    SCIP_Real       max_supp_rel_;         /**< maximum cut support - relative value in the formula */
    int max_noFind_high_;
    int max_noFind_low_;
    int max_nTries_high_;
    int max_nTries_low_;
    int nFails_ = 0;
    double rel_gap_high_;
    double rel_gap_low_;
    SCIP_Real       away_;               /**< minimal fractionality of a basis variable in order to try CPC cut */

    int             maxrounds_;          /**< maximal number of GMI separation rounds per node (-1: unlimited) */
    int             maxroundsroot_;      /**< maximal number of GMI separation rounds in the root node (-1: unlimited) */
    int             maxsepacuts_;        /**< maximal number of GMI cuts separated per separation round */
    int             maxsepacutsroot_;    /**< maximal number of GMI cuts separated per separation round in root node */
    int             maxDepth_ = -1;
    bool descending_supp_ = false;
    int maxsupp_ = -1;
    int freq_ = -1;
    double root_supp_rel_ = 1;
    /** default constructor */
    explicit ConshdlrGO(
            SCIP*       scip
    );

    /** destructor */
    ~ConshdlrGO() override= default;

    /** frees specific constraint data */
    virtual SCIP_DECL_CONSDELETE(scip_delete){
        return SCIP_OKAY;
    };

    /** transforms constraint data into data belonging to the transformed problem */
    virtual SCIP_DECL_CONSTRANS(scip_trans);

    /** separation method of constraint handler for LP solution */
    virtual SCIP_DECL_CONSSEPALP(scip_sepalp);

    /** constraint enforcing method of constraint handler for LP solutions */
    virtual SCIP_DECL_CONSENFOLP(scip_enfolp){
        return SCIP_OKAY;
    }

    /** constraint enforcing method of constraint handler for pseudo solutions */
    virtual SCIP_DECL_CONSENFOPS(scip_enfops){
            return SCIP_OKAY;
    }

    /** feasibility check method of constraint handler for integral solutions */
    virtual SCIP_DECL_CONSCHECK(scip_check){
            return SCIP_OKAY;
    }

    /** variable rounding lock method of constraint handler */
    virtual SCIP_DECL_CONSLOCK(scip_lock){
            return SCIP_OKAY;
    }

    virtual SCIP_DECL_CONSCOPY(scip_copy){
        return SCIP_OKAY;
    }


    /** checks if model contains only variables with integer values in feasible solutions */
    bool isModelInteger(SCIP* scip);

    SCIP_Real getMaxDynamism() const{return maxdynamism_;}
    SCIP_Real getAway() const{return away_;}

private:
    bool            checkedIfInt_;
    long long int   lastnode_;
    int             ncallsnode_;
    SCIP_Real       maxdynamism_;        /**< maximal valid range max(|weights|)/min(|weights|) of cut coefficients */
};


#endif //VRP_ConshdlrCPC_H
