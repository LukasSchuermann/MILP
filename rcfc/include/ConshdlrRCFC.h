//
// Created by lukasschuermann on 27.09.23.
//

#ifndef VRP_ConshdlrRCFC_H
#define VRP_ConshdlrRCFC_H

#include "objscip/objscip.h"
#include "scip/scip.h"
#include "vector"
#include "chrono"

using namespace scip;

enum probingMode {
    PROBING_HEURISTIC = 0,
    PROBING_EXTENDED = 1,
    PROBING_STRONGBRANCHING = 2
};

class ConshdlrRCFC : public ObjConshdlr
{
public:
    long long int lastNode_;

    /* settings */
    bool withCutoff_;
    bool bothDir_;
    int max_noFind_;
    int max_nTries_;
    double validGap_;
    double iter_perc_;
    int iter_init_;
    double lowest_fail_down_;
    double lowest_fail_up_;
    double factor_;
    int probingMode_;
    double away_;
    double time_limit_ = 0;
    int maxDepth_ = -1;
    double fail_reset_ = 1.5;
    double max_nFixVar_rel_;
    int max_nFix_;
    int max_nFix_root_;
    int max_nCalls_;
    int max_nCallsRoot_;
    int nCalls_Node_;
    int freq_ = 1;
    double timePerc_ = 0.1;
    long long int maxTimeIter_ = 0;
    bool noRoot_ = false;
    bool firstRun_ = true;
    bool useTime_ = false;

    /* statistics */
    int nFixed_Root_ = 0;
    int nSuccess_;
    int nFixed_;
    int nFixed_frac_;
    int nCalls_;
    int nCheckedVars_;
    int nCutoffs_;
    long long int nLpIter_;
    long long int propTime_ = 0;

    std::vector<double> lowestFails_down_;
    std::vector<int> nFails_down_;
    std::vector<double> lowestFails_up_;
    std::vector<int> nFails_up_;


    /** default constructor */
    explicit ConshdlrRCFC(
            SCIP*       scip,
            int         mode,
            bool        withCutoff,
            bool        bothDir,
            double      timePerc,
            bool        noRoot,
            bool        useTime
    );

    /** destructor */
    ~ConshdlrRCFC() override= default;

    /** frees specific constraint data */
    virtual SCIP_DECL_CONSDELETE(scip_delete){
        return SCIP_OKAY;
    }

    /** transforms constraint data into data belonging to the transformed problem */
    virtual SCIP_DECL_CONSTRANS(scip_trans);

    /** separation method of constraint handler for LP solution */
    virtual SCIP_DECL_CONSSEPALP(scip_sepalp){
        return SCIP_OKAY;
    }

    /** constraint enforcing method of constraint handler for LP solutions */
    virtual SCIP_DECL_CONSENFOLP(scip_enfolp);

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

    virtual SCIP_DECL_CONSINIT(scip_init){
        lowestFails_down_ = std::vector<double>(SCIPgetNVars(scip), SCIPinfinity(scip));
        nFails_down_ = std::vector<int>(SCIPgetNVars(scip), 0);
        lowestFails_up_ = std::vector<double>(SCIPgetNVars(scip), SCIPinfinity(scip));
        nFails_up_ = std::vector<int>(SCIPgetNVars(scip), 0);
        return SCIP_OKAY;
    }
};


#endif //VRP_ConshdlrRCFC_H
