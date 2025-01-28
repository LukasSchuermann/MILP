//
// Created by lukasschuermann on 18.06.24.
//
#include "iostream"
#include "event_setobjlimit.h"

const auto EVENT = SCIP_EVENTTYPE_FIRSTLPSOLVED;

EventhdlrObjLimit::EventhdlrObjLimit(SCIP *scip, double obj)
        : ObjEventhdlr(scip, "objlimit", "event handler for setting objective limit") {
    primal_bound_ = obj;
}

/** destructor of event handler to free user data (called when SCIP is exiting) */
SCIP_DECL_EVENTFREE(EventhdlrObjLimit::scip_free) {  /*lint --e{715}*/
    return SCIP_OKAY;
}


/** initialization method of event handler (called after problem was transformed) */
SCIP_DECL_EVENTINIT(EventhdlrObjLimit::scip_init) {  /*lint --e{715}*/
    return SCIP_OKAY;
}


/** deinitialization method of event handler (called before transformed problem is freed) */
SCIP_DECL_EVENTEXIT(EventhdlrObjLimit::scip_exit) {  /*lint --e{715}*/
    return SCIP_OKAY;
}


/** solving process initialization method of event handler (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The event handler may use this call to initialize its branch and bound specific data.
 *
 */
SCIP_DECL_EVENTINITSOL(EventhdlrObjLimit::scip_initsol) {
    SCIP_CALL(SCIPcatchEvent(scip, EVENT, eventhdlr, nullptr, nullptr));

    return SCIP_OKAY;
}


/** solving process deinitialization method of event handler (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The event handler should use this call to clean up its branch and bound data.
 */
SCIP_DECL_EVENTEXITSOL(EventhdlrObjLimit::scip_exitsol) {
    return SCIP_OKAY;
}


/** frees specific constraint data */
SCIP_DECL_EVENTDELETE(EventhdlrObjLimit::scip_delete) {  /*lint --e{715}*/
    return SCIP_OKAY;
}


/** execution method of event handler
 *
 *  Processes the event. The method is called every time an event occurs, for which the event handler
 *  is responsible. Event handlers may declare themselves responsible for events by calling the
 *  corresponding SCIPcatch...() method. This method creates an event filter object to point to the
 *  given event handler and event data.
 */
SCIP_DECL_EVENTEXEC(EventhdlrObjLimit::scip_exec) {
    if(SCIPgetDepth(scip) == 0)
    {
        SCIPsetObjlimit(scip, primal_bound_);
    }
    return SCIP_OKAY;
}
