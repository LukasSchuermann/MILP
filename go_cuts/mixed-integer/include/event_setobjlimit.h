//
// Created by lukasschuermann on 18.06.24.
//

#ifndef VRP_EVENT_SETOBJLIMIT_H
#define VRP_EVENT_SETOBJLIMIT_H
#include "scip/scip.h"
#include "objscip/objscip.h"

/** C++ wrapper object for event handlers */
class EventhdlrObjLimit : public scip::ObjEventhdlr {
public:
    double primal_bound_;
    /** default constructor */
    explicit EventhdlrObjLimit(SCIP *scip, double obj);

    /** destructor */
    ~EventhdlrObjLimit() override = default;

    /** destructor of event handler to free user data (called when SCIP is exiting) */
    virtual SCIP_DECL_EVENTFREE(scip_free);

    /** initialization method of event handler (called after problem was transformed) */
    virtual SCIP_DECL_EVENTINIT(scip_init);

    /** deinitialization method of event handler (called before transformed problem is freed) */
    virtual SCIP_DECL_EVENTEXIT(scip_exit);

    /** solving process initialization method of event handler (called when branch and bound process is about to begin)
     *
     *  This method is called when the presolving was finished and the branch and bound process is about to begin.
     *  The event handler may use this call to initialize its branch and bound specific data.
     *
     */
    virtual SCIP_DECL_EVENTINITSOL(scip_initsol);

    /** solving process deinitialization method of event handler (called before branch and bound process data is freed)
     *
     *  This method is called before the branch and bound process is freed.
     *  The event handler should use this call to clean up its branch and bound data.
     */
    virtual SCIP_DECL_EVENTEXITSOL(scip_exitsol);

    /** frees specific constraint data */
    virtual SCIP_DECL_EVENTDELETE(scip_delete);

    /** execution method of event handler
     *
     *  Processes the event. The method is called every time an event occurs, for which the event handler
     *  is responsible. Event handlers may declare themselves resposible for events by calling the
     *  corresponding SCIPcatch...() method. This method creates an event filter object to point to the
     *  given event handler and event data.
     */
    virtual SCIP_DECL_EVENTEXEC(scip_exec);


}; /*lint !e1712*/

#endif //VRP_EVENT_SETOBJLIMIT_H
