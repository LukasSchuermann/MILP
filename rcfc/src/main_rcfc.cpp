/**@file   main_vrp.cpp
 * @brief  Main program file for "column generation for VRP-HTW"
 * @author Lukas Schürmann
 */

#include <stdio.h>
#include <time.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"
#include "objscip/objscip.h"

#include "ConshdlrRCFC.h"
#include "event_setobjlimit.h"

#include "iostream"
#include "fstream"
#include <iomanip>

static
SCIP_RETCODE readArguments(
        int                 argc,               /**< number of shell parameters */
        char**              argv,               /**< array with shell parameters */
        char**              input_file,
        char**              settings_file,
        bool*               useRCFC,
        int*                seed,
        double*             optObj,
        bool*               usedOptObj,
        int*                mode,
        std::string&        outFile,
        bool*               withCutoff,
        bool*               bothDir,
        double*             timePerc,
        bool*               noOpt,
        bool*               noRoot,
        bool*               useTime
)
{
    char usage[SCIP_MAXSTRLEN];
    int status;
    char* locstr;

    assert( argc >= 1 );
    assert( argv != nullptr );
    assert( input_file != nullptr );

    /* init usage text */
    status = snprintf(usage, SCIP_MAXSTRLEN - 1, "usage: %s <path of inputfile> ", argv[0]);
    assert( 0 <= status && status < SCIP_MAXSTRLEN );

    /* init arguments */
    *input_file = nullptr;

    /* mandatory argument: inputfile */
    *input_file = argv[1];
    if ( *input_file == nullptr )
    {
        fprintf(stderr, "No path of data supplied.\n");
        fprintf(stderr, "%s\n", usage);
        return SCIP_ERROR;
    }
    /* check for branching strategy */
    for (int i = 2; i < argc; i++)
    {
        if ( ! strcmp(argv[i], "+useRCFC"))
        {
            *useRCFC = true;
        }else if(! strcmp(argv[i], "+objValue"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "+",1)))
            {
                fprintf(stderr, "Missing obj value. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *optObj = atof(locstr);
            *usedOptObj = true;
            free(locstr);
        }else if ( ! strcmp(argv[i], "+seed"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "+",1)))
            {
                fprintf(stderr, "Missing seed number. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *seed = atoi(locstr);
            if(*seed < 0)
            {
                fprintf(stderr, "Invalid seed! -> Choose seed >= 0!");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            free(locstr);
        }else if ( ! strcmp(argv[i], "+settings"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "+",1)))
            {
                fprintf(stderr, "Missing settings path. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            *settings_file = argv[i];
            if ( *settings_file == nullptr )
            {
                fprintf(stderr, "No path of settings supplied.\n");
                fprintf(stderr, "%s\n", usage);
                return SCIP_ERROR;
            }
        }else if(! strcmp(argv[i], "+outFile"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "+",1)))
            {
                fprintf(stderr, "Missing output file. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            outFile = argv[i];
        }else if ( ! strcmp(argv[i], "+mode"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "+",1)))
            {
                fprintf(stderr, "Missing mode number. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *mode = atoi(locstr);
            if(*mode < 0 || *mode > 2)
            {
                fprintf(stderr, "Invalid mode! -> Choose value from {0, 1, 2}!");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            free(locstr);
        }else if ( ! strcmp(argv[i], "+withCutoff"))
        {
            *withCutoff = true;
        }else if ( ! strcmp(argv[i], "+bothDir"))
        {
            *bothDir = true;
        }else if ( ! strcmp(argv[i], "+timePerc"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "+",1)))
            {
                fprintf(stderr, "Missing max time percentage. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *timePerc = atof(locstr);
            if(*timePerc <= 0 || *timePerc >= 1)
            {
                fprintf(stderr, "Invalid max time percentage! -> Choose value from (0,1)!");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            free(locstr);
        }else if ( ! strcmp(argv[i], "+noOpt")){
            *noOpt = true;
        }else if ( ! strcmp(argv[i], "+noRoot")){
            *noRoot = true;
        }else if ( ! strcmp(argv[i], "+useTime")){
            *useTime = true;
        }
    }

    return SCIP_OKAY;
}

void readSetupFile(
        SCIP*   scip,
        char*   file
){
    std::ifstream data_file;
    data_file.open(file);
    auto* consobj = dynamic_cast<ConshdlrRCFC*>(SCIPfindObjConshdlr(scip, "RCFC"));

    std::string s;
    /* read setup data */
    std::cout << "Prop settings: " << std::endl;
    if(data_file.is_open())
    {
        data_file >> s;
        assert(s == "max_nofind");
        data_file >> s;
        consobj->max_noFind_ = stoi(s);

        data_file >> s;
        assert(s == "max_ntries");
        data_file >> s;
        consobj->max_nTries_ = stoi(s);

        std::cout << "\tmax_nofind: " << consobj->max_noFind_ << ", max_ntries: " << consobj->max_nTries_ << std::endl;

        data_file >> s;
        assert(s == "valid_gap");
        data_file >> s;
        consobj->validGap_ = stof(s);

        data_file >> s;
        assert(s == "sepa_freq");
        data_file >> s;
        SCIPsetIntParam(scip, "constraints/RCFC/sepafreq", stoi(s));
        consobj->freq_ = stoi(s);

        std::cout << "\tvalid_gap: " << consobj->validGap_ << ", freq: " << stoi(s) << std::endl;

        data_file >> s;
        assert(s == "away");
        data_file >> s;
        consobj->away_ = stof(s);

        data_file >> s;
        assert(s == "factor");
        data_file >> s;
        consobj->factor_ = stof(s);

        data_file >> s;
        assert(s == "reset");
        data_file >> s;
        consobj->fail_reset_ = stof(s);

        std::cout << "\taway: " << consobj->away_ << ", factor: " << consobj->factor_ << ", reset: " << consobj->fail_reset_ << std::endl;

        data_file >> s;
        assert(s == "iter_init");
        data_file >> s;
        consobj->iter_init_ = stoi(s);

        data_file >> s;
        assert(s == "iter_perc");
        data_file >> s;
        consobj->iter_perc_ = stof(s);

        std::cout << "\titer_init: " << consobj->iter_init_ << ", iter_perc: " << consobj->iter_perc_ << std::endl;

        data_file >> s;
        assert(s == "maxrounds");
        data_file >> s;
        consobj->max_nCalls_ = stoi(s);

        data_file >> s;
        assert(s == "maxroundsroot");
        data_file >> s;
        consobj->max_nCallsRoot_ = stoi(s);

        data_file >> s;
        assert(s == "maxfix_rel");
        data_file >> s;
        consobj->max_nFixVar_rel_ = stof(s);

        data_file >> s;
        assert(s == "maxfix");
        data_file >> s;
        consobj->max_nFix_ = stoi(s);

        data_file >> s;
        assert(s == "maxfixroot");
        data_file >> s;
        consobj->max_nFix_root_ = stoi(s);

        data_file >> s;
        assert(s == "max_depth");
        data_file >> s;
        consobj->maxDepth_ = stoi(s);

        std::cout << "\tmax_depth: " << consobj->maxDepth_ << ", maxRounds: " << consobj->max_nCalls_ << ", ";
        std::cout << consobj->max_nCallsRoot_ << std::endl;
        std::cout << "\tmaxfix: " << consobj->max_nFix_ << ", rel: " << consobj->max_nFixVar_rel_ << ", root: ";
        std::cout << consobj->max_nFix_root_ << std::endl;
    }
}

static
SCIP_RETCODE setUpScip(
        SCIP**              scip,
        bool                usedOptObj
) {
    /* initialize SCIP */
    SCIP_CALL(SCIPcreate(scip));

    /* include default SCIP plugins */
    SCIP_CALL( SCIPincludeDefaultPlugins(*scip) );

    /* turn off all heuristics if opt objective value was included */
    if(usedOptObj){
        SCIP_CALL( SCIPsetHeuristics(*scip, SCIP_PARAMSETTING_OFF, TRUE) );
    }

    /* change display columns */
    SCIP_CALL( SCIPsetIntParam(*scip,"display/nfrac/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/curcols/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/cuts/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/lpobj/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/primalgap/active",0) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/gap/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/lpavgiterations/active",0) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/lpiterations/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/vars/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/conflicts/active",0) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/strongbranchs/active",0) );
    /* Branching columns */
    SCIP_CALL( SCIPsetIntParam(*scip,"display/nnodes/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/nodesleft/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/maxdepth/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/freq", 200));

    SCIP_CALL( SCIPsetIntParam(*scip,"display/verblevel",4) );

    /* activate reduced costs propagator */
//    SCIP_CALL(SCIPsetBoolParam(*scip, "propagating/redcost/force", true));

    SCIP_CALL( SCIPsetRealParam(*scip, "limits/time", 3600));


    return SCIP_OKAY;
}


/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runSCIPModel(
    int                   argc,               /**< number of shell parameters */
    char**                argv                /**< array with shell parameters */
)
{
    SCIP* scip = nullptr;

    int seed = 0;
    double optObj = 0;
    bool usedOptObj = false;
    bool withRCFC = false;
    bool withCutoff = false;
    bool bothDir = false;
    double timePerc = 0.075;
    int mode = 1;

    char* setup_file = nullptr;
    char* instance = nullptr;
    std::string outFile;
    bool noOpt = false;
    bool noRoot = false;
    bool useTime = false;

    /* read user arguments */
    SCIP_CALL(readArguments(argc, argv, &instance, &setup_file, &withRCFC, &seed, &optObj, &usedOptObj, &mode,
                            outFile, &withCutoff, &bothDir, &timePerc, &noOpt, &noRoot, &useTime));

    std::cout << instance << " " << (setup_file == nullptr ? "no-setup-file" : setup_file) << " useRCRC: ";
    std::cout << withRCFC << " seed: " << seed << " optObj: " << optObj << " mode: " << mode << std::endl;

    if(noOpt){
        usedOptObj = false;
    }
    setUpScip(&scip, usedOptObj);

    if(withRCFC){
        std::cout << "Activated RCFC method on mode " << mode << std::endl;
    }
    /* Seed */
    SCIPsetIntParam(scip, "randomization/randomseedshift", seed);

    /* Load optimal value as cutoff bound */
    if(usedOptObj){
        SCIPincludeObjEventhdlr(scip, new EventhdlrObjLimit(scip, optObj), true);
    }

    SCIP_Result result;
    SCIPreadMps(scip, SCIPfindReader(scip, "mpsreader"), instance, &result, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);

    /* include the constraint handler and its settings */
    if(withRCFC)
    {
        SCIP_CALL(SCIPincludeObjConshdlr(scip, new ConshdlrRCFC(scip, mode, withCutoff, bothDir, timePerc, noRoot, useTime), true));
        if(setup_file != nullptr)
        {
            std::cout << "read setup file: " << setup_file << std::endl;
            readSetupFile(scip, setup_file);
        }
    }

    SCIP_CALL(SCIPsolve(scip));

    if(!outFile.empty()){
        auto* obj = dynamic_cast<ConshdlrRCFC*>(SCIPfindObjConshdlr(scip, "RCFC"));
        SCIP_Conshdlr* RCFC = SCIPfindConshdlr(scip, "RCFC");

        std::cout << "OUTFILE: " << outFile << std::endl;
        std::ofstream sol_file;
        sol_file.open(outFile, std::ios_base::app);
        if(withRCFC)
        {
            sol_file << "Seed;SolTime;nNodes;Gap;LPIter;nSepCalls;sepTime;nFixed;nFixedFrac;propLPIter;nCutOff;nFixedRoot;propTime" << std::endl;
        }else {
            sol_file << "Seed;SolTime;nNodes;Gap;LPIter" << std::endl;
        }
        sol_file << seed << ";" << SCIPgetSolvingTime(scip) << ";" << SCIPgetNNodes(scip) << ";" << SCIPgetGap(scip);
        sol_file << ";" << SCIPgetNLPIterations(scip);
        if(withRCFC)
        {
            sol_file << ";" << SCIPconshdlrGetNEnfoLPCalls(RCFC) << ";" << SCIPconshdlrGetEnfoLPTime(RCFC) << ";";
            sol_file << obj->nFixed_ << ";" << obj->nFixed_frac_ << ";" << obj->nLpIter_ << ";" << obj->nCutoffs_;
            sol_file << ";" << obj->nFixed_Root_ << ";" << obj->propTime_ / 1000;
        }
        sol_file << std::endl;
    }

    /********************
    * Deinitialization *
    ********************/

    SCIP_CALL( SCIPfree(&scip) );

    return SCIP_OKAY;
}

int
main(
    int                        argc,
    char**                     argv
)
{
    SCIP_RETCODE retcode;

    retcode = runSCIPModel(argc, argv);
    if( retcode != SCIP_OKAY )
    {
      SCIPprintError(retcode);
      return -1;
    }

    return 0;
}

