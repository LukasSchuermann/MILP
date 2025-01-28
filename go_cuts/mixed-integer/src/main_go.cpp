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

#include "ConshdlrGO.h"
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
        bool*               useCutting,
        int*                seed,
        double*             optObj,
        bool*               usedOptObj,
        std::string&        outFile,
        double*             scale_zo,
        double*             scale_gen
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
        if ( ! strcmp(argv[i], "+useCutting"))
        {
            *useCutting = true;
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
        }else if(! strcmp(argv[i], "+scaleZO"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "+",1)))
            {
                fprintf(stderr, "Missing scale value. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *scale_zo = atof(locstr);
            free(locstr);
        }else if(! strcmp(argv[i], "+scaleGEN"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "+",1)))
            {
                fprintf(stderr, "Missing scale value. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *scale_gen = atof(locstr);
            free(locstr);
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
    auto* consobj = dynamic_cast<ConshdlrGO *>(SCIPfindObjConshdlr(scip, "CPC"));

    std::string s;
    /* read setup data */
    if(data_file.is_open())
    {
        data_file >> s;
        assert(s == "max_nofind_0");
        data_file >> s;
        consobj->max_noFind_high_ = stoi(s);

        data_file >> s;
        assert(s == "max_nofind_1");
        data_file >> s;
        consobj->max_noFind_low_ = stoi(s);

        data_file >> s;
        assert(s == "max_ntries_0");
        data_file >> s;
        consobj->max_nTries_high_ = stoi(s);

        data_file >> s;
        assert(s == "max_ntries_1");
        data_file >> s;
        consobj->max_nTries_low_ = stoi(s);

        data_file >> s;
        assert(s == "rel_gap_0");
        data_file >> s;
        consobj->rel_gap_high_ = stof(s);

        data_file >> s;
        assert(s == "rel_gap_1");
        data_file >> s;
        consobj->rel_gap_low_ = stof(s);

        data_file >> s;
        assert(s == "max_supp_abs");
        data_file >> s;
        consobj->max_supp_abs_ = stoi(s);

        data_file >> s;
        assert(s == "max_supp_rel");
        data_file >> s;
        consobj->max_supp_rel_ = stof(s);

        data_file >> s;
        assert(s == "sepa_freq");
        data_file >> s;
        SCIPsetIntParam(scip, "constraints/CPC/sepafreq", stoi(s));
        consobj->freq_ = stoi(s);

        data_file >> s;
        assert(s == "away");
        data_file >> s;
        consobj->away_ = stof(s);

        data_file >> s;
        assert(s == "maxrounds");
        data_file >> s;
        consobj->maxrounds_ = stoi(s);

        data_file >> s;
        assert(s == "maxroundsroot");
        data_file >> s;
        consobj->maxroundsroot_ = stoi(s);

        data_file >> s;
        assert(s == "maxsepacuts");
        data_file >> s;
        consobj->maxsepacuts_ = stoi(s);

        data_file >> s;
        assert(s == "maxsepacutsroot");
        data_file >> s;
        consobj->maxsepacutsroot_ = stoi(s);

        data_file >> s;
        assert(s == "root_supp_rel");
        data_file >> s;
        consobj->root_supp_rel_ = stof(s);

        data_file >> s;
        assert(s == "max_depth");
        data_file >> s;
        consobj->maxDepth_ = stoi(s);

        data_file >> s;
        assert(s == "descend_supp");
        data_file >> s;
        if ( s == "1" ){
            consobj->descending_supp_ = true;
        }

        data_file >> s;
        assert(s == "scale_zo");
        data_file >> s;
        consobj->scale_zo_ = stof(s);

        data_file >> s;
        assert(s == "scale_gen");
        data_file >> s;
        consobj->scale_gen_ = stof(s);
        std::cout << "Scales set to: " << consobj->scale_zo_ << " and: " << consobj->scale_gen_ << std::endl;
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

    /** change parameters */

    /* change display columns */
    SCIP_CALL( SCIPsetIntParam(*scip,"display/nfrac/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/curcols/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/cuts/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/lpobj/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/primalgap/active",0) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/gap/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/lpavgiterations/active",0) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/lpiterations/active",0) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/vars/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/conflicts/active",0) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/strongbranchs/active",0) );

    /* Branching */
    SCIP_CALL( SCIPsetIntParam(*scip,"display/nnodes/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/nodesleft/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/maxdepth/active",2) );

//    /* for column generation instances, disable restarts */
    SCIP_CALL( SCIPsetIntParam(*scip,"display/verblevel",4) );

    /* ACTIVATE THIS FOR ROOT NODE EXPERIMENT */
//    SCIP_CALL(SCIPsetLongintParam(*scip, "limits/nodes", 1));
//    SCIP_CALL (SCIPsetRealParam(*scip, "limits/gap", 0.0));
//    SCIP_CALL (SCIPsetIntParam(*scip, "separating/maxcutsroot", 2000000));
//    SCIP_CALL (SCIPsetIntParam(*scip, "separating/maxstallroundsroot", -1));
//    SCIP_CALL (SCIPsetRealParam(*scip, "separating/minefficacyroot", 0.0000001));
//    SCIP_CALL (SCIPsetRealParam(*scip, "separating/minefficacy", 0.0000001));

    SCIP_CALL( SCIPsetIntParam(*scip, "display/freq", 200));
    SCIP_CALL( SCIPsetRealParam(*scip, "limits/time", 3600));


    return SCIP_OKAY;
}


/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runGOcutModel(
    int                   argc,               /**< number of shell parameters */
    char**                argv                /**< array with shell parameters */
)
{
    SCIP* scip = nullptr;

    assert(argc >= 5);
    int seed = 0;
    double optObj = 0;
    bool usedOptObj = false;

    bool withCPC = false;

    char* setup_file = nullptr;
    char* instance = nullptr;
    std::string outFile;
    double scale_zo_ = 4;
    double scale_gen_ = 1;

    SCIP_CALL(readArguments(argc, argv, &instance, &setup_file, &withCPC, &seed, &optObj, &usedOptObj, outFile,
                            &scale_zo_, &scale_gen_));

    std::cout << instance << " " << (setup_file == nullptr ? "no-setup-file" : setup_file) << " ";
    std::cout << withCPC << " " << seed << " " << optObj << std::endl;

    SCIP_CALL(setUpScip(&scip, usedOptObj));

    if(withCPC)
    {
        SCIP_CALL(SCIPincludeObjConshdlr(scip, new ConshdlrGO(scip), true));
        if(setup_file != nullptr)
        {
            std::cout << "read setup file: " << setup_file << std::endl;
            readSetupFile(scip, setup_file);
        }
    }

    /* Seed */
    SCIPsetIntParam(scip, "randomization/randomseedshift", seed);
    /* Include optimal solution value (can also be any primal bound) */
    if(usedOptObj){
        SCIPincludeObjEventhdlr(scip, new EventhdlrObjLimit(scip, optObj), true);
    }

    SCIP_Result result;
    SCIPreadMps(scip, SCIPfindReader(scip, "mpsreader"), instance, &result, nullptr, nullptr, nullptr, nullptr,nullptr, nullptr);


    SCIP_CALL(SCIPsolve(scip));

    if(!outFile.empty()){

        auto* obj = dynamic_cast<ConshdlrGO *>(SCIPfindObjConshdlr(scip, "CPC"));
        SCIP_Conshdlr* CPC = SCIPfindConshdlr(scip, "CPC");

        std::cout << "OUTFILE: " << outFile << std::endl;
        std::ofstream sol_file;
        sol_file.open(outFile, std::ios_base::app);
        if(withCPC)
        {
            sol_file << "Seed;SolTime;nNodes;Gap;LPIter;nSepCalls;sepTime;nCuts;nCutsRoot" << std::endl;
        }else {
            sol_file << "Seed;SolTime;nNodes;Gap;LPIter" << std::endl;
        }
        sol_file << seed << ";" << SCIPgetSolvingTime(scip) << ";" << SCIPgetNNodes(scip) << ";" << SCIPgetGap(scip);
        sol_file << ";" << SCIPgetNLPIterations(scip);
        if(withCPC)
        {
            sol_file << ";" << SCIPconshdlrGetNSepaCalls(CPC) << ";" << SCIPconshdlrGetSepaTime(CPC);
            sol_file << ";" << obj->nCuts_ << ";" << obj->nCutsRoot_;
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

    retcode = runGOcutModel(argc, argv);
    if( retcode != SCIP_OKAY )
    {
      SCIPprintError(retcode);
      return -1;
    }

    return 0;
}

