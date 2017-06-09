/********************************************************************/
/*                                                                  */
/* YAMAS - Yet Another Meta-Analysis Software                       */
/*                                                                  */
/********************************************************************/

// $Rev: 932 $
// $LastChangedDate: 2013-05-28 14:18:57 +0200 (Di, 28 Mai 2013) $
// $Author: meesters $

#define MAIN

#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <stdlib.h> // for exit values
#include <sys/stat.h>
#include <signal.h>
using namespace std;

// include boost's program_options
#include <boost/program_options.hpp>
#include <boost/tuple/tuple.hpp> 
#include <boost/tuple/tuple_io.hpp>
#include <boost/algorithm/string.hpp>
using namespace boost;
namespace po = boost::program_options;

// include OpenMP to be able to set the number of threads globally
#include <omp.h>

// include own software parts
#include "../include/ya_version.hpp"
#include "../include/ya_utils.hpp"
#include "../include/ya_project.hpp"
#include "../include/ya_algorithms.hpp"
#include "../include/ya_interface.hpp"
#include "../include/ya_assoc.hpp"
#include "../include/ya_errors.hpp"
#include "../include/ya_global.hpp"
#include "../include/ya_configwrite.hpp"
#ifdef SANDBOX
#include "../sandbox/include/ya_diagnostic.hpp"
#endif
#ifdef TIMER
#include "../include/timer.hpp"
#endif
using namespace yamas;

typedef void (*signal_handler)(int);

//! internal use, only.
bool file_exists(const std::string &fname) {
    struct stat finfo;
    int fstat;
    fstat = stat(fname.c_str(), &finfo);
    if (fstat == 0) return true;
    return false;
}

void userbreak(int) {
   cerr << endl << endl << "YAMAS has been stopped by you." << endl;
   exit(EXIT_SUCCESS);
}

void genericerror(int) {
   cerr << endl << endl
        << "YAMAS has crashed for an unkown reason." << endl
        << "Perhaps there is a memory problem." << endl
        << "To help improve this program, please send a bug" << endl
        << "report to meesters@imbie.uni-bonn.de." << endl
        << "Please include all data in your report, which are" << endl
        << "necessary to reproduce the crash." << endl;
   exit(EXIT_FAILURE);   
}

void parsingerror(const string msg) {
   cerr << endl << endl
        << "YAMAS was unable to parse your commandline." << endl
        << "The parser says: " << endl
        << msg << endl;
   exit(EXIT_FAILURE);
}

void setupcrashhandlers() {
   signal(SIGINT,  (signal_handler) userbreak);
   signal(SIGSEGV, (signal_handler) genericerror);
   signal(SIGABRT, (signal_handler) genericerror);
}

int main(int args, char* argv[]) {
    cout << "*------------------------------------------------*" << endl
         << "|   YAMAS   |    ";
    cout << tuples::set_open(' ') << tuples::set_close(' ') 
         << tuples::set_delimiter('.') << VERSION;
    // consider lenght of the version tuple
    unsigned short int space = 5;
    if (VERSION.get<0>() >= 10) space--;
    if (VERSION.get<1>() >= 10) space--;
    if (VERSION.get<2>() >= 10) space--;
    for (size_t i = 0; i <= space; i++) cout << " ";
    cout << " |  " << DATE << "  |" << endl;
    cout << "*------------------------------------------------*" << endl;
    cout << "| (c) Christian Meesters                         |" << endl
         << "|     Markus Leber                               |" << endl
         << "|     Christine Herold                           |" << endl
         << "|     Manuel Mattheisen                          |" << endl
         << "|     Dmitriy Drichel                            |" << endl
         // TODO: Fix Andre's name (unicode error)
         << "|     André Lacour                               |" << endl
         << "|     Tim Becker                                 |" << endl
         << "| Gnu Public License v3                          |" << endl
         << "*------------------------------------------------*" << endl;
    
    #ifndef DEBUG
    // before anything else, set up the crash handlers
    setupcrashhandlers();
    #endif // DEBUG

    // defining non-general options or storage items
    string assocfname("");
    #ifdef SANDBOX
    string diagnosticfname("");
    #endif
    string missing("0");
    string fnames("");
    bool   version  = false;
    bool   citation = false;
    #ifdef __linux__
    bool   nonsense = false;
    #endif

    /* WARNING: Do not include duplicate option names for the shell
       and the configuration file.
    */
    // options available only on the commandline
    po::options_description general("General Command Line Options");
    general.add_options()
            ("help,h", "produce help message")
            ("threads,t", po::value<unsigned int>(&ud.nthreads), 
                          "set number of threads to use")
            ("config,c", po::value<string>(&ud.configfname),
                         "specifies name of the configuration file, defaults to 'yamas.cfg'")
            ("logfile,l", po::value<string>(&ud.logfname),
                          "specifies name of the logfile")
            ("algo,a", po::value<string>(&ud.algorithm),
                       "specify algorithm to perform") 
            ("odds_ratio,o", po::value<bool>(&ud.odds_ratio)->zero_tokens(),
                             "if set effect sizes are considered odds-ratios")
            ("verbose,v", po::value<bool>(&ud.verbose)->zero_tokens(), 
                          "verbose or taciturn")
            ("cl", po::value<double>(&ud.cl),
                       "set confidence level required for computation of confidence intervals (default: 0.95)")
            ("pthreshold,p", po::value<double>(&ud.pvaluethreshold), 
                             "select threshold: p-values < threshold will be written out")
            ("r2threshold", po::value<double>(&ud.r2threshold),
                            "select threshold: only SNP pairs with R² >= will be considered")
            #ifdef INCLUDEPRODUCTION
            ("weight_studies", po::value<bool>(&ud.use_study_weights)->zero_tokens(),
                              "if set, studies will be weighted according to the 'study_weigths'"
                              " setting in the configuration file")
            #endif
            ("chromosome", po::value<short int>(&ud.chromosome), 
                           "set specific chromosome to analyse, defaults to all")
            ("assocfilename", po::value<string>(&assocfname),
                              "select a file for which association data will be tablulated")
            ("missing,m", po::value<string>(&missing),
                          "select the character (as string) which marks the missing genotype")
            ("comparison", po::value<bool>(&ud.write_input)->zero_tokens(),
                           "if set, input data will be written out side-by-side with the result")
            ("equal_effects,e", po::value<bool>(&ud.equal_effects)->zero_tokens(),
                           "if set, all effects are assumed to be > 0 (for betas) or > 1 for ORs")
            ("trim", po::value<std::string>(&ud.trim)->implicit_value("NONE"),
                     "if set, only wanted markers will be analysed")
            ("skiplines", po::value<unsigned int>(&ud.skiplines),
                          "number of lines to skip from an inputfile")
            ("stats", po::value<bool>(&ud.stats)->zero_tokens(),
                      "if set, a summary statistic will be written out")
            ("writePLINKconfigfile", po::value<string>(&fnames),
                                     "ask YAMAS to write a configuration file for PLINK input files")
            ("writeINTERSNPconfigfileLogistic", po::value<string>(&fnames),
                                     "ask YAMAS to write a configuration file for INTERSNP logistic input files")
            ("writeINTERSNPconfigfileLinear", po::value<string>(&fnames),
                                     "ask YAMAS to write a configuration file for INTERSNP linear regression input files")                         
            ("writeYAMASconfigfile", po::value<string>(&fnames),
                                     "ask YAMAS to write a configuration file for YAMAS input files")
            ("version", po::value<bool>(&version)->zero_tokens(),
                        "if set, YAMAS will print its version and exit")
            ("citation", po::value<bool>(&citation)->zero_tokens(), "")
    ;
    #ifdef SANDBOX
    po::options_description plotting("Plotting options");
    plotting.add_options()
            ("plot", po::value<bool>(&ud.plot)->zero_tokens(),
                     "if set, all possible plots for the required test will be produced")
            ("selection", po::value<string>(&ud.selection)->implicit_value(""),
                     "select markers for forest plots")
    ;
    #endif

    #ifdef SANDBOX
    po::options_description diagnostic("Options specific for diagnostic tests");
    diagnostic.add_options()
            ("diagnostic_data", po::value<string>(&diagnosticfname),
                                "name of the file containing diagnostic test raw-data")
            ;
    #endif

    #ifdef __linux__
    po::options_description hidden("---");
    hidden.add_options()
            ("shit", po::value<bool>(&nonsense)->zero_tokens(), "");
    #endif
    po::options_description config("Configuration");
    config.add_options()
           ("assocfiles",            po::value<string>(),
                                     "association files")
           #ifdef INCLUDEPRODUCTION
           ("study_weights",         po::value<string>(),
                                     "study weights, e.g. sample size")
           #endif
           ("marker_cols",           po::value<string>(),
                                     "column positions of marker ids")
           ("chromosome_cols",       po::value<string>(),
                                     "column positions of chromosome ids")
           ("position_cols",         po::value<string>(),
                                     "column positions of marker positions")
           ("effect_allele_cols",    po::value<string>(),
                                     "column positions of marker's effect allele")
           ("other_allele_cols",     po::value<string>(),
                                     "column positions of marker's 'other' allele")
           ("effect_cols",           po::value<string>(),
                                     "column positions of marker's OR")
           ("weight_cols",           po::value<string>(),
                                     "column positions of marker's OR-STDERR")
           ("P-value_cols",          po::value<string>(),
                                     "column positions of marker's p-value")
	   ("LD-file",               po::value<string>(),
	                             "name of the file containing LD-data")
           ("Proxyfile",             po::value<string>(),
                                     "name of the file listing allele proxies")
    ;

    // add options 
    po::options_description cmdline_options;
    cmdline_options.add(general);
    #ifdef SANDBOX
    cmdline_options.add(plotting);
    #endif
    #ifdef SANDBOX
    cmdline_options.add(diagnostic);
    #endif
    #ifdef __linux__
    cmdline_options.add(hidden);
    #endif

    po::options_description config_file_options;
    config_file_options.add(config);

    // parse command line
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(args, argv, cmdline_options), vm);
        po::notify(vm);
    }
    // catching unknown options
    catch (exception_detail::clone_impl<exception_detail::error_info_injector<program_options::unknown_option> > &msg) {
        parsingerror(msg.what());
    }
    // catching missing argument
    catch (exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::invalid_command_line_syntax> > &msg) {
        error(msg.what());
    }
    // catch multiple occurences of a command line option
    catch (exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::multiple_occurrences> > &msg) {
        string message = " of command line option(s)";
        parsingerror(msg.what() + message);
    }

    cout << endl; // just one empty line

    if (vm.count("help") || vm.count("h")) {
       cout << general << endl;
       #ifdef SANDBOX
       cout << plotting << endl;
       #endif
       #ifdef SANDBOX
       cout << diagnostic << endl;
       #endif
       exit(EXIT_SUCCESS);
    }
    #ifdef __linux__
    else if (nonsense) {
       cout << "It could be a lot worse ... " << endl;
       if (system("xdg-open http://www.dilbert.com/fast/1994-10-22/")) exit(EXIT_SUCCESS);
       exit(EXIT_FAILURE);
    }
    #endif
    else if (citation) {
       cout << "Please cite:" << endl
            << "BMC Bioinformatics. 2012 Sep 12;13:231. doi: 10.1186/1471-2105-13-231." << endl
            << "Quick, \"imputation-free\" meta-analysis with proxy-SNPs" << endl
            << "Meesters C, Leber M, Herold C, Angisch M, Mattheisen M, Drichel D, Lacour A, Becker T" << endl;
       #ifdef __linux__
       cout << "Please wait for your browser to be opened. The Pubmed entry will be displayed, soon." << endl;
       if (system("xdg-open www.ncbi.nlm.nih.gov/pubmed/22971100")) exit(EXIT_SUCCESS);
       else
          cerr << "Sorry. Browser or operation system not recognized." << endl;
       #endif
       exit(EXIT_SUCCESS);
    }
    else if (vm.count("version")) {
       cout << "Revision number          : " << REVISION << endl;
       #ifdef __GNUC__
       cout << "Compiled with GCC version:" << GCC_VERSION << endl;
       #endif
       exit(EXIT_SUCCESS);
    }
    else if (vm.count("writePLINKconfigfile")) {
       write_PLINKconfig(fnames);
       exit(EXIT_SUCCESS);
    }
    else if (vm.count("writeINTERSNPconfigfileLogistic")) {
       write_INTERSNPconfigLogistic(fnames);
       exit(EXIT_SUCCESS);
    }
    else if (vm.count("writeINTERSNPconfigfileLinear")) {
       write_INTERSNPconfigLinear(fnames);
       exit(EXIT_SUCCESS);
    }
    else if (vm.count("writeYAMASconfigfile")) {
       write_YAMASconfig(fnames);
       exit(EXIT_SUCCESS);
    }

    if (assocfname.size()) { // instead of a meta-analysis a file
                             // name for association testing is given
        if (assocfname.empty()) {
            cerr << "ERROR: no input file specified." << endl
                 << "Please use the '--association' option to " << endl
                 << "specify a tped/tfam base name" << endl;
            exit(EXIT_FAILURE);
        }
        if (missing.size() != 1) {
            cerr << "ERROR: The missing character must be "
                 << "*one* character." << endl;
            exit(EXIT_FAILURE);
        }
        association_test(assocfname, missing);
        return EXIT_SUCCESS;
    }
    #ifdef SANDBOX
    else if (diagnosticfname.size()) { // instead of any standard routine
                                       // a file name for diagnostic data
                                       // is given and those will be evaluated
        if (diagnosticfname.empty()) {
            cerr << "ERROR: no input file specified." << endl
                 << "Please use the '--diagnostic' option to " << endl
                 << "specify a file containing data of a diagnostic test." << endl;
            exit(EXIT_FAILURE);
        }
        if (ud.logfname.size() == 0) {
            // set a standard log file name
            ud.logfname = "yamas_diagnostic";
        }
        if (file_exists(ud.logfname + ".log")) existinglogfile();
        diagnostic_test(diagnosticfname);
        return EXIT_SUCCESS;
    }
    #endif

    // test for presence of the configuration file -- parse if possible
    ifstream configfile(ud.configfname.c_str());
    if (!configfile) {
       cout << "Cowardly refusing to work, unable to open the config file: " 
            << ud.configfname << std::endl;
       exit(EXIT_FAILURE);
    }
    else {
       po::store(po::parse_config_file(configfile, config_file_options), vm);
       po::notify(vm);
    }

    // set parsed variables (if needed)
    omp_set_num_threads(ud.nthreads);

    // extract project setting from the configuration file
    const vector<Project> projects = extractor(vm);

    // issue report
    cout << endl;
    cout << "YAMAS is going to use the following settings globally:" << endl;
    cout << "threads:     "             << ud.nthreads << endl;
    cout << "will work on chrosomome: " << ud.chromosome
	                                << " ['0' = ALL, '-1' = ALL at once]" << endl; 
    cout << "and assumes odds ratios "  << (ud.odds_ratio ? "(YES)" : "(NO)")
         << " or beta estimates "       << (ud.odds_ratio ? "(NO)": "(YES)") << endl;
    #ifdef INCLUDEPRODUCTION
    if (ud.use_study_weights) {
       cout << "YAMAS will try to use weights per study ('Stouffer's method')" << endl;
       for (size_t i = 0; i < projects.size(); i++) {
          if (projects[i].get_weight() == 0) {
             yamas::error("Cannot work with study weights equalling zero.");
          }
       }
    }
    #endif
    cout << "config file: " << ud.configfname << endl << endl;

    // print header information per project, if requested
    if (ud.verbose) {
        cout << "Values found in configuration file:" << endl;
        for (vector<Project>::const_iterator it = projects.begin(); 
            it != projects.end(); it++) {
            if (!(*it).check()) {
                cerr << "Unable to proceed with input for Project " << (*it).id << "." <<endl;
                cerr << "Please check your configuration file." << endl;
                exit(EXIT_FAILURE);
            }
            if (ud.verbose) (*it).print_header();
            cout << endl;
       }
    }

    // now select the choosen algorithm
    map<string, short unsigned int> algomap;
    algomap["pointwise"]       = 1;
    #ifdef INCLUDEPRODUCTION
    algomap["ldblockwise"]     = 2;
    #endif
    algomap["fillwithproxies"] = 3;
    switch(algomap[ud.algorithm]) {
        case 1: // pointwise
            // no logfilename set?
            if (ud.logfname.size() == 0) {
                // set a standard name
                ud.logfname = "yamas_pointwise";
            }
            // before starting to work, check whether the specified log field exists
            if (file_exists(ud.logfname + ".log")) existinglogfile();
	    pointwise_comparison(projects);
            break;
        #ifdef INCLUDEPRODUCTION
	case 2: // per LD block
	    if (ud.logfname.size() == 0) {
		// set a standard log file name
		ud.logfname = "yamas_perLDblock";
	    }
            if (file_exists(ud.logfname + ".log")) existinglogfile();
	    per_ldblock(projects, vm["LD-file"].as<string>(), vm["Proxyfile"].as<string>());
	    break;
        #endif
        case 3: // according to proxy allele
            if (ud.logfname.size() == 0) {
                // set a standard log file name
                ud.logfname = "yamas_proxyanalysis";
            }
            if (file_exists(ud.logfname + ".log")) existinglogfile();
            #ifdef TIMER
            {
            Timer t0(0);
            #endif
	    pointwise_proxy(projects, vm["Proxyfile"].as<string>());
            #ifdef TIMER
            }
            #endif
            break;
        default:
            cerr << "Unable to comply: No or unknown algorithm selected." << endl;
            cerr << "The algorithms YAMAS knows are:" << endl;
            // erase key introduced by the switch statement
            algomap.erase(ud.algorithm);
            for (map<string, short unsigned int>::const_iterator ait = algomap.begin();
                 ait != algomap.end(); ait++) {
                cerr << "- " << (*ait).first << endl;
            }
            cerr << "YAMAS defaults to 'pointwise'." << endl;
            return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
