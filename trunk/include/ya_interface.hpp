/***************************************************
 *                                                 *
 * YAMAS - Yet Another Meta-Analysis               *
 *                                                 *
 * This file contains the YAMAS user interface as  *
 * part of an all-public class. All members can be *
 * made available by importing the ya_global.hpp-  *
 * file and using the EXTERN statement.            *
 *                                                 *
 **************************************************/

// $Rev: 907 $
// $LastChangedDate: 2012-06-04 14:15:13 +0200 (Mo, 04 Jun 2012) $
// $Author: meesters $

#include <vector>
#include <string>

#ifndef YA_INTERFACE
#define YA_INTERFACE

namespace yamas {
    /*! This struct with public members, only, is to provide easy
        accessible flags globally.
    */
    struct UsageData {
         //! number of threads
         unsigned int nthreads;
         //! which chromosome to work on - all, per default
         short int chromosome; 
         //! p-values <= threshold will be written out
         double pvaluethreshold;

         /*! only SNP-pairs with R2 >= threshold will be considered for
             proxy based analysis
         */
         double r2threshold;
         #ifdef INCLUDEPRODUCTION
         // whether or not to use study weights
         bool use_study_weights;
         #endif
         // confidence level
         double cl;
         //! default log file name
         std::string logfname;
         //! verbosity flag
         bool verbose;
         //! flag: shall the input data be logged, too?    
         bool write_input;
         //! algorithm choise
         std::string algorithm;
         //! name of the configuration file
         std::string configfname;
         //! determines whether or not Beta estimates are assumed as effect values
         bool odds_ratio;
         //! determines whether or not the user wants all effects to be in one direction
         bool equal_effects;
         //! if set, only wanted markers will be in the output
         std::string trim;
         //! number of lines to skip in an inputfile
         unsigned int skiplines;
         //! whether or not to write out summary stats
         bool stats;
         //! if set, decides which plot to do
         bool plot;
         //! string of markers for a blobbogram / forest plot
         std::string selection;
         //! diagnostic data file name
         std::string diagnostic_data_fname;
         //! diagnostic data MA file name
         std::string diagnostic_MA_fname;


         UsageData(): nthreads(1),
                      chromosome(0),
                      pvaluethreshold(1.0),
                      r2threshold(0.2),
                      cl(0.95), 
                      logfname(""),
                      verbose(false),
                      write_input(false),
                      algorithm("pointwise"),
                      configfname("yamas.cfg"),
                      odds_ratio(false),
                      equal_effects(false),
                      trim(""),
                      skiplines(0),
                      stats(false),
                      plot(false),
                      selection(""),
                      diagnostic_data_fname(""),
                      diagnostic_MA_fname("")
                      {};
    };
} // end namespace yamas
    
#endif // YA_INTERFACE
