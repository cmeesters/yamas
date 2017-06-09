/********************************************************************/
/*                                                                  */
/* YAMAS - Yet Another Meta-Analysis Software                       */
/*                                                                  */
/* Functions for YAMAS statistical overview support.                */
/*                                                                  */
/********************************************************************/

// $Rev: 807 $
// $LastChangedDate: 2011-08-09 13:15:42 +0200 (Di, 09 Aug 2011) $
// $Author: meesters $

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "../include/ya_io.hpp"
#include "../include/ya_metamarker.hpp"
#include "../include/ya_stats.hpp"

// defining singleton pointer
yamas::Statistics* yamas::Statistics::inst_ = NULL;

yamas::Statistics::Statistics() : fixed_stats_(4, 0),
                                  random_stats_(4, 0) {}

namespace yamas {
    yamas::Statistics* Statistics::get_Statistics() {
        if (inst_ == NULL) {
            inst_ = new Statistics();
        }
        return inst_;
    }

    void Statistics::update(const yamas::ResultMarker &rmarker) {
        double fpvalue = rmarker.get_pvalue();
        double rpvalue = rmarker.get_rpvalue();

        // update for P-values <= 5e-2
        if (fpvalue <= 5e-2) fixed_stats_[0]++;
        if (rpvalue <= 5e-2) random_stats_[0]++;
        // update for P-values <= 1e-4
        if (fpvalue <= 1e-4) fixed_stats_[1]++;
        if (rpvalue <= 1e-4) random_stats_[1]++;
        // update for P-values <= 1e-6
        if (fpvalue <= 1e-6) fixed_stats_[2]++;
        if (rpvalue <= 1e-6) random_stats_[2]++;
        // update for P-values <= 1e-8
        if (fpvalue <= 1e-8) fixed_stats_[3]++;
        if (rpvalue <= 1e-8) random_stats_[3]++;
    }

    void Statistics::writeout(const unsigned long int &total) {
        // short hands:
        std::vector<unsigned long int> f = fixed_stats_;
        std::vector<unsigned long int> r = random_stats_;

        // temporary string stream for optional output redirection
        std::stringstream msg (std::stringstream::in | std::stringstream::out);

        msg << "# Summary:" << std::endl << std::scientific << std::setprecision(4)
            << "# YAMAS found: (fractions)"  << std::endl
            << "# Shell   |  fixed effect P-values  |  random effect P-values" << std::endl
            << "# <= 5e-2 | " << std::setw(16) << double(f[0]) / total << std::setw(9) << "|" 
                              << std::setw(16) << double(r[0]) / total << std::endl
            << "# <= 1e-4 | " << std::setw(16) << double(f[1]) / total << std::setw(9) << "|" 
                              << std::setw(16) << double(r[1]) / total << std::endl
            << "# <= 1e-6 | " << std::setw(16) << double(f[2]) / total << std::setw(9) << "|" 
                              << std::setw(16) << double(r[2]) / total << std::endl
            << "# <= 1e-8 | " << std::setw(16) << double(f[3]) / total << std::setw(9) << "|" 
                              << std::setw(16) << double(r[3]) / total << std::endl;
        std::ofstream out((ud.logfname + ".stats").c_str());
        out << msg.str();
        if (ud.verbose) std::cout << msg.str();
    }
} // end namespace yamas

