/********************************************************************/
/*                                                                  */
/* YAMAS - Yet Another Meta-Analysis Software                       */
/*                                                                  */
/* Functions for YAMAS I/O support.                                 */
/*                                                                  */
/********************************************************************/

// $Rev: 908 $
// $LastChangedDate: 2012-06-13 14:39:38 +0200 (Mi, 13 Jun 2012) $
// $Author: meesters $

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/tuple/tuple.hpp> 
#include <boost/tuple/tuple_io.hpp>

#include <omp.h>

#include "../include/ya_io.hpp"
#include "../include/ya_errors.hpp"
#include "../include/ya_project.hpp"
#include "../include/ya_version.hpp"

namespace yamas {
    void write_header(const std::vector<yamas::Project> &projects, std::ofstream &out,
                      const yamas::UsageData &ud,
	              bool write_status) {
        // short hand for project.size() -- used several times
        unsigned int nprojects = projects.size();
        // first write some general information about the selected
        // parameters and the version of the program
        out << "# *--------------------------------------------*" << std::endl
            << "# | YAMAS - Yet Another Meta-Analysis Software |" << std::endl
            << "# *--------------------------------------------*" << std::endl
            << "#" << std::endl
            << "# Used YAMAS version: " << boost::tuples::set_open(' ') 
                                        << boost::tuples::set_close(' ') 
                                        << boost::tuples::set_delimiter('.') 
                                        << yamas::VERSION << std::endl
            << "# Selected algorithm:      " << ud.algorithm   << std::endl
            << "# to run on chromosome:    " << ud.chromosome << " ['0' = ALL, '-1' = ALL at once]" << std::endl
            << "# data were specified in:  " << ud.configfname << std::endl;
            #ifdef INCLUDEPRODUCTION
            if (ud.use_study_weights ) {
               out << "# Selected study weighting: " << std::endl;
               for (size_t i = 0; i < projects.size(); i++) {
                   out << "# Project " << projects[i].get_id() << ": " << projects[i].get_weight() << std::endl;
               }
            }
            else
               out << "# No particular study weighting was requested." << std::endl;
            #endif
        out << "# ----" << std::endl
            << "# all data are written to: " << ud.logfname << ".log/.err" << std::endl
            << "# ----" << std::endl
            << "# P-value threshold was:   " << ud.pvaluethreshold << std::endl
            << "# confidence level was:    " << ud.cl << std::endl
            << "# " << ud.skiplines << " line(s) were skipped per input file" << std::endl
            << "# and OR were considered: " << (ud.odds_ratio ? "(YES)" : "(NO)") << std::endl
            << (ud.equal_effects ? "# yamas equalized all effects" : "# yamas tried adjusting effects based on marker alleles") << std::endl;
        if (ud.trim == "NONE")       out << "# output is trimmed to Project 1's content" << std::endl;
        else if (ud.trim.size() > 0) out << "# output is trimmed to (the content of): " << ud.trim << std::endl;
        out << "# ----" << std::endl;

        if (ud.algorithm == "fillwithproxies") {
            out << "# r^2 threshold:          "  << ud.r2threshold << std::endl
                << "# ----" << std::endl;
        }
        out << "#" << std::endl;

	if (write_status) {
            out << std::setw(23) << std::left << "#STATUS"
                << std::right 
                << std::setw(12) << "ID";
	}
	else {
            out << std::setw(12) << std::left << "# ID";
	}
        out << std::setw(4)  << "CHR"
            << std::setw(12) << "BP"
            << std::setw(4)  << "EA"
            << std::setw(4)  << "OA"
            << (nprojects <= 20 ? std::setw(20) : std::setw(nprojects + 1))
                             << "DIRECTIONS"
            << std::setw(12) << "EFFECT"
            << std::setw(10) << "SE"
            << std::setw(15) << "CI"
            << std::setw(14) << "P"
            << std::setw(14) << "Q"
            << std::setw(7)  << "I^2"
            << std::setw(12) << "EFFECT(R)" 
            << std::setw(11) << "SE(R)"
            << std::setw(17) << "CI(R)"
            << std::setw(14) << "P(R)";
        #ifdef INCLUDEPRODUCTION
        // conditional output of Stouffer's method's values
        if (ud.use_study_weights) {
            out << std::setw(12) << "Stouffer-Z"
                << std::setw(12) << "P(Stouffer)";
        }
        #endif
        if (ud.write_input) {
            for (std::size_t i = 0; i < nprojects; i++) {
                out << std::setw(11) << "PROJECT_ID"
                    << std::setw(12) << "ID"
                    << std::setw(4)  << "CHR"
                    << std::setw(12) << "BP" 
                    << std::setw(4)  << "EA"
                    << std::setw(4)  << "OA"
                    << std::setw(12) << "EFFECT"
                    << std::setw(11) << "SE"
                    << std::setw(14) << "P";
                if (ud.algorithm == "fillwithproxies")
                    out << std::setw(6) << "R^2"
                        << std::setw(12) << "PROXYHAPLO";
                        //<< std::setw(13) << "STRANDSWITCH";
            }
        }
        out << std::endl;
    }

    //! writes generalized header for error log files
    void write_errorlog_header(const std::vector<Project> &projects, 
                                      std::ofstream &out) {
        out << "# This file contains all markers for which no meta-analysis" << std::endl
            << "# could be carried out."       << std::endl
            << "# Possible reasons include:"   << std::endl
            << "# - less than 2 valid markers" << std::endl
            << "# - P-value above given threshold (not strictly an error," << std::endl
            << "#   but still included to provide some overview)"          << std::endl
            << "# In addition all projects which contain information"      << std::endl
            << "# about a particular marker are listed."                   << std::endl;
        out << std::left 
            << std::setw(12) << "#ID" << std::right
            << std::setw(11) << "PROJECT_ID"
            << std::setw(12) << "ID"
            << std::setw(4)  << "CHR"
            << std::setw(12) << "BP" 
            << std::setw(4)  << "EA"
            << std::setw(4)  << "OA"
            << std::setw(12) << "EFFECT"
            << std::setw(11) << "SE"
            << std::setw(14) << "P" 
            << std::setw(12) << "PROJECT_ID"
            << std::setw(12) << "ID"
            << std::setw(4)  << "CHR"
            << std::setw(12) << "BP" 
            << std::setw(4)  << "EA"
            << std::setw(4)  << "OA"
            << std::setw(12) << "EFFECT"
            << std::setw(11) << "SE"
            << std::setw(14) << "P"
            << std::setw(11) << "etc. ..."
            << std::endl;
    }
} // end namespace yamas

