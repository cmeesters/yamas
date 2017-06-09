/********************************************************************/
/*                                                                  */
/* YAMAS - Yet Another Meta-Analysis Software                       */
/*                                                                  */
/* Functions enabling YAMAS to write default configuration files    */
/* for various input formats.                                       */
/*                                                                  */
/********************************************************************/

// $Rev: 640 $
// $LastChangedDate: 2010-12-21 09:23:03 +0100 (Di, 21 Dez 2010) $
// $Author: meesters $

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <boost/algorithm/string.hpp>

#include "../include/ya_errors.hpp"

namespace yamas {
    void write_yamasconfig(const std::vector<std::string> &infnames,
                           const unsigned int cols[],
		           const std::string &outfname  = "yamas.cfg",
                           const std::string &ldfile    = "ldfile.txt.gz",
                           const std::string &proxyfile = "proxyreference.txt.gz") {
        // number of inputfiles
        unsigned int infiles = infnames.size();
        // column names
        std::vector<std::string> colnames;
        colnames.push_back("marker_cols");
        colnames.push_back("chromosome_cols");
        colnames.push_back("position_cols");
        colnames.push_back("effect_allele_cols");
        colnames.push_back("other_allele_cols");
        colnames.push_back("effect_cols");
        colnames.push_back("weight_cols");
        colnames.push_back("P-value_cols");

        // create outstram
        std::ofstream outfile(outfname.c_str());

        // write assocfilenames
        outfile << "assocfiles  = "; 
        for (size_t i = 0; i < infiles; i++) {
            outfile << infnames[i] << " ";
        }
        outfile << std::endl;

        // write column ids and cols
        for (size_t i = 0; i < colnames.size(); i++) {
            outfile << colnames[i] << " = ";
            for (size_t j = 0; j < infiles; j++) {
                outfile << cols[i] << " ";
            }
            outfile << std::endl;
        }
        //  finally write the two file descriptors
        outfile << "LD-file = "   << ldfile << std::endl;
        outfile << "Proxyfile = " << proxyfile << std::endl;

        outfile.close();
    }

    void write_PLINKconfig(const std::string &infnames) {
        std::vector<std::string> fnames;
        std::cout << "Writing a configuration file ('yamas.cfg') "   << std::endl
                  << "for PLINK. Note that this configuration file " << std::endl
                  << "enables YAMAS to parse association data cal- " << std::endl
                  << "culated with PLINK (version 1.07) using the "  << std::endl
                  << "'--assoc' and '--ci' options, only." << std::endl;
        boost::algorithm::split(fnames, infnames, boost::is_any_of(",; "), 
                                        boost::token_compress_on);
        // defining column integers for PLINK association output format
        unsigned int cols[] = {2, 1, 3, 4, 7, 10, 11, 9};

        yamas::write_yamasconfig(fnames, cols);
    }

    void write_INTERSNPconfigLogistic(const std::string &infnames) {
        std::cout << "Writing a configuration file ('yamas.cfg') "   << std::endl
                  << "for YAMAS input file, generated by INTERSNP"   << std::endl
                  << "with the 'linear' option "                     << std::endl
                  << "(INTERSNP version 1.0.8)."   << std::endl;
        std::vector<std::string> fnames;
        boost::algorithm::split(fnames, infnames, boost::is_any_of(",; "), 
                                        boost::token_compress_on);
        // defining column integers for INTERSNP logistic association output format
        unsigned int cols[] = {3, 2, 4, 6, 7, 28, 29, 12};

        yamas::write_yamasconfig(fnames, cols);
    }
    
    void write_INTERSNPconfigLinear(const std::string &infnames) {
        std::cout << "Writing a configuration file ('yamas.cfg') "   << std::endl
                  << "for YAMAS input file, generated by INTERSNP"   << std::endl
                  << "with the 'linear' option "                     << std::endl
                  << "(INTERSNP version 1.0.8)."   << std::endl;
        std::vector<std::string> fnames;
        boost::algorithm::split(fnames, infnames, boost::is_any_of(",; "), 
                                        boost::token_compress_on);
        // defining column integers for INTERSNP linear regression output format
        unsigned int cols[] = {3, 2, 4, 6, 7, 17, 18, 12};

        yamas::write_yamasconfig(fnames, cols);
    }

    void write_YAMASconfig(const std::string &infnames) {
        std::cout << "Writing a configuration file ('yamas.cfg') "    << std::endl
                  << "for YAMAS-output. Note that this configuration" << std::endl
                  << "file enables YAMAS to parse association data"   << std::endl
                  << "cal culated with YAMAS using the '--assoc'"     << std::endl
                  << "option." << std::endl;
        std::vector<std::string> fnames;
        boost::algorithm::split(fnames, infnames, boost::is_any_of(",; "), 
                                        boost::token_compress_on);
        unsigned int cols[] = {2, 1, 3, 4, 7, 9, 11, 12};

        yamas::write_yamasconfig(fnames, cols);
    }

} // end namespace yamas

