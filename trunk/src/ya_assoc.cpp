/********************************************************************/
/*                                                                  */
/* YAMAS - Yet Another Meta-Analysis Software                       */
/*                                                                  */
/* module for calculating association data                                                       */
/*                                                                  */
/********************************************************************/

// $Rev: 800 $
// $LastChangedDate: 2011-06-20 16:00:23 +0200 (Mo, 20 Jun 2011) $
// $Author: meesters $

#include <vector>
#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <cmath>
#include <limits>
#include <utility> // for std::pair

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/distributions/chi_squared.hpp>
using boost::math::chi_squared;

#include "../include/ya_interface.hpp"
#include "../include/ya_io.hpp"
#include "../include/ya_errors.hpp"
#include "../include/ya_global.hpp"
#include "../include/ya_project.hpp"

EXTERN yamas::UsageData ud;
namespace yamas {
    //! local map for hettype counting
    typedef std::map<std::string, unsigned int> hettypecounts;
    //! short hand to find the min/maj allele
    bool compare(std::pair<std::string, unsigned int> i, std::pair<std::string, unsigned int> j) {
        return i.second < j.second;
    }
    /*! 
        basic association test
        writes out all data from a tped/tfam file tuple
        needed for meta analysis with yamas
        --> requires standard conformant tped/tfam files
    */
    void association_test(const std::string &fname, 
                          const std::string &missing) {
        // TODO: add association test for X chromosome
        // the file name might or might not have an extension,
        // therefore try to extract the extension
        std::size_t i = fname.rfind('.', fname.length());
        std::string basefname;
        // so it has an extension?
        if (i != std::string::npos) basefname = fname.substr(0, i);
        else basefname = fname;
        std::string tpedfname(basefname); tpedfname.append(".tped");
        std::string tfamfname(basefname); tfamfname.append(".tfam");
        
        if (ud.verbose) std::cout << "Performing association test on " << tpedfname
                                  << "/" << tfamfname << "." << std::endl;
        // open input streams
        std::ifstream tpedfile(tpedfname.c_str());
        std::ifstream tfamfile(tfamfname.c_str());
        // open output stream
        std::ofstream outfile(basefname.append(".assoc").c_str());
        // write header
        outfile << "# CHR"
                << std::setw(12) << "SNP"
                << std::setw(11) << "BP"
                << std::setw(3)  << "EA"
                << std::setw(7)  << "F_CA"
                << std::setw(7)  << "F_CO"
                << std::setw(3)  << "RA"
                << std::setw(9)  << "OR"
                << std::setw(12) << "BETA"
                << std::setw(12) << "CHISQ"
                << std::setw(11) << "SE"
                << std::setw(12) << "P"
                << std::setw(12) << "TREND"
                << std::setw(12) << "P_ARMITAGE"
                << std::endl;

        // temporary dummy vector
        std::vector<std::string> fields;
        // temporary dummy string
        std::string line;
        // allele counts
        hettypecounts cacounts, cocounts, cacounts_d, cocounts_d, merged;
        // temporary allele or hettype
        std::string allele1(""), allele2(""), het(""), hommin(""), hommaj("");
        // temporary indeces
        unsigned int index = 0;
        // temporary minor / major allele counts
        std::pair<std::string, unsigned int> minca, majca, minco, majco, minor, major;
        // temporary counts for the Armitage trend test
        double n = 0, n_k = 0, n_not_k = 0, n1 = 0, n2 = 0;

        // values to be calculated
        double odds_ratio = 0, beta = 0, chisq = 0, se = 0, p = 0, f_ca = 0, f_co = 0, trend = 0,
               p_armitage = 0, ca = 0, ea = 0;
        // chisq distribution with different degrees of freedom
        chi_squared dist_df1(1);
        // temp holders for easier four-fold table handling
        double a = 0, b = 0, c = 0, d = 0;

        // read in the tfam file to sort out cases and controls
        std::vector<bool> caco; // case-control list: true == case
	// check whether the tfam-file exists
        if (!tfamfile) yamas::ioerror(tfamfname);

        while (!tfamfile.eof()) {
            std::getline(tfamfile, line);
            // split line
            boost::algorithm::split(fields, line, boost::is_any_of(" "));
            if (fields.size() < 5) break;
            if (fields[5] == "2") caco.push_back(true);
            else                  caco.push_back(false);
        }
        tfamfile.close();
        // calculate number of cases for later use
        unsigned int cases = 0;
        for (std::vector<bool>::const_iterator it = caco.begin();
            it != caco.end(); it++) if (*it) cases++;

        if (ud.verbose) {
            std::cout << std::endl;
            std::cout << "Found " << cases << " cases and " << caco.size() - cases
                      << " controls." << std::endl;
        }
        
        // set dummy controls
        if (!tpedfile) {
            std::cerr << "ERROR: Unable to open: " << tpedfname << "." << std::endl;
            std::cerr << "Please check this file's name and path." << std::endl;
            std::exit(1);
        }

        // walk over the file
        while (!tpedfile.eof()) {
            std::getline(tpedfile, line);
            // split line
            boost::algorithm::split(fields, line, boost::is_any_of(" "));
            if (fields.size() < 5) break;

            outfile << std::setw(5)  << fields[0]  // chromosome
                    << std::setw(12) << fields[1]  // rsid
                    << std::setw(11) << fields[3]; // base pair position
            
            // start counting
            index = 0;
            for (std::size_t i = 4; i < fields.size(); i += 2) {
                allele1 = boost::lexical_cast<std::string>(fields[i]);
                merged[allele1]++;
                allele2 = boost::lexical_cast<std::string>(fields[i + 1]);
                merged[allele2]++;
                // just sort unambigously, lexically
                het = allele1 < allele2 ? allele1 + allele2 : allele2 + allele1;
                if (caco[index]) { // case? 
                    cacounts[allele1]++;
                    cacounts[allele2]++;
                    cacounts_d[het]++;
                }
                else { // control?
                    cocounts[allele1]++;
                    cocounts[allele2]++;
                    cocounts_d[het]++;
                }
                // calculate index in vector of case/control-flags
                index++;
            }

            // now start calculating the counts for the Armitage trend test
            n       = double(caco.size() - cacounts_d[missing + missing] 
                             - cocounts_d[missing + missing]);                     // total
            n_k     = double(cases - cacounts_d[missing + missing]);               // # of cases
            n_not_k = double(caco.size() - cases - cocounts_d[missing + missing]); // # of controls
            // we aren't interested in the missing values
            cacounts.erase(missing);
            cocounts.erase(missing);
            cacounts_d.erase(missing + missing);
            cocounts_d.erase(missing + missing);
            merged.erase(missing); // missing is a single character

            // di-allelic standard case
            if (cacounts.size() == 2 and cocounts.size() == 2) {
                minca = *min_element(cacounts.begin(), cacounts.end(), compare); // a
                majca = *max_element(cacounts.begin(), cacounts.end(), compare); // b
                minco = *min_element(cocounts.begin(), cocounts.end(), compare); // c
                majco = *max_element(cocounts.begin(), cocounts.end(), compare); // d

                // calculate minor allele frequencies for cases and controls
                f_ca = double(minca.second) / double(minca.second + majca.second);
                f_co = double(minco.second) / double(minco.second + majco.second);

                a = double(minca.second);
                b = double(majca.second);
                c = double(minco.second);
                d = double(majco.second);

                odds_ratio = ((a * d) / (b * c));
                beta       = std::log(odds_ratio);
                se         = std::sqrt(1/(a + b + c + d));
                chisq      = (n * std::pow(((a * d) - (b * c)), 2)) / 
                             ((a + b) * (a + c) * (b + d) * (c + d)) * 2; 
                p          = cdf(complement(dist_df1, chisq));
            }
            if (merged.size() == 2) { // marker not monomorph for cases or controls?
                minor = *min_element(merged.begin(), merged.end(), compare);
                major = *max_element(merged.begin(), merged.end(), compare);

                // heterozygous hettype
                if (minor.first[0] < major.first[0]) { 
                    het = minor.first[0];
                    het.push_back(major.first[0]);
                }
                else {
                    het = major.first[0];
                    het.push_back(minor.first[0]);
                }
                hommin = minor.first + minor.first;
                hommaj = major.first + major.first;

                // Armitage trend-test (first calculate those values, 
                // which aren't determined, yet
                n1         = double(cacounts_d[het] + cocounts_d[het]);
                n2         = double(cacounts_d[hommaj] + 
                                    cocounts_d[hommaj]);
                ca         = double(cacounts_d[het]); // heterozygous cases
                ea         = double(cacounts_d[hommaj]);
                try {
                    trend = (n * (std::pow(n * (ca + 2 * ea) - (n_k * (n1 + 2 * n2)), 2))) /
                            (n_k * n_not_k * (n * (n1 + 4 * n2) - std::pow((n1 + 2 * n2), 2)));
                    p_armitage = cdf(complement(dist_df1, trend));
                }
                catch (std::domain_error) {
                    trend = std::numeric_limits<double>::quiet_NaN();
                    p_armitage = std::numeric_limits<double>::quiet_NaN();
                }
            }
            else {
                trend = std::numeric_limits<double>::quiet_NaN();
                p_armitage = std::numeric_limits<double>::quiet_NaN();
            }

            // output
            outfile << std::setw(3)  << minca.first << std::fixed
                    << std::setw(7)  << std::setprecision(4) << f_ca << std::fixed
                    << std::setw(7)  << std::setprecision(4) << f_co
                    << std::setw(3)  << majca.first 
                    << std::setw(9)  << odds_ratio  
                    << std::setw(12) << beta << std::setw(12);
            if (chisq > 100) {
                outfile << std::scientific << chisq << std::fixed;
            }
            else
                outfile << chisq;
            outfile << std::setw(11) << se << std::setw(12);
            if (p < 0.0001) {
                outfile << std::scientific << p << std::fixed;
            }
            else 
                outfile << p;
            if (trend > 100) {
                outfile << std::scientific << std::setw(12) << trend 
                        << std::fixed << std::setw(12);
            }
            else
                outfile << std::setw(12) << trend << std::setw(12);
            if (p_armitage < 0.0001) {
                outfile << std::scientific << p_armitage 
                        << std::fixed << std::endl;
            }
            else
                outfile << p_armitage << std::endl;

            // clear counts
            cacounts.erase(cacounts.begin(), cacounts.end());
            cocounts.erase(cocounts.begin(), cocounts.end());
            cacounts_d.erase(cacounts_d.begin(), cacounts_d.end());
            cocounts_d.erase(cocounts_d.begin(), cocounts_d.end());
            merged.erase(merged.begin(), merged.end());
            het.erase();
            hommin.erase(); hommaj.erase();
        }

        tpedfile.close();
        outfile.close();

        // this as output to stdin
        if (ud.verbose) {
            std::cout << "Analysis output has been written to " << basefname 
                      << "." << std::endl;
            std::cout << "You can use the following as a template for your " << std::endl
                      << "yamas configuration file:" << std::endl; 
            std::cout << "-------------------------------------------------" << std::endl;
            std::cout << "assocfiles            = " << basefname << std::endl
                      << "marker_cols           = 2"             << std::endl
                      << "chromosome_cols       = 1"             << std::endl
                      << "position_cols         = 3"             << std::endl
                      << "effect_allele_cols    = 4"             << std::endl
                      << "other_allele_cols     = 7"             << std::endl
                      << "effect_cols           = 9"             << std::endl
                      << "weight_cols           = 11"            << std::endl
                      << "P-value_cols          = 12"            << std::endl;
            std::cout << "-------------------------------------------------" << std::endl;
            std::cout << "Please merge these values into your final confi- " << std::endl
                      << "figuration file." << std::endl;

        }
        std::exit(0);
    }
} // end namespace yamas

