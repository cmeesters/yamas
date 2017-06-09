/********************************************************************/
/*                                                                  */
/* YAMAS - Yet Another Meta-Analysis Software                       */
/*                                                                  */
/* Class holding project properties                                 */
/*                                                                  */
/********************************************************************/

// $Rev: 925 $
// $LastChangedDate: 2012-09-10 16:19:51 +0200 (Mo, 10 Sep 2012) $
// $Author: meesters $

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <stdlib.h>
#include <assert.h>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <boost/program_options.hpp>
using namespace boost;
namespace po = program_options;

#include <omp.h>

// include own software parts
#include "../include/ya_project.hpp"
#include "../include/ya_metamarker.hpp"
#include "../include/ya_utils.hpp"
#include "../include/ya_io.hpp"
#include "../include/ya_errors.hpp"
#include "../include/ya_interface.hpp"
#include "../include/ya_global.hpp"
#include "../include/ya_lddata.hpp"

EXTERN yamas::UsageData ud;

unsigned int yamas::Project::next_id = 0;

yamas::Project::Project() : fname(""),
                            current_line(0)
                            #ifdef INCLUDEPRODUCTION
                            ,
                            weight(0)
                            #endif
                            {
    id = ++next_id;
    // write column assignment out explicitly
    // in order to ensure having keys - easily 
    // to be extended
    columns["marker"]           = 0;
    columns["chromosome"]       = 0;
    columns["position"]         = 0;
    columns["effect_allele"]    = 0;
    columns["other_allele"]     = 0;
    columns["effect"]           = 0;
    columns["weight"]           = 0;
    columns["P-value"]          = 0;
}

// copy-constructor
yamas::Project::Project(const yamas::Project &p) :
                        id(p.id),
                        columns(p.columns),
                        fname(p.fname),
                        current_line(p.current_line)
                        #ifdef INCLUDEPRODUCTION
                        ,
                        weight(p.weight)
                        #endif
                   {
}

yamas::Project::Project(const std::map<std::string, short unsigned int> &cols, 
                        const std::string &ifname
                        #ifdef INCLUDEPRODUCTION
                        ,
                        const double &iweight
                        #endif
                        ) : current_line(0) {
    id = ++next_id;
    columns = cols;
    fname   = ifname;
    #ifdef INCLUDEPRODUCTION
    weight  = iweight;
    #endif
    std::ifstream instream;
}

yamas::Project::~Project() {}

namespace yamas {
    std::ostream& operator<< (std::ostream& out, const Project &proj) {
        std::string ns("not specified");
        out << "Project " << proj.id << ":" << std::endl;
        out << "Associated data file: " << proj.fname      << std::endl;
        #ifdef INCLUDEPRODUCTION
        out << "Associated study weight: " << proj.weight;
        out << (ud.use_study_weights ? " (will be used)" : " (will not be used)") << std::endl; 
        #endif
        for (std::map<std::string, short unsigned int>::const_iterator it = proj.columns.begin();
             it != proj.columns.end(); it++) {
            out << "\"" << (*it).first << "\" column is " << (*it).second << std::endl;
        }

        return out;
    }

    yamas::MarkerMap yamas::Project::readChromosome(const short int &c) {
        // open input stream depending on the name suffix
        std::ifstream istream; // note: stream must be defined before(!) filtering_stream!!!
        boost::iostreams::filtering_stream<boost::iostreams::input> inputfile;
        if (fname.find(".gz") != std::string::npos) { // dealing with a gzipped file?
            istream.open(fname.c_str());
            if (!istream) yamas::ioerror(fname);
            inputfile.push(boost::iostreams::gzip_decompressor());
            inputfile.push(istream);
        }
        else {
            istream.open(fname.c_str());
            if (!istream) yamas::ioerror(fname);
            inputfile.push(istream);
        }
        // return value: map<string, marker>
        yamas::MarkerMap markertable;
        // temporary vector of Markers (per chromosome and project only a single
        // marker is stored
        std::vector<yamas::MetaMarker> temp_marker;
        // temporary marker
        yamas::MetaMarker marker;
        // temporary dummy line 
        std::string line;
        // temporary dummy vector
        std::vector<std::string> fields;
        // internal representation of the required chromosome
        short int ic = c;

        // check whether the user wanted no chromosome
        bool check_chromosome = true;
        if (columns["chromosome"] == yamas::WRONG - 1 and ud.algorithm != "fillwithproxies" and
                                                          ud.algorithm != "ldblockwise") { 
            check_chromosome = false;
            ic = -1;
        }
        else if (columns["chromosome"] == yamas::WRONG - 1) {
            yamas::error("Unable to parse the 'chromosome' column for '" + fname + "'.");
        }
        // check whether the user is using 1/2 coded allele
        // (== no other allele)
        bool check_other_allele = true;
        if (columns["other_allele"] == yamas::WRONG - 1 and ud.algorithm != "fillwithproxies" and
                                                            ud.algorithm != "ldblockwise") {
            check_other_allele = false;
            yamas::warning("You did not specify 'other_alleles' for '" + fname + "'. No save direction check possible!");
        }
        else if (columns["other_allele"] == yamas::WRONG - 1) {
            yamas::error("You did not specify 'other_alleles' for '" + fname + "'. This parameter is required for the selected algorithm.");
        }
        // check for position cols
        if (columns["position"] == yamas::WRONG - 1 and ud.algorithm != "fillwithproxies" and
                                                        ud.algorithm != "ldblockwise") {
            yamas::warning("You did not specify 'positions' for '" + fname + "'.");
        }
        else if (columns["other_allele"] == yamas::WRONG - 1) {
            yamas::error("You did not specify 'positions' for '" + fname + "'. This parameter is required for the selected algorithm.");
        }

        // check all other col positions
        for (std::map<std::string, short unsigned int>::const_iterator cit = columns.begin(); cit != columns.end(); cit++) {
             // this key already checked?   
             if ((*cit).first == "chromosome" or (*cit).first == "other_allele" or 
                 (*cit).first == "position") 
                 continue;
             if ((*cit).second == yamas::WRONG - 1) {
                 yamas::error("Unable to parse '" + (*cit).first + "' for '" 
                                                  + fname + "'." +
                              " (The parameter name may deviate from the configration file.)");
             }
        }

        // line number within this context
        unsigned long int lnumber = 0;
        // position to be used as key
        unsigned int position = 0;
        // current chromosome
        unsigned short int cc = 0;
	// temporary numerical values
	double ceffect, cweight, cpvalue;
	// temporary character values
	char ceffect_allele, cother_allele;

        // skip lines of the input file, if required
        if (ud.skiplines > 0) {
            for (size_t i = 0; i < ud.skiplines; i++) {\
                std::getline(inputfile, line);
            }
        }
        
        // walk over the file
        while (!inputfile.eof()) {
            std::getline(inputfile, line);
            lnumber++;
            // check how far we are in this project
            if (lnumber < current_line) continue;
            current_line++;
            // skip if first character is '#' or when we have an empty line
            if (not line.size() or line[0] == '#') continue;

            // split the line
            boost::algorithm::split(fields, line, boost::is_any_of("\t,; "), 
                                    boost::token_compress_on);

            // keep this line to prevent nonsense in leading or trailing fields
            fields = yamas::compress(fields);

            // check if beyond the chromosome
            if (check_chromosome) {
                try {
                    cc = boost::lexical_cast<unsigned short int>(fields[columns["chromosome"]]);
                }
                catch (boost::bad_lexical_cast &msg) {
                    std::cerr << "ERROR: Unable to read chromosome column for marker \"" 
                              << fields[columns["marker"]] << "\" in " << fname << "." << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            }
            else {
                cc = 0;
            }
            if (!(ic < 0) && (c > cc)) continue;
            if (!(ic < 0) && (c < cc)) break;

            try {
                position = boost::lexical_cast<unsigned int>(fields[columns["position"]]);
            }
            catch (boost::bad_lexical_cast &msg) {
                std::cerr << "ERROR: Unable to read base pair position for marker \"" 
                          << fields[columns["marker"]] << "\"" << std::endl;
                std::exit(EXIT_FAILURE);
            }

            // The following casts are a bit more tolerant, as we
	    // shouldn't stop during a simple run
            try {
		ceffect = boost::lexical_cast<double>(fields[columns["effect"]]);
        	cweight = boost::lexical_cast<double>(fields[columns["weight"]]);
		cpvalue = boost::lexical_cast<double>(fields[columns["P-value"]]);
	    }
	    catch (boost::bad_lexical_cast &msg) {
		ceffect  =  0;
		cweight  =  0;
		cpvalue  = -1;
	    }

	    try {
		ceffect_allele = boost::lexical_cast<char>(fields[columns["effect_allele"]]);
	    }
	    catch (boost::bad_lexical_cast &msg) {
		ceffect_allele = '0';
	    }

            if (check_other_allele) {
	        try {
		    cother_allele = boost::lexical_cast<char>(fields[columns["other_allele"]]);
	        }
	        catch (boost::bad_lexical_cast &msg) {
	            cother_allele = '0'; 
	        }
            }
            else {
                cother_allele = '0';
            }

            marker = MetaMarker(id, // own project id for reference
                                fields[columns["marker"]], // no casting necessary
                                cc,
                                position,
                                ceffect_allele,
                                cother_allele,
                                ceffect, //boost::lexical_cast<double>(fields[columns["effect"]]),
                                cweight, //boost::lexical_cast<double>(fields[columns["weight"]]),
                                cpvalue //boost::lexical_cast<double>(fields[columns["P-value"]])
                                );
            // in case of 'usual' algos save the marker regardless of the status
            if (ud.algorithm != "fillwithproxies")
                markertable[fields[columns["marker"]]] = marker;
            // otherwise the marker should be valid
            else if (marker.is_valid())
                markertable[fields[columns["marker"]]] = marker;
        }
        //inputfile.close();

        // if we found no suitable marker, we need to reset the lnumber-count:
        current_line--;
        return markertable;
    }

    void yamas::Project::print_header() const {
        // open input stream 
        std::ifstream inputfile(fname.c_str());

        if (!inputfile) yamas::ioerror(fname);

        // temporary dummy line and rsid variables
        std::string line;
        // temporary dummy vector
        std::vector<std::string> fields;
        // identifier names
        std::string ccolnames = "marker,chromosome,position,effect allele,other allele,odds ratio/beta,wheigt,p-value";
        std::vector<std::string> colnames;
        boost::algorithm::split(colnames, ccolnames, boost::is_any_of(","));
        // print out id
        std::cout << "Project " << id << " - " << fname << std::endl;

        // assuming first line we hit is a header line
        std::getline(inputfile, line);
        if (line[0] == '#') {
            // split the line
            boost::algorithm::split(fields, line, boost::is_any_of("\t,; "),
                                    boost::token_compress_on);
            // keep this line to prevent nonsense in leading or trailing fields
            fields = yamas::compress(fields);
            std::cout << std::setw(16) << "Column" << " - " << "You specified : YAMAS read" << std::endl;
            std::cout << std::setw(16) << "------" << "---" << "--------------------------" << std::endl;
            for (std::map<std::string, short unsigned int>::const_iterator it = columns.begin();
                 it != columns.end(); it++) {
                // check whether a valid value is present
                if ((*it).second != yamas::WRONG - 1) {
                    std::cout << std::setw(16)  << (*it).first << " - "
                              << std::setw(13)  << (*it).second + 1 << " : "
                              << fields[(*it).second + 1]
                              << std::endl;
                }
                else {
                    std::cout << std::setw(16) << (*it).first << " - "
                              << std::setw(13) << "None" << " : "
                              << "None"
                              << std::endl;
                }
            }
        }
        else {
            std::cout << "No appropriate header line found." << std::endl;
        }
    }

    bool yamas::Project::check() const {
        for (std::map<std::string, unsigned short int>::const_iterator it = columns.begin(); 
             it != columns.end(); it++) {
            if (isnan((*it).second)) {
                std::cerr << "ERROR: '" << (*it).first << "' column in Project " << id
                          << " not recognized but mandatory." << std::endl;
                return false;
            }
        }
        return true;
    }

/**********************************************************/
/* non-member functions from hereon                       */
/**********************************************************/
    // defining a function to produce a vector of Projects
    std::vector<yamas::Project> extractor (boost::program_options::variables_map &vm) {
        // define a return value
        std::vector<yamas::Project> projects;
        // convert project file names
        std::vector<std::string> assocfilenames;
        yamas::convert(vm["assocfiles"].as<std::string>(), assocfilenames);
        std::string test = vm["assocfiles"   ].as<std::string>();
        #ifdef INCLUDEPRODUCTION
        std::vector<double> project_weights;
        // check whether 'study weights' are available
        if (vm.find("study_weights") != vm.end()) {
            std::vector<std::string> tmp_weights;
            yamas::convert(vm["study_weights"].as<std::string>(), tmp_weights);
            for (size_t i = 0; i < tmp_weights.size(); i++) {
                project_weights.push_back(boost::lexical_cast<double>(tmp_weights[i]));
            }
            if (project_weights.size() != assocfilenames.size())
                yamas::error("No. of filenames does not match No. of given weights.");
        }
        else if (ud.use_study_weights) {
            yamas::error("Requested weighting of studies, but no 'study_weights' argument in configuration file.");
        }
        #endif
        // convert marker ids columns
        std::vector<short unsigned int> marker_cols;
        try {
            yamas::convert(vm["marker_cols"].as<std::string>(), marker_cols);
        }
        catch (boost::bad_lexical_cast &msg) {
	    yamas::error("Unable to read configuration file: wrong marker columns.");
	}
        // convert chromosome columns
        std::vector<short unsigned int> chromosome_cols;
        try {
            yamas::convert(vm["chromosome_cols"].as<std::string>(), chromosome_cols);
        }
        catch (boost::bad_lexical_cast &msg) {
            yamas::error("Unable to read configuration file: wrong column columns.");
        }
        // convert position columns
        std::vector<short unsigned int> position_cols;
        try {
            yamas::convert(vm["position_cols"].as<std::string>(), position_cols);
        }
        catch (boost::bad_lexical_cast &msg) {
            yamas::error("Unable to read configuration file: wrong position columns.");
        }
        // convert effect allele columns
        std::vector<short unsigned int> effect_allele_cols;
        try {
            yamas::convert(vm["effect_allele_cols"].as<std::string>(), effect_allele_cols);
        }
        catch (boost::bad_lexical_cast &msg) {
            yamas::error("Unable to read configuration file: wrong effect allele columns.");
        }
        // convert other allele columns
        std::vector<short unsigned int> other_allele_cols;
        try {
            yamas::convert(vm["other_allele_cols"].as<std::string>(), other_allele_cols);
        }
        catch (boost::bad_lexical_cast &msg) {
	    yamas::error("Unable to read configuration file: wrong 'other' allele columns.");
        }
        // convert OR columns
        std::vector<short unsigned int> effect_cols;
        try {
            yamas::convert(vm["effect_cols"].as<std::string>(), effect_cols);
        }
        catch (boost::bad_lexical_cast &msg) {
	    yamas::error("Unable to read configuration file: wrong effect columns.");
        }
        // convert stderr columns
        std::vector<short unsigned int> weight_cols;
        try {
            yamas::convert(vm["weight_cols"].as<std::string>(), weight_cols);
        }
        catch (boost::bad_lexical_cast &msg) {
	    yamas::error("Unable to read configuration file: wrong weight columns.");
        }
        // convert p-value columns
        std::vector<short unsigned int> pvalue_cols;
        try {
            yamas::convert(vm["P-value_cols"].as<std::string>(), pvalue_cols);
        }
        catch (boost::bad_lexical_cast &msg) {
            yamas::error("Unable to read configuration file: wrong P-value columns.");
        }

        // check that all vectors have the same size as the assoc-file names
        if (assocfilenames.size() != marker_cols.size()) {
            yamas::error("No. of marker columns does not equal no. of file names.");
        }
        if (assocfilenames.size() != chromosome_cols.size()) {
            yamas::error("No. of chromosome columns does not equal no. of file names.");
        }
        if (assocfilenames.size() != position_cols.size()) {
	    yamas::error("No. of position columns does not equal no. of file names.");
        }
        if (assocfilenames.size() != effect_allele_cols.size()) {
            yamas::error("No. of effect allele columns does not equal no. of file names.");
        }
        if (assocfilenames.size() != other_allele_cols.size()) {
            yamas::error("No. of 'other' allele columns does not equal no. of file names.");
        }
        if (assocfilenames.size() != effect_cols.size()) {
            yamas::error("No. of effect columns does not equal no. of file names.");
        }
        if (assocfilenames.size() != weight_cols.size()) {
	    yamas::error("No. of weight columns does not equal no. of file names.");
        }
        if (assocfilenames.size() != pvalue_cols.size()) {
            yamas::error("No. of P-value columns does not equal no. of file names.");
        }

        projects.reserve(assocfilenames.size());
        std::map<std::string, short unsigned int> tmp;
        for (std::size_t i = 0; i < assocfilenames.size(); i++) {
            tmp["marker"]           = marker_cols[i] - 1;
            tmp["chromosome"]       = chromosome_cols[i] - 1;
            tmp["position"]         = position_cols[i] - 1;
            tmp["effect_allele"]    = effect_allele_cols[i] - 1;
            tmp["other_allele"]     = other_allele_cols[i] - 1;
            tmp["effect"]           = effect_cols[i] - 1;
            tmp["weight"]           = weight_cols[i] - 1;
            tmp["P-value"]          = pvalue_cols[i] -1;
            #ifdef INCLUDEPRODUCTION
            projects.push_back(yamas::Project(tmp, assocfilenames[i], project_weights[i]));
            #else
            projects.push_back(yamas::Project(tmp, assocfilenames[i]));
            #endif
            tmp.clear();
        }
        return projects;
    }
} // end namespace yamas


