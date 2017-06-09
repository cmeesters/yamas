/********************************************************************/
/*                                                                  */
/* YAMAS - Yet Another Meta-Analysis Software                       */
/*                                                                  */
/* Classes holding LD-propertis used in  meta-analysis.             */
/*                                                                  */
/********************************************************************/

// $Rev: 916 $
// $LastChangedDate: 2012-06-15 21:58:55 +0200 (Fr, 15 Jun 2012) $
// $Author: meesters $

#include <iostream>
#include <fstream>
#include <iomanip> // manipulator for fixed column width formatting
#include <vector>
#include <string>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

// include own software parts
#include "../include/ya_lddata.hpp"
#include "../include/ya_errors.hpp"
#include "../include/ya_io.hpp"

// copy-constructor
yamas::Proxy::Proxy(const yamas::Proxy &p) : 
                    marker_id_(p.marker_id_),
                    allele_(p.allele_),
                    inphase_allele_(p.inphase_allele_),
                    distance_(p.distance_),
                    r2_(p.r2_),
                    switch_(p.switch_){}

yamas::Proxy::Proxy(const std::string  &marker_id,
		    const char         &allele,
		    const char         &inphase_allele,
		    const unsigned int &distance,
		    const double       &r2,
                    const bool         &switchr) {
    marker_id_        = marker_id;
    allele_           = allele;
    inphase_allele_   = inphase_allele;
    distance_         = distance;
    r2_               = r2;
    switch_           = switch_;
}

yamas::Proxy::~Proxy() {}

namespace yamas {
    std::ostream& operator<< (std::ostream& out, 
                              const Proxy& proxy) {
        // save flags
        std::ios_base::fmtflags flags = out.flags();

        out << std::right    << std::setw(12) << proxy.get_marker_id() 
	    << std::setw(12) << proxy.get_distance()
	    << std::setw(4)  << proxy.get_allele()
	    << std::setw(4)  << proxy.get_inphase_allele()
	    << std::setw(6)  << proxy.get_r2()
            << std::setw(3)  << proxy.get_switch();
        // restore flags 
        out.flags(flags);
        return out;
    }

    yamas::BlocksLT readldblocks(const std::string &fname) {
        yamas::BlocksLT finalblocks; // return value
	// temporary map of block: {begin, end}
	yamas::LDBlocksMap ldblocks;
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

	// string containing the file lines
	std::string line;
	// vector of string fields
	std::vector<std::string> fields;
	// several integers, which we need to store
	unsigned int minpos = 0, maxpos = 0, bin = 0;
        // chromosome
	short int cc = 0, pc = 1;

	if (ud.verbose) std::cout << "Reading LD-block data from '" << fname
		                  << "'." << std::endl;

        // walk over the file
	while (std::getline(inputfile, line)) {
            boost::algorithm::split(fields, line, boost::is_any_of("\t,; "),
	                           boost::token_compress_on);
            // keep this line to prevent nonsense in leading or trailing fields
            fields = yamas::compress(fields);
	    try {
                cc     = boost::lexical_cast<short int>(fields[0]);
	        minpos = boost::lexical_cast<unsigned int>(fields[2]);
	        maxpos = boost::lexical_cast<unsigned int>(fields[3]);
	        bin    = boost::lexical_cast<unsigned int>(fields[1]);
	    }
	    catch (boost::bad_lexical_cast &msg) {
                std::string message = "Unable to parse the LD-block file, parser says: ";
		yamas::error(message + msg.what());
	    }
           
            ldblocks[bin] = boost::tuple<unsigned int, unsigned int>(minpos, maxpos);

            // got to the next chromosome
            if (cc != pc) {
                // store data
                finalblocks[pc] = ldblocks;
                // erase temporary data
                ldblocks.erase(ldblocks.begin(), ldblocks.end());
                // previous = current
                pc = cc;
            }
	}
        // finally 
        finalblocks[cc] = ldblocks;
        ldblocks.erase(ldblocks.begin(), ldblocks.end());

	return finalblocks;
    }

    yamas::ProxyLT readldfile(const std::string &fname, 
                              const yamas::ChrMap &table,
                              const short int &c) {
        // return value
        yamas::ProxyLT proxies;
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

        // r2-value threshold
        const double r2threshold = ud.r2threshold;

        // string containing file lines
	std::string line;
        // line number within this context
        unsigned long int lnumber = 1;
	// vector of string fields
	std::vector<std::string> fields;
        // temporary markers
        std::string marker, marker2;
	// temporary r2
	double r2;
	// temporary proxy assignment
	bool assign;
	// alleles
	char allele = ' ', inphase = ' ', allele2 = ' ', inphase2 = ' ';
        // temporary physical distance
        unsigned int distance = 0;
	// temporary Proxy
	Proxy proxy, proxy2;

        if (not (table.find(c) != table.end())) {
            yamas::warning("Data for the required chromosome not found in the reference file.");
            return proxies;
        }

        long unsigned int begin = (table.find(c))->second.get<0>();
        long unsigned int end   = (table.find(c))->second.get<1>();

        // read in as many lines as to begin with the required chromosome
        if (lnumber < begin) {
            while (std::getline(inputfile, line)) {
                lnumber++;
                if (lnumber == begin) break;
            }
        }

        // walk over the file
        for (std::size_t index = begin; index <= end; index++) {
            std::getline(inputfile, line);
            // split the line
            boost::algorithm::split(fields, line, boost::is_any_of("\t "), 
                                    boost::token_compress_on);
            // keep this line to prevent nonsense in leading or trailing fields
            fields = yamas::compress(fields);
            // skip when we have an empty line
            if (not fields.size() > 1) continue;
            // parse values
	    r2 = boost::lexical_cast<double>(fields[8]);
            // don't proceed, if LD is below threshold
	    if (r2 > r2threshold) {  
                #pragma omp sections
                {
                #pragma omp section
                {
	        marker   = boost::lexical_cast<std::string>(fields[1]);
	        marker2  = boost::lexical_cast<std::string>(fields[4]);
                distance = boost::lexical_cast<unsigned int>(fields[7]);
                }
                #pragma omp section
                {
	        assign   = boost::lexical_cast<bool>(fields[9]);
	        // select allele: if assign == 0, the first one
                allele   = boost::lexical_cast<char>(fields[2]);
                inphase  = assign ? boost::lexical_cast<char>(fields[6]) :
		                    boost::lexical_cast<char>(fields[5]);
                allele2  = boost::lexical_cast<char>(fields[5]);
                inphase2 = assign ? boost::lexical_cast<char>(fields[3]) :
		                    boost::lexical_cast<char>(fields[2]);
                }
                }
                proxy = Proxy(marker2, allele, inphase, distance, r2, assign);
                proxies[marker].push_back(proxy);

                proxy2 = Proxy(marker, allele2, inphase2, distance, r2, assign);
                proxies[marker2].push_back(proxy2);
	    }
        } // end of for loop
        
        // sorting all vectors of proxies in the proxy map
        std::vector<std::string> keys = yamas::get_keys(proxies);
        std::vector<std::string>::const_iterator key;
        #pragma omp parallel for private (key)
        for (key = keys.begin(); key < keys.end(); key++) {
            if (proxies[*key].size() > 1) { // only sort, if 'sortable'
                std::stable_sort(proxies[*key].begin(), proxies[*key].end());
                std::reverse(proxies[*key].begin(), proxies[*key].end());
	    }
        }
        lnumber--;
        return proxies;
    }

    ChrMap proxy_parse(const std::string &fname) {
        // return value
        yamas::ChrMap table;

        // check for the test flag to save time, when testing
        #ifdef TEST
        table[1] = boost::tuple<unsigned long int, unsigned long int>(1, 1000000);
        return table;
        #endif // TEST

        // local instream
	std::ifstream istream;
	// boost variant
	boost::iostreams::filtering_stream<boost::iostreams::input> infile;
        if (fname.find(".gz") != std::string::npos) { // dealing with a gzipped file?
	    istream.open(fname.c_str());
            if (!istream) yamas::ioerror(fname);
            infile.push(boost::iostreams::gzip_decompressor());
            infile.push(istream);
        }
        else {
            istream.open(fname.c_str());
            if (!istream) yamas::ioerror(fname);
            infile.push(istream);
        }

	// string containing lines
	std::string line;
	// vector containing the fields of a line
	std::vector<std::string> fields;
        // local line numbers
	unsigned long int lnumber = 0, begin = 1, end = 0;
	// current chromosome
	short int cc;
	cc = ud.chromosome != 0 ? ud.chromosome : 1;
	// previous chromosome
	short int pc;
	pc = cc;
        // desired chromosome
        short int dc;
        dc = ud.chromosome != 0 ? ud.chromosome : 100;

	// walk over the file
	while(std::getline(infile, line)) {
	    // increment line number under all circumstances
	    lnumber++;
	    // split the line to be able to look for the chromosome field
	    boost::algorithm::split(fields, line, boost::is_any_of("\t,; "),
			            boost::token_compress_on);
            fields = yamas::compress(fields);
	    // apparently necessary, because screening might otherwise fail
	    if (not (fields.size() > 1)) break;
            // current chromosome
            try {
	        cc = boost::lexical_cast<short int>(fields[0]);
	    }
            catch (boost::bad_lexical_cast &msg) {
		std::cerr << "ERROR: Apparently the proxy file is not " 
			  << "correctly formatted" << std::endl;
		std::exit(EXIT_FAILURE);
	    }
	    // did we hit a new chromosome?
	    if (cc != pc) {
		end       = lnumber - 1;
		table[pc] = boost::tuple<unsigned long int, unsigned long int>(begin, end);
		// new begin
		begin     = lnumber;
		// new pc
		pc        = cc;
	    }
            // are we beyond the required chromosome?
            if (cc > dc) break;
	}
	// add the last chromosome data
	end = lnumber - 1;
	table[pc] = boost::tuple<unsigned long int, unsigned long int>(begin, end);

	return table;
    }
} // end namespace yamas

