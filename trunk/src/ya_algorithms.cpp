/********************************************************************/
/*                                                                  */
/* YAMAS - Yet Another Meta-Analysis Software                       */
/*                                                                  */
/* algorithms                                                       */
/*                                                                  */
/********************************************************************/

// $Rev: 929 $
// $LastChangedDate: 2012-09-17 10:00:37 +0200 (Mo, 17 Sep 2012) $
// $Author: meesters $

#include <vector>
#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iterator>

#include <boost/algorithm/string.hpp>

#include <omp.h>

#include "../include/ya_io.hpp"
#include "../include/ya_metamarker.hpp"
#include "../include/ya_project.hpp"
#include "../include/ya_utils.hpp"
#include "../include/ya_interface.hpp"
#include "../include/ya_global.hpp"
#include "../include/ya_lddata.hpp"
#include "../include/ya_stats.hpp"

#ifdef SANDBOX
#include "../sandbox/include/ya_plotting.hpp"
#endif

// turns on memory profiling on some(!) Linux systems
// must be enabled using scons
#ifdef MEMORY
#include <sys/resource.h>
#endif

EXTERN yamas::UsageData ud;

namespace yamas {
    //! list of available chromosomes
    short int carr[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
    std::vector<short int> chromosomes(carr, carr + 26);
    short int chromosome = 0;
    /*!
       performs marker- or point-wise meta-analysis

       input: vector of yamas::Project - each project
       input: std::string file name for the log file
       input: std::bool flag for verbosity - true, if verbose
       input: unsigned int chromosome - if 0 all chromosomes will be evaluated
    */
    void pointwise_comparison(const std::vector<yamas::Project> &projects) {
        if (ud.verbose) std::cout << "Setting up point-wise meta-analysis ..." << std::endl;
        // make non-const vector of projects (copy)
        std::vector<yamas::Project> iprojects;
        for (std::vector<yamas::Project>::const_iterator it = projects.begin(); it < projects.end(); it++) {
            iprojects.push_back(*it);
        }
        unsigned int project_size = projects.size(); // number of projects
        // assign a particular chromosome, if requested
        if (ud.chromosome != 0) {
            chromosomes.clear();
            chromosomes.push_back(ud.chromosome);
        }


        // make up string of 'directions'
        std::string directions;
	directions.reserve(iprojects.size());

        // map of markers and a vector of those maps
        yamas::MarkerMap chrlt;
        yamas::MMarkersLT multilt;
        yamas::MMarkersPos mpos;
        // output result marker
        yamas::ResultMarker resmarker; // will be used as a dummy
        // will hold the total count of analysed markers
        unsigned long int total = 0;

        // create ostreams
        std::ofstream outfile((ud.logfname + ".log").c_str());
        std::ofstream errfile((ud.logfname + ".err").c_str());

        // write out generalized header
        yamas::write_header(iprojects, outfile, ud);
        yamas::write_errorlog_header(iprojects, errfile);
        // get a statistics pointer
        yamas::Statistics* sp = yamas::Statistics::get_Statistics();

        #ifdef SANDBOX
        // enable gathering plotting data
        PValues pvalues;
        std::vector<double> effects;
        std::vector<MarkerEffects> markereffects;
        // transform ud.selection in a vector of strings to be searched in
        std::vector<std::string> selected_markers;
        boost::algorithm::split(selected_markers, ud.selection, boost::is_any_of(",; "),
                                boost::token_compress_on);
        #endif


        // read in chromosome-wise
        for (std::vector<short int>::const_iterator c = chromosomes.begin();
             c != chromosomes.end(); c++) {
            if (ud.verbose and *c >= 0) std::cout << "Reading data for chromosome " << *c << std::endl;
            else if (ud.verbose and *c == -1) std::cout << "Reading all data at once" << std::endl;
            #pragma omp parallel for ordered private(chrlt)
            for (std::vector<yamas::Project>::size_type i = 0; i < iprojects.size(); i++) {
                chrlt = iprojects[i].readChromosome(*c);
                #pragma omp ordered
                {
                if (ud.verbose) std::cout << chrlt.size() << " markers found in project " << iprojects[i].id << std::endl;
                // merge the chromosome wide lookup table with the entire one
                for (yamas::MarkerMap::iterator ci = chrlt.begin(); ci != chrlt.end(); ci++) {
                    multilt[(*ci).first].push_back((*ci).second);
                }
                }
                // get rid of data in chrlt
                chrlt.erase(chrlt.begin(), chrlt.end());
            }
            if (ud.verbose) std::cout << "A total of " 
                                      << multilt.size() 
                                      << " unique markers were found." << std::endl;

            // start calculating
            if (ud.verbose) std::cout << "Sorting ... "; std::cout.flush();

            // transfer rsid:markers-map to position:markers-map
            for (yamas::MMarkersLT::const_iterator mit = multilt.begin(); mit != multilt.end(); mit++) {
                // retrieve first marker's position
                mpos[(*mit).second[0].get_position()] = (*mit).second;
            }
            // erase multilt to save space
            multilt.erase(multilt.begin(), multilt.end());
            // drop the position keys - we only need the sorted vector of vectors
            std::vector<std::vector<yamas::MetaMarker> > markers_only;
            markers_only.reserve(mpos.size());
            for (yamas::MMarkersPos::const_iterator mit = mpos.begin(); mit != mpos.end(); mit++){
                markers_only.push_back((*mit).second);
            }
            mpos.erase(mpos.begin(), mpos.end());

            if (ud.verbose) std::cout << "Done" << std::endl; std::cout.flush();

            // trim output for first marker set, only?
            if (ud.trim.size() > 0) yamas::trim_markers(markers_only);

            if (ud.verbose) std::cout << "Calculating ... ";
            std::cout.flush();

            for (std::vector<std::vector<yamas::MetaMarker> >::iterator mit = markers_only.begin(); 
                 mit < markers_only.end(); mit++) {
                // first: save the rsID to the resmarker
                resmarker.set_rsid((*mit)[0].get_rsid());

                // check all markers for consistency, correct, if necessary
                yamas::check_markers(*mit);

                // start calculations
                #pragma omp sections
                {
                #pragma omp section
                resmarker.fixed_effect(*mit);
                #pragma omp section
                resmarker.random_effect(*mit);
                #ifdef INCLUDEPRODUCTION
                #pragma omp section
                resmarker.stouffers_method(*mit, projects);
                #endif
                }
                resmarker.check_directions(*mit, project_size);
                check_reasonable_pvalue(resmarker);
 
                #ifdef SANDBOX
                if (ud.plot) {
                    #pragma omp sections
                    {
                    #pragma omp section
                    {
                    pvalues.append_pos(*c, resmarker.get_position());
                    pvalues.append_val_fixed(*c, resmarker.get_pvalue());
                    pvalues.append_val_random(*c, resmarker.get_rpvalue());
                    } // section 1
                    #pragma omp section
                    {
                    if (selected_markers.size() > 0) { // any markers left?
                        for (std::vector<std::string>::iterator sit = selected_markers.begin();
                             sit != selected_markers.end(); sit++) {
                           if ((*sit) == resmarker.get_rsid()) { // proceed only, if equal
                              markereffects.push_back(yamas::MarkerEffects(*mit,
                                                                           resmarker));
                              selected_markers.erase(sit); // erase this string, as we already saved the info
                              break;
                           }
                        }

                    }
                    } // section 2
                    }
                }
                #endif

                // write to log file
                // clear previous 'directions' content, such that
                // 'foul' markers cannot mess up the output
                directions.clear();
                directions = resmarker.get_directions();
                
                if ((resmarker.get_pvalue() <= ud.pvaluethreshold or resmarker.get_rpvalue() <= ud.pvaluethreshold) &&
                    // meta-analysis was possible? (at least two valid markers?)
                    directions.size() - (int)count(directions.begin(), directions.end(), '?')
                                      - (int)count(directions.begin(), directions.end(), 'x') >= 2) {
                    // output status (OK, error, or warnings) 
                    resmarker.set_all(*mit, outfile);
                    outfile << resmarker;
                    // write the input markers, only if requested
                    if (ud.write_input) {
                        write_vector(*mit, outfile, ""); // using ya_io template function
                    }
                    outfile << std::endl;
                    // apparently we hid a valid marker
                    total++;                      // total number of valid markers
                    // update statistics
                    if (ud.stats) sp->update(resmarker);
                }
                else {
                    errfile << std::setw(12) << (*mit)[0].get_rsid();
                    write_vector(*mit, errfile, "");
                    errfile << std::endl << std::flush;
                }
		(*mit).erase((*mit).begin(), (*mit).end());
            }

            if (ud.verbose) std::cout << "Done." << std::endl << std::endl;
            // finally, clean up all resources
            markers_only.erase(markers_only.begin(), markers_only.end());
            // iteriert bis hierher 9 mal, warum? ML
        }
        outfile.close();
	errfile.close();


        // The following output should be done regardless of the
        // requested verbosity:
        std::cout << "YAMAS analysed " << total << " valid markers." << std::endl;
        // write out statistics, if requested
        if (ud.stats) {
            std::cout << std::endl;
            sp->writeout(total);
        }

        #ifdef SANDBOX
        // finally plot, if required
        if (ud.plot) {
           if (ud.verbose) std::cout << std::endl << "Producing required plots." << std::endl;
           pvalues.sort();
           #pragma omp parallel sections
           {
           #pragma omp section
           yamas::manhattanplot(pvalues, true);  // fixed P-values
           #pragma omp section
           yamas::manhattanplot(pvalues, false); // random P-values
           #pragma omp section
           yamas::forestplot_gwas(markereffects);  
           }
           if (ud.verbose) 
              std::cout << "Done." << std::endl;
        }
        #endif
    }

    /*!
       performs meta-analysis per LD block

       input: vector of yamas::Project - each project
       input: string - file name containing LD-info
    */
    #ifdef INCLUDEPRODUCTION
    void per_ldblock(const std::vector<yamas::Project> &projects, 
		     const std::string &ldfname,
		     const std::string &proxyfname) {
        // if this compile flag is given, we open a special log file for profiling!
        #ifdef MEMORY
        std::ofstream memlog("yamas_memory.log");
        rusage ru;
        #endif
        if (ud.verbose) std::cout << "Setting up LD-block-wise meta-analysis ..." << std::endl;
        // make non-const vector of projects (copy)
        std::vector<yamas::Project> iprojects;
        for (std::vector<yamas::Project>::const_iterator it = projects.begin(); it < projects.end(); it++) {
            iprojects.push_back(*it);
        }
        unsigned int project_size = projects.size(); // number of projects

        // create a vector of missing project IDs
        std::vector<unsigned int> missings;
        missings.reserve(project_size);

        // assign a particular chromosome, if requested
        if (ud.chromosome != 0) {
            chromosomes.clear();
            chromosomes.push_back(ud.chromosome);
        }

        // make up string of 'directions'
        std::string directions;
        directions.reserve(project_size);

        // map of markers and a vector of those maps
        yamas::MarkerMap chrlt;
        yamas::MMarkersLT multilt;
        yamas::MMarkersPos mpos;
	yamas::LDBlocksMap ldblocks_per_chromosome;
        yamas::BlocksLT ldblocks;

        // output result marker
        yamas::ResultMarker resmarker; // will be used as a dummy
        yamas::MetaMarker marker;      // another dummy
        yamas::MetaMarker proxymarker;
        // temporary vector of markers
        std::vector<yamas::MetaMarker> proxymarkers;
        proxymarkers.reserve(project_size);
        // will hold the total count of analysed markers
        unsigned long int total = 0;

        // create ostreams
        std::ofstream outfile((ud.logfname + ".log").c_str());
        std::ofstream errfile((ud.logfname + ".err").c_str());

        // write out generalized header
        yamas::write_header(iprojects, outfile, ud, false);
        yamas::write_errorlog_header(iprojects, errfile);

        // read in ld blocks -- for the entire genome at once (should easily fit in memory)
        if (ud.verbose) std::cout << "Reading LD blocks from '" << ldfname << "'." << std::endl;
        #ifdef MEMORY
        #pragma omp master
        {
        getrusage(RUSAGE_SELF, &ru);
        memlog << "before reading LD-blockfile: " << ru.ru_maxrss << " kB." << std::endl;
        }
        #endif
        ldblocks = yamas::readldblocks(ldfname);
        if (ud.verbose) std::cout << "Done" << std::endl;
        #ifdef MEMORY
        #pragma omp master
        {
        getrusage(RUSAGE_SELF, &ru);
        memlog << "after reading LD-blockfile: " << ru.ru_maxrss << " kB." << std::endl;
        }
        #endif

	// parse the proxyfile to get a look-up table: {chromosome: (begin, end)}
        yamas::ChrMap table;

        if (ud.verbose) {
            std::cout << "Screening reference file: " << proxyfname << std::endl;
	    std::cout << "(Screening is necessary to speed up the read-in, subsequently.)"
		      << std::endl
                      << "(The screening itself is always master-threaded.)" << std::endl
                      << std::endl;
	}
        #ifdef MEMORY
        #pragma omp master
        {
        getrusage(RUSAGE_SELF, &ru);
        memlog << "before reading proxy file: " << ru.ru_maxrss << " kB." << std::endl;
        }        
        #endif
	table = yamas::proxy_parse(proxyfname);
        #ifdef MEMORY
        #pragma omp master
        {
        getrusage(RUSAGE_SELF, &ru);
        memlog << "after reading proxy file: " << ru.ru_maxrss << " kB." << std::endl;
        }        
        #endif

	yamas::ProxyLT proxies;
        std::vector<Proxy> proxy_vector;
        std::vector<Proxy>::const_iterator pit;
        std::string rsid;
        // process indicator
        bool done;
        // get a statistics pointer
        yamas::Statistics* sp = yamas::Statistics::get_Statistics();

        std::vector<yamas::MetaMarker>::const_iterator mmit, mmit2;
        std::vector<yamas::Proxy>::const_iterator proxy_it;

        // read in chromosome-wise
        for (std::vector<short int>::const_iterator c = chromosomes.begin();
             c != chromosomes.end(); c++) {
	    // read in LD-data
            if (ud.verbose) std::cout << "Reading SNP-Proxy-data from " << proxyfname << std::endl;
	    proxies = readldfile(proxyfname, table, *c);

            if (ud.verbose and *c >= 0) std::cout << "Reading data for chromosome " << *c << std::endl;
            else if (ud.verbose and *c == -1) std::cout << "Reading all data at once" << std::endl;
            #pragma omp parallel for ordered private(chrlt)
            for (std::size_t i = 0; i < iprojects.size(); i++) {
                chrlt = iprojects[i].readChromosome(*c);
                #pragma omp ordered
                {
                if (ud.verbose) std::cout << chrlt.size() << " markers found in project " << iprojects[i].id << std::endl;
                // merge the chromosome wide lookup table with the entire one
                for (yamas::MarkerMap::iterator ci = chrlt.begin(); ci != chrlt.end(); ci++) {
                    multilt[(*ci).first].push_back((*ci).second);
                }
                }
                // get rid of data in chrlt
                chrlt.erase(chrlt.begin(), chrlt.end());
            }
            if (ud.verbose) std::cout << "A total of " 
                                      << multilt.size()
                                      << " unique markers were found." << std::endl;

            // start calculating
            if (ud.verbose) std::cout << "Sorting ... "; std::cout.flush();

            // transfer rsid:markers-map to position:markers-map
            for (yamas::MMarkersLT::const_iterator mit = multilt.begin(); 
                mit != multilt.end(); mit++) {
                // retrieve first marker's position
                mpos[(*mit).second[0].get_position()] = (*mit).second;
            }
	    
            // drop the position keys - we only need the sorted vector of vectors
            std::vector<std::vector<yamas::MetaMarker> > markers_only;
	    markers_only.reserve(mpos.size());

            for (yamas::MMarkersPos::const_iterator mit = mpos.begin();
	         mit != mpos.end(); mit++) {
                markers_only.push_back((*mit).second);
            }
            mpos.erase(mpos.begin(), mpos.end());

            // trim output for first marker set, only?
            if (ud.trim.size() > 0) yamas::trim_markers(markers_only);

            if (ud.verbose) std::cout << "Done" << std::endl
		                      << "Calculating ... "; std::cout.flush();

            // copy LD-blocks for this chromosome
            ldblocks_per_chromosome = ldblocks[*c];

	    // associate markers with ldblock
            std::vector<yamas::MetaMarker> markers;
            markers.reserve(project_size);

            for (std::map<unsigned int, yamas::LDBlock>::const_iterator bit = ldblocks_per_chromosome.begin(); 
                 bit != ldblocks_per_chromosome.end(); bit++) {
                markers.clear();
                // ensure dummy markers in the markers vector
                for (size_t i = 0; i < project_size; i++) {
                    markers.push_back(yamas::MetaMarker());
                }
                for (size_t i = 0; i < markers_only.size(); i++) {
                    // check whether the current marker is within the ld-block
	            // do this check only for the first marker in the vector
	            // (they all should have the same position ...)
		    // if not, we don't care for that marker(s)
                    if (!yamas::isin(markers_only[i][0], (*bit).second)) continue;
                    // rsid of the marker we are dealing with
                    rsid = markers_only[i][0].get_rsid(); // [0] because validity does not matter for the time being
                    // is a projects marker missing?
                    yamas::reset(missings, project_size);
                    for (mmit = multilt[rsid].begin(); mmit != multilt[rsid].end(); mmit++) {
                        missings[(*mmit).get_project_id() - 1]++; // the -1 is because we need the position
                                                                  // and YAMAS counts from 1, vector from 0
                    }
                    // now insert the proxy marker for each missing marker per project

                    // for this we need the appropriate vector of proxies, first
                    proxy_vector = proxies[rsid];
                    // and we iterate over all project IDs with missing markers
                    for (size_t mindex = 0; mindex < missings.size(); ++mindex) {
                        if (missings[mindex]) continue; // did not actually hit a missing marker?
                        // do we actually possess proxies for this missing marker?
                        if (not proxy_vector.size()) continue;
                        // the we iterate over the proxy vector to see, whether or not we have a substitute
                        for (proxy_it = proxy_vector.begin(); proxy_it != proxy_vector.end(); proxy_it++) {
                           // key present?
                           if (multilt.find(rsid) != multilt.end()) {}
                               // retrieve vector of markers with that key
                               proxymarkers = multilt[(*proxy_it).get_marker_id()];
                               // now we ask whether this marker is actually present for the required project
                               for (mmit2 = proxymarkers.begin(); mmit2 != proxymarkers.end(); mmit2++) {
                                  // skip marker, if that one is not valid within the project
                                  if (not (*mmit2).is_valid()) continue;
                                  proxymarker = (*mmit2);
                                  // do we hit the right project?
                                  if (proxymarker.get_project_id() == (mindex + 1)) {
                                      // now, we just set the allele and check later:
                                      proxymarker.set_proxy_allele((*proxy_it).get_allele());
                                      proxymarker.set_inphase_allele((*proxy_it).get_inphase_allele());
                                      // assign the proxy marker's switch status
                                      proxymarker.set_switch_status((*proxy_it).get_switch());
                                      // remember the r2!
                                      proxymarker.set_r2((*proxy_it).get_r2());
                                      // then insert that very marker into the vector of markers 
                                      // in place of the missing one
                                      markers_only[i].push_back(marker);
                                      done = true;
                                      break;
                                  }
                               }
                        }
                        if (done) break;
                    }
                    done = false;
                    
                    // for each project: find the marker with the smallest p-value
                    for (size_t j = 0 ; j < project_size; j++) {
                        // marker we are dealing with
                        marker = markers_only[i][j];
                        // this expression holds, because dummy markers have p-values = 1
                        if (marker.is_valid() and
		            markers[j].get_pvalue() > marker.get_pvalue()) {
			    markers[j] = marker;
                        }
                    }
		}

                // first: save the rsID to the resmarmker -- we are
                // using the block number instead!!!
                resmarker.set_rsid(yamas::to_string((*bit).first));
                // start calculations
                resmarker.fixed_effect(markers);
                resmarker.random_effect(markers);
                #ifdef INCLUDEPRODUCTION
                resmarker.stouffers_method(markers, projects);
                #endif
                resmarker.check_directions(markers, project_size);
                check_reasonable_pvalue(resmarker);

                // write to log file
                // clear previous 'directions' content, such that
                // 'foul' markers cannot mess up the output
                directions.clear();
                directions = resmarker.get_directions();

                if ((resmarker.get_pvalue() <= ud.pvaluethreshold or resmarker.get_rpvalue() <= ud.pvaluethreshold)  && 
                    // meta-analysis was possible? (more then two valid markers?)
                    directions.size() - (int)count(directions.begin(), directions.end(), '?') 
                                      - (int)count(directions.begin(), directions.end(), 'x') >= 2) {
                    resmarker.set_all(markers, outfile, false);
                    outfile << resmarker;
                    // write the input markers, only if requested
                    if (ud.write_input) {
		        write_vector(markers, outfile, ""); // using ya_io template function
                    }
                    outfile << std::endl;
                    // apparently we hid a valid marker
                    total++;                      // total number of valid markers
                    // update statistics
                    if (ud.stats) sp->update(resmarker);
                }
                else {
                    errfile << std::setw(12) << markers[0].get_rsid();
                    write_vector(markers, errfile, "");
                    errfile << std::endl << std::flush;
                }
	        // erase content to avoid mangling between datasets  
                markers.erase(markers.begin(), markers.end());
	    }
            // markers are only valid for one chromosome - so erase them and save space
            markers_only.erase(markers_only.begin(), markers_only.end());
            // erase LD-Blocks per Chromosome
            ldblocks_per_chromosome.erase(ldblocks_per_chromosome.begin(), ldblocks_per_chromosome.end());
            // erase proxies per Chromosome
            proxies.erase(proxies.begin(), proxies.end());
            #ifdef MEMORY
            #pragma omp master
            {
            getrusage(RUSAGE_SELF, &ru);
            memlog << "after deleting proxies and LD-table for chromosome " << *c << " : " << ru.ru_maxrss << " kB." << std::endl;
            }        
            #endif
            //TODO: check whether no information is lost at this point
            // erase multilt per Chromosome
            multilt.erase(multilt.begin(), multilt.end());
	}
	outfile.close();
	errfile.close();

	// erase ld blocks
	ldblocks.erase(ldblocks.begin(), ldblocks.end());

        // The following output should be done regardless of the
        // requested verbosity:
        std::cout << "Done" << std::endl 
		  << "YAMAS analysed " << total << " valid markers." << std::endl;
        // write out statistics, if requested
        if (ud.stats) {
            std::cout << std::endl;
            sp->writeout(total);
        }
    }
    #endif // INCLUDEPRODUCTION

    /*!
       performs point-wise meta-analysis, falls back on proxy-alleles,
       if necessary

       input: vector of yamas::Project - each project
       input: file name of the proxyfile
    */
    void pointwise_proxy(const std::vector<yamas::Project> &projects, 
		         const std::string &proxyfname) {
        // if this compile flag is given, we open a special log file for profiling!
        #ifdef MEMORY
        std::ofstream memlog("yamas_memory.log");
        rusage ru;
        #endif
        // trim output for first marker set, only?
        if (ud.trim.size() > 0) {
            std::cerr << "Trimming is not possible with this algorithm." << std::endl
                      << "It would not make sense, either." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (ud.verbose) std::cout << "Setting pointwise meta-analysis with proxy alleles ..." << std::endl;
        // make non-const vector of projects (copy)
        std::vector<yamas::Project> iprojects;
        for (std::vector<yamas::Project>::const_iterator it = projects.begin(); it < projects.end(); it++) {
            iprojects.push_back(*it);
        }
        unsigned int project_size = projects.size(); // number of projects

        // create a vector of missing project IDs
        std::vector<unsigned int> missings;
        missings.reserve(project_size);

        // assign a particular chromosome, if requested
        if (ud.chromosome != 0) {
            chromosomes.clear();
            chromosomes.push_back(ud.chromosome);
        }

        // make up string of 'directions'
        std::string directions;
	directions.reserve(project_size);

        // a temporary marker
        yamas::MetaMarker marker;
        // temporary vector of markers
        std::vector<yamas::MetaMarker> proxymarkers, newmarkers;
        proxymarkers.reserve(project_size);
        newmarkers.reserve(project_size);
        // a flag
        bool done = false;

        // map of markers and a vector of those maps
        yamas::MarkerMap chrlt;
        yamas::MMarkersLT multilt;
        yamas::MMarkersPos mpos;
	yamas::ProxyLT proxies;
        std::vector<Proxy> proxy_vector;

        // output result marker
        yamas::ResultMarker resmarker; // will be used as a dummy
        // will hold the total count of analysed markers
        unsigned long int total = 0;
        // a temporary position marker
        unsigned long int position = 0;

        // create ostreams
        std::ofstream outfile((ud.logfname + ".log").c_str());
        std::ofstream errfile((ud.logfname + ".err").c_str());

        // write out generalized header
        yamas::write_header(iprojects, outfile, ud);
        yamas::write_errorlog_header(iprojects, errfile);

	// parse the proxyfile to get a look-up table: {chromosome: (begin, end)}
        yamas::ChrMap table;

        if (ud.verbose) {
            std::cout << "Screening reference file: " << proxyfname << std::endl;
	    std::cout << "(Screening is necessary to speed up the read-in, subsequently.)"
		      << std::endl
                      << "(The screening itself is always master-threaded.)" << std::endl
                      << std::endl;
	}
        #ifdef MEMORY
        #pragma omp master
        {
        getrusage(RUSAGE_SELF, &ru);
        memlog << "before reading proxy file: " << ru.ru_maxrss << " kB." << std::endl;
        }        
        #endif

	table = yamas::proxy_parse(proxyfname);
        
        #ifdef MEMORY
        #pragma omp master
        {
        getrusage(RUSAGE_SELF, &ru);
        memlog << "after reading proxy file: " << ru.ru_maxrss << " kB." << std::endl;
        }        
        #endif

        // get a statistics pointer
        yamas::Statistics* sp = yamas::Statistics::get_Statistics();

        // read in chromosome-wise
        for (std::vector<short int>::const_iterator c = chromosomes.begin();
             c != chromosomes.end(); c++) {
            if (ud.verbose and *c >= 0) std::cout << "Reading data for chromosome " << *c << std::endl;
            else if (ud.verbose and *c == -1) std::cout << "Reading all data at once" << std::endl;
	    // read in LD-data
            if (ud.verbose) std::cout << "Reading LD-data from " << proxyfname << std::endl;
	    proxies = readldfile(proxyfname, table, *c);
            if (ud.verbose) {
                std::cout << "Found " << proxies.size() << " proxies." << std::endl
                          << std::endl
                          << "Reading project data files:" << std::endl;
            }

            #pragma omp parallel for ordered private(chrlt)
            for (std::size_t i = 0; i < iprojects.size(); i++) {
                chrlt = iprojects[i].readChromosome(*c);
                #pragma omp ordered
                {
                if (ud.verbose) std::cout << chrlt.size() << " markers found in project " << iprojects[i].id << std::endl;
                // merge the chromosome wide lookup table with the entire one
                for (yamas::MarkerMap::iterator ci = chrlt.begin(); ci != chrlt.end(); ci++) {
                    multilt[(*ci).first].push_back((*ci).second);
                }
                }
                // get rid of data in chrlt
                chrlt.erase(chrlt.begin(), chrlt.end());
            }
            if (ud.verbose) std::cout << "A total of " 
                                      << multilt.size()
                                      << " unique markers were found." << std::endl;
            #ifdef MEMORY
            #pragma omp master
            {
            getrusage(RUSAGE_SELF, &ru);
            memlog << "after reading data for chromosome " << *c << ": " << ru.ru_maxrss << " kB." << std::endl;
            }            
            #endif

            // start calculating
            if (ud.verbose) std::cout << "Associating with proxies ... " << std::endl;

            // declare iterators to be prepared for multithreading
            yamas::MMarkersLT::const_iterator mit;
            std::vector<std::string> keys;
            keys.reserve(multilt.size());
            for (mit = multilt.begin(); mit != multilt.end(); ++mit) {
                keys.push_back((*mit).first);
            }
            std::vector<yamas::MetaMarker>::iterator tmpit;
            std::vector<yamas::MetaMarker>::const_iterator tmpit2, tmpit3, mmit, mmit2;
            std::vector<yamas::Proxy>::const_iterator proxy_it;

            std::vector<std::string>::const_iterator key;

            // Note: Cannot be performed in parallel, due to multiple iterator accesses on
            //       the same memory adress!
	    for (key = keys.begin(); key < keys.end(); ++key) {
                // Is this iterator still matching the map? -- if not abort loop.
                // This may happen, in threaded code due to master thread advancements - or is it?
                if (not multilt.count(*key) > 0) continue;
                // save the position for later use
                position = multilt[*key][0].get_position();
                // regardless of missings or not -- save the new markers to a position-wise lookup table
                for (tmpit = multilt[*key].begin(); tmpit != multilt[*key].end(); tmpit++) {    
                    mpos[position].push_back((*tmpit));
                }
                newmarkers.clear();
		// for a particular marker ID: are there less markers than projects?
		if (multilt[*key].size() < project_size) {
                    // save markers we already have
                    for (tmpit2 = multilt[*key].begin(); tmpit2 != multilt[*key].end(); ++tmpit2) {
                        newmarkers.push_back((*tmpit2));
                    }
	            // which project's marker is missing?
                    // first see which project markers are present:
                    yamas::reset(missings, project_size);
                    for (mmit = multilt[*key].begin(); mmit != multilt[*key].end(); mmit++) {
                        missings[(*mmit).get_project_id() - 1]++; // the -1 is because we need the position
                                                                  // and YAMAS counts from 1, vectors from 0
                    }

                    // now insert the proxy marker for each missing marker per project

                    // for this we need the appropriate vector of proxies, first
                    proxy_vector = proxies[*key];
                    // and we remember that we aren't done for that project and marker
                    done = false;
                    // and we iterate over all project IDs with missing markers
                    for (size_t mindex = 0; mindex < missings.size(); ++mindex) {
                        if (missings[mindex]) continue; // did not actually hit a missing marker?
                        // do we actually possess proxies for this missing marker?
                        if (not proxy_vector.size()) continue;
                        // the we iterate over the proxy vector to see, whether or not we have a substitute
                        for (proxy_it = proxy_vector.begin(); proxy_it != proxy_vector.end(); proxy_it++) {
                           // key present?
                           if (multilt.find(*key) != multilt.end()) {
                               // retrieve vector of markers with that key
                               proxymarkers = multilt[(*proxy_it).get_marker_id()];
                               // now we ask whether this marker is actually present for the required project
                               for (mmit2 = proxymarkers.begin(); mmit2 != proxymarkers.end(); mmit2++) {
                                  // skip marker, if that one is not valid within the project
                                  if (not (*mmit2).is_valid()) continue;
                                  marker = (*mmit2);
                                  // do we hit the right project?
                                  if (marker.get_project_id() == (mindex + 1)) {
                                      // now, we just set the allele and check later:
                                      marker.set_proxy_allele((*proxy_it).get_allele());
                                      marker.set_inphase_allele((*proxy_it).get_inphase_allele());
                                      // assign the proxy marker's switch status
                                      marker.set_switch_status((*proxy_it).get_switch());
                                      // remember the r2!
                                      marker.set_r2((*proxy_it).get_r2());
                                      // then insert that very marker into the vector of markers 
                                      // in place of the missing one
                                      newmarkers.push_back(marker);
				      done = true;
                                      break;
                                  }
                               }
                           }
                           if (done) break;
                        }
                        done = false;
                    }
		}
                else {
                    for (tmpit3 = multilt[*key].begin();
                         tmpit3 != multilt[*key].end(); ++tmpit3) {
                        newmarkers.push_back((*tmpit3));
                    }
                }
                // retrieve first marker's position
                // transfer rsid:markers-map to position:markers-map
                mpos[position] = newmarkers;
            }

            if (ud.verbose) std::cout << "Preparing meta-analysis" << std::endl;
	    
            // erase multilt to save space
            multilt.erase(multilt.begin(), multilt.end());
            // drop the position keys - we only need the sorted vector of vectors
            std::vector<std::vector<yamas::MetaMarker> > markers_only;
	        markers_only.reserve(mpos.size());

            for (yamas::MMarkersPos::const_iterator mit = mpos.begin();
	         mit != mpos.end(); mit++) {
                markers_only.push_back((*mit).second);
            }
            mpos.erase(mpos.begin(), mpos.end());

            if (ud.verbose) std::cout << "Starting meta-analysis" << std::endl;
            std::vector<std::vector<yamas::MetaMarker> >::iterator mitn;
            #pragma omp parallel for ordered private(resmarker,directions,mitn)
            for (mitn = markers_only.begin(); mitn < markers_only.end(); mitn++) {
                // first: save the rsID to the resmarker
                resmarker.set_rsid((*mitn)[0].get_rsid());
                // check all markers for consistency, correct, if necessary
                yamas::check_markers(*mitn); 
                // start calculations
                resmarker.fixed_effect(*mitn);
                resmarker.random_effect(*mitn);
                #ifdef INCLUDEPRODUCTION
                resmarker.stouffers_method(*mitn, projects);
                #endif
                resmarker.check_directions(*mitn, project_size);
                check_reasonable_pvalue(resmarker);

                // write to log file
                // clear previous 'directions' content, such that
                // 'foul' markers cannot mess up the output
                directions.clear();
                directions = resmarker.get_directions();

                #pragma omp ordered
                {
                if ((resmarker.get_pvalue() <= ud.pvaluethreshold or resmarker.get_rpvalue() <= ud.pvaluethreshold) &&
                    // meta-analysis was possible? (more then two valid markers?)
                    directions.size() - (int)count(directions.begin(), directions.end(), '?') 
                                      - (int)count(directions.begin(), directions.end(), 'x') >= 2) {
                    resmarker.set_all(*mitn, outfile);
                    outfile << resmarker;
                    // write the input markers, only if requested
                    if (ud.write_input) {
                        write_vector(*mitn, outfile, ""); // using ya_io template function
                    }
                    outfile << std::endl;
                    // apparently we hid a valid marker
                    total++;                      // total number of valid markers
                    // update statistics
                    if (ud.stats) sp->update(resmarker);
                }
                else {
                    errfile << std::setw(12) << (*mitn)[0].get_rsid();
                    write_vector(*mitn, errfile, "");
                    errfile << std::endl << std::flush;
                }
		(*mitn).erase((*mitn).begin(), (*mitn).end());
                }
            }

            // finally, clean up all resources
            markers_only.erase(markers_only.begin(), markers_only.end());
	    proxies.erase(proxies.begin(), proxies.end());

            // output that we are done for this chromosome
            if (ud.verbose) std::cout << "Done for this chromosome." << std::endl << std::endl;
            #ifdef MEMORY
            #pragma omp master
            {
            getrusage(RUSAGE_SELF, &ru);
            memlog << "after MA for chromosome " << *c << ": " << ru.ru_maxrss << " kB." << std::endl;
            }
            #endif
	}
	outfile.close();
	errfile.close();

        // The following output should be done regardless of the
        // requested verbosity:
        std::cout << "Done" << std::endl 
		  << "YAMAS analysed " << total << " valid markers." << std::endl;
        // write out statistics, if requested
        if (ud.stats) {
            std::cout << std::endl;
            sp->writeout(total);
        }
    } // end of function
} // end: namespace - yamas

