/********************************************************************/
/*                                                                  */
/* YAMAS - Yet Another Meta-Analysis Software                       */
/*                                                                  */
/* Class holding project properties                                 */
/*                                                                  */
/********************************************************************/

// $Rev: 907 $
// $LastChangedDate: 2012-06-04 14:15:13 +0200 (Mo, 04 Jun 2012) $
// $Author: meesters $

#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <cmath>

#include <boost/program_options.hpp>

#include "ya_io.hpp"
#include "ya_lddata.hpp"

#ifndef YA_PROJECT
#define YA_PROJECT

namespace yamas {
    class MetaMarker; // foward declaration

    //! a map with key: rsid and value: MetaMarker object
    typedef std::map<std::string, MetaMarker> MarkerMap;
    //! a map with key: rsid and value: vector of MetaMarker objects sharing the key
    typedef std::map<std::string, std::vector<MetaMarker> > MMarkersLT;
    //! a map with key: position and value: vector of MetaMarker objects sharing the position
    typedef std::map<unsigned long int, std::vector<MetaMarker> > MMarkersPos;

    /*! \class Project
        yamas::Project will hold all project related (meta)information
    */
    class Project {
         friend std::ostream& operator<< (std::ostream& out, const Project& proj);
         protected:
             static unsigned int next_id;
         public:
             unsigned int id;
             Project();
             Project(const Project &p);
             #ifdef INCLUDEPRODUCTION
             Project(const std::map<std::string, short unsigned int> &cols,
                     const std::string &ifname,
                     const double &iweight);
             #else
             Project(const std::map<std::string, short unsigned int> &cols,
                     const std::string &ifname);
             #endif
             ~Project();
             Project& operator=(const Project &p) {
                  // no need to check for self-reference in this case
                  id           = p.id;
                  columns      = p.columns;
                  fname        = p.fname;
                  current_line = p.current_line;
                  #ifdef INCLUDEPRODUCTION
                  weight       = p.weight;
                  #endif
                  
                  return *this;
             };

             // member functions
             MarkerMap readChromosome(const short int &chromosome);
             void print_header() const;
             bool check() const;
             inline unsigned int get_id() const {return id;};
             inline void set_cl(const unsigned long int &cl) {current_line = cl;};
             inline unsigned long int get_cl() {return current_line;};
             #ifdef INCLUDEPRODUCTION
             inline double get_weight() const {return weight;};
             #endif

         private:
             /*! the column's keys are (in order):
                 marker, chr, position, effect_allele, other_allele
                 effect, wheight, P-value
             */
             std::map<std::string, short unsigned int> columns;
             //! file name of the input file
             std::string fname;
             //! keeps track of the current line in the associated input file
             unsigned long int current_line;
             #ifdef INCLUDEPRODUCTION
             //! the weight associated with this study / project
             double weight;
             #endif
    };    

    //! initializes a vector of Projects from a boost::program_options::variables_map
    std::vector<Project> extractor(boost::program_options::variables_map &vm);
} // end namespace yamas
    
#endif // YA_PROJECT
