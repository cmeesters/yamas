// $Rev: 906 $
// $LastChangedDate: 2012-06-04 10:04:45 +0200 (Mo, 04 Jun 2012) $
// $Author: meesters $

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>

#include "ya_interface.hpp"

#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>

#ifndef YA_IO
#define YA_IO

/*! This header file only defines basic I/O routines.
    More sophisticated ones, e.g. for parsing LD-information,
    please refer to the designated headerfiles.
*/

namespace yamas {
    // forward class declaration
    class Project;
    //! writes a generalized header for output files
    void write_header(const std::vector<yamas::Project> &projects, std::ofstream &out,
                      const yamas::UsageData &ud,
                      bool write_status = true);

    //! writes generalized header for error log files
    void write_errorlog_header(const std::vector<Project> &projects, 
                               std::ofstream &out); 
} // end namespace yamas
    
#endif // YA_IO


