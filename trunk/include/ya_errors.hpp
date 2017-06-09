/********************************************************************/
/*                                                                  */
/* YAMAS - Yet Another Meta-Analysis Software                       */
/*                                                                  */
/* Error handling functions                                         */
/*                                                                  */
/********************************************************************/

// $Rev: 541 $
// $LastChangedDate: 2010-11-22 08:06:51 +0100 (Mo, 22 Nov 2010) $
// $Author: meesters $

#include <iostream>
#include <string>
#include <stdlib.h>

#include "./ya_global.hpp"

#ifndef YA_ERRORS
#define YA_ERRORS

namespace yamas {
    //! generic error function -> will issue a retrieved message and abort the program
    inline void error(const std::string &msg) {
	std::cerr << std::endl << "ERROR: " << msg << std::endl;
	std::exit(EXIT_FAILURE);
    }

    //! generic warning function -> will write retrieved warning message to std::cerr
    inline void warning(const std::string &msg) {
	std::cerr << std::endl << "WARNING: " << msg << std::endl;
    }

    //! Unknown file / path name -> will issue a message and abort the program
    inline void ioerror(const std::string &fname) {
	std::cerr << std::endl << "ERROR: Unable to open '" << fname << "'." << std::endl;
	std::cerr << "Please check this file's name and path."  << std::endl;
	std::exit(EXIT_FAILURE);
    }

    //! will be issued with existing log files
    inline void existinglogfile() {
         std::cerr << std::endl << "Log file with name '" << ud.logfname  
		   << ".log' already exists." << std::endl;
         std::cerr << "Please specify a different log file name and try again." << std::endl;
         exit(EXIT_FAILURE);
    }
} // end namespace yamas

#endif // YA_ERRORS
