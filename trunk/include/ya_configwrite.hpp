// $Rev: 609 $
// $LastChangedDate: 2010-12-16 06:57:41 +0100 (Do, 16 Dez 2010) $
// $Author: meesters $

#include <iostream>
#include <fstream>
#include <string>

#ifndef YA_CONFIGWRITER
#define YA_CONFIGWRITER

namespace yamas {
    /*! function to write a simple header file
        Attention: This function assumes that all input
                   files adhere to the same format.
    */
    void write_yamasconfig(const std::vector<std::string> &infnames,
                           const unsigned int cols[],
		           const std::string &outfname  = "yamas.cfg",
                           const std::string &ldfile    = "ldfile.txt.gz",
                           const std::string &proxyfile = "proxyreference.txt.gz");

    //! default writer for PLINK
    void write_PLINKconfig(const std::string &infnames);

    //! default writer for INTERSNP logistic
    void write_INTERSNPconfigLogistic(const std::string &infnames);

    //! default writer for INTERSNP linear
    void write_INTERSNPconfigLinear(const std::string &infnames);

    //! default writer for YAMAS
    void write_YAMASconfig(const std::string &infnames);

} // end namespace yamas
    
#endif // YA_CONFIGWRITER
