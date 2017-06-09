// $Rev: 354 $
// $LastChangedDate: 2010-06-22 16:55:31 +0200 (Di, 22 Jun 2010) $
// $Author: meesters $

#ifndef YA_ASSOC
#define YA_ASSOC

namespace yamas {
    struct UsageData; // forward declaration

    /*! will tabulate association data
        input should be the path and filename to
        a tped file
    */
    void association_test(const std::string &fname,
                          const std::string &missing);
} // end namespace yamas
    
#endif // YA_ASSOC
