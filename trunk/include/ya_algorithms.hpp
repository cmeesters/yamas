// $Rev: 873 $
// $LastChangedDate: 2012-01-26 00:43:46 +0100 (Do, 26 Jan 2012) $
// $Author: meesters $

#ifndef YA_ALGORITHMS
#define YA_ALGORITHMS

namespace yamas {
    class Project;   // forward declaration
    struct UsageData; // forward declaration

    void pointwise_comparison(const std::vector<Project> &projects);
    void per_ldblock(const std::vector<Project> &projects,
		     const std::string &ldfname,
                     const std::string &proxyfname);
    void pointwise_proxy(const std::vector<Project> &projects,
		         const std::string &proxyfname);
} // end namespace yamas
    
#endif // YA_ALGORITHMS
