/********************************************************************/
/*                                                                  */
/* YAMAS - Yet Another Meta-Analysis Software                       */
/*                                                                  */
/* Statistical overview function                                    */
/*                                                                  */
/********************************************************************/

// $Rev: 858 $
// $LastChangedDate: 2011-12-21 10:13:29 +0100 (Mi, 21 Dez 2011) $
// $Author: meesters $

#include <iostream>
#include <fstream>
#include <vector>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>

#ifndef YA_STATS
#define YA_STATS

namespace yamas {
    // forward declaration
    class ResultMarker;

    /*! a statistic overview class
        -- implemented as a singleton --
        This class will gather and report about the
        resulting P-values of the meta-analysis.

        'P-value shells' are: 5e-2, 1e-4, 1e-6, 1e-8, and smaller
    */
    class Statistics {
        public:
            /*! retrieve a pointer to a Statistics instance
               invoke like: Statistics* pointer = Statistics::get_Statistics();
            */
            static Statistics* get_Statistics();
            //! update the statistics
            void update(const yamas::ResultMarker &rmarker);
            //! write out the statistical summary 
            void writeout(const unsigned long int &total);
        protected:
            // summary for fixed effects
            std::vector<unsigned long int> fixed_stats_;
            // summary for random effects;
            std::vector<unsigned long int> random_stats_;
        private:
            static Statistics* inst_; // the singleton instance
            Statistics();
            Statistics(const Statistics&);
            Statistics& operator=(const Statistics&);
    };
} // end namespace yamas

#endif // YA_STATS
