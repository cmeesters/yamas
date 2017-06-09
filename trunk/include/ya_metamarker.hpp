/********************************************************************/
/*                                                                  */
/* YAMAS - Yet Another Meta-Analysis Software                       */
/*                                                                  */
/* Class holding meta marker properties                             */
/*                                                                  */
/********************************************************************/

// $Rev: 929 $
// $LastChangedDate: 2012-09-17 10:00:37 +0200 (Mo, 17 Sep 2012) $
// $Author: meesters $

#include <iomanip>
#include <cmath>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/tuple/tuple.hpp>

#include "./ya_interface.hpp"
#include "./ya_global.hpp"
#include "./ya_utils.hpp"
//#include "./ya_project.hpp"

#ifndef YA_METAMARKER
#define YA_METAMARKER

namespace yamas {
    // forward declaration
    class UsageData;
    class Proxy;
    class Project;

    typedef boost::tuple<double, double> ConfInterv;

    /*! \brief Class for holding GWAS output information 

        \class MetaMarker
        Will hold per Marker information.
        Each project will have severall MetaMarkers (limited by memory).
    */
    class MetaMarker {
        friend class ResultMarker;
        friend std::ostream& operator<< (std::ostream& out, const MetaMarker& marker);
        public:
            MetaMarker();
            MetaMarker(const MetaMarker &m);
            MetaMarker(const unsigned short int &id         ,//= 0,
                       const std::string &rsid              ,//= "",
                       const unsigned short int &chromosome ,//= 0,
                       const unsigned int       &position   ,//= 0,
                       const char &effect_allele            ,//= ' ',
                       const char &other_allele             ,//= ' ',
                       const double &effect                 ,//= 0.0,
                       const double &weight                 ,//= 0.0,
                       const double &pvalue                 );//= 0.0,

            virtual ~MetaMarker();
            virtual inline const MetaMarker& operator=(const MetaMarker &m) {
                rsid_             = m.rsid_;
                id                = m.id;
                chromosome_       = m.chromosome_;
                position_         = m.position_;
                effect_allele_    = m.effect_allele_;
                other_allele_     = m.other_allele_;
                effect_           = m.effect_;
                weight_           = m.weight_;
                pvalue_           = m.pvalue_;
                valid_            = m.valid_;
                r2_               = m.r2_;
                proxy_allele_     = m.proxy_allele_;
                inphase_allele_   = m.inphase_allele_;
                return *this;
            };

            //! comparison operator< for sorting according to p-values
            bool operator< (const MetaMarker &other) 
                 {return this->pvalue_ < other.pvalue_;};
            //! relational operator> for sorting according to p-values
            bool operator> (const MetaMarker &other) 
                 {return this->pvalue_ > other.pvalue_;};

            // define getters

            //! retrieve id of the project the marker belongs to
            inline unsigned short int get_project_id() const {return id;};
            //! retrieve marker id (rs number) of the marker
            inline std::string get_rsid()       const {return rsid_;};
            //! retrieve chromosome_ number (numeric, only) of a marker
            inline unsigned short int get_chromosome() const {return chromosome_;};
            //! retrieve position_ (in bp) of a marker
            inline unsigned int get_position()  const {return position_;};
            inline char get_effect_allele()     const {return effect_allele_;};
            inline char get_other_allele()      const {return other_allele_;};
            inline char get_proxy_allele()      const {return proxy_allele_;}; // can be NULL
            inline char get_inphase_allele()    const {return inphase_allele_;}; // can be NULL
            inline bool get_switch()            const {return switch_;};
            inline double get_effect()          const {return effect_;};
            inline double get_weight()          const {return weight_;};
            inline double get_pvalue()          const {return pvalue_;};
            inline double get_r2()              const {return r2_;};
            //! return a confidence interval
            ConfInterv get_confidence_interval() const;
            //! determine whether the valid_ flag is true
            inline bool   is_valid()            const {return valid_;};
            //! determine whether the marker is a proxy_marker
            inline bool   is_proxy_marker()     const {return proxy_allele_ != ' ';};
            // define setters
            inline void set_rsid(std::string id)      {rsid_             = id;};
            inline void set_chromosome(unsigned short int c) {chromosome_ = c;};
            inline void set_position(unsigned int p)  {position_         = p;};
            inline void set_effect_allele(char e)     {effect_allele_    = e;};
            inline void set_other_allele(char o)      {other_allele_     = o;};
            inline void set_proxy_allele(char p)      {proxy_allele_     = p;};
            inline void set_inphase_allele(char p)    {inphase_allele_   = p;};
            inline void set_switch_status(bool s)     {switch_           = s;};
            inline void set_effect(double e)          {effect_           = e;};
            inline void set_weight(double w)          {weight_           = w;
                                                       if (weight_ <= 0.0) valid_ = false;};
            inline void set_pvalue(double p)          {pvalue_           = p;
                                                       if (pvalue_ <= 0.0) valid_ = false;};
            inline void set_valid(bool v)             {valid_ = v;};
            inline void set_r2(double r2)             {r2_ = r2;};
            /*! 
               switch effect
               - *= -1 in case of betas
               - 
            */
            inline void switch_effect();
	    //! test whether effect allele is set
	    inline bool effect_allele_set()     const {return effect_allele_ != ' ';};

         private:
            unsigned short int id;

         protected:
            std::string rsid_;
	    unsigned short int chromosome_;
            unsigned int position_;
            char effect_allele_, other_allele_;
            double effect_, weight_, pvalue_;
            bool valid_; // indicating whether dealing with a 
                         // marker or not
         private:
            double r2_;
            char proxy_allele_; // in case a common marker gets the status
                                // of a proxy marker
            char inphase_allele_; // in case a common marker gets the status
                                // of a proxy marker
            bool switch_;        // proxy switched?
    };

    //! returns number of valid markers within a vector of MetaMarkers
    inline unsigned int count_valid_(const std::vector<MetaMarker> &markers) {
        unsigned int count = 0;
        for (std::vector<MetaMarker>::const_iterator mit = markers.begin();
             mit != markers.end(); mit++) {
            if (mit->is_valid()) count += 1;
        }
        return count;
    }

    /*! will check whether or not all rsids within a vector of markers
        are equal
    */
    inline bool equal_ids(const std::vector<MetaMarker> &markers) {
        // get first rsid
        std::string rsid = markers[0].get_rsid();
        for (size_t i = 1; i < markers.size(); i++) {
            if (markers[i].get_rsid() != rsid) return false;
        }
        return true;
    }

    /*! ResultMarker member variables should be defined on the fly:
        While initialization sets all to 0, "", or ' ', all values
        need to be calculated and cannot be read in. Ensure that a
        ResultMarker is indeed initialized before(!) usage, which means:
        Don't use pointers, only.
    */
    class ResultMarker : public MetaMarker {
        /*! not inhereted from MetaMarker, because the private MetaMarker project id
            should not be written out.
        */
        friend std::ostream& operator<< (std::ostream& out, const ResultMarker& marker);
        
        private:
            //! random effect_ - effect_
            double reffect_;
            //! random effect_ - weight_
            double rweight_;
            //! random effect_ - p-value
            double rpvalue_;
            // confidence interval
            ConfInterv fci_, rci_;
            //! Q - represents the total variance
            double q_;
            //! I2-statistics
            double i2_;
            #ifdef INCLUDEPRODUCTION
            //! Stouffer's Z-Value
            double sz_;
            //! P-values from Stouffer's method
            double spvalue_;
            #endif
            //! indices for direction and valid_ity of input markers
            std::string directions_;
        public:
            //TODO: add stouffers method's member variable to all relevant functions
            //! default constructor
            ResultMarker();
            //! copy constructor
            ResultMarker(const ResultMarker &r);
            //! explicit assignment constructor
            ResultMarker(// project id is irrelevant
                       const std::string &mrsid,
                       const unsigned short int &chr,
                       const unsigned int &pos,
                       const char &e_allele,
                       const char &o_allele,
                       const double &e,
                       const double &w,
                       const ConfInterv &fci,
                       const double &p,
                       const double &re,
                       const double &rw,
                       const ConfInterv &rci,
                       const double &rp,
                       const double &qvalue,
                       const double &i2value,
                       #ifdef INCLUDEPRODUCTION
                       const double &sz,
                       const double &spvalue,
                       #endif
                       const std::string &directions);
                       
            //! assignment operator
            const ResultMarker& operator=(const ResultMarker &r) {
                chromosome_       = r.chromosome_;
                position_         = r.position_;
                effect_allele_    = r.effect_allele_;
                other_allele_     = r.other_allele_;
                effect_           = r.effect_;
                weight_           = r.weight_;
                fci_              = r.fci_;
                pvalue_           = r.pvalue_;
                reffect_          = r.reffect_;
                rweight_          = r.rweight_;
                rci_              = r.rci_;
                rpvalue_          = r.pvalue_;
                q_                = r.q_;
                i2_               = r.i2_;
                #ifdef INCLUDEPRODUCTION
                sz_               = r.sz_;
                spvalue_          = r.spvalue_;
                #endif
                directions_       = r.directions_;
                return *this;
            }

            inline double get_reffect()            const {return reffect_;};
            inline double get_rweight()            const {return rweight_;};
            inline double get_rpvalue()            const {return rpvalue_;};
            inline double get_qvalue()             const {return q_;};
            //! returns confidence interval of the effect
            inline ConfInterv get_fci()            const {return fci_;}
            //! returns confidence interval of the OR
            inline ConfInterv get_fci_or()         const {return ConfInterv(std::exp(fci_.get<0>()), std::exp(fci_.get<1>()));}
            //! returns confidence interval of the effect
            inline ConfInterv get_rci()            const {return rci_;}
            //! returns confidence interval of the OR
            inline ConfInterv get_rci_or()         const {return ConfInterv(std::exp(rci_.get<0>()), std::exp(rci_.get<1>()));}
            inline double get_i2()                 const {return i2_;};
            #ifdef INCLUDEPRODUCTION
            inline double get_stouffers_z()        const {return sz_;};
            inline double get_stouffers_pvalue()   const {return spvalue_;};
            #endif
            inline std::string get_directions()    const {return directions_;};

            inline void check_directions(const std::vector<MetaMarker> &markers,
                                         const unsigned int project_size) {
                directions_.clear();
                directions_.reserve(project_size);
                bool flag = false;
                for (size_t i = 1; i <= project_size; i++) {
                    for (std::vector<MetaMarker>::const_iterator mit = markers.begin();
                         mit != markers.end(); mit++) {
                        if (mit->id == i) {
                            if (mit->valid_ and mit->effect_ > 0)
                                directions_.append("+");
                            else if (mit->valid_)
                                directions_.append("-");
                            else
                                directions_.append("x");
                            flag = true;
                            break;
                        }
                    }
                    // missing marker?
                    if (! flag) directions_.append("?");
                    flag = false;
                }
            }

            /*!
              will set several member attributes and write whether or not errors are
              detected to the input outstream
              
              some 'mismatch checking' is performed

              -> members 'rsid_', 'chromosome_', 'position_', 'effect_allele_', and 'reference_allele_'
                 are set
            */
            void set_all(const std::vector<MetaMarker> &markers, 
			 std::ostream &outfile,
			 bool warnings = true);

            /*!
              will calculate and set effect_, weight_, p-value for a
              ResultMarker considering fixed effect_s

              \f$E = \frac{\sum_{i=1}^{k} w_i \cdot E_i}{\sum_{i=1}^{k} w_i}\f$ with \f$w = \frac{1}{\mathrm{standard~error}^2}~.\f$
              Finally a \f$z\f$-value is computed, \f$z = \frac{\bar{E}}{SE_{\bar{E}}}~,\f$ which is used to calculate a p-value with a two-tailed test:
              \f$p = 2 \cdot \left(1 - \Theta(| z |) \right)~,\f$where \f$\Theta(| z |)\f$ is the standard normal cumulative density distribution function.
            */
            void fixed_effect(const std::vector<MetaMarker> &markers);

            /*!
              will calculate and set effect_, weight_, p-value for a
              Result Marker considering fixed effect_s
            */
            void random_effect(const std::vector<MetaMarker> &markers);

            /*! MA according to Stouffer's method
            */
            #ifdef INCLUDEPRODUCTION
            void stouffers_method(const std::vector<MetaMarker> &markers,
                                  const std::vector<Project>    &projects);
            #endif
    };

    /*! trims the input vector for valid ones for the first project
        Attention: The function's behaviour depends on ud.trim.
    */
    void trim_markers(std::vector<std::vector<MetaMarker> >&markers); 

    //! will check for the occurence of unreasonable P-Values and issue a warning.
    void check_reasonable_pvalue(const ResultMarker &marker);
   
    //! will check for genotype equality and correct effect_ direction if necessary,
    void check_markers(std::vector<MetaMarker> &markers);
} // end namespace yamas

#endif // YA_METAMARKER
