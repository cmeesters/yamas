/********************************************************************/
/*                                                                  */
/* YAMAS - Yet Another Meta-Analysis Software                       */
/*                                                                  */
/* Classes holding Markers used in meta-analysis.                   */
/* Also: The actual analysis is carried out within member           */
/*       functions of the result marker class.                      */
/*                                                                  */
/********************************************************************/

// $Rev: 940 $
// $LastChangedDate: 2014-08-30 12:57:41 +0200 (Sa, 30 Aug 2014) $
// $Author: meesters $

#include <iostream>
#include <fstream>
#include <iomanip> // manipulator for fixed column width formatting
#include <vector>
#include <string>
#include <cmath>
#include <math.h> // for isnan()
#include <limits>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
using boost::math::normal;
using boost::math::chi_squared;

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/lexical_cast.hpp>
using namespace boost;

// include own software parts
#include "../include/ya_errors.hpp"
#include "../include/ya_lddata.hpp"
#include "../include/ya_metamarker.hpp"
#include "../include/ya_utils.hpp"
#include "../include/ya_interface.hpp"
#include "../include/ya_global.hpp"

#ifdef INCLUDEPRODUCTION
#include "../include/ya_project.hpp"
#endif

EXTERN yamas::UsageData ud;

yamas::MetaMarker::MetaMarker() : id(0),
                                  rsid_(""),
                                  chromosome_(0),
                                  position_(0),
                                  effect_allele_(' '),
                                  other_allele_(' '),
                                  effect_(0),
                                  weight_(0),
                                  pvalue_(1),
                                  r2_(0),
                                  proxy_allele_(' '),
                                  inphase_allele_(' '),
                                  switch_(false) {}

// copy-constructor
yamas::MetaMarker::MetaMarker(const yamas::MetaMarker &m) :
                              id(m.id),
                              rsid_(m.rsid_),
                              chromosome_(m.chromosome_),
                              position_(m.position_),
                              effect_allele_(m.effect_allele_),
                              other_allele_(m.other_allele_),
                              effect_(m.effect_),
                              weight_(m.weight_),
                              pvalue_(m.pvalue_),
                              valid_(m.valid_),
                              r2_(m.r2_),
                              proxy_allele_(m.proxy_allele_),
                              inphase_allele_(m.inphase_allele_),
                              switch_(m.switch_) {}

yamas::MetaMarker::MetaMarker(const unsigned short int &pid,
                              const std::string &rsid,
                              const unsigned short int &chromosome,
                              const unsigned int       &position,
                              const char &effect_allele,
                              const char &other_allele,
                              const double &effect,
                              const double &weight,
                              const double &pvalue) {
    id                = pid;
    rsid_             = rsid;
    chromosome_       = chromosome;
    position_         = position;
    effect_allele_    = effect_allele;
    other_allele_     = other_allele;
    effect_           = effect;
    weight_           = weight;
    pvalue_           = pvalue;
    r2_               = 0.0; // r2 should not be set upon construction
    proxy_allele_     = ' '; // proxy_allele_ cannot be set upon construction!
    inphase_allele_   = ' '; // inphase_allele_ cannot be set upon construction!
    switch_           = false;
    if (pvalue_ < 0 || pvalue >= 1 ||  weight_ <= 0.0 ||
        position_ <= 0 ) valid_ = false;
    else valid_ = true;
}

yamas::MetaMarker::~MetaMarker() {}

yamas::ResultMarker::ResultMarker() : reffect_(0),
                                      rweight_(0),
                                      rpvalue_(1.0),
                                      fci_(ConfInterv(0, 0)),
                                      rci_(ConfInterv(0, 0)),
                                      q_(0),
                                      i2_(0),
                                      directions_("") {}

yamas::ResultMarker::ResultMarker(const yamas::ResultMarker &r) :
                                  reffect_(r.reffect_),
                                  rweight_(r.rweight_),
                                  rpvalue_(r.pvalue_),
                                  fci_(r.fci_),
                                  rci_(r.rci_),
                                  q_(r.q_),
                                  i2_(r.i2_),
                                  #ifdef INCLUDEPRODUCTION
                                  sz_(r.sz_),
                                  spvalue_(r.spvalue_),
                                  #endif
                                  directions_(r.directions_) {
                                  this->id             = r.id;
                                  this->rsid_          = r.rsid_;
                                  this->chromosome_    = r.chromosome_;
                                  this->position_      = r.position_;
                                  this->effect_allele_ = r.effect_allele_;
                                  this->other_allele_  = r.other_allele_;
                                  this->effect_        = r.effect_;
                                  this->weight_        = r.weight_;
                                  this->pvalue_        = r.pvalue_;
                                  this->valid_         = true; // can a result marker be false?
                                  //this->r2_             = r.r2_;             // not needed
                                  //this->proxy_allele_   = r.proxy_allele_;   // not needed
                                  //this->inphase_allele_ = r.inphase_allele_; // not needed
                                  }

yamas::ResultMarker::ResultMarker(const std::string &mrsid,
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
                                  const std::string &d) {
    rsid_             = mrsid;
    chromosome_       = chr;
    position_         = pos;
    effect_allele_    = e_allele;
    other_allele_     = o_allele;
    effect_           = e;
    weight_           = w;
    fci_              = fci;
    pvalue_           = p;
    reffect_          = re;
    rweight_          = rw;
    rci_              = rci;
    rpvalue_          = rp;
    q_                = qvalue;
    i2_               = i2value;
    #ifdef INCLUDEPRODUCTION
    sz_               = sz;
    spvalue_          = spvalue;
    #endif
    directions_       = d;
}

namespace yamas {
    std::ostream& operator<< (std::ostream& out,
                              const MetaMarker& marker) {
    // temporary variable for doubles
        double tmp;
        // temporary variable for alleles
        char allele;
        char inphase;
        // save flags
        std::ios_base::fmtflags flags = out.flags();

        out << std::setw(11) << marker.get_project_id()
            << std::setw(12) << marker.get_rsid()
            << std::setw(4)  << marker.get_chromosome()
            << std::setw(12) << marker.get_position()
            << std::setw(4)  << marker.get_effect_allele()
            << std::setw(4)  << marker.get_other_allele();
        if (ud.odds_ratio)
            out << std::right << std::setw(12) << std::exp(marker.get_effect());
        else
            out << std::setw(12) << marker.get_effect();
        out << std::fixed << std::setw(11) << marker.get_weight();
        tmp = marker.get_pvalue();
        if (tmp < 0.0001) {
            out << std::right << std::scientific // set to scientific notation
                << std::setw(14)  << tmp
                << std::fixed;      // set to 'normal' notation
        }
        else {
            out << std::right << std::setw(14) << std::setprecision(4) << tmp;
        }
        if (ud.algorithm == "fillwithproxies") {
            allele  = marker.get_proxy_allele();
            inphase = marker.get_inphase_allele();
            allele == ' ' ?
            out << std::setw(6)  << std::setprecision(3) << marker.get_r2()
                << std::setw(12) << "No~Proxy"
                //<< std::setw(13) << "No~Proxy" // Former STRANDSWITCH column
                          :
            out << std::setw(6)  << std::setprecision(3) << marker.get_r2()
                << std::setw(10) << allele << "-" << inphase;
                //<< std::setw(13) << std::boolalpha << (bool)marker.get_switch(); // Former STRANDSWITCH column
        }
        // restore flags
        out.flags(flags);
        return out;
    }

    std::ostream& operator<< (std::ostream& out,
                              const ResultMarker& marker) {
        double tmp;
        unsigned int nprojects = marker.get_directions().size();
        // save flags
        std::ios_base::fmtflags flags = out.flags();
        // consider first column only for algos where status warnings
        // are issued
        if (ud.algorithm != "pointwise") out << std::right;
        out << std::setw(12) << marker.get_rsid()
            << std::setw(4)  << marker.get_chromosome()
            << std::setw(12) << marker.get_position()
            << std::setw(4)  << marker.get_effect_allele()
            << std::setw(4)  << marker.get_other_allele()
            << (nprojects <= 20 ? std::setw(20) : std::setw(nprojects + 1))
                             << marker.get_directions()
            << std::setw(12) << std::setprecision(5);
            if (ud.odds_ratio)
                out << std::exp(marker.get_effect());
            else
                out << marker.get_effect();
        out << std::setw(10) << std::setprecision(5) << marker.get_weight() << std::fixed;
        // confidence interval output for fixed effect
        out << boost::tuples::set_open('[') << boost::tuples::set_close(']')
            << boost::tuples::set_delimiter(',') << std::setw(15)
            << std::setprecision(2);
        if (ud.odds_ratio) out << marker.get_fci_or();
        else               out << marker.get_fci();

        out << std::fixed << std::setw(14) << std::right;
        tmp = marker.get_pvalue();
        if (tmp < 0.001) {
            out << std::right << std::scientific // set to scientific notation
                << tmp
                << std::fixed;      // set to 'normal' notation
        }
        else {
            out << std::setprecision(4) << tmp;
        }

        // output for heterogeneous effect indicators
        out << std::setw(14) << std::right;
        tmp = marker.get_qvalue();
        if (tmp < 0.001) {
            out << std::scientific // set to scientific notation
                << tmp
                << std::fixed;    // set to 'normal' notation
        }
        else {
            out << std::setprecision(4) << tmp;
        }
        out << std::setw(7) << std::setprecision(1) << marker.get_i2();
        // restore flags
        out.flags(flags);

        out << std::setw(12) << std::setprecision(5);

        if (ud.odds_ratio)
            out << std::exp(marker.get_reffect());
        else
            out << marker.get_reffect();
        out << std::setw(11) << marker.get_rweight();

        // confidence interval output for fixed effect
        out << boost::tuples::set_open('[') << boost::tuples::set_close(']')
            << boost::tuples::set_delimiter(',') << std::setw(17)
            << std::setprecision(2) << std::fixed;
        if (ud.odds_ratio) out << marker.get_rci_or();
        else               out << marker.get_rci();

        tmp = marker.get_rpvalue();
        if (tmp < 0.001) {
            out << std::right    << std::scientific // set to scientific notation
                << std::setw(14) << tmp;
        }
        else {
            out << std::right << std::setprecision(4) << std::setw(14)  << tmp;
        }

        #ifdef INCLUDEPRODUCTION
        // conditional output of values according to Stouffer's method
        if (ud.use_study_weights) {
            out << std::setw(12)  << marker.get_stouffers_z();
            tmp = marker.get_stouffers_pvalue();
            if (tmp < 0.001) {
                out << std::right << std::scientific << std::setw(12) << tmp;
            }
            else {
                out << std::right << std::setprecision(4) << std::setw(12) << tmp;
            }
        }
        #endif

        // restore flags
        out.flags(flags);
        return out;
    }

    /* definition see hpp */
    void yamas::MetaMarker::switch_effect() {
         // as the effect already is transformed to 'beta':
         effect_ *= -1;
         // This swapping is only performed prior to marker output,
         // hence it does not harm the evaluation of markers.

         // std::swap involves generic casting and is a potential
         // performance bottleneck!
         char temp;
         temp = effect_allele_;
         effect_allele_ = other_allele_;
         other_allele_  = temp;
    }

    yamas::ConfInterv yamas::MetaMarker::get_confidence_interval() const {
         yamas::ConfInterv ci;

         normal nd; // normal distribution
         double limit = quantile(complement(nd, (1 - ud.cl) / 2)) * weight_;
         if (!isnan(limit) && ud.odds_ratio) {
             // lower limit of the confidence interval
             ci.get<0>() = exp(effect_ - limit);
             // upper limit of the confidence interval
             ci.get<1>() = exp(effect_ + limit);
         }
         else if (!isnan(limit)) {
             // lower limit of the confidence interval
             ci.get<0>() = effect_ - limit;
             // upper limit of the confidence interval
             ci.get<1>() = effect_ + limit;
         }
         else {
             ci.get<0>() = std::numeric_limits<double>::quiet_NaN();
             ci.get<1>() = std::numeric_limits<double>::quiet_NaN();
         }
         return ci;
    }

    /*!
       will set several member attributes and write whether or not errors are
       detected to the input outstream

       some 'mismatch checking' is performed

       -> members 'rsid_', 'chromosome_', 'position_', 'effect_allele_', and 'reference_allele_'
          are set
    */
    void yamas::ResultMarker::set_all(const std::vector<MetaMarker> &markers,
                                      std::ostream &outfile,
                                      bool warnings) {
        // save flags
        std::ios_base::fmtflags flags = outfile.flags();
        std::vector<MetaMarker>::const_iterator mit;
        chromosome_       = markers[0].chromosome_;
        position_         = markers[0].position_;
        effect_allele_    = markers[0].effect_allele_;
        other_allele_     = markers[0].other_allele_;
        // own flag - did we already meet an error?
        bool eflag = false;
        // counter for valid_ markers
        unsigned int count = 0;
        // effect and reference allele of the first marker
        char eff = ' ', ref = ' ', eff2 = ' ', ref2 = ' ';
        // flag to know whether the previous alles have been set
        bool alflag = false;
        for (mit = markers.begin(); mit != markers.end(); mit++) {
            if (!alflag && mit->valid_) {
               eff = mit->get_effect_allele();
               ref = mit->get_other_allele();
               alflag = true;
            }
            if (mit->valid_) count++;
        }
        if (count >= 2 and warnings) {
            outfile << std::left << std::setw(23);
            for (mit = markers.begin() + 1; mit != markers.end(); mit++) {
                eff2 = mit->get_effect_allele();
                ref2 = mit->get_other_allele();
                // a particular marker has an allele not matching the
                // 1st marker at all?
                if (!((eff2 == eff) || (eff2 == ref)) && 
                    !((ref2 == ref) || (ref2 == eff))) {
                    outfile << "ERROR:~Not~all~alleles~matching.";
                    eflag = true;
                    break;
                }
                else if (chromosome_ != mit->chromosome_) {
                    outfile << "WARNING:~chr.~mismatch";
                    eflag = true;
                    break;
                }
            }
            if (!eflag) {
               if (this->pvalue_ < 1e-50) {
                  outfile << "WARNING:~Insane~P~Value";
                  eflag = true;
               }
            }

            if (!eflag) {
                for (mit = markers.begin() + 1; mit != markers.end(); mit++) {
                    /*
                      This constant is a bit arbitrary. It assumes an effective
                      sample size of ~1000. If the actual number deviates only a
                      more sophisticated alternative would help (simulation?).
                    */
                    if (ud.odds_ratio
                            and mit->valid_ // only useful for included markers
                            and std::fabs(mit->effect_ - std::log(0.001)) < std::numeric_limits<double>::min()
                            and ud.odds_ratio) {
                        outfile << "WARNING:~1~OR~was~zero?";
                        eflag = true;
                        break;
                    }
                }
            }
            if (!eflag) {
                for (mit = markers.begin(); mit != markers.end(); mit++) {
                    if ((ud.algorithm == "fillwithproxies" or
                         ud.algorithm == "ldblockwise") and mit->r2_ == 1.0) {
                        outfile << "INFO:~proxy~with~r^2=1";
                        eflag = true;
                        break;
                     }
                }
            }
            if (!eflag) {
                for (mit = markers.begin() + 1; mit != markers.end(); mit++) {
                    // write out only one warning, even if more could be issued
                    if (position_ != mit->position_ and ud.algorithm != "fillwithproxies"
                                                    and ud.algorithm != "ldblockwise") {
                        outfile << "WARNING:~pos.~mismatch";
                        eflag = true;
                        break;
                    }
                    else if (ud.algorithm == "fillwithproxies" and not equal_ids(markers)) {
                        outfile << "INFO:~entered~proxy";
                        eflag = true;
                        break;
                    }
                }
            }
            if (not eflag) outfile << "OK";

            // restore flags
            outfile.flags(flags);
        }
    }

    void yamas::ResultMarker::fixed_effect(const std::vector<yamas::MetaMarker> &markers) {
        std::vector<yamas::MetaMarker>::const_iterator mit;
        double w       = 0; // sum of inverse weights as variances
        double wse     = 0; // sum of weighted effects
        unsigned int l = 0; // number of valid markers
        double limit   = 0; // boundary of ci
        normal nd;          // normal distribution
        for (mit = markers.begin(); mit != markers.end(); mit++) {
            // check for validity of each marker
            if (mit->valid_) {
                wse += mit->effect_ / (mit->weight_ * mit->weight_);
                w   +=            1 / (mit->weight_ * mit->weight_);
                l++;
            }
        }
        effect_  = wse / w;
        weight_  = std::sqrt(1 / w);

        limit = quantile(complement(nd, (1 - ud.cl) / 2)) * weight_; // fix Andre / std::sqrt(l-1);
        if (!isnan(limit)) {
            // lower limit of the confidence interval
            fci_.get<0>() = effect_ - limit;
            // upper limit of the confidence interval
            fci_.get<1>() = effect_ + limit;
        }
        else {
            fci_.get<0>() = std::numeric_limits<double>::quiet_NaN();
            fci_.get<1>() = std::numeric_limits<double>::quiet_NaN();
        }
        // if zero division takes place, 'inf' will be set
        // which is enough for an error indication
        if (!isnan(std::fabs(effect_)))
            pvalue_ = cdf(complement(nd, std::fabs(effect_ / weight_))) * 2;
        else
            pvalue_ = std::numeric_limits<double>::quiet_NaN();
    }

    void yamas::ResultMarker::random_effect(const std::vector<MetaMarker> &markers) {
        std::vector<yamas::MetaMarker>::const_iterator mit;
        double c        = 0; // scaling factor
        double w        = 0; // weights
        double sw       = 0; // sum of weights
        double swe      = 0; // sum of weight * effect
        double swe2     = 0; // sum of weight * effect^2
        double sw2      = 0; // sum of the squared weights
        unsigned int df = 0; // degrees of freedom within the study -
                             // cannot be a member variable, because
                             // the degrees of freedom depend on
                             // the number of valid markers for each test
        unsigned int k  = 0; // number of valid markers
        double wr       = 0; // sum of 1 / (within study variance + r2)
        double we       = 0; // sum of wr * effect
        double tausqr   = 0; // correction factor for heterogenous effects / weights
        double limit    = 0; // boundary of ci
        double Q        = 0; // Q statistic
        normal nd;           // standard normal distribution

        // start with testing if a random effect should be accounted for
        // calculate Cochrane's Q
        for (mit = markers.begin(); mit != markers.end(); mit++) {
            if (mit->valid_) {
                w    =  1 / std::pow(mit->weight_, 2);
                sw   += w;
                sw2  += std::pow(w, 2);
                swe  += w * mit->effect_;
                swe2 += w * std::pow(mit->effect_, 2);
                k++;
            }
        }

        Q  = swe2 - (std::pow(swe, 2) / sw);
        df = k - 1; // degrees of freedom = N - 1
        // if df == 0 -> c is also == 0 and all values set to nan or inf
        c = sw - (sw2 / sw); // calculate scaling factor
        // calculate between study variance

        if (Q > df) {
            tausqr = (Q - df) / c;
            i2_    = 100.0 * ((Q - df) / Q);
        }
        else {
            i2_    = 0;
        }

        // calculate random effect and random se
        for (mit = markers.begin(); mit != markers.end(); mit++) {
            if (mit->valid_) {
                we += (1 / (std::pow(mit->weight_, 2) + tausqr)) * mit->effect_;
                wr +=  1 / (std::pow(mit->weight_, 2) + tausqr);
            }
        }
        reffect_  = we / wr;
        rweight_  = 1 / std::sqrt(wr);

        limit = quantile(complement(nd, (1 - ud.cl)/2)) * rweight_; // fix Andre / std::sqrt(k-1);
        if (!isnan(limit)) {
            // lower limit of the confidence interval
            rci_.get<0>() = reffect_ - limit;
            // upper limit of the confidence interval
            rci_.get<1>() = reffect_ + limit;
        }
        else {
            rci_.get<0>() = std::numeric_limits<double>::quiet_NaN();
            rci_.get<1>() = std::numeric_limits<double>::quiet_NaN();
        }

        // re-use wr
        wr = std::sqrt(wr);
        if (!isnan(reffect_) && wr > 0)
            rpvalue_ = cdf(complement(nd, std::fabs(reffect_ / rweight_))) * 2;
        else rpvalue_ = std::numeric_limits<double>::quiet_NaN();

        // in order to calculate the P-value for Q we need a chi2 statistic
        if (Q > 0 and df > 0 and not yamas::isinf(Q)) {
            chi_squared chi2(df);
            q_ = cdf(complement(chi2, Q));
        }
        else {
            q_ = 1;
        }
    }

    #ifdef INCLUDEPRODUCTION
    void yamas::ResultMarker::stouffers_method(const std::vector<MetaMarker> &markers,
                                               const std::vector<yamas::Project> &projects) {
        normal nd;               // standard normal distribution
        double zi        = 0;    // the ith z-value
        unsigned int pid = 0;    // project id
        double weights   = 0;    // sum of all project weights
        double weight    = 0;    // individual project weight
        unsigned int n   = 0;    // number of valid markers
        for (std::vector<MetaMarker>::const_iterator mit = markers.begin();
             mit != markers.end(); mit++) {
           if (mit->is_valid()) {
             pid      = mit->get_project_id() - 1;
             weight   = projects[pid].get_weight();
             zi      += cdf(complement(nd, 1 - mit->get_pvalue())) * weight;
             weights += weight;
             n++;
           }
        }
        if (n >= 2) {
           sz_      = zi / sqrt(weights);
           chi_squared chi2(n); // df = n
           spvalue_ = cdf(complement(chi2, sz_));
        }
        else {
           sz_ = 0;
           spvalue_ = std::numeric_limits<double>::quiet_NaN();
        }
    }
    #endif

    void check_reasonable_pvalue(const ResultMarker &marker) {
       static unsigned int times_called = 0;
       static bool issued               = false;
       if ((!issued) && (times_called < 100) && !ud.odds_ratio && (marker.get_pvalue() < 1e-50)) {
          yamas::warning("P-value for marker '" + marker.get_rsid() + "' (and perhaps for other markers) is unreasonably low. " +
                         "Did you forget to put the -o flag?");
          issued = true;
       }
       times_called++;
    }

    void check_markers(std::vector<MetaMarker> &markers) {
        // transfer effect into beta estimated, if necessary
        double beta = 0.0, effect = 0.0;
        // temporary index
        size_t index = 0;
        // characters for the effect and other allele
        char eff   = ' ', ref  = ' ';
        // same for different marker
        char eff2  = ' ', ref2 = ' ';
        // proxy allele, if needed
        char prox  = ' ';
        // proxy inphase-allele, if needed
        char iprox  = ' ';
        // normal distribution to check for correct P-values
        normal nd;
        // convert to beta estimates, if demanded
        if (ud.odds_ratio) {
           for (size_t i = 0; i < markers.size(); i++) {
               effect = markers[i].get_effect();
               if (effect < 0) yamas::error("Odds ratio for marker '" + markers[i].get_rsid() + "' in project" +
                                            boost::lexical_cast<std::string>(markers[i].get_project_id()) + " was < 0");
               // catch OR's which are zero:
               if (std::fabs(effect) < std::numeric_limits<double>::min() and markers[i].is_valid()) {
                   effect += 0.001; // do not change to std::numeric_limits<double>::min()
                                    // as this might be different on different systems
               }
               // beta = log(OR)
               beta = std::log(effect);
               if (!isnan(beta)) // sucess?
                   markers[i].set_effect(beta);
               else
               markers[i].set_valid(false);
           }
        }

        // does the user want to be all effects interpreted in one direction?
        if (ud.equal_effects) {
            for (std::size_t i = 0; i < markers.size(); i++) {
                beta = markers[i].get_effect();
                if (beta < 0)
                    markers[i].set_effect(beta * -1.0);
            }
        }
        // if not check the alleles
        else {
            for (index = 0; index < markers.size(); index++) {
                // we just take the first (which cannot be a proxy marker!) as reference
                if (markers[index].is_valid()) {
                    eff = markers[index].get_effect_allele();
                    ref = markers[index].get_other_allele();
                    break;
                }
            }
            index++;

            if (ud.algorithm != "pointwise") {
                for (std::size_t i = index; i < markers.size(); i++) {
                    eff2 = markers[i].get_effect_allele();
                    ref2 = markers[i].get_other_allele();
                    if (markers[i].is_valid()) {
                        if (markers[i].is_proxy_marker()) {
                            prox = markers[i].get_proxy_allele();
                            iprox = markers[i].get_inphase_allele();

                            // This if/else-block assumes:
                            // substituted marker shows the same allele order as
                            // the study.

                            if (prox == ref || (complement(prox) == ref && eff != complement(ref)))
                            // && in the above clause is because of A/T polymorphisms
                                markers[i].switch_effect();
                            // TODO: check for A/T polymorphisms
                            if (iprox == eff2) // because of A/T poly
                                {
                                }
                            else if (iprox == ref2 || iprox == complement(ref2))
                                markers[i].switch_effect();
                            // NOTE: In case a study will show different alleles
                            // than the first (e.g. A/C with respect to G/A) a
                            // warning/error-message will be written to the
                            // specified error-log file.
                        }
                        else { // treatment as in the code below
                               // NOTE: This code duplication saves checking the proxy case in the non-proxy algorithms
                               //       which might result in a minor speed-up
                               // case 1 (eff == eff2 && ref == ref2) implicitly ignored and not switched
                            if (eff != eff2) {
                                // ignore if complements agree -- case 3
                                if ( eff==complement(ref)) //BUGFIX C/G A/T
                                {
					markers[i].switch_effect();
                                }
                                else if (eff2 == complement(eff) and ref2 == complement(ref)) continue;
                                // no other allele given? - just check effect alelles
                                else if (ref2 == '0' && eff != eff2) markers[i].switch_effect();
                                // case "swapped" -- case 2 and 4 //note: A/T and C/G poly correct under the assumption of identical strand
                                else if ((eff2 == ref && ref2 == eff) or (eff2 == complement(ref) && ref2 == complement(eff)))
                                    markers[i].switch_effect();
                                // total allele mismatch? -- own case -- ignore if alleles cannot be identical!
                                else markers[i].set_valid(false);
                            }
                        }
                    }
                }
            }
            else {
                for (std::size_t i = index; i < markers.size(); i++) {
                    if (markers[i].is_valid()) {
                        eff2 = markers[i].get_effect_allele();
                        ref2 = markers[i].get_other_allele();
                        // case 1 (eff == eff2 && ref == ref2) implicitly ignored and not switched
                        if (eff != eff2) {
                            if ( eff==complement(ref)) //BUGFIX C/G A/T
                            {
                                markers[i].switch_effect();
                            }
                            else if (eff2 == complement(eff) and ref2 == complement(ref)) continue;
                            // no other allele given? - just check effect alelles
                            else if (ref2 == '0' && eff != eff2) markers[i].switch_effect();
                            // case "swapped" -- case 2 and 4 //note: A/T and C/G poly correct under the assumption of identical strand
                            else if ((eff2 == ref && ref2 == eff) or (eff2 == complement(ref) && ref2 == complement(eff)))
                                markers[i].switch_effect();
                            // total allele mismatch? -- own case -- ignore if alleles cannot be identical!
                            else markers[i].set_valid(false);
                        }
                    }
                }
            }
        }
    }

    //! trims the input vector for valid ones for the first project
    void trim_markers(std::vector<std::vector<MetaMarker> >&markers) {
        if (ud.verbose) std::cout << "Trimming ... "; std::cout.flush();

        // boolean flag - file present or not?
        bool infileflag = false;

        // list of allowed marker IDs
        std::vector<std::string> allowed_ids;
        // interpret trim as a file of IDs
        if (ud.trim.size() > 0 and ud.trim != "NONE") {
            // open input stream depending on the file suffix
            std::ifstream istream; // must be declared before any filtering stream
            boost::iostreams::filtering_stream<boost::iostreams::input> inputfile;
            if (ud.trim.find(".gz") != std::string::npos) { // dealing with a gzipped file?
                istream.open(ud.trim.c_str());
                if (istream) {
                    infileflag = true;
                    inputfile.push(boost::iostreams::gzip_decompressor());
                    inputfile.push(istream);
                }
            }
            else {
                istream.open(ud.trim.c_str());
                if (istream) {
                    infileflag = true;
                    inputfile.push(istream);
                }
            }
            // only proceed, if we are indeed dealing with a file
            if (infileflag) {
                std::string line;
                while (!inputfile.eof()) {
                    std::getline(inputfile, line);
                    // skip if first character is '#' or when we have an empty line
                    if (not line.size() or line[0] == '#') continue;
                    allowed_ids.push_back(line);
                }
            }
        }
        // interpret trim as comma separated list of IDs on the command line
        if (ud.trim.size() > 0 and ud.trim != "NONE" and !infileflag) {
            boost::algorithm::split(allowed_ids, ud.trim, boost::is_any_of(","),
                                    boost::token_compress_on);
        }
        // interpret as not set -> read in first data file and take those markers
        else if (ud.trim == "NONE") {
            for (size_t i = 0; i < markers.size(); i++) {
                if (markers[i][0].get_project_id() == 1 and markers[i][0].is_valid()) {
                    allowed_ids.push_back(markers[i][0].get_rsid());
                }
            }
        }

        // final list of markers
        std::vector<std::vector<MetaMarker> > trimmed;
        #pragma omp parallel for
        for (size_t i = 0; i < markers.size(); i++) {
             if (std::find(allowed_ids.begin(), allowed_ids.end(),
                           markers[i][0].get_rsid()) != allowed_ids.end()) {
                 #pragma omp critical
                 {
                 trimmed.push_back(markers[i]);
                 }
             }
        }
        // delete old content
        markers.erase(markers.begin(), markers.end());
        // copy new to old
        markers.assign(trimmed.begin(), trimmed.end());
        // clean up
        trimmed.erase(trimmed.begin(), trimmed.end());

        if (ud.verbose) {
           unsigned int nmarkers = markers.size();
            std::cout << "Done" << std::endl
                      << "After trimming " << nmarkers
                      << ((nmarkers > 1) ? " markers are left." : " marker is left.") << std::endl;
        }
    }
} // end namespace yamas

