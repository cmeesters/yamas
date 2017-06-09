/********************************************************************/
/*                                                                  */
/* YAMAS - Yet Another Meta-Analysis Software                       */
/*                                                                  */
/* Functions dealing with LD properties                             */
/*                                                                  */
/********************************************************************/

// $Rev: 906 $
// $LastChangedDate: 2012-06-04 10:04:45 +0200 (Mo, 04 Jun 2012) $
// $Author: meesters $

#include <iomanip> // stream manipulators
#include <vector>
#include <map>
#include <string>
#include <cmath>

#include <boost/tuple/tuple_io.hpp>

#include "./ya_metamarker.hpp"
#include "./ya_io.hpp"

#ifndef YA_LDDATA
#define YA_LDDATA

namespace yamas {
    //! precision threshold
    const double EPSILON = 0.009;

    //! LD-block
    typedef boost::tuple<unsigned int, unsigned int> LDBlock;

    /*! lookup for LD-blocks
        The key is the block id
    */
    typedef std::map<unsigned int, LDBlock> LDBlocksMap;

    //! lookup for LD-blocks per chromosome
    typedef std::map<unsigned short int, LDBlocksMap> BlocksLT;

    //! map of {block id: vector<markers>}
    typedef std::map<unsigned int, std::vector<yamas::MetaMarker> > LDMarkers;

    //! checker function to determine, whether a marker is in a LD-block
    inline bool isin(const yamas::MetaMarker &marker, const LDBlock &block) {
        return block.get<0>() <= marker.get_position() and marker.get_position() <= block.get<1>();
    }

    //! map: {chromosome: (begin, end)}
    typedef std::map<short int, 
	    boost::tuple<long unsigned int, long unsigned int> > ChrMap;

    /*! function to parse proxy marker file 
        return value is a map, where each key is a chromosome and
	the value is a line number tuple of begin and end lines of the
	chromosome.
    */
    ChrMap proxy_parse(const std::string &fname);

    /*! function to read in LD-block data
    */
    BlocksLT readldblocks(const std::string &fname);

    /*! \brief class which constitutes a proxy marker

        \class A proxy marker is defined by its ID, an allele and its
               in-phase allele, a distance (in base pairs) to a its
               proxy marker, along with the r2 value. In addition, information
               about a strand switch is stored.
    */
    class Proxy {
	friend std::ostream& operator<<(std::ostream& out, const Proxy &proxy);
        friend bool operator<(const Proxy &rhs, const Proxy &lhs);
	public:
            Proxy(const Proxy &p);
	    Proxy(const std::string  &marker_id_ = "",
                  const char         &allele_    = ' ',
                  const char         &inphase_allele_ = ' ',
		  const unsigned int &distance_  = 0,
		  const double       &r2_        = 0.0,
                  const bool         &switch_    = false);
	    ~Proxy();

            // no setters as only getters are needed
	    inline std::string  get_marker_id() const {return marker_id_;};
	    inline char         get_allele()    const {return allele_;};
	    inline char         get_inphase_allele()    const {return inphase_allele_;};
	    inline unsigned int get_distance()  const {return distance_;};
	    inline double       get_r2()        const {return r2_;};
            //! was there a switch defined for that proxy?
            inline bool         get_switch()    const {return switch_;};

            bool operator<(const Proxy &other) {
		// difference of r² values
                double diff = this->r2_ - other.r2_;
		// if r² values are below the resolution threshold, compare distances
                if (std::fabs(diff) < EPSILON) {
                    return (this->distance_ > other.distance_);
		}
		return (diff < 0);
	    }; 

            bool operator<=(const Proxy &other) {
                // difference of r² values
                double diff = this->r2_ - other.r2_;
		// if r² values are below the resolution threshold, compare distances
                if (std::fabs(diff) < EPSILON) {
                    return (this->distance_ >= other.distance_);
		}
		return (diff <= 0);
            };

            bool operator>=(const Proxy &other) {
                // difference of r² values
                double diff = this->r2_ - other.r2_;
		// if r² values are below the resolution threshold, compare distances
                if (std::fabs(diff) < EPSILON) {
                    return (this->distance_ <= other.distance_);
		}
		return (diff >= 0);
            };

	    bool operator>(const Proxy &other) {
                // difference of r² values
                double diff = this->r2_ - other.r2_;
		// if r² values are below the resolution threshold, compare distances
                if (std::fabs(diff) < EPSILON) {
                    return (this->distance_ < other.distance_);
		}
		return (diff > 0);
	    };

	    const Proxy& operator=(const Proxy &p) {
                marker_id_ = p.marker_id_;
                allele_    = p.allele_;
                inphase_allele_ = p.inphase_allele_;
		distance_  = p.distance_;
		r2_        = p.r2_;
		return *this;
	    };

        private:
            std::string  marker_id_;
	    char         allele_;
            char         inphase_allele_;
	    unsigned int distance_;
	    double       r2_;
            bool         switch_;
    };
    
    /*! this operator is necessary for the std::sort algorithm
    */
    inline bool operator<(const Proxy &lhs, const Proxy &rhs) {
        // difference of r² values
        double diff = lhs.r2_ - rhs.r2_;
        // if r² values are below the resolution threshold, compare distances
        if (std::fabs(diff) < EPSILON) {
            return (lhs.distance_ > rhs.distance_);
        }
        return (diff < 0);
    }

    //! map of {markerid : vector of Proxies}
    typedef std::map<std::string, std::vector<Proxy> > ProxyLT;

    /*! returns a vector of string keys from a std::map
    */
    inline std::vector<std::string> get_keys(const ProxyLT &in) {
        std::vector<std::string> out;
        ProxyLT::const_iterator it;
        for (it = in.begin(); it != in.end(); it++) {
            out.push_back(it->first);
        }
        return out;
    }

    //TODO: enter description
    yamas::ProxyLT readldfile(const std::string &fname, 
                              const yamas::ChrMap &table, 
                              const short int &c);

} // end namespace yamas

#endif // YA_LDDATA
