#ifndef __CONFIGURATION_H__
#define __CONFIGURATION_H__

#include <cstring>
#include <valarray>
#include <ext/hash_map>
#include <utility>

using std::valarray;
using std::pair;

/*******************************************************************************
 * class to store sampled configurations
 ******************************************************************************/

class Configuration {

public:
   /* public member functions */
   Configuration(int distinct_tfs);
   ~Configuration();
   void add(int tf, int pos);
   void clear(void);
   size_t finalize_hash() const;
   friend std::ostream& operator<<(std::ostream &s, Configuration &c);
   double prob_expr_given_config(std::valarray<double> &lambda, int L);
   const pair<int,int>& operator[](unsigned int x) const;

   /* public member variables */
   std::vector< pair<int,int> > info;
   std::valarray<int> tf_counts;
   unsigned int csize;
   struct is_equal { /* comparison function for two configuration objs */
      bool operator()(const Configuration &c1, const Configuration &c2) const {
         if (c1.csize != c2.csize) return 0;
         for (int i = 0; i < (int)c1.csize; i++) {
            if (c1.info[i] != c2.info[i]) return false;
         }
         return true;
      }
   };

private:
   /* private member variables */
   uint32_t hash;
};

/* hash specialization */
namespace __gnu_cxx {
   template<> struct hash<Configuration> {
      size_t operator()(const Configuration &c) const {
         return c.finalize_hash();
      }
   };
}

/* in each storage bin is a pair containing an int, for the number of times 
 * sampled and a double for the weight of this config */
typedef __gnu_cxx::hash_map<Configuration, std::pair<int,double>, __gnu_cxx::hash<Configuration>, 
   Configuration::is_equal> ConfigHash;

#endif /* __CONFIGURATION_H__ */

