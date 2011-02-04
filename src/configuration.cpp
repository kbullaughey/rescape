/******************************************************************************/
/* Implementation of class Configuration */
/******************************************************************************/

#include "configuration.h"

#include <iostream>

using std::ostream;

Configuration::Configuration(int distinct_tfs) : info(), 
      tf_counts((int)0, distinct_tfs) { 
   hash = 0;
   csize = 0;
}
Configuration::~Configuration() { }

/* add a TF to the configuration */
/* here the TF is 1-indexed */
void Configuration::add(int tf, int p) { 
   pair<int,int> i;
   i.first = tf;
   i.second = p;
   info.push_back(i);
   csize++;
   tf_counts[tf]++;

   /* update our hash with the new tf and pos */
   uint32_t tmp;
   hash += (uint16_t)tf;
   tmp = ((uint16_t)p << 11) ^ hash;
   hash = (hash << 16) ^ tmp;
   hash += hash >> 11;
}

/* clear out the configuration object, to reuse it */
void Configuration::clear(void) {
   hash = 0;
   csize = 0;
   info.clear();
   for (int i=0; i < (int)tf_counts.size(); i++) {
      tf_counts[i] = 0;
   }
}

/* perform the last bit to avalanch the bits */
size_t Configuration::finalize_hash() const {
   size_t res = hash;
   res ^= res << 3;
   res += res >> 5;
   res ^= res << 4;
   res += res >> 17;
   res ^= res << 25;
   res += res >> 6;
   return res;
}

/* read-only index operator */
const pair<int,int>& Configuration::operator[](unsigned int x) const {
   return info[x];
}

/* print a configuration */
ostream& operator<<(ostream &s, Configuration &c) {
   for (int i = 0; i < (int)c.csize; i++) {
      if (i > 0) s << ",";
      s << c.info[i].first << ":" << c.info[i].second;
   }
   return s;
}

/* END */
