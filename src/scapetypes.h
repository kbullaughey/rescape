#ifndef __SCAPETYPES_H__
#define __SCAPETYPES_H__

#include <string>
#include <string.h>
#include <iostream>
#include <valarray>

#define IS_BINS 100000 /* number of bins for the imporatnce sample hash */

/* cstring compare object for maps and hashes */
struct eqstr {
  bool operator()(const char* s1, const char* s2) const {
    return strcmp(s1, s2) == 0;
  }
};

std::ostream& operator<<(std::ostream &s, std::valarray<double> &x);

#endif /* __SCAPETYPES_H__ */
