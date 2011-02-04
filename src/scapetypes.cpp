#include "scapetypes.h"

#include <iostream>
#include <valarray>

using std::ostream;
using std::valarray;

ostream& operator<<(ostream &s, valarray<double> &x) {
   for (std::size_t i = 0; i < x.size(); i++) {
      if (i!=0) s << ","; 
      s << x[i];
   }
   return s;
}
