#include "array.h"

#include <valarray>
#include <cstdarg>
#include <iostream>
#include <ostream>
#include <ios>

using std::valarray;
using std::cerr;
using std::endl;
using std::slice_array;
using std::slice;
using std::ostream;

Array::Array(int n, int len) : dims(n), mem(len) {
   for (int i=0; i < (int)mem.size(); i++) 
      mem[i] = 0;
}

int prod(std::valarray<int> x) {
   int res = 1;
   for (int i=0; i < (int)x.size(); i++) {
      res *= x[i];
   }
   return res;
}

double& Array::pos(Nop nop, ...) {
   va_list ap;
   va_start(ap, nop);
   int dimprod = 1;
   int cell = 0;
   for (int i = 0; i < (int)dims.size(); i++) {
      int x = va_arg(ap, int);
      cell += dimprod*x;
      dimprod *= dims[i];
   }
   va_end(ap);
   return mem[cell];
}

void Array::set_dims(Nop nop, ...) {
   va_list ap;
   va_start(ap, nop);
   for (int i = 0; i < (int)dims.size(); i++) {
      int d = va_arg(ap, int);
      dims[i] = d;
   }
   va_end(ap);
   return;
}

double& Array::posi(int i) {
   return mem[i];
}

/* A print stream operator for Array objects */
ostream& operator<<(ostream &s, const Array &a) {
   int p = prod(a.dims);
   std::ios_base::fmtflags old = s.flags();
   s.setf(std::ios_base::fixed);
   s.precision(9);
   s << "dims has product p=" << p << endl << "    ";
   for (int i=0; i <(int)a.mem.size(); i++) {
      s << a.mem[i] << " ";
      if ((i+1) % 6 == 0)
         s << endl << "    ";
   }
   s.flags(old);
   return s;
}


/* END */



