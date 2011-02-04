#ifndef __ARRAY_H__
#define __ARRAY_H__

#include <valarray>
#include <ostream>

using std::slice;

enum Nop {nop};

int prod(std::valarray<int> x);

class Array {
public:
   Array(int n, int len);
   void set_dims(Nop nop, ...);

   ~Array() { }
   double& pos(Nop nop, ...);
   double& posi(int i);
   int size(void) { return mem.size(); }
   friend std::ostream& operator<<(std::ostream &s, const Array &a);

   std::valarray<int> dims;
   std::valarray<double> mem;
private:
};

#endif /* __ARRAY_H__ */
