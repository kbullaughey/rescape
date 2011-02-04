#ifndef __SIM_RAND_H
#define __SIM_RAND_H

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "sequence.h"

#include <valarray>

using std::valarray;

class Landscape;

/* Class encapsulating the gsl random number geneartor */
class RandomNumberGenerator {

public:
   RandomNumberGenerator(unsigned long seed);
   ~RandomNumberGenerator(void);
   double rnorm(double mean, double sd) { 
      return gsl_ran_gaussian(rng, sd)+mean; }
   int rpois(double lambda) { return gsl_ran_poisson (rng, lambda); }
   void multinomial(int N, const valarray<double> p, unsigned int result[]) {
      gsl_ran_multinomial(rng, p.size(), N, &p[0], result); }
   double runif(void) { return gsl_ran_flat(rng, 0, 1); }
   int uniform_int(int a, int b) {
      return (int)gsl_ran_flat(rng, (double)a, (double)b); }
   void set_seed(unsigned long seed) { gsl_rng_set(rng, seed); }
   void uniform_dna(Sequence &s, unsigned int len, Landscape *ls);
   unsigned int rnb(double n, double p) {
      return gsl_ran_negative_binomial(rng, p, n);
   }

private:
   gsl_rng *rng;
};

#endif /* __SIM_RAND_H */
