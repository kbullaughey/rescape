#include "sim_rand.h"
#include "landscape.h"

/* Construct a random number generator starting with seed */
RandomNumberGenerator::RandomNumberGenerator(unsigned long seed) {
   rng = gsl_rng_alloc(gsl_rng_taus2);
   gsl_rng_set(rng, seed);
}

/* Destroy the random number generator */
RandomNumberGenerator::~RandomNumberGenerator() {
   gsl_rng_free(rng);
}

/* generate a random DNA sequence of length len */
void RandomNumberGenerator::uniform_dna(Sequence &s, unsigned int len,
      Landscape *ls) {
   for (unsigned int i=0; i<len; i++) {
      s.replace(i, (char)uniform_int(0, s.alphabet_size()));
   }
}

/* END */



