#ifndef __CHAIN_H__
#define __CHAIN_H__

#include "array.h"

#include <vector>
#include <valarray>

class SegalModelImportanceSampler;
class LandscapeSegalModelFitted;

/* Segal parameter object that implements Parameter */
class SegalModelParameter {

public:
   /* public member functions */
   SegalModelParameter(std::valarray<double> &lm, 
      std::valarray<double> &a, Array &g);
   ~SegalModelParameter();
   friend std::ostream& operator<<(std::ostream &s, SegalModelParameter &p);
   
   /* public member variables */
   std::valarray<double> lambda; /* lambda of TFs (ability to recruit POL) */
   std::valarray<double> alpha; /* vector to rescale expression levels */
   Array gamma; /* MxMxD symmetrical matrix specifying cooperativity params */

private:
};

/* MCMC chain class, designed for parallel tempering */
class Chain {

public:
   /* public member functions */
   Chain(double t, SegalModelParameter &smp, LandscapeSegalModelFitted *lsp);
   ~Chain();
   friend std::ostream& operator<<(std::ostream &s, Chain &ch);

   /* public member variables */
   SegalModelParameter theta;
   double energy;
   double temp;

protected:
   LandscapeSegalModelFitted *ls;
};

/* class for a set of parallel tempered chains */
class ParallelTemperedChains {

public:
   /* public member functions */
   ParallelTemperedChains(LandscapeSegalModelFitted *lsp);
   ~ParallelTemperedChains(void);
   void new_chain(double temp, SegalModelParameter &start_param);
   unsigned int num_chains(void);
   void run_chains(int num_samples);
   double energy_function(SegalModelParameter &theta);
   double propose(SegalModelParameter &orig, SegalModelParameter &proposed);
   void init_energies(void);

   /* public member variables */
   std::vector<SegalModelImportanceSampler*> samplers;

private:
   /* private member functions */
   void advance(void);

   /* private member variables */
   std::vector<Chain> chains;
   LandscapeSegalModelFitted *ls;
   int current_step;
};

#endif /* __CHAIN_H__ */
