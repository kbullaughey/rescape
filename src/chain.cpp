#include "chain.h"
#include "array.h"
#include "scapetypes.h"
#include "landscape.h"

#include <valarray>
#include <iostream>
#include <vector>
#include <omp.h>

using std::valarray;
using std::cout;
using std::ostream;
using std::endl;
using std::cerr;
using std::vector;

/******************************************************************************
 * SegalModelParameter implements Parameter interface
 *****************************************************************************/

/* constructor */
SegalModelParameter::SegalModelParameter(valarray<double> &lm, 
      valarray<double> &a, Array &g) : lambda(lm), alpha(a), gamma(g) { 
}

/* destructor */
SegalModelParameter::~SegalModelParameter() { }

/* print the state of the chain */
ostream& operator<<(ostream &s, SegalModelParameter &theta) {
   int M = theta.alpha.size();
   valarray<double> g(theta.gamma.mem[slice(0, (int)(M*M), 1)]);
   g -= 1;
   s << "lambda=" << theta.lambda << " alpha=" << theta.alpha 
      << " gamma=" << g;
   return s;
}

/******************************************************************************
 * Implementation of an MCMC chain class
 *****************************************************************************/

/* constructor */
Chain::Chain(double t, SegalModelParameter &smp, 
      LandscapeSegalModelFitted *lsp) : theta(smp) { 
   temp = t;
   ls = lsp;
   energy = -1;
}

/* destructor */
Chain::~Chain() { }

/* print the state of the chain */
ostream& operator<<(ostream &s, Chain &ch) {
   s << "temp=" << ch.temp << " energy=" << ch.energy << " " << ch.theta;
   return s;
}

/******************************************************************************
 * class for a set of parallel tempered chains
 *****************************************************************************/

/* constructor */
ParallelTemperedChains::ParallelTemperedChains(LandscapeSegalModelFitted *lsp) { 
   ls = lsp;
   current_step = 0;
}

/* destructor */
ParallelTemperedChains::~ParallelTemperedChains(void) { }

/* add a new chain to the set with the given starting parameter */
void ParallelTemperedChains::new_chain(double temp, 
      SegalModelParameter &start_param) {

   Chain ch(temp, start_param, ls);
   cerr << current_step << ": new chain with temp=" << temp << " has energy=" << ch.energy 
      << endl;
   chains.push_back(ch);
}

/* return the number of chains currently in the set */
unsigned int ParallelTemperedChains::num_chains(void) {
   return chains.size();
}

/* run the chains */
void ParallelTemperedChains::run_chains(int num_samples) {

   for (int c = 0; c < (int)chains.size(); c++) {
      cout << "sample: 0 chains[" << c << "]: " << 
         chains[c] << endl;
   }
   /* sample from all the chains */
   for (int i = 1; i <= num_samples; i++) {
      advance();
      for (int c = 0; c < (int)chains.size(); c++) {
         cout << "sample: " << i << " chains[" << c << "]: " << 
            chains[c] << endl;
      }
   }
}

/* advance the chains */
void ParallelTemperedChains::advance(void) {
   current_step++;
   cerr << current_step << ": advancing chains" << endl;
   double u = ls->rng->runif();
   /* currently the relative proportion of swap steps to parallel steps is 
    * hard coded */
   if (u < 0.3) {
      /* swap step */
      int i = ls->rng->uniform_int(0, chains.size()-1);
      cerr << current_step << ": maybe swapping chains " << i << " and " << i+1 << endl;
      cerr << current_step << ": chains[" << i << "]: " << chains[i] << endl;
      cerr << current_step << ": chains[" << i+1 << "]: " << chains[i+1] << endl;
      double pi_i_x_i = exp((-1)*chains[i].energy/chains[i].temp);
      double pi_i_x_ip1 = exp((-1)*chains[i+1].energy/chains[i].temp);
      double pi_ip1_x_i = exp((-1)*chains[i].energy/chains[i+1].temp);
      double pi_ip1_x_ip1 = exp((-1)*chains[i+1].energy/chains[i+1].temp);
      double ratio = pi_i_x_ip1 * pi_ip1_x_i / (pi_i_x_i * pi_ip1_x_ip1);
      double u2 = ls->rng->runif();
      cerr << current_step << ": pi_i_x_i=" << pi_i_x_i << " pi_i_x_ip1=" << pi_i_x_ip1
         << " pi_ip1_x_i=" << pi_ip1_x_i 
         << " pi_ip1_x_ip1=" << pi_ip1_x_ip1 << endl;
      cerr << current_step << ": drew u2=" << u2 << endl;
      if (u2 < ratio) {
         /* swap */
         cerr << current_step << ": swapping" << endl;
         SegalModelParameter tmp_smp(chains[i].theta);
         double tmp_energy = chains[i].energy;
         chains[i].theta = chains[i+1].theta;
         chains[i].energy = chains[i+1].energy;
         chains[i+1].theta = tmp_smp;
         chains[i+1].energy = tmp_energy;
         cerr << current_step << ": chains[" << i << "]: " << chains[i] << endl;
         cerr << current_step << ": chains[" << i+1 << "]: " << chains[i+1] << endl;
      }
   } else {
      /* parallel step */
      cerr << current_step << ": parallel step" << endl;
      for (int c = 0; c < (int)chains.size(); c++) {
         cerr << current_step << ": maybe updating chains[" << c << "]: " << chains[c] << endl;
         SegalModelParameter proposal(chains[c].theta);
         double ratio = propose(chains[c].theta, proposal);
         double new_energy = energy_function(proposal);
         ratio *= exp((-1)*new_energy/chains[c].temp) / 
            exp((-1)*chains[c].energy/chains[c].temp);
         double u2 = ls->rng->runif();
         cerr << current_step << ": sampled u2=" << u2 << " ratio=" << ratio << endl;
         if (u2 < ratio) {
            /* advance */
            chains[c].theta = proposal;
            chains[c].energy = new_energy;
            cerr << current_step << ": updated chains[" << c << "]: " << chains[c] << endl;
         }
      }
   }
}

/* take orig (x) and fill proposed with the proposal (x') conditional on orig
 * Return the ratio: q(x' -> x)/q(x -> x') */
double ParallelTemperedChains::propose(SegalModelParameter &orig, 
      SegalModelParameter &proposed) {
   proposed = orig;
   double u;
//   proposed.lambda[0] += ls->rng->rnorm(0, ls->sd_lambda);
   return 1; // DEBUG
   /* we only propose chainging each parameter with probability 0.3 */
   for (int i=0; i < (int)proposed.lambda.size(); i++) {
      u = ls->rng->runif();
      if (u < 0.3)
         proposed.lambda[i] += ls->rng->rnorm(0, ls->sd_lambda);
   }
   for (int i=0; i < (int)proposed.alpha.size(); i++) {
      u = ls->rng->runif();
      if (u < 0.3) 
         proposed.alpha[i] *= exp(ls->rng->rnorm(0, ls->sd_alpha));
   }
   Array gamma_tmp(2, ls->M*ls->M);
   gamma_tmp.dims[0] = ls->M;
   gamma_tmp.dims[1] = ls->M;
   gamma_tmp.mem = orig.gamma.mem[slice(0, ls->M*ls->M, 1)];

   /* we fill the upper diagnal and then copy it to the lower */
   for (int j=0; j < (int)gamma_tmp.dims[1]; j++) {
      for (int i=0; i <= j; i++) {
         if (ls->rng->runif() < 0.3)
            gamma_tmp.pos(nop, i, j) *= exp(ls->rng->rnorm(0, ls->sd_gamma));
      }
   }
   for (int j=0; j < (int)gamma_tmp.dims[1]; j++) {
      for (int i=j+1; i < (int)gamma_tmp.dims[0]; i++) 
         gamma_tmp.pos(nop, i, j) = gamma_tmp.pos(nop, j, i);
   }
   gamma_tmp.mem -= 1;
   ls->set_gamma(&gamma_tmp.mem, proposed.gamma);

   /* currently we always have a symmetric distribution */
   return 1;
}

/* calculate the likelihood for theta */
double ParallelTemperedChains::energy_function(SegalModelParameter &theta) {

   double sumabs = 0;
   double expr, res;
   /* everything we need is already stored in the samplers */

   /* this loop is parallelized using Open MP */
   omp_set_num_threads(ls->slots);
#pragma omp parallel for private(expr)
   for (int i = 0; i < (int)samplers.size(); i++) {
      expr = samplers[i]->importance_prob_expr(theta.lambda, theta.alpha, theta.gamma);
      cerr << current_step << ": energy_function: expr=" << expr << " obs=" 
         << samplers[i]->observed_expression << endl;
      res = fabs(expr - samplers[i]->observed_expression);
#pragma omp atomic
      sumabs += res;
   }
   cerr << current_step << ": energy_function: returning " << sumabs << endl;
   return sumabs;
}

/* set the energy for each chain */
void ParallelTemperedChains::init_energies(void) {
   for (int c=0; c < (int)chains.size(); c++) 
      chains[c].energy = energy_function(chains[c].theta);
}

/* END */





