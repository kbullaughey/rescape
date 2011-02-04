#include "landscape.h"
#include "error_handling.h"
#include "scapetypes.h"

#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include <ostream>
#include <istream>
#include <vector>
#include <valarray>
#include <utility>
#include <algorithm>
#include <map>

using std::string;
using std::ifstream;
using std::istream;
using std::ostream;
using std::endl;
using std::cerr;
using std::cout;
using std::vector;
using std::valarray;
using std::pair;
using std::slice;
using std::min;
using std::map;

/* storage for static variable used by FitFunc */
map<string, FitCalc> FitFunc::funclist;

/*****************************************************************************
 * Abstract Segal class implementation
 ****************************************************************************/

/* constructor for segal model abstract class */
LandscapeSegalModel::LandscapeSegalModel(unsigned int num_tfs, 
      string &pwm_file, vector< valarray<double>* > &conds, 
      unsigned int d_max, RandomNumberGenerator *rngp, double logistic_par, int min_tf_sep) : 
      Landscape(rngp), r(num_tfs), pwms(num_tfs), lambda((double)0, num_tfs+1),
      gamma(3, num_tfs*num_tfs*d_max), alpha((int)1, num_tfs) {
   M = num_tfs;
   D = d_max;
   load_pwms(pwm_file);
   conditions = conds;
   logistic_parameter = logistic_par;
   minimum_tf_separation = min_tf_sep;
   return;
}

/* destructor */
LandscapeSegalModel::~LandscapeSegalModel() {
   for (vector<PWM*>::iterator i = pwms.begin(); i != pwms.end(); i++)
      if (*i != 0) delete *i;
   if (conditions.size() != 0) {
      for (vector< valarray<double>* >::iterator i = conditions.begin();
            i != conditions.end(); i++)
         delete *i;
   }
}

/* calculate the degree of cooperativity */
double LandscapeSegalModel::coopfunc(double g, int d) {
   return 1+g*exp((-1)*pow(d,2)/80);
}

/* set up a gamma vector of all zero */
void LandscapeSegalModel::set_gamma(valarray<double> *gm, Array &local_gamma) {
   local_gamma.dims[0] = M;
   local_gamma.dims[1] = M;
   local_gamma.dims[2] = D;
   if (gm != 0) {
      for (int i=0; i < (int)M; i++) {
         for (int j=0; j < (int)M; j++) {
            for (int k=0; k < (int)D; k++) 
               local_gamma.pos(nop, i, j, k) = coopfunc((*gm)[((int)M)*i+j], k);
         }
      }
   } else {
      for (int i=0; i < (int)M; i++) {
         for (int j=0; j < (int)M; j++) {
            for (int k=0; k < (int)D; k++) 
               local_gamma.pos(nop, i, j, k) = coopfunc(0, k);
         }
      }
   }
   return;
}

/* load the PWMs frim a file */
void LandscapeSegalModel::load_pwms(const string &fn) {
   ifstream in(fn.c_str(), std::fstream::in);
   if (!in.good()) throw SimError("failed to open "+fn);
   for (unsigned int i=0; i < M; i++) {
      pwms[i] = new PWM(alphabet);
      in >> *(pwms[i]);
      if (in.bad()) throw SimError("failed to read from "+fn);
      r[i] = pwms[i]->w;
      if (r[i] <= 0) throw SimError("invalid TF width, check pwm file");
   }
   in.close();
}

/* calculate the probability of expression given a configuration */
double LandscapeSegalModel::pec(const Configuration &c, valarray<double> &s, int L) {
   double sum_contribution = s[0];
   for (int i=0; i < (int)c.tf_counts.size(); i++) {
      sum_contribution += c.tf_counts[i] * s[i+1];
   }
   sum_contribution = sum_contribution/L*1000;
   return 1/(1+exp((-1)*sum_contribution/logistic_parameter));
}

/* calculate cumulative lookup tables */
void LandscapeSegalModel::calc_prob_tables(valarray<double> *rel_tf_expr, 
      valarray<double> &alpha, Array &gm, const Sequence &seq, 
      ProbTable &ret) {

   int L = (int)seq.length();

   for (int x=0; x < (int)L; x++) {
      for (int i=0; i < (int)M; i++) {
         double pwm_tmp = p_pwm(i, x, seq);
         for (int j=0; j <= (int)M; j++) {
            for (int d=0; d < (int)D; d++) {
               if (x+1 - r[i] < 0) {
                  ret.cpa.pos(nop, i,j,d,x) = 0; 
                  continue;
               }
               if (x+1 - r[i] == 0 && j == 0 && d == 0) {
                  ret.cpa.pos(nop, i,j,d,x) = alpha[i] * 
                     (*rel_tf_expr)[i] * pwm_tmp;
                  continue;
               }
               if (j > 0) { /* the previous TF is less than D away */
                  if (x+1 - r[i] - d - r[j-1] < 0) {
                     ret.cpa.pos(nop, i,j,d,x) = 0; continue;
                  }
                  if (d < minimum_tf_separation) {
                     /* we forbid TFs this close together */
                     ret.cpa.pos(nop, i,j,d,x) = 0; continue;
                  }
                  ret.cpa.pos(nop, i,j,d,x) = 
                     ret.cpt.pos(nop, j-1,x-r[i]-d) * 
                     gm.pos(nop, j-1, i, d) * alpha[i] * 
                     (*rel_tf_expr)[i] * pwm_tmp;
                  continue;
               } else { /* previous TF is D or further away or no other TF */
                  if (d > 0) { /* we don't use this slot */
                     ret.cpa.pos(nop, i,0,d,x) = 0; continue;
                  }
                  if (x+1 - (int)r[i] - (int)D <= 0) { /* no other TF */
                     ret.cpa.pos(nop, i,0,0,x) = alpha[i] * 
                        (*rel_tf_expr)[i] * pwm_tmp;
                     continue;
                  }
                  ret.cpa.pos(nop, i,0,0,x) = (1+ret.cpr.pos(nop, x-r[i]-D)) * 
                     alpha[i] * (*rel_tf_expr)[i] * pwm_tmp;
                  continue;
               }
            }
         }
         valarray<double> ans1(ret.cpa.mem[slice(x*M*(M+1)*D+i, (M+1)*D, M)]);
         ret.cpt.pos(nop, i,x) = ans1.sum();
      }
      if (x == 0) {
         valarray<double> ans2(ret.cpt.mem[slice(0, M, 1)]);
         ret.cpr.posi(0) = ans2.sum();
      } else {
         valarray<double> ans3(ret.cpt.mem[slice(x*M, M, 1)]);
         ret.cpr.posi(x) = ret.cpr.posi(x-1) + ans3.sum();
      }
   }

   /* make the margin tables */
   for (int k1 = 0; k1 < (int)M; k1++) {
      for (int k2 = 0; k2 < (int)L; k2++) {
         double csum = 0;
         ret.cpm.pos(nop, 0, k1, k2) = 0;
         for (int k3 = 0; k3 < (int)((M+1)*D); k3++) {
            csum += ret.cpa.mem[k2*M*(M+1)*D+k1+M*k3];
            ret.cpm.mem[((M+1)*D+1)*(k2*M+k1)+1+k3] = csum;
         }
      }
   }

   /* now make s_array cumulative, we need the first two values to be 0, 1
    * and shifting the result of the cumulative sum over by 2 */
   double tmp1 = 1, tmp2, tmp3;
   int i;
   tmp2 = ret.cpa.posi(0);
   tmp3 = ret.cpa.posi(1);
   ret.cpa.posi(0) = 0;
   for (i=2; i<prod(ret.cpa.dims); i++) {
      ret.cpa.posi(i-1) = tmp1;
      tmp1 += tmp2;
      tmp2 = tmp3;
      tmp3 = ret.cpa.posi(i);
   }
   ret.cpa.posi(i-1) = tmp1;
   tmp1 += tmp2;
   ret.cpa.posi(i) = tmp1;
   tmp1 += tmp3;
   ret.cpa.posi(i+1) = tmp1;
   
   /* now we add one to s_row to account for the empty configuration */
   ret.cpr.mem += 1;
}

/* calculate probability of tf binding to the sequence with right edge at pos */
double LandscapeSegalModel::p_pwm(unsigned int tf, int x, const Sequence &seq) {
   int a = (int)alphabet.length();
   if (x-r[tf]+1 < 0 || x >= (int)seq.length())  
      return 0;
   Sequence part(seq, x+1-r[tf], r[tf]);
   /* I assume the affinity of a TF is the sum of the affinity to both strands
    * at this location */
   double res1 = 1, res2 = 1;
   for (int k=0; k<r[tf]; k++) {
      res1 *= (*(pwms[tf]->pwm))[k*a+part.code(k)];
   }
   for (int k=0; k<r[tf]; k++) {
      res2 *= (*(pwms[tf]->pwm))[(r[tf] - k - 1)*a + (a - part.code(k) - 1)];
   }
   res1 = (res1 + res2);
   res1 *=  pow(a, r[tf]);
   return res1;
}

/*****************************************************************************
 * Fitted model class implementation
 ****************************************************************************/

/* currently for our importance distribution, the expression scalars, alpha, 
 * start off as all 1s and there are no cooperativity interactions */

/* constructor */
LandscapeSegalModelFitted::LandscapeSegalModelFitted(
      unsigned int num_tfs, string &pwm_file, 
      vector< valarray<double>* > &conds, unsigned int d_max, 
      unsigned int num_steps, RandomNumberGenerator *rngp,
      vector< valarray<double>* > &obs_expr,
      double sdlambda, double sdalpha, double sdgamma,
      vector<string> &seqs, int nslots, double logistic_par, int min_tf_sep) : 
      LandscapeSegalModel(num_tfs, pwm_file, conds, d_max, rngp, logistic_par, min_tf_sep)
      /*, alpha_lattice_points(2, (int)pow((double)es_lattice->size(), 
         (double)num_tfs)) */,
      parallel_tempered(this) {
   set_gamma(NULL, gamma);
   obs_expr_profiles = obs_expr;
   sd_lambda = sdlambda;
   sd_alpha = sdalpha;
   sd_gamma = sdgamma;

   /* check to make sure each profile is the right length */
   for (vector< valarray<double>* >::iterator i = obs_expr_profiles.begin();
         i != obs_expr_profiles.end(); i++) {
      if ((*i)->size() != conditions.size()) {
         throw SimError("length of observed expression profile doesn't match"
            " number of conditions");
      }
   }

   /* we need to pick starting lambda. So we set each to either 0.5, or -0.5, 
    * with lambda[0] = 0 */
   double choices[2] = {-0.5, 0.5};
   for (int i = 0; i < (int)lambda.size(); i++) {
      lambda[i] = choices[rng->uniform_int(0, 2)];
      cerr << "picked lambda[" << i << "] = " << lambda[i] << endl;
   }

   /* we use this to control the number of threads when using openmp */
   slots = nslots;
}

/* destructor */
LandscapeSegalModelFitted::~LandscapeSegalModelFitted() {
   /* we assume that when we're done with the landscape object we don't need
    * the expression profiles */
   if (obs_expr_profiles.size() != 0) {
      for (vector< valarray<double>* >::iterator i = obs_expr_profiles.begin();
            i != obs_expr_profiles.end(); i++)
         delete *i;
   }
}

/* create the parallel tempering MCMC chains */
void LandscapeSegalModelFitted::create_chains(int n, valarray<double> &chain_temps) {

   /* we create storage space for our priors. */
   valarray<double> sample_lambda((double)0, (int)M+1);
   valarray<double> sample_alpha((double)0, (int)M);
   Array gamma_aux(2, (int)(M*M));
   gamma_aux.dims[0] = (int)M;
   gamma_aux.dims[1] = (int)M;
   Array sample_gamma(3, (int)(M*M*D));

   /* I use two arbitrary functions that have approximately the shape I want */
//   double T_min = exp((-1)*(double)n/3)*0.6;
//   double T_max = log(n)/1.5 + 0.5;
   while ((int)parallel_tempered.num_chains() < n) {

      /* For now, the priors are hard coded here 
       * The prior for lambda[i] is N(0, 0.4)
       * The prior for alpha[i] is LogNormal(0, 0.5)
       * The prior for gamma[i] is zero with probability 0.8 and 
       *    LogNormal(0, 0.3)-1 with probability 0.2 */
      for (int i=0; i < (int)sample_lambda.size(); i++) 
         sample_lambda[i] = rng->rnorm(0, 0.4);
      for (int i=0; i < (int)sample_alpha.size(); i++) 
         sample_alpha[i] = exp(rng->rnorm(0, 0.5));

      /* DEBUG: here we temporarily fix these, ignoring the prior for testing purposes */
      sample_lambda[0] = -3;
      sample_lambda[1] = 3;
      sample_lambda[2] = -2;
      sample_lambda[3] = -1;
      sample_alpha[0] = 3;
      sample_alpha[1] = 1;
      sample_alpha[2] = 1;

      /* we fill the upper diagnal and then copy it to the lower */
      for (int j=0; j < (int)gamma_aux.dims[1]; j++) {
         for (int i=0; i <= j; i++) {
/* DEBUG
            if (rng->runif() < 0.8) {
               gamma_aux.pos(nop, i, j) = 0;
            } else {
               gamma_aux.pos(nop, i, j) = exp(rng->rnorm(0, 0.3))-1;
            }
*/
gamma_aux.pos(nop, i, j) = 0;
         }
      }
      for (int j=0; j < (int)gamma_aux.dims[1]; j++) {
         for (int i=j+1; i < (int)gamma_aux.dims[0]; i++) 
            gamma_aux.pos(nop, i, j) = gamma_aux.pos(nop, j, i);
      }
      set_gamma(&gamma_aux.mem, sample_gamma);
   
      SegalModelParameter smp(sample_lambda, sample_alpha, sample_gamma);
      cerr << "chain " << parallel_tempered.num_chains() << 
         " sample_lambda=" << smp.lambda << endl; 
      cerr << "chain " << parallel_tempered.num_chains() << 
         " sample_alpha=" << smp.alpha << endl; 
      cerr << "chain " << parallel_tempered.num_chains() << 
         " gamma_aux=" << gamma_aux.mem << endl; 

   
      /* for now, the chains are all started at the same point, smp */
//      parallel_tempered.new_chain((int)parallel_tempered.num_chains() * 
//         (T_max-T_min)/(n-1) + T_min, smp);
      parallel_tempered.new_chain(
         chain_temps[(int)parallel_tempered.num_chains()], smp);
   }
}

/* for each combination of condition, sequence, and alpha latice point,
 * compute a lookup table to allow effecient iid sampling from the IS distr */
void LandscapeSegalModelFitted::precompute_lookup_tables(vector<string> &seqs, 
      unsigned int steps) {

   /* all our sequence-condition combinations start off with the same 
    * importance distribution */
   SegalModelParameter smp(lambda, alpha, gamma);

   /* loop through sequence-condition point combinations, creating
    * importance samplers */
   for (int s = 0; s < (int)seqs.size(); s++) {
      Sequence seq(seqs[s], alphabet);
      for (int c = 0; c < (int)conditions.size(); c++) {
         SegalModelImportanceSampler *is = 
            new SegalModelImportanceSampler(seq, smp, *conditions[c], 
               (*obs_expr_profiles[s])[c], this, (int)steps);
         parallel_tempered.samplers.push_back(is);
      }
   }
   return;
}

/* not implemented yet */
double LandscapeSegalModelFitted::get_fitness(const Sequence&) {
   SimError("fitness not implemented for SegalModelFitted");
   return 0;
}

/*****************************************************************************
 * Specified model class implementation
 ****************************************************************************/

/* constructor for segal model landscapes */
LandscapeSegalModelSpecified::LandscapeSegalModelSpecified(
      unsigned int num_tfs, string &pwm_file, valarray<double> *li,
      valarray<double> *expr_s,
      vector< valarray<double>* > &conds, valarray<double> *gm, 
      FitFunc &fitfunc, unsigned int d_max, RandomNumberGenerator *rngp,
      double logistic_par, int min_tf_sep) : 
      LandscapeSegalModel(num_tfs, pwm_file, conds, d_max, rngp, logistic_par, min_tf_sep),
      eff(fitfunc) {
   lambda = *li;
   if (expr_s != 0) {
      alpha = *expr_s;
   }
   set_gamma(gm, gamma);
   return;
}

/* destructor */
LandscapeSegalModelSpecified::~LandscapeSegalModelSpecified() {
}

/* calculate the fitness of the given sequence */
double LandscapeSegalModelSpecified::get_fitness(const Sequence &seq) {
   valarray<double> expr(conditions.size());
   get_expression(seq, expr);
   return eff.calc(expr);
}

/* fill pe with the expression levels for seq */
void LandscapeSegalModelSpecified::get_expression(const Sequence &seq, 
      valarray<double> &pe) {

   if (pe.size() != conditions.size()) 
   throw SimError("num. of desired expression levels must match num. cond.");

   /* loop through the conditions first pre-computing tables then est. expr. */
   ProbTable ptables(M, D, seq.length());
   for (unsigned int i=0; i<conditions.size(); i++) {
      calc_prob_tables_wrapper(i, seq, ptables);
      pe[i] = prob_expr(ptables, seq, i);
   }
   return;
}

/* wrapper for calculating cumulative lookup tables */
void LandscapeSegalModelSpecified::calc_prob_tables_wrapper(unsigned int condition, 
      const Sequence &seq, ProbTable &ret) {
   calc_prob_tables(conditions[condition], alpha, gamma, seq, ret);
   return;
}

/* fill moc with mean TF occupancy in condition, cond */
void LandscapeSegalModelSpecified::get_tf_occupancy(const Sequence &seq,
      int cond, Array &oc) {
   if (cond >= (int)conditions.size()) 
      throw SimError("condition number exceeds number of conditions");
   int L = (int)seq.length();

   ProbTable ptables(M, D, L);
   calc_prob_tables_wrapper(cond, seq, ptables);

   int up, min_r, cell, left_tf, right_tf, d, tbl_size, pos;
   int i=nsteps();
   int mmp1d = M*(M+1)*D;
   int mmp1 = M*(M+1);
   int mp1d = (M+1)*D;
   min_r = r.min();
   
   while (i > 0) {
      left_tf = 0;
      up = L;
      d = 0;
      while (up >= min_r) {
         if (left_tf == 0) {
            tbl_size = mmp1d * up + 1;
            cell = sample_cdf_wrapper(tbl_size, ptables.cpr.posi(up-1), 
               ptables.cpa.mem, 0); 
            if (cell == 0) { break; }
            right_tf = ((cell-1) % M)+1;
            pos = ((cell-1) / mmp1d)+1;
            left_tf = ((cell-1) % mmp1) / M;
            d = ((cell-1) % mmp1d) / mmp1;
         } else {
            cell = sample_cdf_wrapper(mp1d, ptables.cpt.pos(nop, left_tf-1, up-1),
               ptables.cpm.mem, ((up-1)*M + left_tf-1) * (mp1d+1));
            right_tf = left_tf;
            pos = up;
            left_tf = cell % (M+1);
            d = cell / (M+1);
         }
         if (left_tf == 0) up = pos - r[right_tf-1] - D;
         else up = pos - r[right_tf-1] - d;
         for (int j = 0; j < r[right_tf-1]; j++) {
            if ((right_tf-1) < 0 || right_tf > (int)M) {
               cerr << "Error: right_tf=" << right_tf << endl;
               break;
            } else if (pos-j < 0 || pos-j > L) {
               cerr << "Error: pos=" << pos << ", j=" << j << endl;
               break;
            }  else {
               oc.pos(nop, right_tf-1, pos-j-1)++;
            }
         }
      }
      i--;
   }
   for (i = 0; i < (int)oc.size(); i++) 
      oc.posi(i) /= nsteps();
   return;
}

/******************************************************************************/
/* generic sampling class */
/******************************************************************************/

/* constructor */
Sampler::Sampler(RandomNumberGenerator *rngp) { 
   rng = rngp; 
}

/* sample from the cummulative distribution of cell probabilities 
 * which - right edge of sub-table used
 * cpt - full cumulative probability table, all cells
 * cpr - marginal cumulative probability row
 */
int Sampler::sample_cdf(int tbl_size, double tbl_sum,
      valarray<double> &tbl, int offset) {
   double u;
   int cell, top, bottom;

   u = rng->runif();
   u *= tbl_sum; 
   top = tbl_size+offset;
   bottom = offset;

   /* find the cell containing u */
   do {
      cell = (top+bottom)/2;
      if (u <= tbl[cell]) top = cell;
      else if (u > tbl[cell+1]) bottom = cell;
      else { bottom = cell; break; }
   } while (bottom+1 != top);

   return bottom-offset;
}

/******************************************************************************/
/* IID sampling class */
/******************************************************************************/

/* constructor */
SegalModelIID::SegalModelIID(unsigned int num_steps, RandomNumberGenerator *rng) : 
      Sampler(rng) {
   steps = num_steps;
}

/* calculate the probability of expression using the precomputed tables */
double SegalModelIID::prob_expr(int L, int M, int D, valarray<int> &r, 
      valarray<double> &lambda, ProbTable &ptables, double logistic_parameter,
      bool print_configs) {
   int up, min_r, cell, left_tf, right_tf, d, tbl_size, pos;
   int i=steps;
   double pe_tmp = 0;
   double sum_lambda;

   min_r = r.min();
   int mmp1d = M*(M+1)*D;
   int mmp1 = M*(M+1);
   int mp1d = (M+1)*D;
   
   while (i > 0) {
      sum_lambda = lambda[0];
      left_tf = 0;
      up = L;
      d = 0;
      while (up >= min_r) {
         if (left_tf == 0) {
            tbl_size = mmp1d * up + 1;
            cell = sample_cdf(tbl_size, ptables.cpr.posi(up-1), 
               ptables.cpa.mem, 0); 
            if (cell == 0) { break; }
            right_tf = ((cell-1) % M)+1;
            pos = ((cell-1) / mmp1d)+1;
            left_tf = ((cell-1) % mmp1) / M;
            d = ((cell-1) % mmp1d) / mmp1;
         } else {
            cell = sample_cdf(mp1d, ptables.cpt.pos(nop, left_tf-1, up-1),
               ptables.cpm.mem, ((up-1)*M + left_tf-1) * (mp1d+1));
            right_tf = left_tf;
            pos = up;
            left_tf = cell % (M+1);
            d = cell / (M+1);
         }
         if (left_tf == 0) up = pos - r[right_tf-1] - D;
         else up = pos - r[right_tf-1] - d;
         sum_lambda += lambda[right_tf];
         if (print_configs) 
            cout << right_tf << ":" << pos << ", ";
      }
      sum_lambda = sum_lambda/L*1000;
      if (print_configs) 
         cout << endl;
      pe_tmp += (1/(1+exp((-1)*sum_lambda/logistic_parameter)));
      i--;
   }
   return pe_tmp/steps;
}

/******************************************************************************/
/* IS sampling class */
/******************************************************************************/

/* constructor */
SegalModelIS::SegalModelIS(unsigned int num_steps, RandomNumberGenerator *rngp) : 
      Sampler(rngp) {
   steps = num_steps;
}

/* for the IS model, we need to pre-sample */
void SegalModelIS::generate_importance_sample(int L, int M, int D, 
      valarray<int> &r, ProbTable &ptables, ConfigHash &importance_sample) {

   int up, min_r, cell, left_tf, right_tf, d, tbl_size, pos;
   int i=steps;

   min_r = r.min();
   int mmp1d = M*(M+1)*D;
   int mmp1 = M*(M+1);
   int mp1d = (M+1)*D;
   Configuration c((int)M);
   
   /* sample the configurations */
   while (i > 0) {
      left_tf = 0;
      up = L;
      d = 0;
      while (up >= min_r) {
         if (left_tf == 0) {
            tbl_size = mmp1d * up + 1;
            cell = sample_cdf(tbl_size, ptables.cpr.posi(up-1), 
               ptables.cpa.mem, 0); 
            if (cell == 0) { break; }
            right_tf = ((cell-1) % M)+1;
            pos = ((cell-1) / mmp1d)+1;
            left_tf = ((cell-1) % mmp1) / M;
            d = ((cell-1) % mmp1d) / mmp1;
         } else {
            cell = sample_cdf(mp1d, ptables.cpt.pos(nop, left_tf-1, up-1),
               ptables.cpm.mem, ((up-1)*M + left_tf-1) * (mp1d+1));
            right_tf = left_tf;
            pos = up;
            left_tf = cell % (M+1);
            d = cell / (M+1);
         }
         if (left_tf == 0) up = pos - r[right_tf-1] - D;
         else up = pos - r[right_tf-1] - d;
         c.add(right_tf-1, pos);
      }
      importance_sample[c].first++;
      c.clear();
      i--;
   }

   return;
}

/******************************************************************************/
/* IID sampling version of specified model */
/******************************************************************************/

/* constructor for importance sampling version */
LandscapeSegalModelSpecifiedIID::LandscapeSegalModelSpecifiedIID(unsigned int num_tfs,
      string &pwm_file, valarray<double> *li, valarray<double> *expr_s,
      vector< valarray<double>* > &conds, 
      valarray<double> *gm, unsigned int num_steps,
      FitFunc &fitfunc, unsigned int d_max, RandomNumberGenerator *rngp, 
      double logistic_par, int min_tf_sep) :
      LandscapeSegalModelSpecified(num_tfs, pwm_file, li, expr_s, conds, gm,
      fitfunc, d_max, rngp, logistic_par, min_tf_sep), SegalModelIID(num_steps, rngp) {
   return;
}

/* destructor */
LandscapeSegalModelSpecifiedIID::~LandscapeSegalModelSpecifiedIID(void) { }

/* call IID version */
double LandscapeSegalModelSpecifiedIID::prob_expr(ProbTable &ptables, 
      const Sequence &seq, int condition, bool print_configs) {
   int L = (int)seq.length();
   return SegalModelIID::prob_expr(L, M, D, r, lambda, ptables, logistic_parameter, 
      print_configs);
}

/* fill wsum with the sum of all configuration weights, computed with 
 * dynamic programming */
void LandscapeSegalModelSpecifiedIID::get_wsum(const Sequence &seq, 
      valarray<double> &ws) {

   if (ws.size() != conditions.size()) 
      throw SimError("num. of wsum calculations must match num. cond.");

   /* loop through the conditions first pre-computing tables then est. expr. */
   ProbTable ptables(M, D, seq.length());
   for (unsigned int i=0; i<conditions.size(); i++) {
      calc_prob_tables_wrapper(i, seq, ptables);
      ws[i] = ptables.cpr.posi(seq.length()-1);
   }
   return;
}

/******************************************************************************/
/* IS sampling version of specified model */
/******************************************************************************/

/* constructor for importance sampling version */
LandscapeSegalModelSpecifiedIS::LandscapeSegalModelSpecifiedIS(unsigned int num_tfs,
      string &pwm_file, valarray<double> *li, valarray<double> *expr_s,
      vector< valarray<double>* > &conds, 
      valarray<double> *gm, unsigned int num_steps,
      FitFunc &fitfunc, unsigned int d_max, RandomNumberGenerator *rngp,
      valarray<double> *tar_gm, valarray<double> *tar_a, double logistic_par,
      int min_tf_sep) : LandscapeSegalModelSpecified(num_tfs, pwm_file, li, 
      expr_s, conds, gm, fitfunc, d_max, rngp, logistic_par, min_tf_sep), 
      SegalModelIS(num_steps, rngp), target_alpha(*tar_a), 
      target_gamma(3, num_tfs*num_tfs*d_max) {
   set_gamma(tar_gm, target_gamma);
   num_importance_samples = num_steps;
   return;
}

/* destructor */
LandscapeSegalModelSpecifiedIS::~LandscapeSegalModelSpecifiedIS(void) { }

/* calculate the probability of expression using the precomputed tables, 
 * asjusting by importance weights */
double LandscapeSegalModelSpecifiedIS::prob_expr(ProbTable &ptables, 
      const Sequence &seq, int condition, bool print_configs) {

   int L = (int)seq.length();

   /* vector, one per condition, of hashes of configs containing how many 
    * of each were sampled */
   ConfigHash importance_sample(IS_BINS); 
   generate_importance_sample(L, (int)M, (int)D, r, ptables, 
      importance_sample);
   calc_importance_sample_energy_weights(importance_sample, seq,
      *conditions[condition]);

   double expr_estimator = 0;
   double sum_importance_weights = 0;
   double importance_weight;
   
   /* loop over all the unique configs in the importance sample */
   for (ConfigHash::iterator i = importance_sample.begin(); 
         i != importance_sample.end(); i++) {

      /* the importance weight is the ratio of the config's energy weight 
       * under the target distribution divided by the config's energy weight
       * under the importance distribution */
      importance_weight = wc(i->first, target_alpha, seq, target_gamma, *conditions[condition]);
      importance_weight /= i->second.second;

      /* each term added to the sum is the product of the importance weight
       * and the probability of expression given the configuration and 
       * the number of times the configuration was seen in the importance 
       * sample */ 
      expr_estimator += importance_weight * pec((*i).first, lambda, L) *
         i->second.first;
      sum_importance_weights += importance_weight * i->second.first;
   }

   /* this estimator might be biased, but generally has a much lower error 
    * according to Liu JS. Monte Carlo Stragegies in Scientific Computing */
   return expr_estimator / sum_importance_weights;
}

/* calculate the energy weights for each unique configuration in the sample */
void LandscapeSegalModelSpecifiedIS::calc_importance_sample_energy_weights(
      ConfigHash &importance_sample, const Sequence &seq, 
      valarray<double> &tf_expr) {
   for (ConfigHash::iterator i = importance_sample.begin(); 
         i != importance_sample.end(); i++) {
      i->second.second = wc(i->first, alpha, seq, gamma, 
         tf_expr);
   }
   return;
}

/* calculate the weight of a configuration given the sequence and params */
double LandscapeSegalModelSpecifiedIS::wc(const Configuration &c, const valarray<double> &sc,
      const Sequence &seq, Array &g, const valarray<double> &expr) {
   double energy_weight = 1;
   int tf, pos, prev_tf, prev_pos;
   for (int i = 0; i < (int)c.csize; i++) {
      tf = c[i].first;
      pos = c[i].second-1;
      energy_weight *= expr[tf] * sc[tf] * p_pwm(tf, pos, seq);
      if (i!=0) {
         prev_tf = c[i-1].first;
         prev_pos = c[i-1].second-1;
         if (prev_pos-pos-r[prev_tf] < (int)D) {
            energy_weight *= g.pos(nop, prev_tf, tf, prev_pos-pos-r[prev_tf]);
         }
      }
   }
   return energy_weight;
}


/******************************************************************************/
/* Implementation of class PWM */
/******************************************************************************/

PWM::PWM(string &s) : alphabet(s) {
   w = 0; h = 0; pwm = 0;
}

/* input operator for PWM objects */
istream& operator>>(istream &in, PWM &x) {
   x.h = x.alphabet.length();
   string pwm_header;
   getline(in, pwm_header);
   if (pwm_header.compare(0, 4, "pwm ") != 0) 
      throw SimError("malformatted pwm. Should be 'pwm ', got "+pwm_header);
   pwm_header.erase(0, 4);
   char *end;
   x.w = strtoul(pwm_header.c_str(), &end, 10);
   if (end == pwm_header.c_str())
      throw SimError("malformed pwm header");
   if (x.w <= 0)
      throw SimError("number of pwm columns must be >= 0");
   x.pwm = new valarray<double>(x.h * x.w);
   for (unsigned int k=0; k<x.h; k++) {
      for (unsigned int j=0; j<x.w; j++) {
         in >> (*x.pwm)[j*x.h+k];
         if(in.fail()) throw SimError("failure reading pwm file");
      }
   }
   in >> std::ws; /* get any remaining white space and new lines */
   return in;
}

/* dump a PWM object */
ostream& operator<<(ostream &s, PWM &x) {
   s << "pwm " << x.w << endl;
   for (unsigned int k=0; k<x.h; k++) {
      for (unsigned int j=0; j<x.w; j++) {
         if (j!=0) s << " ";
         s << (*x.pwm)[j*x.h+k];
      }
      s << endl;
   }
   return s;
}

/* constructor for ProbTable */
ProbTable::ProbTable(unsigned int M, unsigned int D, unsigned int L) :
      cpa(4, (int)(M*(M+1)*D*L+2)), cpt(2, (int)(M*L)), cpr(1, (int)L), 
      cpm(3, (int)((D*(M+1)+1)*M*L)) {

   /* s_array is made cumulative at end. is the total weight of all configs that
    * have a TF, i, as the right most TF and whose right most edge is at most
    * position x going right and have a TF, j, to the left with d bp in between.
    * The first position is 0, followed by 1, followed by the cumulative 
    * probabilities of TFs binding at each position. This array has 2 extra
    * slots at the end, hence the hack here, to allow for later switch into
    * a cumulative distribution with an extra state, (the empty config) added */
   cpa.set_dims(nop, (int)M, (int)(M+1), (int)D, (int)L);

   /* s_table is non cumulative. It is the total probability of all configs
    * that have a right-most TF, i, at position x */
   cpt.set_dims(nop, (int)M, (int)L);

   /* s_row is cumulative. It is the total sum of all configurations for which
    * the right most TF (whichever TF that is) has its right-most edge over 
    * position x. Notice neither s_table or s_row contain the weight of the 
    * empty configuration, which is 1, and which is added on later */
   cpr.set_dims(nop, (int)L);

   /* s_mar is the set of marginal tables. Each marginal table is cumulative
    * and starts with zero followed by the cumulative distribution. It is the 
    * linearized array formed by the following R operation:
    *    apply(s_array, c(1,4), function(m) cumsum(c(0, m))) */
   cpm.set_dims(nop, (int)(D*(M+1)+1), (int)M, (int)L);
}

/******************************************************************************
 * class for accessing an importance sample for calculating expr
 *****************************************************************************/

/* constructor */
SegalModelImportanceSampler::SegalModelImportanceSampler(Sequence &seq, 
      SegalModelParameter &smp, valarray<double> &condition,
      double obs_expr, LandscapeSegalModelFitted *lsp, 
      int num_steps) : SegalModelIS(num_steps, lsp->rng), enhancer(seq), 
      theta(smp), tf_expr(condition), ptables(lsp->M, lsp->D, seq.length()), 
      importance_sample(IS_BINS), 
      pwm_lookup_table((double)0, lsp->M*seq.length()) {
   observed_expression = obs_expr;
   ls = lsp;
   ls->calc_prob_tables(&condition, theta.alpha, theta.gamma, enhancer, ptables);
   generate_importance_sample((int)enhancer.length(), (int)ls->M, (int)ls->D, 
      ls->r, ptables, importance_sample);
   compute_pwm_lookup_table();
   calc_importance_sample_energy_weights();
   return;
}

/* destructor */
SegalModelImportanceSampler::~SegalModelImportanceSampler() { }

/* calculate the energy weights for each unique configuration in the sample */
void SegalModelImportanceSampler::calc_importance_sample_energy_weights() {
   for (ConfigHash::iterator i = importance_sample.begin(); 
         i != importance_sample.end(); i++) {
      i->second.second = wc(i->first, theta.alpha, theta.gamma);
   }
   return;
}

/* use the importance weights to calculate expression */
double SegalModelImportanceSampler::importance_prob_expr(valarray<double> &lam,
       valarray<double> &tar_a, Array &tar_g) {

   int L = (int)enhancer.length();
   double expr_estimator = 0;
   double sum_importance_weights = 0;
   double importance_weight;
   
   /* loop over all the unique configs in the importance sample */
   for (ConfigHash::iterator i = importance_sample.begin(); 
         i != importance_sample.end(); i++) {

      /* the importance weight is the ratio of the config's energy weight 
       * under the target distribution divided by the config's energy weight
       * under the importance distribution */
      importance_weight = wc(i->first, tar_a, tar_g);
      importance_weight /= i->second.second;

      /* each term added to the sum is the product of the importance weight
       * and the probability of expression given the configuration and 
       * the number of times the configuration was seen in the importance 
       * sample */ 
      expr_estimator += importance_weight * ls->pec((*i).first, lam, L) *
         i->second.first;
//      expr_estimator += importance_weight * i->second.first;
      sum_importance_weights += importance_weight * i->second.first;
   }

   /* this estimator might be biased, but generally has a much lower error 
    * according to Liu JS. Monte Carlo Stragegies in Scientific Computing */
   return expr_estimator / sum_importance_weights;
}

/* compute lookup tables for the pwm of each TF at each position in each 
 * sequence so we don't need to do this on the fly */
void SegalModelImportanceSampler::compute_pwm_lookup_table(void) {

   /* the lookup table can be immagined as a 2D array with the number of
    * columns equal to the sequence lengths and the number of
    * rows as the number of TFs, but it's really one long vector that is
    * the matrix traversed going down columns then to the next column. */
   int M = (int)ls->M;
   for (int x=0; x < (int)enhancer.length(); x++) {
      for (int tf=0; tf < M; tf++) 
         pwm_lookup_table[x*M+tf] = ls->p_pwm(tf, x, enhancer);
   }
   return;
}

/* calculate the weight of a configuration given the sequence and params
 * but using the pwm lookup table */
double SegalModelImportanceSampler::wc(const Configuration &c, const valarray<double> &sc,
      Array &g) {
   double energy_weight = 1;
   int tf, pos, prev_tf = 0, prev_pos = 0;
   int tmp;
   int M = (int)ls->M;
   for (int i = 0; i < (int)c.csize; i++) {
      tf = c[i].first;
      pos = c[i].second-1;
      energy_weight *= tf_expr[tf] * sc[tf] * pwm_lookup_table[pos*M+tf];
      if (i!=0) {
         tmp = prev_pos-pos-ls->r[prev_tf];
         if (tmp < (int)ls->D) {
            energy_weight *= g.pos(nop, prev_tf, tf, tmp);
         }
      }
      prev_tf = tf;
      prev_pos = pos;
   }
   return energy_weight;
}


/* END */




