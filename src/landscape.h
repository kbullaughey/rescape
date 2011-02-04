#ifndef __LANDSCAPE_H__
#define __LANDSCAPE_H__

#include "scapetypes.h"
#include "stats.h"
#include "sim_rand.h"
#include "array.h"
#include "sequence.h"
#include "configuration.h"
#include "chain.h"

#include <ext/hash_map>
#include <vector>
#include <string>
#include <math.h>
#include <valarray>
#include <utility>
#include <set>
#include <map>

#define DEFAULT_HASH_SIZE 1000000

typedef __gnu_cxx::hash_map<const char*, double, 
   __gnu_cxx::hash<const char*>, eqstr> StringHash; 

/* Valid landscape models */
enum Model {block, segal, neutral, invalid};

/* Valid sampling methods */
enum SamplingMethod {iid, importance};

/*******************************************************************************
 * genetic sampler class to share functions between IID and IS
 ******************************************************************************/

class Sampler {

public:
   Sampler(RandomNumberGenerator *rngp);
   ~Sampler() { }
   int sample_cdf(int tbl_size, double tbl_sum, std::valarray<double> &tbl,
      int offset);

private:
   RandomNumberGenerator *rng;
};

/******************************************************************************/
/* Abstract class Landscape */
/******************************************************************************/
class Landscape {

public: 
   Landscape(RandomNumberGenerator *rngp);
   virtual ~Landscape() { }
   virtual double get_fitness(const Sequence&) = 0;
   static Model get_model(char *m);
   char pick_mutation(char off);
   Sequence sample(unsigned int L);
   void insert_neighbors(std::set<Sequence> &s); 

   /* public member variables */
   std::string alphabet;
   unsigned int cur_rep; /* number samples so far from landscape */
   RandomNumberGenerator *rng;
};

/******************************************************************************/
/* Landscape neutral model class */
/******************************************************************************/
class LandscapeNeutralModel : public Landscape {
public:
   LandscapeNeutralModel(RandomNumberGenerator *rngp);
   ~LandscapeNeutralModel() { }
   double get_fitness(const Sequence&);
};

/******************************************************************************/
/* Landscape block model class */
/******************************************************************************/
class LandscapeBlockModel : public Landscape {

public:
   LandscapeBlockModel(unsigned int B, unsigned int bw, double mean, double sd, 
      RandomNumberGenerator *rng, unsigned int hash_size = DEFAULT_HASH_SIZE);
   ~LandscapeBlockModel();
   double get_fitness(const Sequence&);

private:
   unsigned int blocks; /* number of blocks in this landscape */
   unsigned int blk_width; /* the number of nucleotides per block */
   std::vector<StringHash*> caches;
   double mu; /* the mean of a block's Gaussian distribution */
   double sigma; /* the standard deviation of a block's Gaussian distribution */
   unsigned int L;

   /* private function */
   double sample_fitness(void) { return fmax(0, rng->rnorm(mu, sigma)); }
};

/******************************************************************************/
/* class for representing PWMs */
/******************************************************************************/
class PWM {
public:
   /* public member functions */
   PWM(std::string &);
   ~PWM() { }
   friend std::istream& operator>>(std::istream &, PWM &);
   friend std::ostream& operator<<(std::ostream &, PWM &);

   /* public member variables */
   std::valarray<double> *pwm;
   unsigned int w, h;
   std::string alphabet;
};

/******************************************************************************/
/* class for representing fitness functions */
/******************************************************************************/
typedef double (*FitCalc)(std::valarray<double> &, std::valarray<double> &);
class FitFunc {
public:
   /* public member function */
   FitFunc(std::string &choice);
   FitFunc();
   ~FitFunc() { /*if (c!=0) delete c;*/ }
   void set_coef(std::valarray<double> &c);
   double calc(std::valarray<double> &x) { return (*calc_var)(x, *c); }
   double (*calc_var)(std::valarray<double> &x, std::valarray<double> &c);

private:
   static std::map<std::string, FitCalc> funclist;
   std::valarray<double> *c; /* vector of coefficients for the function */
};

/******************************************************************************/
/* Landscape on Segal's model */
/******************************************************************************/

/* class to store probability tables in for the segal model */
class ProbTable {
public:
   ProbTable(unsigned int M, unsigned int D, unsigned int L);
   ~ProbTable() { }
   Array cpa;
   Array cpt;
   Array cpr;
   Array cpm;
};

class LandscapeSegalModel : public Landscape {

public:
   /* public member functions */
   LandscapeSegalModel(unsigned int num_tfs, 
      std::string &pwm_file, std::vector< std::valarray<double>* > &conds, 
      unsigned int d_max, RandomNumberGenerator *rngp, double logistic_par, int min_tf_sep);
   ~LandscapeSegalModel();
   virtual void get_expression(const Sequence &seq, std::valarray<double> &pe) = 0;
   double pec(const Configuration &c, valarray<double> &s, int L);
   double p_pwm(unsigned int tf, int pos, const Sequence &seq);
   void calc_prob_tables(std::valarray<double> *rel_tf_expr, 
      std::valarray<double> &alpha, Array &gm, const Sequence &seq, 
      ProbTable &ret);
   void set_gamma(valarray<double> *gm, Array &local_gamma);

   /* public member variables */
   unsigned int M; /* number of TFs in the model */
   unsigned int D; /* maximum distance between TFs to include gamma term */
   std::valarray<int> r; /* width of each TF binding function */
   std::vector< std::valarray<double>* > conditions; /* TF levels */
   double logistic_parameter; /* d in 1/(1+expr(-x/d)) */
   int minimum_tf_separation; /* the minimum dist. between TFs in the same config */

protected:
   /* private member functions */
   void load_pwms(const std::string &fn);
   double coopfunc(double g, int distance);

   /* private member variables */
   std::vector<PWM*> pwms; /* each entry is a pwm matrix */
   std::valarray<double> lambda; /* lambda of TFs (ability to recruit POL) */
   Array gamma; /* MxMxD matrix specifying cooperativity parameters */
   std::valarray<double> alpha; /* vector to rescale expression levels */
};

/*******************************************************************************
 * Fitted model class
 ******************************************************************************/

class LandscapeSegalModelFitted : public LandscapeSegalModel {
public:
   /* public member functions */
   LandscapeSegalModelFitted(unsigned int num_tfs, std::string &pwm_file, 
      std::vector< std::valarray<double>* > &conds, unsigned int d_max, 
      unsigned int num_steps, RandomNumberGenerator *rngp,
      std::vector< std::valarray<double>* > &obs_expr,
      double sdlambda, double sdalpha, double sdgamma,
      std::vector<std::string> &seqs, int nslots, double logistic_par, int min_tf_sep);
   ~LandscapeSegalModelFitted();
   void get_expression(const Sequence &seq, std::valarray<double> &pe) {}
   double get_fitness(const Sequence&);
   void create_chains(int n, std::valarray<double> &chain_temps);
   void precompute_lookup_tables(std::vector<std::string> &seqs, 
      unsigned int steps);

   /* public member variables */
   ParallelTemperedChains parallel_tempered;
   double sd_lambda;
   double sd_alpha;
   double sd_gamma;
   int slots;

private:
   /* private member variables */
   std::vector< std::valarray<double>* > obs_expr_profiles;
};

/*******************************************************************************
 * Specified model class
 ******************************************************************************/

class LandscapeSegalModelSpecified : public LandscapeSegalModel {

public:
   /* public member functions */
   LandscapeSegalModelSpecified(unsigned int num_tfs,
      std::string &pwm_file, std::valarray<double> *li, 
      std::valarray<double> *expr_s,
      std::vector< std::valarray<double>* > &conds, 
      std::valarray<double> *gm, FitFunc &fitfunc, unsigned int d_max,
      RandomNumberGenerator *rngp, double logistic_par, int min_tf_sep);
   ~LandscapeSegalModelSpecified();
   double get_fitness(const Sequence&);
   void get_expression(const Sequence &seq, std::valarray<double> &pe);
   void get_tf_occupancy(const Sequence &seq, int cond, Array &oc);
   virtual int nsteps(void) = 0;
   virtual int sample_cdf_wrapper(int tbl_size, double tbl_sum, 
      std::valarray<double> &tbl, int offset) = 0;
   void calc_prob_tables_wrapper(unsigned int condition, const Sequence &seq,
      ProbTable &ret);

private:
   /* private member functions */
   virtual double prob_expr(ProbTable &ptables, const Sequence &seq, int condition,
      bool print_configs=false) = 0;

   /* private member variables */
   FitFunc eff; /* expression to fitness function */
};

/*******************************************************************************
 * IS class
 ******************************************************************************/

class SegalModelIS : public Sampler {

public:
   SegalModelIS(unsigned int num_steps, RandomNumberGenerator *rngp);
   ~SegalModelIS() {}
   void generate_importance_sample(int L, int M, int D, 
      std::valarray<int> &r, ProbTable &ptables, ConfigHash &importance_sample);

protected:
   unsigned int steps; /* number of samples to draw when estimating expr.*/
};

/*******************************************************************************
 * IID class
 ******************************************************************************/

class SegalModelIID : public Sampler {

public:
   SegalModelIID(unsigned int num_steps, RandomNumberGenerator *rngp);
   ~SegalModelIID() {}

protected:
   double prob_expr(int L, int M, int D, std::valarray<int> &r, 
      std::valarray<double> &lambda, ProbTable &ptables, 
      double logistic_parameter, bool print_configs=false);
   unsigned int steps; /* number of iid samples to draw when estimating expr.*/
};

/******************************************************************************
 * Full Segal model Specified IS class
 ******************************************************************************/

class LandscapeSegalModelSpecifiedIS : public LandscapeSegalModelSpecified, 
   public SegalModelIS {

public:
   /* public member functions */
   LandscapeSegalModelSpecifiedIS(unsigned int num_tfs,
      std::string &pwm_file, std::valarray<double> *li, 
      std::valarray<double> *expr_sc,
      std::vector< std::valarray<double>* > &conds, 
      std::valarray<double> *gm, unsigned int num_steps,
      FitFunc &fitfunc, unsigned int d_max, RandomNumberGenerator *rngp,
      std::valarray<double> *tar_gm, std::valarray<double> *tar_a,
      double logistic_par, int min_tf_sep);
   ~LandscapeSegalModelSpecifiedIS();
   int nsteps(void) { return steps; }
   int sample_cdf_wrapper(int tbl_size, double tbl_sum, 
         std::valarray<double> &tbl, int offset) { 
      return sample_cdf(tbl_size, tbl_sum, tbl, offset); }

private:
   /* private member functions */
   double prob_expr(ProbTable &ptables, const Sequence &seq, int condition, 
      bool print_configs=false);
   double pec_importance(Configuration &c);
   void calc_importance_sample_energy_weights(ConfigHash &importance_sample, 
      const Sequence &seq, std::valarray<double> &tf_expr);
   double wc(const Configuration &c, const std::valarray<double> &sc,
      const Sequence &seq, Array &g, const std::valarray<double> &expr);

   /* private member variables */
   std::valarray<double> target_alpha; /* vector to rescale expression levels */
   Array target_gamma; /* cooperativity parameters for target distr */
   int num_importance_samples; /* total number of importance samples drawn */
};

/*******************************************************************************
 * Full Segal model Specified IID class
 *******************************************************************************/

class LandscapeSegalModelSpecifiedIID : public LandscapeSegalModelSpecified, 
      protected SegalModelIID {

public:
   /* public member functions */
   LandscapeSegalModelSpecifiedIID(unsigned int num_tfs,
      std::string &pwm_file, std::valarray<double> *li, 
      std::valarray<double> *expr_sc,
      std::vector< std::valarray<double>* > &conds, 
      std::valarray<double> *gm, unsigned int num_steps,
      FitFunc &fitfunc, unsigned int d_max, RandomNumberGenerator *rngp,
      double logistic_par, int min_tf_sep);
   ~LandscapeSegalModelSpecifiedIID();
   int nsteps(void) { return steps; }
   int sample_cdf_wrapper(int tbl_size, double tbl_sum, 
         std::valarray<double> &tbl, int offset) { 
      return sample_cdf(tbl_size, tbl_sum, tbl, offset); }
   void get_wsum(const Sequence &seq, std::valarray<double> &ws);
   double prob_expr(ProbTable &ptables, const Sequence &seq, int condition, 
      bool print_configs=false);
};

/*******************************************************************************
 * Class to further flush out importance sampling functionality 
 *******************************************************************************/

class SegalModelImportanceSampler : public SegalModelIS {
public:
   SegalModelImportanceSampler(Sequence &seq, SegalModelParameter &smp, 
      std::valarray<double> &condition, double obs_expr,
      LandscapeSegalModelFitted *lsp, int num_steps);
   ~SegalModelImportanceSampler();
   void calc_importance_sample_energy_weights(void);
   double importance_prob_expr(std::valarray<double> &lam, 
      std::valarray<double> &tar_a, Array &tar_g);
   void compute_pwm_lookup_table(void);
   double wc(const Configuration &c, const std::valarray<double> &sc, Array &g);

   Sequence enhancer;
   SegalModelParameter theta;
   std::valarray<double> tf_expr;
   double observed_expression;
   ProbTable ptables;
   ConfigHash importance_sample;
   LandscapeSegalModelFitted *ls;
   std::valarray<double> pwm_lookup_table;
};

#endif /* __LANDSCAPE_H__ */





