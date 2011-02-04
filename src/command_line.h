#ifndef __COMMAND_LINE_H
#define __COMMAND_LINE_H

#include "landscape.h"
#include "sequence.h"

#include <vector>

class Args {
public:
   /* public member functions */
   Args();
   ~Args();
   void populate(int, char**);
   friend std::ostream& operator<<(std::ostream &s, const Args &a);
   
   /* general public member variables */
   int do_popsim;
   int debug_level;
   Model model;
   unsigned int seq_length;
   unsigned int rand_seed;
   std::string cmd;
   std::vector<PopStat> popstats;
   std::vector<ScapeStat> scapestats;
   std::vector<std::string> sequences;

   /* block model-related */
   unsigned int num_blocks;
   double distr_mean;
   double distr_sd;

   /* segal model-related */
   unsigned int num_tfs;
   std::string pwm_filename;
   std::valarray<double> *segal_lambda; /* ability to recruit POL */
   std::valarray<double> *segal_alpha; /* TF-specific scalar for TF levels */
   std::vector< std::valarray<double>* > segal_conditions; /* TF levels */
   std::vector< std::valarray<double>* > segal_obs_expr; /* observed expr */
   unsigned int steps;
   std::valarray<double> *fitness_coefficients;
   FitFunc segal_expr_fit_func;
   std::valarray<double> *segal_gamma; /* cooperativity matrix */
   SamplingMethod sampling_method; /* either iid or importance */
   double segal_logistic_parameter; /* d in 1/1+exp(-x/d), controlling the transition */
   int segal_min_tf_sep; /* minimum separation between TFs in the same configuration */

   /* segal importance sampling-related */
   std::valarray<double> *IS_segal_gamma;
   std::valarray<double> *IS_segal_alpha;

   /* popsim-related */
   unsigned int generations;
   unsigned int pop_size;
   double rec_rate;
   double mu_rate;
   std::valarray<double> *deletion_params;
   std::valarray<double> *duplication_params;

   /* landscape characterization-related */
   unsigned int max_reps;
   unsigned int neighbor_steps;

   /* autofit related */
   int autofit;
   int num_chains;
   double segal_sd_alpha;
   double segal_sd_lambda;
   double segal_sd_gamma;
   int segal_chain_steps;
   std::valarray<double> *segal_chain_temps;

   /* parallelization */
   int num_slots;

private:
   /* private functions */
   bool has_option(char *);
   void fix_negatives(char *str);

   /* flags used to figure out which options we need */
   bool need_fit_func;
};

#endif /* __COMMAND_LINE_H */
