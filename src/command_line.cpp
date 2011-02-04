#include "command_line.h"
#include "error_handling.h"
#include "stats.h"

#include <getopt.h>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>
#include <valarray>

/* identifiers for arguments, below 255 is for ascii characters */
#define ARG_MEAN      300
#define ARG_SD        301
#define RAND_SEED     302
#define PSTAT         303
#define LSTAT         304
#define SEQ_LEN       305
#define NUM_BLOCKS    306
#define NUM_TFS       307
#define PWM           308
#define CONDITION     309
#define LAMBDA        310
#define STEPS         311
#define NEIGHBORS     312
#define FITNESS_FUNC  313
#define COEFFICIENTS  314
#define SEQUENCE      315
#define GAMMA         316
#define SAMPLING      317
#define IS_GAMMA      319
#define IS_ALPHA      320
#define SCALARS       321
#define OBS_EXPR      322
#define NUM_CHAINS    323
#define SD_LAMBDA     324
#define SD_ALPHA      325
#define SD_GAMMA      326
#define CHAIN_STEPS   327
#define CHAIN_TEMPS   328
#define DEL_PARAMS    329
#define DUP_PARAMS    330
#define LOGISTIC_PAR  331
#define TF_SEP        332

using std::cerr;
using std::endl;
using std::ostream;
using std::string;
using std::vector;
using std::valarray;

void strsplit(const string &s, vector<string> &res, char sep);
valarray<double>* valdouble_from_string(const string &);
valarray<double>* valdouble_from_string(const char *);

/* set default options */
Args::Args() {
   do_popsim = 0;
   debug_level = 0;
   model = block;
   seq_length = 0;
   rand_seed = 0;

   num_blocks = 4;
   distr_mean = 1.0;
   distr_sd = 0.02;

   num_tfs = 0;
   segal_lambda = 0;
   segal_gamma = 0;
   segal_alpha = 0;
   steps = 20000;
   fitness_coefficients = 0;
   sampling_method = iid;
   IS_segal_gamma = 0;
   IS_segal_alpha = 0;
   segal_logistic_parameter = 1;
   segal_min_tf_sep = 0;

   generations = 10;
   pop_size = 50000;
   rec_rate = 0;
   mu_rate = 2.0*pow(10.0, -8.0);
   deletion_params = 0;
   duplication_params = 0;

   max_reps = 0;
   neighbor_steps = 0;

   /* flags used to sort out which options are required in the presence 
    * of other options */
   need_fit_func = false;

   /* autofit parameters */
   autofit = 0;
   num_chains = 0;
   segal_sd_lambda = 0;
   segal_sd_gamma = 0;
   segal_sd_alpha = 0;
   segal_chain_steps = 0;
   segal_chain_temps = 0;

   /* parallelization */
   num_slots = 1;
}

Args::~Args() {
   /* currently I'm leaving this mostly un-implemented */
   if (fitness_coefficients)
      delete fitness_coefficients;
}

/* if the string begins with N, change this back to '-' */
void Args::fix_negatives(char *str) {
   if (str[0] == 'N') str[0] = '-';
   int pos;
   string s(str);
   pos = s.find("=");
   if (pos != (int)string::npos) 
      if ((int)s.length() > pos+1 && s[pos+1] == 'N') str[pos+1] = '-';
}

/* parse the argv list and populate the relevant member variables */
void Args::populate(int argc, char *argv[]) {
   int c;

   /* read in information from the environment */
   char *num_slots_str = getenv("NSLOTS");   
   char *end;
   if (num_slots_str != NULL) {
      num_slots = (int)strtol(num_slots_str, &end, 10);
      if (num_slots_str == end)
         num_slots = 1;
   }

   /* re-create the original command line as best we can */
   string sep("");
   for (int i=0; i<argc; i++) {
      cmd += sep+argv[i];
      sep = " ";
   }

   /* the following is a kluge to get negative arguments to work. Basically 
    * if -X occurs at the beginning of the string, I replace -X with NX 
    * where X is in [0-9]. Then later, before sending the optarg string on 
    * to a function that expects negative numbers I change it back to -X. */
   int pos, pos2;
   for (int i=0; i<argc; i++) {
      string s(argv[i]);
      pos = s.find_first_of("0123456789.");
      if (pos == 1 && *(argv[i]) == '-') {
         *(argv[i]) = 'N';
      } else {
         /* perhaps this arg is of the form --option=-X */
         pos = s.find("=");
         if (pos != (int)string::npos) {
            pos2 = s.find_first_of("0123456789.", pos);
            if (pos2-pos == 2 && (argv[i])[pos2-1] == '-') {
               (argv[i])[pos2-1] = 'N';
            }  
         }
      }
   }

   opterr = 1;

   while (1) {
      static struct option long_options[] = {
         /* These options set a flag. */
         {"popsim", no_argument, &do_popsim, 1},
         /* These options don't set a flag. */
         {"model",  required_argument, 0, 'm'},
         {"pstat",  required_argument, 0, PSTAT},
         {"lstat",  required_argument, 0, LSTAT},
         {"debug",  required_argument, 0, 'd'},
         {"seed",  required_argument, 0, RAND_SEED},
         {"seq",  required_argument, 0, SEQUENCE},
         /* population simulation options */
         {"generations",  required_argument, 0, 'g'},
         {"mutation",  required_argument, 0, 'u'},
         {"popsize",  required_argument, 0, 'N'},
         /* landscape sampling/traversal options */
         {"neighbors",  required_argument, 0, NEIGHBORS},
         /* options for models */
         {"mean",  required_argument, 0, ARG_MEAN},
         {"length",  required_argument, 0, SEQ_LEN},
         {"blocks",  required_argument, 0, NUM_BLOCKS},
         {"tfs",  required_argument, 0, NUM_TFS},
         {"pwm",  required_argument, 0, PWM},
         {"condition",  required_argument, 0, CONDITION},
         {"lambda",  required_argument, 0, LAMBDA},
         {"alpha",  required_argument, 0, SCALARS},
         {"steps",  required_argument, 0, STEPS},
         {"coefficients",  required_argument, 0, COEFFICIENTS},
         {"ffunc",  required_argument, 0, FITNESS_FUNC},
         {"gamma",  required_argument, 0, GAMMA},
         {"tfsep",  required_argument, 0, TF_SEP},
         {"autofit",  no_argument, &autofit, 1},
         {"sampling",  required_argument, 0, SAMPLING},
         {"IS_gamma",  required_argument, 0, IS_GAMMA},
         {"IS_alpha",  required_argument, 0, IS_ALPHA},
         {"expression",  required_argument, 0, OBS_EXPR},
         {"chains",  required_argument, 0, NUM_CHAINS},
         {"psd_lambda",  required_argument, 0, SD_LAMBDA},
         {"psd_alpha",  required_argument, 0, SD_ALPHA},
         {"psd_gamma",  required_argument, 0, SD_GAMMA},
         {"chain_steps",  required_argument, 0, CHAIN_STEPS},
         {"chain_temps",  required_argument, 0, CHAIN_TEMPS},
         {"deletion",  required_argument, 0, DEL_PARAMS},
         {"duplication",  required_argument, 0, DUP_PARAMS},
         {"logistic_param",  required_argument, 0, LOGISTIC_PAR},
         {0, 0, 0, 0}
      };
      /* getopt_long stores the option index here. */
      int option_index = 0;
  
      c = getopt_long(argc, argv, "m:d:g:u:N:", long_options, &option_index);

      char *end;
      /* Detect the end of the options. */
      if (c == -1) break;
      switch (c) {
         case '?': {
            throw SimUsageError("extraneous option");
            break;
         }
         case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0) break;

            /* every non-flag option should not have c==0 */
            throw SimUsageError("unrecognized option");
            break;
  
         case 'm':
            if (!has_option(optarg))
               throw SimUsageError("must specify model");
            model = Landscape::get_model(optarg);
            break;
  
         case 'd':
            if (!has_option(optarg))
               throw SimUsageError("must specify debug level");
            debug_level = strtol(optarg, &end, 10);
            if (optarg == end) throw SimUsageError("non-numeric debug level");
            break;
  
         case 'g':
            if (!has_option(optarg))
               throw SimUsageError("must specify number of generations");
            generations = strtoul(optarg, &end, 10);
            if (optarg == end) 
               throw SimUsageError("non-numeric number of generations");
            break;
  
         case 'N':
            if (!has_option(optarg))
               throw SimUsageError("must specify population size");
            pop_size = strtoul(optarg, &end, 10);
            if (optarg == end) 
               throw SimUsageError("non-numeric population size");
            break;
  
         case 'u':
            if (!has_option(optarg))
               throw SimUsageError("must specify mutation rate");
            mu_rate = strtod(optarg, &end);
            if (optarg == end) 
               throw SimUsageError("non-numeric mutation rate");
            break;
  
         case ARG_MEAN:
            if (!has_option(optarg))
               throw SimUsageError("must specify mean");
            fix_negatives(optarg);
            distr_mean = strtod(optarg, &end);
            if (optarg == end) 
               throw SimUsageError("non-numeric mean");
            break;
  
         case ARG_SD:
            if (!has_option(optarg))
               throw SimUsageError("must specify standard deviation");
            distr_sd = strtod(optarg, &end);
            if (optarg == end) 
               throw SimUsageError("non-numeric standard deviation");
            break;
  
         case RAND_SEED:
            if (!has_option(optarg))
               throw SimUsageError("must specify random number seed");
            rand_seed = strtoul(optarg, &end, 10);
            if (optarg == end) 
               throw SimUsageError("non-numeric random number seed");
            break;
  
         case TF_SEP:
            if (!has_option(optarg))
               throw SimUsageError("must specify minimum TF separation");
            segal_min_tf_sep = strtoul(optarg, &end, 10);
            if (optarg == end) 
               throw SimUsageError("non-numeric min TF sep");
            break;
  
         case PSTAT: {
            if (!has_option(optarg))
               throw SimUsageError("must specify statistic");
            PopStat pstat(optarg);
            /* we only need to queue queued stats */
            if (pstat.stype == queued) 
               popstats.push_back(pstat);
            break;
         }
         case LSTAT: {
            if (!has_option(optarg))
               throw SimUsageError("must specify statistic");
            ScapeStat lstat(optarg);
            if (lstat.name.compare("fitness") == 0) need_fit_func = true;
            scapestats.push_back(lstat);
            break;
         }
         case NUM_BLOCKS:
            if (!has_option(optarg))
               throw SimUsageError("must specify number of blocks");
            num_blocks = strtoul(optarg, &end, 10);
            if (optarg == end) 
               throw SimUsageError("number of blocks must be non-negative int");
            break;

         case SEQ_LEN: {
            if (!has_option(optarg))
               throw SimUsageError("must specify sequence length");
            unsigned int len_tmp = strtoul(optarg, &end, 10);
            if (seq_length > 0 && seq_length != len_tmp) 
               throw SimUsageError("seq len and len param don't match");
            seq_length = len_tmp;
            if (optarg == end) 
               throw SimUsageError("sequence length must be non-negative int");
            break;
         }
         case NUM_TFS:
            if (!has_option(optarg))
               throw SimUsageError("must specify number of tfs");
            num_tfs = strtoul(optarg, &end, 10);
            if (optarg == end) 
               throw SimUsageError("number of tfs must be non-negative int");
            break;

         case PWM:
            if (!has_option(optarg))
               throw SimUsageError("must specify pwm filename");
            pwm_filename = string(optarg);
            break;

         case LAMBDA:
            if (!has_option(optarg))
               throw SimUsageError("must specify lambda");
            fix_negatives(optarg);
            segal_lambda = valdouble_from_string(optarg);
            break;

         case SCALARS:
            if (!has_option(optarg))
               throw SimUsageError("must specify TF expr. scalar vec. alpha");
            fix_negatives(optarg);
            segal_alpha = valdouble_from_string(optarg);
            break;

         case CONDITION:
            if (!has_option(optarg))
               throw SimUsageError("must specify condition");
            fix_negatives(optarg);
            segal_conditions.push_back(valdouble_from_string(optarg));
            break;

         case OBS_EXPR:
            if (!has_option(optarg))
               throw SimUsageError("must specify observed enhancer expression");
            fix_negatives(optarg);
            segal_obs_expr.push_back(valdouble_from_string(optarg));
            break;

         case STEPS:
            if (!has_option(optarg))
               throw SimUsageError("must specify number steps");
            steps = strtoul(optarg, &end, 10);
            if (optarg == end) 
               throw SimUsageError("number of steps must be non-negative int");
            break;

         case NEIGHBORS:
            if (!has_option(optarg))
               throw SimUsageError("need steps in neighborhood");
            neighbor_steps = strtoul(optarg, &end, 10);
            if (optarg == end) 
               throw SimUsageError("number of steps must be non-negative int");
            break;

         case COEFFICIENTS:
            if (!has_option(optarg))
               throw SimUsageError("must specify fitness function coef.");
            fix_negatives(optarg);
            fitness_coefficients = valdouble_from_string(optarg);
            break;

         case FITNESS_FUNC: {
            if (!has_option(optarg))
               throw SimUsageError("must specify fitness function choice");
            string s(optarg);  
            segal_expr_fit_func = FitFunc(s);
            break;
         }
         case SEQUENCE: {
            if (!has_option(optarg))
               throw SimUsageError("must specify sequence");
            string s(optarg);
            sequences.push_back(s);
            if (seq_length > 0 && seq_length != s.length()) 
               throw SimUsageError("seq len and len param don't match");
            if (seq_length == 0) seq_length = s.length();
            break;
         }
         case GAMMA:
            if (!has_option(optarg))
               throw SimUsageError("must specify cooperativity matrix");
            fix_negatives(optarg);
            segal_gamma = valdouble_from_string(optarg);
            break;
         case IS_GAMMA:
            if (!has_option(optarg))
               throw SimUsageError("must specify IS cooperativity matrix");
            fix_negatives(optarg);
            IS_segal_gamma = valdouble_from_string(optarg);
            break;
         case IS_ALPHA:
            if (!has_option(optarg))
               throw SimUsageError("must specify IS expr. scalar vec. alpha");
            fix_negatives(optarg);
            IS_segal_alpha = valdouble_from_string(optarg);
            cerr << "setting IS_segal_expr_scalers. len=" << IS_segal_alpha->size() << endl;
            break;
         case SAMPLING: {
            if (!has_option(optarg))
               throw SimUsageError("must specify sampling method");
            string s(optarg);
            if (s == "iid") {
               sampling_method = iid;
            } else if (s == "importance") {
               sampling_method = importance;
            } else {
               throw SimUsageError("invalid sampling type");
            }
            break;
         }
         case NUM_CHAINS:
            if (!has_option(optarg))
               throw SimUsageError("must specify number of MCMC chains");
            num_chains = strtol(optarg, &end, 10);
            if (optarg == end || num_chains <= 0) 
               throw SimUsageError("number of chains must be positive int");
            break;
         case CHAIN_STEPS:
            if (!has_option(optarg))
               throw SimUsageError("must specify number of chain steps");
            segal_chain_steps = strtol(optarg, &end, 10);
            if (optarg == end || segal_chain_steps <= 0) 
               throw SimUsageError("num. of chain steps must be positive int");
            break;
         case CHAIN_TEMPS:
            if (!has_option(optarg))
               throw SimUsageError("must specify chain temperatures");
            fix_negatives(optarg);
            segal_chain_temps = valdouble_from_string(optarg);
            break;
         case SD_LAMBDA:
            if (!has_option(optarg))
               throw SimUsageError("must specify s.d. of lambda proposal");
            segal_sd_lambda = strtod(optarg, &end);
            if (optarg == end || segal_sd_lambda <= 0) 
               throw SimUsageError("lambda proposal s.d. must be positive");
            break;
         case SD_ALPHA:
            if (!has_option(optarg))
               throw SimUsageError("must specify s.d. of alpha proposal");
            segal_sd_alpha = strtod(optarg, &end);
            if (optarg == end || segal_sd_alpha <= 0) 
               throw SimUsageError("alpha proposal s.d. must be positive");
            break;
         case SD_GAMMA:
            if (!has_option(optarg))
               throw SimUsageError("must specify s.d. of gamma proposal");
            segal_sd_gamma = strtod(optarg, &end);
            if (optarg == end || segal_sd_gamma <= 0) 
               throw SimUsageError("gamma proposal s.d. must be positive");
            break;
         case DEL_PARAMS:
            if (!has_option(optarg))
               throw SimUsageError("must specify deletion parameters");
            fix_negatives(optarg);
            deletion_params = valdouble_from_string(optarg);
            break;
         case DUP_PARAMS:
            if (!has_option(optarg))
               throw SimUsageError("must specify duplication parameters");
            fix_negatives(optarg);
            duplication_params = valdouble_from_string(optarg);
            break;
         case LOGISTIC_PAR:
            if (!has_option(optarg))
               throw SimUsageError("must specify segal logistic parameter");
            segal_logistic_parameter = strtoul(optarg, &end, 10);
            if (optarg == end || segal_logistic_parameter <= 0) 
               throw SimUsageError("segal logistic parameter must be non-negative int");
            break;

         default:
            char o[10];
            snprintf(o, 10, "%d", c);
            string s(o);
            throw SimUsageError("invalid command option" + s);
      }
   }

   /* throw exception if there are any extra options */
   if (optind < argc) {
      throw SimUsageError("unrecognized additional command options");
   }

   /* detect any conflicts among options specified */

   if (!autofit) {
      if ((do_popsim && popstats.size() == 0) ||
            (!do_popsim && scapestats.size() == 0))
         throw SimUsageError("must specify one or more statistics to print");
   }
   if (do_popsim) {
      need_fit_func = true;
      if (deletion_params != 0 && (int)deletion_params->size() != 3)
         throw SimUsageError("must have 3 deletion parameters");
      if (duplication_params != 0 && (int)duplication_params->size() != 3)
         throw SimUsageError("must have 3 duplication parameters");
   }

   if (model == block && (seq_length % num_blocks) != 0) 
      throw SimUsageError("sequence length must be multiple of num_blocks");

   /* search the scape stats to find the maximum number of reps */
   for (vector<ScapeStat>::const_iterator i = scapestats.begin();
         i != scapestats.end(); i++) 
      if (i->reps > max_reps) max_reps = i->reps;

   if (model == segal) {
      if (pwm_filename.length() == 0) 
         throw SimUsageError("pwm file must be specified");
      if (num_tfs == 0) 
         throw SimUsageError("must specify number of tfs for Segal model");
      if (segal_min_tf_sep < 0) 
         throw SimUsageError("minimum TF separation must be non-negative");
      if (steps == 0) throw SimUsageError("steps must be > 0");

      /* make sure lambda and conditions are all the right length */
      if (segal_conditions.size() == 0) 
         throw SimUsageError("conditions must be specified");
      for (vector< valarray<double>* >::iterator i = segal_conditions.begin();
            i != segal_conditions.end(); i++) {
         if ((*i)->size() != num_tfs)
            throw SimUsageError("must be one TF level per TF for condition");
      }

      if (autofit && do_popsim)
         throw SimUsageError("we can only autofit when not doing a popsim");

      /* specify parameters only if not autofitting */
      if (autofit) {
         if (segal_lambda != 0 || segal_alpha != 0 || segal_gamma != 0 )
            throw SimUsageError("segal parameters given but told to autofit");
         if (sampling_method != importance)
            throw SimUsageError("autofit only implemented for IS");
         if (segal_obs_expr.size() != sequences.size()) 
            throw SimUsageError("number of expression profiles must match "
               "number of enhancers");
         if (num_chains <= 0) 
            throw SimUsageError("number of chains must be positive");
         if (segal_chain_steps <= 0) 
            throw SimUsageError("number of chain steps must be positive");
         if (segal_sd_lambda <= 0)
            throw SimUsageError("need to specify lambda proposal s.d.");
         if (segal_sd_alpha <= 0)
            throw SimUsageError("need to specify alpha proposal s.d.");
         if (segal_sd_gamma <= 0)
            throw SimUsageError("need to specify gamma proposal s.d.");
         if (segal_chain_temps == 0 || 
               (int)segal_chain_temps->size() != num_chains)
            throw SimUsageError("chain temps must be num_chains in length");
/*         if (segal_alpha_lattice == 0)
            throw SimUsageError("must specify alpha lattice");
*/
      } else {
         /* if we're not autofitting, we need to specify the params */
         if (segal_lambda == 0 || segal_lambda->size() != num_tfs+1)
            throw SimUsageError("must be 1+num_tfs lambda");
         if (segal_alpha != 0 && segal_alpha->size() != num_tfs)
            throw SimUsageError("segal alpha must have num_tfs values");
      }
      if (segal_gamma != 0 && segal_gamma->size() != num_tfs*num_tfs)
         throw SimUsageError("cooperativity matrix must be num_tfs^2 entries");

      /* perhaps we need a fitness function */
      if (need_fit_func) {
         if (fitness_coefficients == 0) 
            throw SimUsageError("must specify fitness function coef.");
         if (segal_expr_fit_func.calc_var == 0)
            throw SimUsageError("must specify fitness function choice");
         segal_expr_fit_func.set_coef(*fitness_coefficients);
      }

      /* there are some parameters contingent on the sampling method, but 
       * currently we only specify them if we're not autofitting */
      if (!autofit && sampling_method == importance) {
         if (IS_segal_gamma == 0)
            throw SimUsageError("must specify IS cooperativity matrix, gamma");
         if (IS_segal_alpha == 0)
            throw SimUsageError("must specify IS alpha");
         if (IS_segal_gamma->size() != num_tfs*num_tfs)
            throw SimUsageError("IS coop matrix must be num_tfs^2 in length");
         if (IS_segal_alpha->size() != num_tfs)
            throw SimUsageError("IS alpha must be num_tfs in length");
      }
      
   } else {
      if (segal_conditions.size() > 0) 
         throw SimUsageError("conditions only for segal model");
   }
   
   if (model == neutral) {
      if (seq_length == 0) 
         throw SimUsageError("sequence length must be greater than zero");
   }
   return;
}

/* make sure optarg doesn't point to the next option */
bool Args::has_option(char *optarg) {
   return !((optarg == 0) || (strlen(optarg) > 1 && optarg[0] == '-'));
}

/* print the arguments to the screen */
ostream& operator<<(ostream &s, const Args &a) {
   s << "cmd: " << a.cmd << endl;
   s << "params: debug=" << a.debug_level << " model=" << a.model 
      << " popsim=" << a.do_popsim << " pop_size=" << a.pop_size 
      << " seq_length=" << a.seq_length << " mu_rate=" << a.mu_rate 
      << " recrate=" << a.rec_rate << " generations=" << a.generations 
      << " distr_mean=" << a.distr_mean << " distr_sd=" << a.distr_sd 
      << " seed=" << a.rand_seed << " max_reps=" << a.max_reps
      << " num_blocks=" << a.num_blocks << " num_tfs=" << a.num_tfs
      << " pwm_filename=" << a.pwm_filename 
      << " logistic_param=" << a.segal_logistic_parameter;

   /* print lambda */
   if (a.segal_lambda == 0) {
      s << " segal_lambda=none";
   } else {
      s << " segal_lambda=";
      for (unsigned int i=0; i<a.segal_lambda->size(); i++) {
         if (i!=0) s << ",";
         s << (*a.segal_lambda)[i];
      }
   }

   /* print cooperativity matrix */
   if (a.segal_gamma == 0) {
      s << " segal_gamma=none";
   } else {
      s << " segal_gamma=";
      for (unsigned int i=0; i<a.segal_gamma->size(); i++) {
         if (i!=0) s << ",";
         s << (*a.segal_gamma)[i];
      }
   }

   /* print IS cooperativity matrix */
   if (a.IS_segal_gamma == 0) {
      s << " IS_segal_gamma=none";
   } else {
      s << " IS_segal_gamma=";
      for (unsigned int i=0; i<a.IS_segal_gamma->size(); i++) {
         if (i!=0) s << ",";
         s << (*a.IS_segal_gamma)[i];
      }
   }

   /* print IS alpha */
   if (a.IS_segal_alpha == 0) {
      s << " IS_segal_alpha=none";
   } else {
      s << " IS_segal_alpha=";
      for (unsigned int i=0; i<a.IS_segal_alpha->size(); i++) {
         if (i!=0) s << ",";
         s << (*a.IS_segal_alpha)[i];
      }
   }

   /* print conditions arguments */
   if (a.segal_conditions.size() == 0) {
      s << " segal_conditions=none";
   } else {
      s << " segal_conditions=";
      string sep("");
      for (vector< valarray<double>* >::const_iterator i = 
            a.segal_conditions.begin(); i != a.segal_conditions.end(); i++) {
         s << sep+"(";
         sep = ",";
         for (unsigned int j=0; j < (*i)->size(); j++) {
            if (j!=0) s << ",";
            s << (**i)[j];
         }
         s << ")";
      }
   }

   /* print deletion parameters */
   if (a.deletion_params == 0) {
      s << " deletion_params=none";
   } else {
      s << " deletion_params=";
      for (unsigned int i=0; i<a.deletion_params->size(); i++) {
         if (i!=0) s << ",";
         s << (*a.deletion_params)[i];
      }
   }

   /* print duplication parameters */
   if (a.duplication_params == 0) {
      s << " duplication_params=none";
   } else {
      s << " duplication_params=";
      for (unsigned int i=0; i<a.duplication_params->size(); i++) {
         if (i!=0) s << ",";
         s << (*a.duplication_params)[i];
      }
   }

   return s;
}

/* split a string by sep filling the vector */
void strsplit(const string &s, vector<string> &res, char sep) {
   int newpos = 0, oldpos = 0;
   while (1) {
      newpos = s.find(sep, oldpos);
      if (newpos == (int)string::npos) break;
      if (newpos-oldpos <= 0) throw SimError("nothing between separators");
      res.push_back(s.substr(oldpos, newpos-oldpos));
      oldpos = newpos+1;
   }
   res.push_back(s.substr(oldpos, s.length()-oldpos));
}

/* allocate a valarray for doubles, filling it from the comma-separated str */
valarray<double>* valdouble_from_string(const string &s) {
   vector<string> pieces;
   strsplit(s, pieces, ',');
   if (pieces.size() == 0) throw SimError("error parsing comma-separated list");
   valarray<double> *ret = new valarray<double>(pieces.size());
   unsigned int k = 0;
   for (vector<string>::iterator i = pieces.begin(); i != pieces.end(); i++) {
      char *end;
      (*ret)[k++] = strtod(i->c_str(), &end);
      if (i->c_str() == end) throw SimError("error in comma-separated list");
   }
   return ret;
}

/* same as above but for char* */
valarray<double>* valdouble_from_string(const char *s) {
   string tmp(s);
   return valdouble_from_string(tmp);
}

/* END */



