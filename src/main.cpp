/* my own headers */
#include "landscape.h"
#include "error_handling.h"
#include "command_line.h"
#include "popsim.h"
#include "stats.h"
#include "sequence.h"

/* system headers */
#include <iostream>
#include <string>

using namespace std;

void usage(void);

int main(int argc, char *argv[]) { try {
   Args ar;
   if (argc == 1) {
      usage();
      return 0;
   }

   /* read in the command line arguments and print them out */
   ar.populate(argc, argv);
   cout << "args: " << ar << endl;

   /* create the random number generator */
   RandomNumberGenerator *rng = 0;
   rng = new RandomNumberGenerator(ar.rand_seed);

   /* create the landscape */
   Landscape *l = 0;
   switch (ar.model) {
      case block:
         l = new LandscapeBlockModel(ar.num_blocks, 
            ar.seq_length/ar.num_blocks, ar.distr_mean, ar.distr_sd, rng);
         break;
      case segal:
         if (ar.autofit) {
            cerr << "autofitting" << endl;
            l = new LandscapeSegalModelFitted(ar.num_tfs, ar.pwm_filename,
               ar.segal_conditions, (unsigned int)19, ar.steps, rng, 
               ar.segal_obs_expr, ar.segal_sd_lambda, ar.segal_sd_alpha, 
               ar.segal_sd_gamma, ar.sequences, ar.num_slots, ar.segal_logistic_parameter, ar.segal_min_tf_sep);
         } else {
            if (ar.sampling_method == importance) {
               cerr << "using importance sampling method" << endl;
               /* we use the IS_* parameters as the importance distribution 
                * from which we'll sample. And then we pass the segal_* parameters
                * as the target distribution parameters that is then used to 
                * calculate the importance weights */
               l = new LandscapeSegalModelSpecifiedIS(ar.num_tfs, ar.pwm_filename,
                  ar.segal_lambda, ar.IS_segal_alpha, ar.segal_conditions,
                  ar.IS_segal_gamma, ar.steps, ar.segal_expr_fit_func, 19, rng, 
                  ar.segal_gamma, ar.segal_alpha, ar.segal_logistic_parameter, ar.segal_min_tf_sep);
            } else { /* regular iid sampling */
               l = new LandscapeSegalModelSpecifiedIID(ar.num_tfs, ar.pwm_filename,
                  ar.segal_lambda, ar.segal_alpha, ar.segal_conditions, 
                  ar.segal_gamma, ar.steps, ar.segal_expr_fit_func, 19, rng,
                  ar.segal_logistic_parameter, ar.segal_min_tf_sep);
            }
         }
         break;
      case neutral:
         l = new LandscapeNeutralModel(rng);
         break;
      default:
         throw SimError("invalid landscape type");
   }

   if (ar.autofit) {
      /* We're trying to fit fit the model to the provided expression data */

      LandscapeSegalModelFitted *lpt = (LandscapeSegalModelFitted*)l;
      lpt->create_chains(ar.num_chains, *ar.segal_chain_temps);
      cerr << "auto fit created " << lpt->parallel_tempered.num_chains() << 
         " chains" << endl;
      lpt->precompute_lookup_tables(ar.sequences, ar.steps);
      cerr << "initing energies" << endl;
      lpt->parallel_tempered.init_energies();
      cerr << "running chains" << endl;
      lpt->parallel_tempered.run_chains(ar.segal_chain_steps);
   
   } else if (ar.do_popsim) {
      /* we're doing a population simulation */

      PopState ps(ar.pop_size, ar.mu_rate, ar.rec_rate, l, ar.deletion_params,
         ar.duplication_params);

      /* give the popstats access to the popstate and set stats marked print
       * at end, to print at the end. */
      for (vector<PopStat>::iterator i = ar.popstats.begin(); 
            i != ar.popstats.end(); i++) {
         i->set_popstate(&ps);
         if (i->modulo == PRINT_STAT_AT_END)
            i->modulo = ar.generations;
      }

      /* create an allele unless we have ones from the command line */
      if (ar.sequences.size() < 1) {
         string empty(ar.seq_length, l->alphabet[0]);
         Sequence s(empty, l->alphabet);
         l->rng->uniform_dna(s, ar.seq_length, l);
         ps.alleles[s] = new Allele(s, ps.N, ps.generation, l);
      }
      for (vector<string>::iterator i = ar.sequences.begin(); 
            i != ar.sequences.end(); i++) {
         Sequence s(*i, l->alphabet);
         ps.alleles[s] = new Allele(s, ps.N/ar.sequences.size(), 
            ps.generation, l);
      }

      /* run the simulation */
      ps.run(ar.popstats, ar.generations);

   } else {
      /* characterization of landscape part */

      /* make sure we're not asked for popstats if not doing a popsim */
      if (ar.popstats.size() > 0) {
         cerr << "pop stats require --popsim option" << endl;
         usage();
      }

      /* give the scapetats access to the landscape */
      for (vector<ScapeStat>::iterator i = ar.scapestats.begin(); 
            i != ar.scapestats.end(); i++) 
         i->set_landscape(l);

      /* loop, sampling nodes in the landscape uniformly until we do enough */
      set<Sequence> draws;
      for (vector<string>::iterator i = ar.sequences.begin(); 
            i != ar.sequences.end(); i++) {
         Sequence s(*i, l->alphabet);
         draws.insert(s);
      }
      sort(ar.scapestats.rbegin(), ar.scapestats.rend());
      for (l->cur_rep = 0; l->cur_rep < ar.max_reps; l->cur_rep++) {
         if (draws.size() == 0) {
            Sequence draw = l->sample(ar.seq_length);
            draws.insert(draw);
         }
         unsigned int ns = ar.neighbor_steps;
         while (ns-- > 0) {
            /* add neighbors */
            l->insert_neighbors(draws);
         }
         for (set<Sequence>::iterator dr=draws.begin(); dr!=draws.end(); dr++) {
            for (vector<ScapeStat>::iterator i = ar.scapestats.begin();
                  i != ar.scapestats.end(); i++) {
               if (l->cur_rep < i->reps) {
                  i->calc_and_show(cout, *dr);
               } else {
                  /* since we've sorted the stats by reps required, when we
                   * encounter a stat that's already done enough reps we know
                   * the rest have also done enough reps */
                  break;
               }
            }
         }
         draws.clear();
      }
   }

   if (rng != 0) delete rng;
   if (l!=0) delete l;

/* deal with various possible errors */
} catch (SimUsageError e) {
   cerr << endl << "detected usage error: " << e.detail << endl << endl;
   usage();
   return 1;
} catch(SimError &e) {
   cerr << "uncaught exception: " << e.detail << endl;

/* exit */
} return 0; }

/* print out usage information */
void usage(void) {
   cerr << 
"usage: regevscape [options]\n"
"\n"
"  general options:\n"
"     -m/--model (block|segal|neutral)\n" 
"     --popsim                   - run a population simulation\n"
"     --seed <int>               - random number generator seed\n"
"     --seq <string>             - a sequence to queue\n"
"\n"
"  model-specific options: \n"
"     block model:\n"
"        --mean <float>          - mean of within-block distr.\n"
"        --sd <float>            - sd of within-block distr.\n"
"        --blocks <int>          - number of blocks\n"
"        --length <int>          - length of sequence\n"
"\n"
"     segal model, general:\n"
"        --length <int>          - length of sequence\n"
"        --tfs <int>             - number of transcription factors\n"
"        --pwm <filename>        - filename containing PWM definitions\n"
"        --condition=c1,c2,...   - list of floats specifying TF levels\n"
"                                  multiple of conditions may be specified\n"
"        --steps <int>           - number of iid samples for estimating expr.\n"
"\n"
"     segal model automatic parameter estimation:\n"
"        --autofit               - when specified, will try and fit params\n"
"        --chains <int>          - number of parallel tempering chains\n"
"        --chain_steps <int>     - number of steps to run the chain\n"
"        --expression=e1,e2,...  - observed expression. One entry per seq.\n"
"        --psd_lambda <float>    - proposal variance for lambda updates\n"
"        --psd_alpha <float>     - proposal variance for alpha updates\n"
"        --psd_gamma <float>     - proposal variance for gamma updates\n"
//"        --alpha_lattice=a1,a2,  - grid of alpha (expr scalar) values for IS\n"
"\n"
"     segal model paremeters\n"
"        --lambda=l0,l1,l2,...   - list of floats for baseline and TF lambda\n"
"        --gamma=g11,g21,g12,g22 - cooperativity matrix by column\n"
"        --alpha=s1,s2,...       - expression scaling parameters\n"
"        --logistic_param=d      - use 1/(1+exp(-x/d))\n"
"        --tfsep <int>           - minimum spacing between TFs in a config\n"
"\n"
"     segal expression-fitness function\n"
"        --ffunc <choice>        - here are the fitness function choices.\n"
"                                  --coef must bi given with the correct\n"
"                                  number of coefficients.\n" 
"           linear : f = <x> dot <c> where x is the expr. vec., c coef vec.\n"
"           nop : f = 1\n"
"           f1 : f = a * (x1-b)^2 + c * (x2-d)^2 + e\n"
"           f2 : f = a^(-b*(x1-c)^2) * d^(-e*(x2-f)^2)\n"
"           f3 : f = prod(a^(-b*(xi-ci)^2))^(1/C) : i in (1..C) for C cond.\n"
"           f4 : f = exp((-1)*sum((xi-bi)^2/a)) : i in (1...C) for C cond.\n"
"        --coef=a,b,c,..         - list of coefficients for chosen fit. func.\n"
"\n"
"     segal model sampling parameters:\n"
"        --sampling=<iid|importance>  - which sampling method to employ\n"
"        --IS_gamma=g11,g21,g12,g22   - IS distr. cooperativity matrix\n"
"        --IS_alpha=s1,s2,...         - list of expr. scalars, one per TF\n"
"\n"
"  population simulation options:\n"
"     -g/--generations <int>     - number of gen to run\n"
"     -u/--mutation <float>      - per bp mutation rate\n"
"     -N/--popsize <int>         - effective population size\n"
"     --deletion=u,n,p           - rate, and negative binomial dup params\n" 
"     --duplication=u,n,p        - rate, and negative binomial dup params\n" 
"\n"
"  population statistics options (only valid in the presence\n"
"           of --popsim. May use multiple pstat options. Of the form\n"
"           --pstat=statname to print at the end or --pstat=statname,2\n"
"           if one wants statname printed every 2 generations. Default is\n"
"           each generation).\n"
"     --pstat=most_freq_seq      - most frequent allele\n"
"     --pstat=most_freq_noseq    - most frequent allele minus the sequence\n"
"     --pstat=mean_fitness       - population mean fitness\n"
"     --pstat=all_alleles        - details of each allele, sorted by # copies\n"
"     --pstat=allele_counter     - number of alleles queried so far\n"
//"     --pstat=site_frequencies   - positions and frequencies of seg. sites\n"
"\n"
"  population real-time statistics:\n"
"     --pstat=mutational_effects - print fitness info for new mutants\n"
"     --pstat=allele_loss        - print the generation and allele id for\n"
"                                  alleles when they're lost\n"
"\n"
"  landscape sampling/traversal options:\n"
"     --neighbors <int>          - for each sample, include neighbors within\n"
"                                  the specified number of steps\n"
"\n"
"  landscape statistic options (only valid in the absence of --popsim.\n"
"           May use multiple lstat options. In all cases, reps is the number\n"
"           of nodes for which the statistic will be printed, all uniformly\n"
"           sampled from the landscape. If multiple statistics are requested\n"
"           the same starting points will be used up until the number of reps\n"
"           for each particular statistic is reached. All parameters are\n"
"           always required for each statistic.)\n"
"     --lstat=fitness,reps          - node fitness\n"
//"     --lstat=n_step_mean,reps,n,s  - mean fitness (s samples) n steps away\n"
//"     --lstat=banding,reps,x        - band statistics for tried size x\n"
"     --lstat=sequence,reps         - print out sequence per rep\n"
"\n"
"  landscape statistics specifically for Segal model\n"
"     --lstat=expression,reps       - print vector of expression levels\n"
"     --lstat=tf_occupancy,reps     - print list of TF occupancy levels for\n"
"                                     each condition and TF combination\n"
"     --lstat=wsum,reps             - the sum of all configuration weights\n"
"     --lstat=configs,1             - print out configurations\n"
"\n";
   return;
}
