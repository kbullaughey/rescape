#include "stats.h"
#include "error_handling.h"
#include "popsim.h"
#include "sequence.h"
#include "landscape.h"

#include <string>
#include <ostream>
#include <map>
#include <iostream>

using std::string;
using std::ostream;
using std::map;
using std::endl;
using std::cerr;

/******************************************************************************/
/* Stat functions */
/******************************************************************************/

Stat::Stat() {
}

/******************************************************************************/
/* PopStat member functions */
/******************************************************************************/

map<string, PopStatPut> PopStat::valid_stats;

PopStat::PopStat(const char *cstr) : Stat() {
   /* when first PopStat object is created, populate the valid stats list */
   if (valid_stats.size() == 0) populate_valid_stats();

   string s(cstr);
   int pos = s.find(',');
   popstate = 0;

   if (pos == (int)string::npos) {
      modulo = PRINT_STAT_AT_END;
      name = s;
   } else {
      string gen(s, pos+1);
      if (gen.length() == 0) 
         throw SimUsageError("zero length modulo-generation parameter");
      char *end;
      modulo = (unsigned int) strtoul(gen.c_str(), &end, 10);
      if (gen.c_str() == end)
         throw SimUsageError("invalid modulo-generation parameter");
      if (pos <= 1)
         throw SimUsageError("zero length modulo-generation parameter");
      name = s.substr(0, pos);
   }

   map<string,PopStatPut>::iterator i = valid_stats.find(name);
   if (i == valid_stats.end()) 
      throw SimUsageError("invalid population statistic: " + name);
   if (i->second == 0) {
      /* this means real-time flag needs to be set */
      stype = realtime;
      map<string,bool>::iterator j = PopState::real_time_flags.find(name);
      if (j == PopState::real_time_flags.end())
         throw SimUsageError("invalid popsim real-time statistic: " +name);
      j->second = true;
   } else {
      put = i->second;
      stype = queued;
   }
}

/* list of function declarations for pop stat functions */
void pstat_test(ostream &s, const PopStat &x);
void pstat_most_freq_seq(ostream &s, const PopStat &x);
void pstat_most_freq_noseq(ostream &s, const PopStat &x);
void pstat_mean_fitness(ostream &s, const PopStat &x);
void pstat_all_alleles(ostream &s, const PopStat &x);
void pstat_allele_counter(ostream &s, const PopStat &x);
void pstat_site_frequencies(ostream &s, const PopStat &x);

/* populate the pop stats with types and functions */
void PopStat::populate_valid_stats(void) {
   valid_stats[string("most_freq_seq")] = pstat_most_freq_seq;
   valid_stats[string("most_freq_noseq")] = pstat_most_freq_noseq;
   valid_stats[string("mean_fitness")] = pstat_mean_fitness;
   valid_stats[string("all_alleles")] = pstat_all_alleles;
   valid_stats[string("allele_counter")] = pstat_allele_counter;
   valid_stats[string("site_frequencies")] = pstat_site_frequencies;
   valid_stats[string("mutational_effects")] = 0;
   valid_stats[string("allele_loss")] = 0;

   /* population real-time stat flags, by default, not enabled */
   PopState::real_time_flags[string("mutational_effects")] = false;
   PopState::real_time_flags[string("allele_loss")] = false;
}

ostream& operator<<(ostream &s, const PopStat &x) {
   s << "gen: " << x.popstate->generation << " ";
   x.put(s, x);
   s << endl;
   return s;
}

/* functions that output pop statistics */

/* a test function */
void pstat_test(ostream &s, const PopStat &x) {
   s << "null";
}

/* output the most frequent allele, its fitness and number of copies */
void pstat_most_freq_seq(ostream &s, const PopStat &x) {
   AlleleList::const_iterator max, i;
   const PopState *p = x.popstate;
   for (max = i = p->alleles.begin();i!=p->alleles.end();i++) {
      if (max->second->copies < i->second->copies)
         max = i;
   }
   s << "pstat_most_freq_seq: " << *max->second;
}

/* output the most frequent allele, its fitness and number of copies */
void pstat_most_freq_noseq(ostream &s, const PopStat &x) {
   AlleleList::const_iterator max, i;
   const PopState *p = x.popstate;
   for (max = i = p->alleles.begin();i!=p->alleles.end();i++) {
      if (max->second->copies < i->second->copies)
         max = i;
   }
   s << "pstat_most_freq_noseq: " << max->second->allele_id << " " 
      << max->second->fitness << " " << max->second->copies;
}

/* display the mean fitness of the population */
void pstat_mean_fitness(ostream &s, const PopStat &x) {
   const PopState *p = x.popstate;
   valarray<double> fitnesses(p->alleles.size());
   valarray<double> frequencies(p->alleles.size());
   unsigned int k = 0;
   for (AlleleList::const_iterator i=p->alleles.begin(); 
         i!=p->alleles.end(); i++) {
      fitnesses[k] = i->second->fitness;
      frequencies[k] = (double)i->second->copies / (double)p->N;
      k++;
   }
   double mean_fitness = (fitnesses * frequencies).sum();
   s << "pstat_mean_fitness: " << mean_fitness;
}

/* print details of all alleles */
void pstat_all_alleles(ostream &s, const PopStat &x) {
   s << "pstat_all_alleles: ";
   for (AlleleList::const_iterator i=x.popstate->alleles.begin(); 
         i!=x.popstate->alleles.end(); i++) {
      s << endl << *(i->second);
   }
   s << endl;
}

/* print out the total number of alleles queried so far */
void pstat_allele_counter(ostream &s, const PopStat &x) {
   s << "pstat_allele_counter: " << Allele::alleles_queried();
}

/* print out the frequencies of segregating sites */
void pstat_site_frequencies(ostream &s, const PopStat &x) {
   s << "pstat_site_frequencies: no longer implemented";
   return;
#if 0
   s << "pstat_site_frequencies: ";
   bool last_was_segregating = false;
   for (unsigned int z = 0; z <= x.popstate->scape->L; z++) {
      char base = 0;
      bool this_is_segregating = false;
      for (AlleleList::const_iterator i=x.popstate->alleles.begin(); 
            i!=x.popstate->alleles.end(); i++) {
         /* don't do this when we're past the end */
         if (z < x.popstate->scape->L) {
            /* determine if this nucleotide segregates, we use this next iter */
            if (base == 0) {
               base = i->first.letter(z);
            } else if (base != i->first.letter(z)) {
               this_is_segregating = true;
            }
         }

         if (last_was_segregating) {
            s << endl << "site: " << z-1 << " allele: " << i->second->allele_id 
               << " base: " << i->first.letter(z-1) << " copies: " 
               << i->second->copies;
         }
      }
      last_was_segregating = this_is_segregating;
   }
   s << endl;
#endif
}

/******************************************************************************/
/* ScapeStat member functions */
/******************************************************************************/

map<string,ScapeStatPut> ScapeStat::valid_stats;

ScapeStat::ScapeStat(const char *cstr) : Stat() {
   /* when first ScapeStat object is created, populate the valid stats list */
   if (valid_stats.size() == 0) populate_valid_stats();

   string s(cstr);
   int pos = s.find(',');

   if (pos == (int)string::npos) 
      throw SimUsageError("lstat " + s + " must specify params (like reps)");

   string reps_str;
   name = s.substr(0, pos);

   map<string,ScapeStatPut>::iterator i = valid_stats.find(name);
   if (i == valid_stats.end()) 
      throw SimUsageError("invalid landscape statistic: " + name);
   put = i->second;

   string remainder(s, pos+1);
   if (remainder.length() == 0) 
      throw SimUsageError("zero length parameter string");

   /* look for any more commas */
   pos = remainder.find(',');
   if (pos == (int)string::npos) {
      char *end;
      reps = (unsigned int) strtoul(remainder.c_str(), &end, 10);
      if (remainder.c_str() == end)
         throw SimUsageError("invalid reps parameter for lstat " + name);
   } else {
      /* more params */
      if (pos <= 1) throw SimUsageError("zero length string reps parameter");
      throw SimError("multi-parameter lstats not implemented yet");
   }
}

/* list of function declarations for scape stats */
void lstat_fitness(ostream &s, const ScapeStat &x, const Sequence &seq);
void lstat_test(ostream &s, const ScapeStat &x, const Sequence &seq);
void lstat_expression(ostream &s, const ScapeStat &x, const Sequence &seq);
void lstat_sequence(ostream &s, const ScapeStat &x, const Sequence &seq);
void lstat_tf_occupancy(ostream &s, const ScapeStat &x, const Sequence &seq);
void lstat_wsum(ostream &s, const ScapeStat &x, const Sequence &seq);
void lstat_configs(ostream &s, const ScapeStat &x, const Sequence &seq);

void ScapeStat::populate_valid_stats(void) {
   valid_stats[string("fitness")] = lstat_fitness;
   valid_stats[string("test")] = lstat_test;
   valid_stats[string("expression")] = lstat_expression;
   valid_stats[string("sequence")] = lstat_sequence;
   valid_stats[string("tf_occupancy")] = lstat_tf_occupancy;
   valid_stats[string("wsum")] = lstat_wsum;
   valid_stats[string("configs")] = lstat_configs;
}

void ScapeStat::calc_and_show(ostream &s, const Sequence &seq) {
   s << "rep: " << scape->cur_rep << " ";
   put(s, *this, seq);
   s << endl;
}

/* ordering of ScapeStat objects is based on how many times they should be 
 * printed */
bool operator<(const ScapeStat &x, const ScapeStat &y) {
   return x.reps < y.reps;
}

/* functions that output landscape statistics */

/* print out the fitness for this sequence */
void lstat_fitness(ostream &s, const ScapeStat &x, const Sequence &seq) {
   s << "lstat_fitness: " <<x.scape->get_fitness(seq);
}

/* nop function */
void lstat_test(ostream &s, const ScapeStat &x, const Sequence &seq) {
   s << "test function";
}

/* print out the expression vector for this sequence */
void lstat_expression(ostream &s, const ScapeStat &x, const Sequence &seq) {
   LandscapeSegalModel* ls = dynamic_cast<LandscapeSegalModel*>(x.scape);
   if (ls == NULL) throw SimError("lstat_expression only for Segal model");

   valarray<double> expr(ls->conditions.size());
   ls->get_expression(seq, expr);

   s << "lstat_expression:";
   for (unsigned int i=0; i<ls->conditions.size(); i++)
      s << " " << expr[i];
}

/* print out the mean tf occupancy for each condition */
void lstat_tf_occupancy(ostream &s, const ScapeStat &x, const Sequence &seq) {
   LandscapeSegalModelSpecified* ls = 
      dynamic_cast<LandscapeSegalModelSpecified*>(x.scape);
   if (ls == NULL) throw SimError("lstat_tf_occupancy only for Segal model");

   s << "lstat_tf_occupancy: ";
   for (int cond = 0; cond < (int)ls->conditions.size(); cond++) {
      Array oc(2, ls->M * seq.length());
      oc.dims[0] = (int)ls->M;
      oc.dims[1] = (int)seq.length();
      ls->get_tf_occupancy(seq, cond, oc);
      for (int i=0; i < (int)ls->M; i++) {
         s << "condition: " << cond << " tf: " << i << " oc: ";
         for (int x=0; x < (int)seq.length(); x++)
            s << " " << oc.pos(nop, i, x);
         if (cond+1 != (int)ls->conditions.size() || i+1 != (int)ls->M) 
            s << endl << "rep: X lstat_tf_occupancy: ";
      }
   }
}

/* just print out the sequence */
void lstat_sequence(ostream &s, const ScapeStat &x, const Sequence &seq) {
   s << seq;
}

/* print out the sum of all configuration weights for each condition */
void lstat_wsum(ostream &s, const ScapeStat &x, const Sequence &seq) {
   LandscapeSegalModelSpecifiedIID* ls = 
      dynamic_cast<LandscapeSegalModelSpecifiedIID*>(x.scape);
   if (ls == NULL) throw SimError("lstat_wsum only for SegalModelSpecifiedIID");
   
   valarray<double> wsum((int)ls->conditions.size());
   ls->get_wsum(seq, wsum);
   s << "wsum:";
   for (int cond = 0; cond < (int)ls->conditions.size(); cond++) 
      s << " " << wsum[cond];
   s << endl;
}

/* print out each sampled configuration */
void lstat_configs(ostream &s, const ScapeStat &x, const Sequence &seq) {
   LandscapeSegalModelSpecifiedIID* ls = 
      dynamic_cast<LandscapeSegalModelSpecifiedIID*>(x.scape);
   if (ls == NULL) throw SimError("lstat_configs only for SegalModelSpecifiedIID");

   ProbTable ptables(ls->M, ls->D, seq.length());
   s << "lstat_configs: ";
   for (unsigned int i=0; i<ls->conditions.size(); i++) {
      if (i > 0)
         s << "rep: X lstat_configs: ";
      s << "condition: " << i << " printing configs" << endl;
      ls->calc_prob_tables_wrapper(i, seq, ptables);
      ls->prob_expr(ptables, seq, i, true);
   }
}

/* END */
