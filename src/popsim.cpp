#include "popsim.h"
#include "error_handling.h"
#include "sequence.h"

#include <string>
#include <math.h>
#include <ostream>
#include <map>
#include <iostream>
#include <valarray>
#include <vector>
#include <utility>

using std::vector;
using std::string;
using std::endl;
using std::ostream;
using std::map;
using std::cerr;
using std::cout;
using std::valarray;
using std::pair;

/*****************************************************************************/
/* member functions for class Allele */
/*****************************************************************************/

unsigned int Allele::next_allele_id = 0;
std::map<std::string, bool> PopState::real_time_flags;

Allele::Allele(Sequence s, unsigned int num_copies, unsigned int gen, 
      Landscape *ls) : seq(s) {
   landscape = ls;
   fitness = ls->get_fitness(s);
   allele_id = next_allele_id++;
   copies = num_copies;
   arisal_gen = gen;
   mutations = 0;
}

Allele::~Allele() {
   if (this!=0)
      delete this;
}

/* A print stream operator for Allele objects */
ostream& operator<<(ostream &s, const Allele &a) {
   s << a.allele_id << " " << a.mutations << " " << a.seq << " " 
      << a.fitness << " " << a.copies;
   return s;
}

/* ordering of Alleles is based on number of copies */
bool operator<(const Allele &x, const Allele &y) {
   return x.copies < y.copies;
}


/*****************************************************************************/
/* member functions for lass PopState */
/*****************************************************************************/

PopState::PopState(unsigned int popsize, double mu_rate, double rec_rate, 
      Landscape *ls, valarray<double> *del_params, valarray<double> *dup_params)
       : alleles() {
   scape = ls;
   N = popsize;
   mu = mu_rate;
   r = rec_rate;
   scape = ls;
   generation = 0;
   
   /* set up the deletion and duplication parameters */
   if (del_params == 0) {
      deletion_prob = 0;
      deletion_neg_bin_n = 0;
      deletion_neg_bin_p = 0;
   } else {
      deletion_prob = (*del_params)[0];
      deletion_neg_bin_n = (*del_params)[1];
      deletion_neg_bin_p = (*del_params)[2];
   }
   if (dup_params == 0) {
      duplication_prob = 0;
      duplication_neg_bin_n = 0;
      duplication_neg_bin_p = 0;
   } else {
      duplication_prob = (*dup_params)[0];
      duplication_neg_bin_n = (*dup_params)[1];
      duplication_neg_bin_p = (*dup_params)[2];
   }
}

/* A print stream operator for PopState objects */
ostream& operator<<(ostream &s, const PopState &ps) {
   s << "N=" << ps.N << ", mu=" << ps.mu << ", r=" << ps.r;
   s << ", alleles=" << ps.alleles.size();
   s << ", generation=" << ps.generation << endl;
   
   for (AlleleList::const_iterator i = ps.alleles.begin(); 
         i != ps.alleles.end(); i++) {
      s << "allele: " << *(i->second);
   }
   return s;
}

/* Advance the simulation one generation using the prefix increment opeartor */
PopState& PopState::operator++() {
   advance();
   return *this;
}

void PopState::mutate(Mutation m, int tot) {
   int pos = scape->rng->uniform_int(0, tot);
   pair<Allele *, int> res = find_mutation_allele_and_position(pos);
   Allele *a = res.first;
   int x = res.second;
   int len;
   Sequence seq(a->get_seq());
   string type;

   /* do mutation type-specific stuff */
   switch (m) {
      case point:{
         type = "point";
         len = 1;
         char replacement = scape->pick_mutation(seq.code(x));
         seq.replace(x, replacement);
         break;
      }
      case deletion: {
         len = (int)scape->rng->rnb(deletion_neg_bin_n, deletion_neg_bin_p);
         type =  "deletion";

         /* throw out mutations that go beyond the end of the sequence or are 
          * zero length */
         if (x+len > (int)seq.length() || len == 0) return; 

         seq.delete_part(x, len);
         break;
      }
      case duplication: {
         len = (int)scape->rng->rnb(duplication_neg_bin_n, 
            duplication_neg_bin_p);
         type = "duplication";

         /* throw out mutations that go beyond the end of the sequence or are 
          * zero length */
         if (x+len > (int)seq.length() || len == 0) return; 

         seq.duplicate_part(x, len);
         break;
      }
      default:
         throw SimError("invalid mutation type");
   }

   /* see if the allele already exists */
   AlleleList::iterator i = alleles.find(seq);
   bool isnew = false;
   if (i == alleles.end()) {
      alleles[seq] = new Allele(seq, 1, generation, scape);
      alleles[seq]->mutations = a->mutations+1;
      isnew = true;
   } else {
      i->second->copies++;
   }

   /* see if we should print mutation information */
   if (real_time_flags[string("mutational_effects")]) {
      Sequence bg = a->get_seq();
      string from, to;
      from = bg.subseq(from, x, len);
      switch (m) {
         case point:
            to = seq.subseq(to, x, len);
            break;
         case deletion: 
            to.clear();
            break;
         case duplication: 
            to = bg.subseq(to, x, len);
            to = to + to;
            break;
         default:
            throw SimError("invalid mutation type");
      }
      cout << "gen: " << generation << " pstat_mutational_effects: " 
         << "background: " << bg 
         << " old_id: " << a->allele_id 
         << " new_id: " << alleles[seq]->allele_id
         << " copies: " << a->copies 
         << " type: " << type
         << " site: " << x
         << " len: " << len
         << " from: '" << from << "'"
         << " to: '" << to << "'"
         << " bfit: " << a->fitness
         << " mfit: " << alleles[seq]->fitness
         << " new: " << seq
         << " isnew: " << isnew
         << endl;
   }

   if (a->copies <= 0) throw SimError("too few alleles");
   if (a->copies == 1) {
      AlleleList::iterator k = alleles.find(a->get_seq());
      if (k == alleles.end()) throw SimError("where did the allele go?");
      if (real_time_flags[string("allele_loss")]) {
         cout << "gen: " << generation << " pstat_allele_loss: " 
            << k->second->allele_id << endl;
      }
      alleles.erase(k);
   } else {
      a->copies--;
   }
}

/* Advance the the simulation */
void PopState::advance(void) {
   drift();

   int tot = total_allele_sequence();
   /* add mutations */
   int point_mutations = scape->rng->rpois(mu*tot);
   int deletions = scape->rng->rpois(deletion_prob*tot);
   int duplications = scape->rng->rpois(duplication_prob*tot);
   while (point_mutations-- > 0) {
      mutate(point, tot);
   }
   while (deletions-- > 0) {
      mutate(deletion, tot);
   }
   while (duplications-- > 0) {
      mutate(duplication, tot);
   }
   generation++;
}

/* return the concatenated length of all alleles */
int PopState::total_allele_sequence(void) {
   int tot = 0;
   for (AlleleList::iterator i=alleles.begin(); i!=alleles.end(); i++) {
      tot += (i->second->copies)*(i->second->length());
   }
   return tot;
}

/* take an offset into the total concatenated set of allele sequences and 
 * determine which allele and position the offset corresponds to */
pair<Allele*,int> PopState::find_mutation_allele_and_position(int pos) {
   pair<Allele*,int> res;
   int running_sum = 0;
   int prev_sum;
   for (AlleleList::iterator i=alleles.begin(); i!=alleles.end(); i++) {
      prev_sum = running_sum;
      running_sum += ((int)i->second->copies)*(i->second->length());
      if (pos < running_sum) {
         res.first = i->second;
         res.second = (pos-prev_sum) % i->second->length();
         return res;
      }
   }
   /* if we got here we didn't find the individual */
   throw SimError("Failed to find allele");
}


/* Adjust the allele counts according to one generation of drift */
void PopState::drift(void) {
   valarray<double> fitnesses(alleles.size());
   valarray<double> frequencies(alleles.size());
   unsigned int k = 0;
   unsigned int total_copies = 0;
   for (AlleleList::const_iterator i=alleles.begin(); i!=alleles.end(); i++) {
      fitnesses[k] = i->second->fitness;
      frequencies[k] = (double)i->second->copies / (double)N;
      total_copies += i->second->copies;
      k++;
   }
   if (total_copies != N) throw SimError("population size changed!");
   double mean_fitness = (fitnesses * frequencies).sum();
   valarray<double> p = fitnesses * frequencies / mean_fitness;
   
   /* draw the next generation's counts from a multinomial distribution */
   valarray<unsigned int> next_gen_counts(alleles.size());
   scape->rng->multinomial(N, p, &next_gen_counts[0]);

   /* update the allele counts with the new numbers */
   k = 0;
   AlleleList::iterator to_erase;
   for (AlleleList::iterator i = alleles.begin(); i != alleles.end();) {
      if (next_gen_counts[k] == 0) {
         to_erase = i;
         i++;
         if (real_time_flags[string("allele_loss")]) {
            cout << "gen: " << generation << " pstat_allele_loss: " 
               << to_erase->second->allele_id << endl;
         }
         alleles.erase(to_erase);
      } else {
         i->second->copies = next_gen_counts[k];
         i++;
      }
      k++;
   }
}

/* run the population simulation */
void PopState::run(const vector<PopStat> &popstats, unsigned int gens) {

   while (generation < gens) {
      /* cycle through the stats to see if any should be displayed */
      for (vector<PopStat>::const_iterator i = popstats.begin(); 
            i != popstats.end(); i++)
         if (generation % i->modulo == 0) cout << *i;
      ++(*this); /* advance to the next generation */
   }

   /* cycle through the stats to see if any should be displayed */
   for (vector<PopStat>::const_iterator i = popstats.begin(); 
         i != popstats.end(); i++)
      if (generation % i->modulo == 0) cout << *i;

   return;
}

/* END */
