#include "landscape.h"
#include "error_handling.h"
#include "sim_rand.h"
#include "array.h" /* this gives us prod() */

#include <iostream>
#include <string>
#include <set>
#include <map>
#include <string.h>

using std::string;
using std::set;
using std::cerr;
using std::endl;
using std::map;

/******************************************************************************/
/* Implementation of Landscape class members */
/******************************************************************************/

/* constructor */
Landscape::Landscape(RandomNumberGenerator *rngp) : alphabet("agct") {
   rng = rngp;
   cur_rep = 0;
}

/* function that does a model lookup based on a c string */
Model Landscape::get_model(char *m) {
   if (strncmp(m, "block", 5) == 0) 
      return block;
   else if (strncmp(m, "segal", 5) == 0)
      return segal;
   else if (strncmp(m, "neutral", 7) == 0)
      return neutral;
   else 
      throw SimUsageError("invalid model");
   return invalid;
}

/* given an offset into the alphabet, pick a different one */
/* note that this returns the position in the alphabet string, not the char */
char Landscape::pick_mutation(char off) {
   int steps = rng->uniform_int(1, alphabet.length());
   char ret = (off + steps) % alphabet.length();
   return ret;
}

/* generate the next random sample from the landscape */
Sequence Landscape::sample(unsigned int L) {
   Sequence seq(alphabet, L);
   rng->uniform_dna(seq, L, this);
   return seq;
}

/* insert all one-step away neighbors of all the sequences in s into s */
void Landscape::insert_neighbors(set<Sequence> &s) {
   /* make a copy */
   set<Sequence> old = s;
   for (set<Sequence>::iterator i = old.begin(); i != old.end(); i++) {
      for (unsigned int x = 0; x < i->length(); x++) {
         for (unsigned int j=0; j< alphabet.length(); j++) {
            /* make a copy of i */
            Sequence neighbor(*i);
            if (neighbor.letter(x) == alphabet[j]) continue;
            neighbor.replace(x, j);
            s.insert(neighbor);
         }
      }
   }
   return;
}

/******************************************************************************/
/* functions for neutral Model */
/******************************************************************************/
LandscapeNeutralModel::LandscapeNeutralModel(RandomNumberGenerator *rng) : 
      Landscape(rng) { }

double LandscapeNeutralModel::get_fitness(const Sequence &seq) {
   return 1.0;
}

/******************************************************************************/
/* functions for Block Model */
/******************************************************************************/

/* set up the cache and parameters of this landscape model 
 * hash_size and seed are optional
 */
LandscapeBlockModel::LandscapeBlockModel(unsigned int B, unsigned int bw, 
      double mean, double sd, RandomNumberGenerator *rng, unsigned int hash_size) 
      : Landscape(rng), caches(B) {
   L = B*bw;
   blocks = B;
   blk_width = bw;
   mu = mean;
   sigma = sd;
   for (unsigned int i=0; i < blocks; i++) 
      caches[i] = new StringHash(hash_size);
}

/* destructor */
LandscapeBlockModel::~LandscapeBlockModel() {
   for (unsigned int i=0; i < blocks; i++) {
      delete caches[i];
   }
}

/* query the fitness corresponding to the nucleotide sequence, s */
double LandscapeBlockModel::get_fitness(const Sequence &s) {
   double fitness, total_fitness = 1.0;
   if (s.length() != L) throw SimError("sequence wrong length");
   for (unsigned int i=0; i <blocks; i++) {
      Sequence thisBlock(s, i*blk_width, blk_width);
      StringHash::iterator j = caches[i]->find(thisBlock.c_str());
      if (j == caches[i]->end()) {
         char *buf = (char*)malloc(sizeof(char)*thisBlock.length());
         if (buf == 0) throw SimError("not enough memory");
         strncpy(buf, thisBlock.c_str(), thisBlock.length());
         fitness = (*caches[i])[buf] = sample_fitness();
      } else {
         fitness = j->second;
      }
      fitness = sample_fitness();
      total_fitness += fitness;
   }
   return total_fitness;
}

/******************************************************************************/
/* Implementation of FitFunc class members */
/******************************************************************************/

double fitness_func_f1(valarray<double> &x, valarray<double> &c);
double fitness_func_linear(valarray<double> &x, valarray<double> &c);
double fitness_func_nop(valarray<double> &x, valarray<double> &c);
double fitness_func_f2(valarray<double> &x, valarray<double> &c);
double fitness_func_f3(valarray<double> &x, valarray<double> &c);
double fitness_func_f4(valarray<double> &x, valarray<double> &c);

/* default constructor */
FitFunc::FitFunc(void) {
   calc_var = 0;
   c = 0;
}

/* real constructor */
FitFunc::FitFunc(string &choice) {
   calc_var = 0;
   c = 0;
   /* if there are no fitness functions loaded, load the set of possible ones */
   if (funclist.size() == 0) {
      funclist[string("linear")] = fitness_func_linear;
      funclist[string("nop")] = fitness_func_linear;
      funclist[string("f1")] = fitness_func_f1;
      funclist[string("f2")] = fitness_func_f2;
      funclist[string("f3")] = fitness_func_f3;
      funclist[string("f4")] = fitness_func_f4;
   }
   map<string,FitCalc>::iterator i = funclist.find(choice);
   if (i == funclist.end()) throw SimError("invalid fitness function choice");
   calc_var = i->second;
}

/* set the coefficients for this function */
void FitFunc::set_coef(valarray<double> &coef) {
   c = new valarray<double>(coef.size());
   *c = coef;
   return;
}

/* fitness functions */

double fitness_func_linear(valarray<double> &x, valarray<double> &c) {
   if (x.size() != c.size() || x.size() == 0)
      throw SimError("expression-coefficient dimmension mismatch");
   return (x*c).sum();
}

double fitness_func_nop(valarray<double> &x, valarray<double> &c) {
   return 1;
}
double fitness_func_f1(valarray<double> &x, valarray<double> &c) {
   if (x.size() != 2 || c.size() != 5)
      throw SimError("expression-coefficient dimmension mismatch");
   return c[0]*pow(x[0]-c[1],2)+c[2]*pow(x[1]-c[3],2)+c[4];
}

double fitness_func_f2(valarray<double> &x, valarray<double> &c) {
   if (x.size() != 2 || c.size() != 6)
      throw SimError("expression-coefficient dimmension mismatch");
   return pow(c[0], c[1]*pow(x[0]-c[2],2)) * pow(c[3], c[4]*pow(x[1]-c[5],2));
}

double fitness_func_f3(valarray<double> &x, valarray<double> &c) {
   if (x.size()+2 != c.size())
      throw SimError("expression-coefficient dimmension mismatch");
   double res = 1;
   for (int i = 0; i < (int)x.size(); i++) {
      res *= pow(c[0], c[1]*pow(x[i]-c[i+2],2));
   }
   return pow(res, 1/((double)x.size()));
}

double fitness_func_f4(valarray<double> &x, valarray<double> &c) {
   if (x.size()+1 != c.size())
      throw SimError("expression-coefficient dimmension mismatch");
   double res = 0;
   for (int i = 0; i < (int)x.size(); i++) {
      res -= (x[i]-c[i+1])*(x[i]-c[i+1]);
   }
   res /= c[0];
   return exp(res);
}


/* END */



