#ifndef __POPSIM_H__
#define __POPSIM_H__

#include "landscape.h"
#include "stats.h"
#include "command_line.h"
#include "sequence.h"

#include <string>
#include <map>
#include <ostream>
#include <utility>
#include <valarray>

#define DEFAULT_K 40 /* how many individual classes to put in the queue */

enum Mutation {point, deletion, duplication};

class Allele {
public:
   /* member functions */
   Allele(Sequence s, unsigned int num_copies, unsigned int gen,
      Landscape *ls);
   ~Allele();
   Sequence& get_seq(void) { return seq; }
   static unsigned int alleles_queried(void) { return next_allele_id-1; }
   int length(void) { return seq.length(); }

   /* operators */
   friend std::ostream& operator<<(std::ostream &s, const Allele &a);
   friend bool operator<(const Allele&, const Allele&);

   double fitness;
   unsigned int copies;
   unsigned int allele_id;
   unsigned int arisal_gen;
   unsigned int mutations;
private:
   /* member variables */
   Landscape *landscape;
   Sequence seq;
   static unsigned int next_allele_id;
};

typedef std::map<const Sequence, Allele *> AlleleList;

class PopState {
public:
   /* public member functions */
   PopState(unsigned int popsize, double mu_rate, double rec_rate, 
      Landscape *ls, std::valarray<double> *del_params, 
      std::valarray<double> *dup_params);
   std::pair<Allele*,int> find_mutation_allele_and_position(int pos);
   void run(const std::vector<PopStat> &s, unsigned int generations);

   /* operators */
   friend std::ostream& operator<<(std::ostream &s, const PopState &ps);
   PopState& operator++();   /* prefix */

   /* public member variables */
   AlleleList alleles;
   unsigned int generation;
   unsigned int N;
   Landscape *scape;
   static std::map<std::string, bool> real_time_flags;

private:
   double mu;
   double r;
   double deletion_prob;
   double deletion_neg_bin_n;
   double deletion_neg_bin_p;
   double duplication_prob;
   double duplication_neg_bin_n;
   double duplication_neg_bin_p;

   /* private member functions */
   void advance(void);
   void drift(void);
   int total_allele_sequence();
   void mutate(Mutation m, int tot);
};

#endif /* __POPSIM_H__ */

/* END */



