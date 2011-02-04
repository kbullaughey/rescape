#ifndef __STATS_H__
#define __STATS_H__

#include "sequence.h"

#include <string>
#include <ostream>
#include <map>

/* abstract stat class */
class Stat {
public:
   Stat(); 
   virtual ~Stat() { }
   friend std::ostream& operator<<(std::ostream&, const Stat&);
   std::string name;
};

#define PRINT_STAT_AT_END 0

/* PopStat needs to know about PopState, but it may not get included 
 * until later */
class PopState;
class PopStat;

typedef void (*PopStatPut)(std::ostream &s, const PopStat &x);
enum StatType {queued, realtime};

/* stat class for popsim statistics */
class PopStat : public Stat {
public:
   /* public member functions */
   PopStat(const char *);
   ~PopStat() { }
   friend std::ostream& operator<<(std::ostream&, const PopStat&);
   void set_popstate(PopState *p) { popstate=p; }

   /* public member variables */
   unsigned int modulo; /* period in generatioins between printings */
   const PopState *popstate;
   StatType stype;

private:
   /* private member functions */
   static void populate_valid_stats(void);

   /* private member variables */
   static std::map<std::string, PopStatPut> valid_stats;
   PopStatPut put;
};

class Landscape;
class ScapeStat;
typedef void (*ScapeStatPut)(std::ostream &s, const ScapeStat &x, 
   const Sequence &seq);

/* stat class for landscape statistics */
class ScapeStat : public Stat {
public:
   /* public member functions */
   ScapeStat(const char *); 
   ~ScapeStat() { }
   void set_landscape(Landscape *l) { scape=l; }
   void calc_and_show(std::ostream &s, const Sequence &seq);
   friend bool operator<(const ScapeStat&, const ScapeStat&);

   /* public member variables */
   unsigned int reps;
   Landscape *scape;

private:
   static void populate_valid_stats(void);

   /* private member variables */
   static std::map<std::string, ScapeStatPut> valid_stats;
   ScapeStatPut put;
};


#endif /* __STATS_H__ */
