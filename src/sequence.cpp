#include "sequence.h"
#include "error_handling.h"

#include <string>
#include <iostream>
#include <sstream>

using std::string;
using std::ostream;
using std::cerr;
using std::endl;
using std::ostringstream;

Sequence::Sequence(string &alpha, unsigned int len) : seq(len, (char)0), 
   alphabet(alpha) { }

Sequence::Sequence(string &s, string &alpha) : seq(s.length(), (char)0), 
      alphabet(alpha) {
   int which;
   for (unsigned int k = 0; k < s.length(); k++) {
      which = alphabet.find(s[k]);
      if (which == (int)string::npos) { 
         ostringstream tmp;
         tmp << "letter #" << (int)s[k] << " not found in alphabet, k=" << k;
         throw SimError(tmp.str());
      }
      seq[k] = which;
   }
}

Sequence::Sequence(const Sequence &s, unsigned int from, unsigned int len) 
   : seq(string(s.seq, from, len)), alphabet(s.alphabet) { }

ostream& operator<<(ostream &s, const Sequence &seq) {
   string out;
   for (unsigned int k = 0; k < seq.length(); k++) {
      out += seq.alphabet[seq.seq[k]];
   }
   s << out;
   return s;
}

/* which is a valid offset into alphabet, not the letter itself */
void Sequence::replace(unsigned int pos, char which) {
   if (pos >= seq.length()) throw SimError("sequence position out of range");
   seq[pos] = which;
}

/* delete part of the sequence */
void Sequence::delete_part(unsigned int pos, unsigned int len) {
   seq.erase(pos, len);
}

/* duplicate part of the sequence */
void Sequence::duplicate_part(unsigned int pos, unsigned int len) {
   seq.insert(pos, seq, pos, len);
}

string& Sequence::subseq(string &s, unsigned int pos, unsigned int len) {
   s.clear();
   for (unsigned int k = 0; k < len; k++) {
      s += alphabet[seq[pos+k]];
   }
   return s;
}


/* END */



