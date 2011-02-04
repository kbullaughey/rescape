#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include <string>
#include <ostream>

class Sequence {
public:
   /* public function */
   Sequence(std::string &s, std::string &alphabet);
   Sequence(std::string &alphabet, unsigned int len);
   Sequence(const Sequence &s, unsigned int from, unsigned int len);
   ~Sequence() { }
   unsigned int alphabet_size() { return alphabet.length(); }
   const char* c_str() { return seq.c_str(); }
   unsigned int length() const { return seq.size(); }

   /* public operators */
   
   void duplicate_part(unsigned int pos, unsigned int len);
   void delete_part(unsigned int pos, unsigned int len);
   char letter(unsigned int pos) const { return alphabet[seq[pos]]; }
   char code(unsigned int pos) const { return seq[pos]; }
   char code_to_letter(unsigned int code) const { return alphabet[code]; }
   void replace(unsigned int pos, char x);
   std::string& subseq(std::string &s, unsigned int pos, unsigned int len);
   friend std::ostream& operator<<(std::ostream &s, const Sequence &seq);
   friend bool operator<(const Sequence &s1, const Sequence &s2) {
      return s1.seq < s2.seq; }
   friend bool operator!=(const Sequence &s1, const Sequence &s2) {
      return s1.seq != s2.seq || s1.alphabet != s2.alphabet; }
protected:
   std::string seq;
   std::string alphabet;
};

#endif /* __SEQUENCE_H__ */
