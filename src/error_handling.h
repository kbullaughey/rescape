#ifndef __ERROR_HANDLING_H__
#define __ERROR_HANDLING_H__

#include <string>

class SimError {
public:
   SimError(const char *s = "unknown error");
   SimError(const std::string &s);
   std::string detail;
};

class SimUsageError : public SimError {
public:
   SimUsageError(const char *s);
   SimUsageError(const std::string &s);
};

#endif /* __ERROR_HANDLING_H__ */
