#include "error_handling.h"

#include <string>

SimError::SimError(const char *s) {
   detail = s;
}

SimError::SimError(const std::string &s) {
   detail = s;
}

SimUsageError::SimUsageError(const char *s) 
      : SimError(s) {
}

SimUsageError::SimUsageError(const std::string &s) 
      : SimError(s) {
}

/* END */
