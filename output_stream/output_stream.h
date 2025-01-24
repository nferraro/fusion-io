#ifndef OUTPUT_STREAM
#define OUTPUT_STREAM

#include <iostream>

namespace output {
    // Control variable for output
    inline bool quiet = true;

    // Custom output stream derived from std::ostream
    class ConditionalCerr : public std::ostream {
    private:
        bool* condition;

    public:
        // Constructor to initialize with a condition variable
        explicit ConditionalCerr(bool* conditionVar)
            : condition(conditionVar) {}

        // Overloaded << operator for output, checks the condition
        template <typename T>
        ConditionalCerr& operator<<(const T& value) {
            if (!*condition) {
                std::cerr<<(value);
            }
            return *this;
        }

        // Overloaded << for manipulators like std::endl
        ConditionalCerr& operator<<(std::ostream& (*manip)(std::ostream&)) {
            if (!*condition) {
                std::cerr<<manip; // Apply the manipulator to std::cerr
            }
            return *this;
        }
    };

    // Global instance of ConditionalCerr
    inline ConditionalCerr qerr(&quiet);
} // namespace ConditionalOutput

#endif
