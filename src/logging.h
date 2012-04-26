/*=============================================================================

  Matt Rasmussen

  Logging functions

=============================================================================*/


#ifndef ARGHMM_LOGGING_H
#define ARGHMM_LOGGING_H

// c/c++ includes
#include <stdio.h>
#include <sys/time.h>

namespace arghmm {


// logging
enum {
    LOG_QUIET=0,
    LOG_LOW=1,
    LOG_MEDIUM=2,
    LOG_HIGH=3
};

void printError(const char *fmt, ...);
void printLog(int level, const char *fmt, ...);
bool openLogFile(const char *filename);
void openLogFile(FILE *stream);
void closeLogFile();
FILE *getLogFile();
void setLogLevel(int level);
bool isLogLevel(int level);


// timing
class Timer
{
public:
    Timer(bool begin=true)
    {
        if (begin)
            start();
    }

    void start()
    {
        gettimeofday(&start_time, NULL);
    }

    float time()
    {
        timeval result, stop;
        gettimeofday(&stop, NULL);
        timersub(&stop, &start_time, &result);

        return result.tv_sec + result.tv_usec/1000000.0;
    }

    timeval start_time;
};



} // namespace spidir

#endif // ARGHMM_COMMON_H
