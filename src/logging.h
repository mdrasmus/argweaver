/*=============================================================================

  Matt Rasmussen

  Logging functions

=============================================================================*/


#ifndef ARGHMM_LOGGING_H
#define ARGHMM_LOGGING_H

// c/c++ includes
#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
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


inline float timeval2seconds(const timeval &t)
{
    return t.tv_sec + t.tv_usec/1000000.0;
}


// timing
class Timer
{
public:
    Timer(bool begin=true)
    {
        if (begin)
            start();
    }

    // clear total time, same as start(true)
    void clear()
    {
        timerclear(&total_time);
        gettimeofday(&start_time, NULL);
    }

    // start timer, optional reset total time
    void start(bool clear=true)
    {
        if (clear)
            timerclear(&total_time);
        gettimeofday(&start_time, NULL);
    }

    // pause timer
    float pause()
    {
        timeval result, stop;
        gettimeofday(&stop, NULL);
        timersub(&stop, &start_time, &result);
        timeradd(&result, &total_time, &total_time);
        return timeval2seconds(total_time);
    }

    // resume timer, same as start(false)
    void resume()
    {
        start(false);
    }
    
    // return total time so far
    float time() const
    {
        timeval result, stop;
        gettimeofday(&stop, NULL);
        timersub(&stop, &start_time, &result);
        timeradd(&result, &total_time, &result);
        return timeval2seconds(result);
    }

    timeval start_time;
    timeval total_time;
};


class Logger 
{
public:
     Logger(FILE *stream, int level):
        logstream(stream),
        loglevel(level)
    {}

    
    void printLog(int level, const char *fmt, ...)
    {
        if (level <= loglevel) {
            va_list ap;   
            va_start(ap, fmt);
            vfprintf(logstream, fmt, ap);
            va_end(ap);
            fflush(logstream);
        }
    }

    void printLog(int level, const char *fmt, va_list ap)
    {
        if (level <= loglevel) {
            vfprintf(logstream, fmt, ap);
            fflush(logstream);
        }
    }

    void setLogLevel(int level)
    {
        loglevel = level;
    }

    bool isLogLevel(int level) const
    {
        return level <= loglevel;
    }

    void printTimerLog(const Timer &timer, int level, const char *fmt, 
                       va_list ap);
    void printTimerLog(const Timer &timer, int level, const char *fmt, ...);


    //============================
    // log stream

    bool openLogFile(const char *filename)
    {
        FILE *stream = fopen(filename, "w");
        
        if (stream != NULL) {
            logstream = stream;
            return true;
        } else {
            return false;
        }
    }

    void openLogFile(FILE *stream)
    {
        logstream = stream;
    }

    void closeLogFile()
    {
        fclose(logstream);
    }

    FILE *getLogFile() const
    {
        return logstream;
    }


protected:

    FILE *logstream;
    int loglevel;
};


//=============================================================================
// global logging functions

extern Logger g_logger;

inline bool openLogFile(const char *filename)
{ return g_logger.openLogFile(filename); }

inline void openLogFile(FILE *stream)
{ return g_logger.openLogFile(stream); }

inline void closeLogFile()
{ g_logger.closeLogFile(); }

inline FILE *getLogFile()
{ return g_logger.getLogFile(); }

inline void setLogLevel(int level)
{ return g_logger.setLogLevel(level); }

inline bool isLogLevel(int level)
{ return g_logger.isLogLevel(level);
}


// global function API
void printLog(int level, const char *fmt, ...);
void printTimerLog(const Timer &timer, int level, const char *fmt, ...);

void printError(const char *fmt, va_list ap);
void printError(const char *fmt, ...);






} // namespace arghmm

#endif // ARGHMM_COMMON_H
