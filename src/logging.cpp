/*=============================================================================

  Matt Rasmussen

  Logging functions

=============================================================================*/

// arghmm headers
#include "logging.h"

namespace arghmm {

//=============================================================================
// Errors and Logging

Logger g_logger(stderr, LOG_QUIET);


void Logger::printTimerLog(const Timer &timer, int level, const char *fmt, ...)
{
    if (level <= loglevel) {
        va_list ap;
        va_start(ap, fmt);
        printTimerLog(timer, level, fmt, ap);
        va_end(ap);
    }
}


void Logger::printTimerLog(const Timer &timer, int level, const char *fmt, 
                           va_list ap)
{
    if (level <= loglevel) {
        // print message
        vfprintf(logstream, fmt, ap);

        // print time
        float time = timer.time();
        if (time < 1e-3)
            fprintf(logstream, " %5.1f us", time * 1e6);
        else if (time < 1.0)
            fprintf(logstream, " %5.1f ms", time * 1e3);
        else if (time < 60.0)
            fprintf(logstream, " %5.1f s", time);
        else if (time < 3600.0)
            fprintf(logstream, " %5.1f m", time / 60.0);
        else
            fprintf(logstream, " %5.1f h", time / 3600.0);
        
        fprintf(logstream, "\n");
        fflush(logstream);
    }
}


void printLog(int level, const char *fmt, ...)
{
    if (g_logger.isLogLevel(level)) {
        va_list ap;
        va_start(ap, fmt);
        g_logger.printLog(level, fmt, ap);
        va_end(ap);
    }
}


void printTimerLog(const Timer &timer, int level, const char *fmt, ...)
{
    if (g_logger.isLogLevel(level)) {
        va_list ap;
        va_start(ap, fmt);
        g_logger.printTimerLog(timer, level, fmt, ap);
        va_end(ap);
    }
}


void printError(const char *fmt, va_list ap)
{
    fprintf(stderr, "error: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
}


void printError(const char *fmt, ...)
{
    va_list ap;   
    va_start(ap, fmt);
    fprintf(stderr, "error: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
}


//=============================================================================
// C interface

extern "C" {

void setLogLevel(int level)
{ return g_logger.setLogLevel(level); }

} // extern "C"
    
} // namespace spidir
