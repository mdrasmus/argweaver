/*=============================================================================

  Matt Rasmussen

  Logging functions

=============================================================================*/

// c++ headers
#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>


// arghmm headers
#include "logging.h"

namespace arghmm {

//=============================================================================
// Errors and Logging

// stream for logging
static FILE *g_logstream = stderr;
static int g_loglevel = LOG_QUIET;


void printError(const char *fmt, ...)
{
    va_list ap;   
    va_start(ap, fmt);
   
    fprintf(stderr, "error: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
   
    va_end(ap);
}


void printLog(int level, const char *fmt, ...)
{
    if (level <= g_loglevel) {
        va_list ap;   
        va_start(ap, fmt);
        vfprintf(g_logstream, fmt, ap);
        va_end(ap);
	fflush(g_logstream);
    }
}


bool openLogFile(const char *filename)
{
    FILE *stream = fopen(filename, "w");
    
    if (stream != NULL) {
        openLogFile(stream);
        return true;
    } else {
        return false;
    }
}

void openLogFile(FILE *stream)
{
    g_logstream = stream;
}


void setLogLevel(int level)
{
    g_loglevel = level;
}

bool isLogLevel(int level)
{
    return level <= g_loglevel;
}

void closeLogFile()
{
    fclose(g_logstream);
}

FILE *getLogFile()
{
    return g_logstream;
}



} // namespace spidir
