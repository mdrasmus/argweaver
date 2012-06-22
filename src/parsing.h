#ifndef ARGHMM_PARSING_H
#define ARGHMM_PARSING_H


// headers c++ 
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>

namespace arghmm {

using namespace std;

int fgetline(char **lineptr, int *linesize, FILE *stream);
char *fgetline(FILE *stream);

// Returns true if c is in chars
inline bool inChars(char c, const char *chars)
{
   if (!chars)
      return false;
   for (;*chars; chars++)
      if (c == *chars) return true;
   return false;
}


// remove possible newlines '\n', '\r' from end of line
inline bool chomp(char *str)
{
   int len = strlen(str);
   if (len > 0 && str[len-1] == '\n') {
      str[len-1] = '\0';
      if (len > 1 && str[len-2] == '\r')
          str[len-2] = '\0';
      return true;
   } else
      return false;
}

// Returns string with whitespace removed from beginning and end of word
inline string trim(const char *word)
{
    char *start = (char*) word;
    while (*start == ' ' ||
           *start == '\t' ||
           *start == '\n' ||
           *start == '\r')
        start++;

    char *end = (char*) &word[strlen(word) - 1];
    while (end > start && (*end == ' ' ||
                           *end == '\t' ||
                           *end == '\n' ||
                           *end == '\r'))
        end--;

    // pure white space case
    if (end < start)
        return string("");

    return string(start, end - start + 1);
}


vector<string> split(const char *str, const char *delim, 
                     bool multiDelim = true);


} // namespace arghmm

#endif // ARGHMM_PARSING_H

