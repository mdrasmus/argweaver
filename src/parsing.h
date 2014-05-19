#ifndef ARGWEAVER_PARSING_H
#define ARGWEAVER_PARSING_H


// headers c++
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>

namespace argweaver {

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
inline char *trim(char *word)
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
    end[1] = '\0';

    // pure white space case
    if (end < start)
        return &end[1];

    return start;
}


// Returns string with whitespace removed from beginning and end of word
inline char *trim(string &word)
{
    char *str = (char*) word.c_str();
    char *start = str;
    while (*start == ' ' ||
           *start == '\t' ||
           *start == '\n' ||
           *start == '\r')
        start++;

    char *end = (char*) &str[strlen(str) - 1];
    while (end > start && (*end == ' ' ||
                           *end == '\t' ||
                           *end == '\n' ||
                           *end == '\r'))
        end--;
    word[end - start + 1] = '\0';

    // pure white space case
    if (end < start)
        return &end[1];

    return start;
}


void split(const char *str, const char *delim, vector<string> &tokens);

char *concat_strs(char **strs, int nstrs);

string quote_arg(string text);


} // namespace argweaver

#endif // ARGWEAVER_PARSING_H

