
#include <string.h>

#include "parsing.h"

namespace argweaver {

const int DEFAULT_LINE_SIZE = 4 * 1024;


// similar to libc readline except line is allocated with new/delete
int fgetline(char **lineptr, int *linesize, FILE *stream)
{
    char *line;
    int readsize = 0;

    if (*lineptr == NULL) {
        // allocate line buffer if it is not given
        line = new char [*linesize];
    } else{
        // use existing line buffer
        line = *lineptr;
    }

    char *start = line;
    int capacity = *linesize;

    while (true) {
        // place null char in last position
        char *end = &start[capacity - 2];
        *end = '\0';

        // read line
        if (!fgets(start, capacity, stream)) {
            // ensure eol is present
            start[0] = '\0';
        }

        // check to see if line is longer than buffer
        if (*end != '\0' && *end != '\n') {
            readsize += capacity - 1;

            // resize buffer
            int newsize = *linesize * 2;
            char *tmp = new char [newsize];

            // failed to get memory
            if (!tmp)
                return -1;

            // copy over line read so far
            memcpy(tmp, line, *linesize);

            // update sizes
            capacity = newsize - *linesize;
            *linesize = newsize;

            // update pointers
            delete [] line;
            line = tmp;
            start = &line[readsize];

            continue;
        }
        break;
    }

    // line is shorter than buffer, return it

    // determine how much was read
    readsize += strlen(start);

    *lineptr = line;
    return readsize;
}


// read a line from a file
// return NULL pointer on eof
char *fgetline(FILE *stream)
{
    int linesize = DEFAULT_LINE_SIZE;
    char *line = NULL;
    fgetline(&line, &linesize, stream);

    // detect eof
    if (line[0] == '\0' && feof(stream)) {
        delete [] line;
        return NULL;
    }

    return line;
}


void split(const char *str, const char *delim, vector<string> &tokens)
{
    tokens.clear();
    int i = 0;
    bool end = false;

    while (!end) {
        // walk to end of next token
        int j = i;
        for (; str[j] && !inChars(str[j], delim); j++);
        end = !str[j];

        // save token
        tokens.push_back(string(&str[i], j-i));
        i = j + 1;
    }
}



// concatenate multiple strings into one newly allocated string
char *concat_strs(char **strs, int nstrs)
{
    int total_len = 0;
    int lens[nstrs];

    // count total length
    for (int i=0; i<nstrs; i++) {
        lens[i] = strlen(strs[i]);
        total_len += lens[i];
    }

    char *str = new char [total_len + 1];

    // copy strings
    char *p = str;
    for (int i=0; i<nstrs; i++) {
        memcpy(p, strs[i], lens[i]);
        p += lens[i];
    }
    str[total_len] = '\0';

    return str;
}


// Quote an argument for use in a command shell
// Argument is wrapped in single quotes and every single quote is escaped.
string quote_arg(string text)
{
    int j = 0;
    char text2[text.size() * 4 + 3];
    text2[j++] = '\'';

    for (unsigned int i=0; i<text.size(); i++) {
        if (text[i] == '\'') {
            text2[j++] = '\'';
            text2[j++] = '\\';
            text2[j++] = '\'';
            text2[j++] = '\'';
        } else
            text2[j++] = text[i];
    }

    // terminate string
    text2[j++] = '\'';
    text2[j++] = '\0';

    return string(text2);
}


} // namespace argweaver
