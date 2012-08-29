#ifndef ARGHMM_CONFIG_PARAM_H
#define ARGHMM_CONFIG_PARAM_H


#include <string.h>
#include <string>
#include <vector>

#include "logging.h"
#include "parsing.h"


namespace arghmm {

using namespace std;


template <class T>
class RegionValue {
public:
    RegionValue(string chrom, int start, int end, T value) :
        chrom(chrom), start(start), end(end), value(value)
    {}

    string chrom;
    int start;
    int end;
    T value;
};

template <class T>
class Track : public vector<RegionValue<T> > {
public:

    void append(string chrom, int start, int end, T value) {
        push_back(RegionValue<T>(chrom, start, end, value));
    }

    bool read_track_line(const char *line);
};


template <class T>
bool read_track(FILE *infile, Track<T> *track)
{    
    char *line = NULL;
    int linesize = 1024;
    int lineno = 0;
    while (fgetline(&line, &linesize, infile)) {
        lineno++;
        if (!track->read_track_line(line)) {
            printError("could not read track line %d", lineno);
            delete [] line;
            return false;
        }
    }
    
    delete [] line;
    
    return true;
}


template <class T>
bool read_track(const char *filename, Track<T> *track)
{
    FILE *infile;
    if ((infile = fopen(filename, "r")) == NULL) {
        printError("cannot read file '%s'", filename);
        return false;
    }
    
    bool result = read_track(infile, track);
    
    fclose(infile);
    return result;
}



} // namespace arghmm

#endif // ARGHMM_TRACK
