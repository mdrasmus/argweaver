#ifndef ARGHMM_CONFIG_PARAM_H
#define ARGHMM_CONFIG_PARAM_H


#include <string.h>
#include <string>
#include <vector>

#include "compress.h"
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


template <>
bool Track<double>::read_track_line(const char *line)
{
    const int CHROM_NAME_LENGTH = 1000;
    char chrom[CHROM_NAME_LENGTH+1];
    int start, end;
    double value;

    if (sscanf("%1000s\t%d\t%d\t%lf\n", chrom, &start, &end, &value) != 4)
        return false;
    append(chrom, start, end, value);
    return true;
}


template <>
bool Track<int>::read_track_line(const char *line)
{
    const int CHROM_NAME_LENGTH = 1000;
    char chrom[CHROM_NAME_LENGTH+1];
    int start, end;
    int value;

    if (sscanf("%1000s\t%d\t%d\t%d\n", chrom, &start, &end, &value) != 4)
        return false;
    append(chrom, start, end, value);
    return true;
}


template <class T>
bool read_track(const char *filename, Track<T> *track)
{
    CompressStream stream(filename);
    if (!stream.stream) {
        printError("cannot read file '%s'", filename);
        return false;
    }
    
    char *line = NULL;
    int linesize = 1024;
    int lineno = 0;
    while (fgetline(&line, &linesize, stream.stream)) {
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



} // namespace arghmm

#endif // ARGHMM_TRACK
