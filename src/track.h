#ifndef ARGHMM_CONFIG_PARAM_H
#define ARGHMM_CONFIG_PARAM_H


#include <string.h>
#include <string>
#include <vector>

#include "logging.h"
#include "parsing.h"


namespace arghmm {

using namespace std;


class Region
{
public:
    Region(string chrom="", int start=0, int end=0) :
        chrom(chrom), start(start), end(end)
    {}

    void set(const string &_chrom, int _start, int _end) {
        chrom = _chrom;
        start = _start;
        end = _end;
    }

    int length() const {
        return end - start;
    }

    string chrom;
    int start;
    int end;
};

template <class T>
class RegionValue {
public:
    RegionValue() :
        chrom(""), start(0), end(0)
    {}

    RegionValue(string chrom, int start, int end, T value) :
        chrom(chrom), start(start), end(end), value(value)
    {}

    int length() const {
        return end - start;
    }

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


    bool read_track_line(const char *line)
    {
        string chrom;
        int start, end;
        T value;
        
        if (!read_track_line(line, chrom, start, end, value))
            return false;
        append(chrom, start, end, value);
        return true;
    }

    int start_coord() const {
        if (Track<T>::size() == 0)
            return -1;
        else
            return Track<T>::front().start;
    }

    int end_coord() const {
        if (Track<T>::size() == 0)
            return -1;
        else
            return Track<T>::back().start;
    }
};


template <class T>
bool read_track_line(const char *line, RegionValue<T> &region);


template <class T>
class TrackReader {
public:
    TrackReader() : 
        has_error(false),
        line(NULL),
        linesize(1024)
    {}
    ~TrackReader() {
        if (line)
            delete [] line;
    }

    bool open(const char *filename) {
        FILE *infile;
        if ((infile = fopen(filename, "r")) == NULL) {
            printError("cannot read file '%s'", filename);
            return false;
        }
        
        has_error = false;
        
        return true;
    }
    
    bool open(FILE *_infile) {
        infile = _infile;
        has_error = false;
        return true;
    }
    
    bool next(RegionValue<T> &region) {
        if (fgetline(&line, &linesize, infile)) {
            if (!read_track_line(line, region)) {
                // error reading line
                has_error = true;
                return false;
            }
        } else {
            // no more lines in file
            return false;
        }
        
        // line has been successfully read
        return true;
    }

    bool error() const {
        return has_error;
    }


protected:
    bool has_error;
    FILE *infile;
    char *line;
    int linesize;
};





template <class T>
bool read_track(FILE *infile, Track<T> *track)
{    
    char *line = NULL;
    int linesize = 1024;
    int lineno = 0;

    RegionValue<T> region; 

    while (fgetline(&line, &linesize, infile)) {
        lineno++;
        if (!read_track_line(line, region)) {
            printError("could not read track line %d", lineno);
            delete [] line;
            return false;
        }
        track->push_back(region);
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
