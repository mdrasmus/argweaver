#ifndef ARGWEAVER_TRACK_H
#define ARGWEAVER_TRACK_H

// c++ includes
#include <string.h>
#include <string>
#include <vector>

// arghmm includes
#include "logging.h"
#include "parsing.h"


namespace argweaver {

using namespace std;

// A region within a chromosome
// start is inclusive, end is exclusive
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


// A region within a chromosome associated with a value
// start is inclusive, end is exclusive
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


typedef char NullValue;
typedef RegionValue<char> RegionNullValue;


// A track is a series of regions each associated with a value
template <class T>
class Track : public vector<RegionValue<T> > {
public:


    // Returns start coordinate if regions are available
    // Returns -1 otherwise
    int start_coord() const {
        if (Track<T>::size() == 0)
            return -1;
        else
            return Track<T>::front().start;
    }

    // Returns end coordinate if regions are available
    // Returns -1 otherwise
    int end_coord() const {
        if (Track<T>::size() == 0)
            return -1;
        else
            return Track<T>::back().end;
    }

    // Returns value of region containing position pos
    T find(int pos, const T &default_value) const {
        for (unsigned int i=0; i<Track<T>::size(); i++) {
            const RegionValue<T> &region = Track<T>::at(i);
            if (region.start <= pos && pos < region.end)
                return region.value;
        }
        // region not found
        return default_value;
    }

    // Returns index of region containing position pos
    int index(int pos) const {
        for (unsigned int i=0; i<Track<T>::size(); i++) {
            const RegionValue<T> &region = Track<T>::at(i);
            if (region.start <= pos && pos < region.end)
                return i;
        }
        // region not found
        return -1;
    }

    // Adds one region to the track
    void append(string chrom, int start, int end, T value) {
        this->push_back(RegionValue<T>(chrom, start, end, value));
    }

    // Reads one region from a map file and adds it to the track
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
};


typedef Track<NullValue> TrackNullValue;

// Reads one region from a map file
template <class T>
bool read_track_line(const char *line, RegionValue<T> &region);


// A reader for reading a track from a map file
template <class T>
class TrackReader {
public:
    TrackReader() :
        has_error(false),
        line(NULL),
        linesize(1024),
        lineno(0)
    {}
    ~TrackReader() {
        if (line)
            delete [] line;
    }

    // open a map file by filename
    bool open(const char *filename) {
        FILE *infile;
        if ((infile = fopen(filename, "r")) == NULL) {
            printError("cannot read file '%s'", filename);
            return false;
        }

        lineno = 0;
        has_error = false;

        return true;
    }

    // open a map file by stream
    bool open(FILE *_infile) {
        infile = _infile;
        has_error = false;
        return true;
    }

    // Fetches the next RegionValue from a map file
    // Returns true if region read, false otherwise
    bool next(RegionValue<T> &region) {
        lineno++;
        if (fgetline(&line, &linesize, infile)) {
            // ignore track lines
            //if (!strncmp(line, "track", 5))
            //    continue;

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

    // Returns true if error encountered
    bool error() const {
        return has_error;
    }

    // Returns current line number in map file
    int line_number() const {
        return lineno;
    }


protected:
    bool has_error;
    FILE *infile;
    char *line;
    int linesize;
    int lineno;
};


// Reads a track from a map stream
template <class T>
bool read_track(FILE *infile, Track<T> *track)
{
    TrackReader<T> reader;
    RegionValue<T> region;
    reader.open(infile);

    while (reader.next(region)) {
        track->push_back(region);
    }
    if (reader.error()) {
        printError("could not read track line %d", reader.line_number());
        return false;
    }

    return true;
}


// Reads a track from a map stream
// NOTE: only regions within chrom:start-end are kept
template <class T>
bool read_track_filter(FILE *infile, Track<T> *track,
                       string chrom, int start, int end)
{
    TrackReader<T> reader;
    RegionValue<T> region;
    reader.open(infile);

    while (reader.next(region)) {
        // only keep regions that overlap desired region
        if (region.chrom == chrom &&
            region.end > start && region.start < end) {
            // trim region
            if (region.start < start)
                region.start = start;
            if (region.end > end)
                region.end = end;

            track->push_back(region);
        }
    }
    if (reader.error()) {
        printError("could not read track line %d", reader.line_number());
        return false;
    }

    return true;
}

template <class T>
bool read_track_filter(FILE *infile, Track<T> *track, const Region &region)
{
    return read_track_filter<T>(infile, track,
                                region.chrom, region.start, region.end);
}


// Reads a track from a map file
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

// Reads a track from a map file
// NOTE: only regions with chrom:start-end are kept
template <class T>
bool read_track_filter(const char *filename, Track<T> *track,
                       string chrom, int start, int end)
{
    FILE *infile;
    if ((infile = fopen(filename, "r")) == NULL) {
        printError("cannot read file '%s'", filename);
        return false;
    }

    bool result = read_track_filter(infile, track, chrom, start, end);

    fclose(infile);
    return result;
}



} // namespace argweaver

#endif // ARGWEAVER_TRACK
