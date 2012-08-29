

#include "track.h"


namespace arghmm {

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

} // namespace arghmm


