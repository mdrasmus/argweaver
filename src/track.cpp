

#include "track.h"


namespace arghmm {


template <>
bool read_track_line<double>(const char *line, 
                     string &chrom, int &start, int &end, double &value)
{
    const int CHROM_NAME_LENGTH = 1000;
    char chrom2[CHROM_NAME_LENGTH+1];

    if (sscanf("%1000s\t%d\t%d\t%lf\n", chrom2, &start, &end, &value) != 4)
        return false;
    chrom = string(chrom2);
    return true;
}


template <>
bool read_track_line<int>(const char *line, 
                     string &chrom, int &start, int &end, int &value)
{
    const int CHROM_NAME_LENGTH = 1000;
    char chrom2[CHROM_NAME_LENGTH+1];

    if (sscanf("%1000s\t%d\t%d\t%d\n", chrom2, &start, &end, &value) != 4)
        return false;
    chrom = string(chrom2);
    return true;
}



} // namespace arghmm


