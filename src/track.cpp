

#include "track.h"


namespace arghmm {


template <>
bool read_track_line<double>(const char *line, RegionValue<double> &region)
{
    const int CHROM_NAME_LENGTH = 1000;
    char chrom2[CHROM_NAME_LENGTH+1];

    if (sscanf(line, "%1000s\t%d\t%d\t%lf", 
               chrom2, &region.start, &region.end, &region.value) != 4)
        return false;
    region.chrom = string(chrom2);
    return true;
}


template <>
bool read_track_line<int>(const char *line, RegionValue<int> &region)
{
    const int CHROM_NAME_LENGTH = 1000;
    char chrom2[CHROM_NAME_LENGTH+1];

    if (sscanf(line, "%1000s\t%d\t%d\t%d", 
               chrom2, &region.start, &region.end, &region.value) != 4)
        return false;
    region.chrom = string(chrom2);
    return true;
}


} // namespace arghmm


