#ifndef ARGWEAVER_TABIX_H
#define ARGWEAVER_TABIX_H

#include <stdio.h>
#include <string>

#include "logging.h"

namespace argweaver {

using namespace std;

FILE *read_tabix(const char *filename, const char *region,
                 const char *tabix_dir);
int close_tabix(FILE *stream);

class TabixStream
{
public:
    TabixStream(const char *filename, const char *region,
                const char *tabix_dir=NULL)
    {
        stream = read_tabix(filename, region, tabix_dir);
        if (stream == NULL) {
            printError("Error opening %s, region=%s\n",
                        filename, region == NULL ? "NULL" : region);
        }
    }

    TabixStream(string filename, const char *region, string tabix_dir) {
        stream = read_tabix(filename.c_str(), region,
                            tabix_dir.empty() ? NULL : tabix_dir.c_str());
        if (stream == NULL) {
            printError("Error opening %s, region=%s\n",
                       filename.c_str(), region == NULL ? "NULL" : region);
        }
    }

    ~TabixStream()
    {
        close();
    }

    void close()
    {
        if (stream) {
            close_tabix(stream);
            stream = NULL;
        }
    }
    FILE *stream;
};

} //namespace argweaver
#endif // ARGWEAVER_TABIX_H
