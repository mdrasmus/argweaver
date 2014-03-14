#ifndef ARGWEAVER_TABIX_H
#define ARGWEAVER_TABIX_H

#include <string.h>
#include <stdio.h>

namespace argweaver {

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
	  fprintf(stderr, "Error opening %s, region=%s\n",
		  filename, region == NULL ? "NULL" : region);
	  exit(1);
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
