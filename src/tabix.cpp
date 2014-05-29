#include <unistd.h>
#include <string>
#include <stdlib.h>

#include "tabix.h"
#include "parsing.h"
#include "compress.h"

namespace argweaver {

using namespace std;


FILE *read_tabix(const char *filename, const char *region,
                 const char *tabix_dir) {
    //    bool exists = !access(filename, F_OK);
    //    char index_file[strlen(filename)+5];
    FILE *pipe;
    /*    if (!exists) {
        fprintf(stderr, "%s not found\n", filename);
        return NULL;
        }*/

    if (region == NULL) {
        return read_compress(filename);
    }

    /*    sprintf(index_file, "%s.tbi", filename);
          exists = !access(index_file, F_OK);
    if (!exists) {
        fprintf(stderr, "%s does not have index file. Try running tabix on your file first.\n",
                filename);
        return NULL;
        }*/
    string cmd = "tabix -h " + quote_arg(filename) + " " +  region;
    if (tabix_dir != NULL)
        cmd = string(tabix_dir) + "/" + cmd;
    pipe = popen(cmd.c_str(), "r");
    if (pipe == NULL) {
        fprintf(stderr, "Error opening %s with tabix. Is tabix installed and in your PATH?\n",
                filename);
        exit(1);
    }
    return pipe;
}


int close_tabix(FILE *stream) {
    return pclose(stream);
}

}
