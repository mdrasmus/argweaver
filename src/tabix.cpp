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
    FILE *pipe;

    if (region == NULL) {
        return read_compress(filename);
    }

    string cmd = "tabix -h " + quote_arg(filename) + " " +  region;
    if (tabix_dir != NULL)
        cmd = string(tabix_dir) + "/" + cmd;
    pipe = popen(cmd.c_str(), "r");
    if (pipe == NULL) {
        fprintf(stderr, "Error opening %s with tabix. Is tabix installed and"
                " in your PATH?\n",
                filename);
        exit(1);
    }
    return pipe;
}


int close_tabix(FILE *stream) {
    return pclose(stream);
}

}
