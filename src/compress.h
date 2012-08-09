#ifndef ARGHMM_COMPRESS_H
#define ARGHMM_COMPRESS_H


#include <stdio.h>

namespace arghmm {


// default zip/unzip commands
#define ZIP_COMMAND "gzip -"
#define UNZIP_COMMAND "gunzip -"


FILE *read_compress(const char *filename, const char *command=UNZIP_COMMAND);

FILE *write_compress(const char *filename, const char *command=ZIP_COMMAND);

FILE *open_compress(const char *filename, const char *mode, 
                    const char *command=ZIP_COMMAND);

int close_compress(FILE *stream);


} // namespace arghmm

#endif // ARGHMM_COMPRESS_H
