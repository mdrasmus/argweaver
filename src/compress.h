#ifndef ARGHMM_COMPRESS_H
#define ARGHMM_COMPRESS_H


#include <stdio.h>

namespace arghmm {

FILE *open_compress(const char *filename, const char *mode, 
                    const char *command="gzip -");

int close_compress(FILE *stream);


} // namespace arghmm

#endif // ARGHMM_COMPRESS_H
