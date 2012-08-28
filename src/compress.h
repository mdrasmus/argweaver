#ifndef ARGHMM_COMPRESS_H
#define ARGHMM_COMPRESS_H

#include <string.h>
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


class CompressStream
{
public:
    CompressStream(const char *filename)
    {
        int len = strlen(filename);
        compress = false;
        
        if (len > 3 && strcmp(&filename[len - 3], ".gz") == 0) {
            compress = true;
            stream = read_compress(filename);
        } else {
            stream = fopen(filename, "r");
        }
    }

    CompressStream(FILE *stream, bool compress=true) :
        compress(compress),
        stream(stream)
    {}

    ~CompressStream()
    {
        if (stream) {
            if (compress)
                close_compress(stream);
            else
                fclose(stream);
        }
    }

    bool compress;
    FILE *stream;
};


} // namespace arghmm

#endif // ARGHMM_COMPRESS_H
