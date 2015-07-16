#ifndef ARGWEAVER_COMPRESS_H
#define ARGWEAVER_COMPRESS_H

#include <string.h>
#include <stdio.h>

namespace argweaver {


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
    CompressStream(const char *filename, const char *mode="r",
                   const char *command=NULL)
    {
        int len = strlen(filename);
        compress = false;

        if (len > 3 && strcmp(&filename[len - 3], ".gz") == 0) {
            compress = true;
            if (mode[0] == 'r')
                stream = read_compress(filename, command);
            else
                stream = write_compress(filename, command);
        } else {
            stream = fopen(filename, mode);
        }
    }

    CompressStream(FILE *stream, bool compress=true) :
        compress(compress),
        stream(stream)
    {}

    ~CompressStream()
    {
        close();
    }

    void close()
    {
        if (stream) {
            if (compress)
                close_compress(stream);
            else
                fclose(stream);
            stream = NULL;
        }
    }

    bool compress;
    FILE *stream;
};


} // namespace argweaver

#endif // ARGWEAVER_COMPRESS_H
