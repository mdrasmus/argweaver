
#include <unistd.h>
#include <string>

#include "compress.h"
#include "parsing.h"


namespace argweaver {

using namespace std;

FILE *read_compress(const char *filename, const char *command)
{
    bool exists = !access(filename, F_OK);
    if (!exists)
        return NULL;
    const char *command2 = (command ? command : UNZIP_COMMAND);
    string cmd = string(command2) + " < " + quote_arg(filename);
    return popen(cmd.c_str(), "r");
}


FILE *write_compress(const char *filename, const char *command)
{
    // TODO: add check to prevent write error
    const char *command2 = (command ? command : ZIP_COMMAND);
    string cmd = string(command2) + " > " + quote_arg(filename);
    return popen(cmd.c_str(), "w");
}


FILE *open_compress(const char *filename, const char *mode,
                    const char *command)
{
    string cmd;
    if (mode[0] == 'w')
        cmd = string(command) + " > " + quote_arg(filename);
    else if (mode[0] == 'r')
        cmd = string(command) + " < " + quote_arg(filename);

    return popen(cmd.c_str(), mode);
}


int close_compress(FILE *stream)
{
    return pclose(stream);
}


} // namespace argweaver


