
#include "compress.h"
#include <string>

namespace arghmm {

using namespace std;

string quote_arg(string text)
{
    int j = 0; 
    char text2[text.size() * 4 + 3];
    text2[j++] = '\'';

    for (unsigned int i=0; i<text.size(); i++) {
        if (text[i] == '\'') {
            text2[j++] = '\'';
            text2[j++] = '\\';
            text2[j++] = '\'';
            text2[j++] = '\'';
        } else
            text2[j++] = text[i];
    }

    // terminate string
    text2[j++] = '\'';
    text2[j++] = '\0';
    
    return string(text2);
}


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


} // namespace arghmm


