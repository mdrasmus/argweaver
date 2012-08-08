
#include "compress.h"
#include <string>

namespace arghmm {

using namespace std;

string quote_arg(string text)
{
    char text2[text.size() * 2 + 1];

    int j = 0;
    for (unsigned int i=0; i<text.size(); i++) {
        if (text[i] == '\'')
            text2[j++] = '\\';
        text2[j++] = text[i];
    }

    // terminate string
    text2[j++] = '\0';

    return string(text2);
}


FILE *open_compress(const char *filename, const char *mode, 
                    const char *command)
{
    string cmd = string(command) + " > '" + quote_arg(filename) + "'";
    return popen(cmd.c_str(), mode);
}


int close_compress(FILE *stream)
{
    return pclose(stream);
}


} // namespace arghmm


