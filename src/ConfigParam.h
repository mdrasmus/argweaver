/*=============================================================================

  Matt Rasmussen
  Copyright 2010-2011

  Configuration option parsing code

=============================================================================*/

#ifndef SPIDIR_CONFIG_PARAM_H
#define SPIDIR_CONFIG_PARAM_H

#include <assert.h>
#include <libgen.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>


using namespace std;

namespace arghmm {

enum {
    OPTION_ARG,
    OPTION_COMMENT
};


class ConfigParamBase
{
public:
    ConfigParamBase(string shortarg, string longarg, string argstr,
                    string help="", int debug=0) :
        kind(OPTION_ARG),
        shortarg(shortarg),
        longarg(longarg),
        argstr(argstr),
        help(help),
        debug(debug)
    {
    }

    virtual ~ConfigParamBase()
    {}

    virtual int parse(int argc, const char **argv)
    {
        return -1;
    }

    int kind;
    string shortarg;
    string longarg;
    string argstr;
    string help;
    int debug;
};


class ConfigParamComment : public ConfigParamBase
{
public:
    ConfigParamComment(string msg, int debug=0) :
        ConfigParamBase("", "", "", "", debug),
        msg(msg)
    {
        kind = OPTION_COMMENT;
    }

    virtual ~ConfigParamComment()
    {}

    virtual int parse(int argc, const char **argv)
    {
        return -1;
    }

    string msg;
};


template <class T>
class ConfigParam : public ConfigParamBase
{
public:
    ConfigParam(string shortarg, string longarg, string argstr,
                T *value, string help, int debug=0) :
        ConfigParamBase(shortarg, longarg, argstr, help, debug),
        value(value),
        hasDefault(false)
    {
    }

    ConfigParam(string shortarg, string longarg, string argstr,
                T *value, T defaultValue, string help, int debug=0) :
        ConfigParamBase(shortarg, longarg, argstr, help, debug),
        value(value),
        defaultValue(defaultValue),
        hasDefault(true)
    {
        *value = defaultValue;
    }

    virtual int parse(int argc, const char **argv)
    {
        if (argc > 0) {
            *value = T(argv[0]);
            return 1;
        } else {
            fprintf(stderr, "error: argument value expected\n");
            return -1;
        }
    }

    T *value;
    T defaultValue;
    bool hasDefault;
};


class ConfigSwitch : public ConfigParamBase
{
public:
    ConfigSwitch(string shortarg, string longarg,
                bool *value, string help, int debug=0) :
        ConfigParamBase(shortarg, longarg, "", help, debug),
        value(value)
    {
        *value = false;
    }

    virtual int parse(int argc, const char **argv)
    {
        *value = true;
        return 0;
    }

    bool *value;
};


template <>
int ConfigParam<int>::parse(int argc, const char **argv)
{
    if (argc > 0) {
        if (sscanf(argv[0], "%d", value) != 1) {
            fprintf(stderr, "error: int expected '%s'\n", argv[0]);
            return -1;
        }
        return 1;
    } else {
        fprintf(stderr, "error: argument value expected\n");
        return -1;
    }
}


template <>
int ConfigParam<float>::parse(int argc, const char **argv)
{
    if (argc > 0) {
        if (sscanf(argv[0], "%f", value) != 1) {
            fprintf(stderr, "error: float expected '%s'\n", argv[0]);
            return -1;
        }
        return 1;
    } else {
        fprintf(stderr, "error: argument value expected\n");
        return -1;
    }
}


template <>
int ConfigParam<double>::parse(int argc, const char **argv)
{
    if (argc > 0) {
        if (sscanf(argv[0], "%lf", value) != 1) {
            fprintf(stderr, "error: float expected '%s'\n", argv[0]);
            return -1;
        }
        return 1;
    } else {
        fprintf(stderr, "error: argument value expected\n");
        return -1;
    }
}


//=============================================================================

class ConfigParser
{
public:
    ConfigParser()
    {}

    ~ConfigParser()
    {
        clear();
    }


    bool parse(int argc, const char **argv)
    {
        assert(argc > 0);

        char *argdup = strdup(argv[0]);
        prog = basename(argdup);
        free(argdup);

        if (argc < 2)
            return false;

        int i;
        for (i=1; i<argc; i++) {
            bool parsed = false;

            // detect stop parsing options
            if (!strcmp(argv[i], "--")) {
                i++;
                break;
            }

            for (unsigned int j=0; j<rules.size(); j++) {
                if (strcmp(argv[i], "") &&
                    (argv[i] == rules[j]->shortarg ||
                     argv[i] == rules[j]->longarg))
                {
                    int consume = rules[j]->parse(argc - i-1, &argv[i+1]);
                    if (consume == -1)
                        return false;
                    i += consume;
                    parsed = true;
                    break;
                }
            }

            // report error if no arguments are matched
            if (!parsed) {
                if (argv[i][0] == '-') {
                    fprintf(stderr, "unknown option '%s'\n", argv[i]);
                    return false;
                } else {
                    // no more options
                    break;
                }
            }
        }


        // remaining arguments are "rest" arguments
        for (; i<argc; i++)
            rest.push_back(argv[i]);

        return true;
    }

    void printHelp(FILE *stream=stderr, int debug=0)
    {
        fprintf(stream, "Usage: %s [OPTION]\n\n", prog.c_str());

        for (unsigned int i=0; i<rules.size(); i++) {

            // skip rules that are for debug only if debug mode is not enabled
            if (rules[i]->debug > debug)
                continue;

            if (rules[i]->kind == OPTION_ARG) {
                fprintf(stream, "  ");
                if (rules[i]->shortarg != "")
                    fprintf(stream, "%s,", rules[i]->shortarg.c_str());
                fprintf(stream, "%s  %s\n    %s\n\n",
                        rules[i]->longarg.c_str(),
                        rules[i]->argstr.c_str(), rules[i]->help.c_str());
            } else if (rules[i]->kind == OPTION_COMMENT)
                fprintf(stream, "%s\n",
                        ((ConfigParamComment*) rules[i])->msg.c_str());
            else
                assert(0);
        }
    }


    // add rule
    void add(ConfigParamBase *rule)
    {
        rules.push_back(rule);
    }


    // remove all rules
    void clear()
    {
        for (unsigned int i=0; i<rules.size(); i++)
            delete rules[i];
    }

    string prog;
    vector<ConfigParamBase*> rules;
    vector<string> rest;
};


} // namespace spidir

#endif // SPIDIR_CONFIG_PARAM_H
