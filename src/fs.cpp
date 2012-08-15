
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <string.h>



// ensure a directory exists
static bool ensure_mkdir(const char *path, mode_t mode)
{
    struct stat st;

    if (stat(path, &st) != 0) {
        // directory does not exist. EEXIST for race condition
        if (mkdir(path, mode) != 0 && errno != EEXIST)
            return false;
    } else if (!S_ISDIR(st.st_mode))
        return false;

    return true;
}


// make a path of directories
bool makedirs(const char *path, mode_t mode)
{
    char *pp;
    char *sp;
    bool status = true;
    char *copypath = strdup(path);
    
    pp = copypath;
    while (status && (sp = strchr(pp, '/')) != 0) {
        if (sp != pp) {
            // Neither root nor double slash in path */
            *sp = '\0';
            status = ensure_mkdir(copypath, mode);
            *sp = '/';
        }
        pp = sp + 1;
    }
    if (status)
        status = ensure_mkdir(path, mode);

    free(copypath);

    return status;
}

