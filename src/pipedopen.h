#ifndef PIPEDOPEN_H
#define PIPEDOPEN_H

#include <stdio.h>
#include <stdbool.h>
#include "spawn.h"

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct PFILE_tag {
        FILE *file;
        child_process_t process;
        int fd;
    } PFILE;
    
    PFILE *pipedopen(const char *filename, bool writemode);
    int pipedclose(PFILE *file);

#ifdef __cplusplus
}
#endif

#endif /* PIPEDOPEN_H */
