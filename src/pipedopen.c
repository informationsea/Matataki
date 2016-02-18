#include "pipedopen.h"

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define DEBUG(fmt,...) (fprintf(stderr, __FILE__ ": %3d: " fmt "\n" ,__LINE__, ## __VA_ARGS__))

PFILE *pipedopen(const char *filename, bool writemode)
{
    const char *pipeprog;
    const char *mode = writemode ? "w" : "r";
    const char *option;


    // setup compressor/decompressor
    if (strcmp(filename + strlen(filename) - 3, ".gz") == 0) {
        pipeprog = "gzip";
    } else if (strcmp(filename + strlen(filename) - 4, ".bz2") == 0) {
        pipeprog = "bzip2";
    } else if (strcmp(filename + strlen(filename) - 3, ".xz") == 0) {
        pipeprog = "xz";
    } else if (strcmp(filename + strlen(filename) - 5, ".lzop") == 0) {
        pipeprog = "lzop";
    } else {
        FILE *fp = fopen(filename, mode);
        if (fp == NULL) return NULL;
        PFILE *f = (PFILE *)malloc(sizeof(PFILE));
        if (f == NULL) return NULL;
        f->file = fp;
        f->fd = -1;
        return f;
    }

    if (writemode) {
        option = "-c";
    } else {
        option = "-dc";
    }

    // open file
    int fd;
    if (writemode) {
        fd = open(filename, O_WRONLY|O_CREAT, 0644);
    } else {
        fd = open(filename, O_RDONLY);
    }

    if (fd < 0) {
        return NULL;
    }

    
    PFILE *f = (PFILE *)malloc(sizeof(PFILE));
    bzero(f, sizeof(PFILE));
    if (f == NULL) return NULL; // failed to allocate
    f->fd = fd;
    

    char *args[] = {strdup(pipeprog), strdup(option), NULL};
    child_process_t p;

    if (writemode)
        p = spawn(pipeprog, args, SPAWN_PIPE_STDIN|SPAWN_REDIRECT_STDOUT, true, -1, fd, -1);
    else
        p = spawn(pipeprog, args, SPAWN_REDIRECT_STDIN|SPAWN_PIPE_STDOUT, true, fd, -1, -1);
    
    if (p.pid == 0) {
        free(f);
        return NULL; // failed to spawn
    }

    FILE *fp;
    if (writemode)
        fp = fdopen(p.fd_stdin, mode);
    else
        fp = fdopen(p.fd_stdout, mode);
        
    if (fp == NULL) {
        free(f);
        return NULL; // failed to open file
    }
    
    f->file = fp;
    return f;
}


int pipedclose(PFILE *file)
{
    fclose(file->file);
    if (file->process.pid != 0)
        waitpid(file->process.pid, NULL, 0);
    if (file->fd >= 0)
        close(file->fd);
    free(file);
    return 0;
}

