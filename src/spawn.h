#ifndef SPWAN_H
#define SPWAN_H

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <string.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct child_process_tag
{
    int fd_stdin;
    int fd_stdout;
    int fd_stderr;
    pid_t pid;
} child_process_t;

#define SPAWN_PIPE_STDIN  (1 << 0)
#define SPAWN_PIPE_STDOUT (1 << 1)
#define SPAWN_PIPE_STDERR (1 << 2)
#define SPAWN_REDIRECT_STDIN  (1 << 3)
#define SPAWN_REDIRECT_STDOUT (1 << 4)
#define SPAWN_REDIRECT_STDERR (1 << 5)

child_process_t spawn_basic(const char *file, char *argv[], bool usepath);
child_process_t spawn_pipe(const char *file, char *argv[], int pipe_flags, bool usepath);
child_process_t spawn(const char *file, char *argv[], int pipe_flags, bool usepath, int stdin_fd, int stdout_fd, int stderr_fd);

#ifdef __cplusplus
}
#endif

#endif /* SPWAN_H */
