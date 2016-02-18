#include "spawn.h"
#include <stdlib.h>
#include <stdio.h>

child_process_t spawn(const char *file, char *argv[], int pipe_flags, bool usepath, int stdin_fd, int stdout_fd, int stderr_fd)
{
    int stdin_pipe[2];
    int stdout_pipe[2];
    int stderr_pipe[2];
    child_process_t process;
    bzero(&process, sizeof(process));
    process.fd_stdin = -1;
    process.fd_stdout = -1;
    process.fd_stderr = -1;

    if (pipe_flags & SPAWN_PIPE_STDIN) {
        pipe(stdin_pipe);
    }
    
    if (pipe_flags & SPAWN_PIPE_STDOUT) {
        pipe(stdout_pipe);
    }
    if (pipe_flags & SPAWN_PIPE_STDERR) {
        pipe(stderr_pipe);
    }
    
    process.pid = fork();
    if (process.pid < 0) // failed to fork
        return process;

    if (process.pid > 0) { // parent process
        if (pipe_flags & SPAWN_PIPE_STDIN) {
            process.fd_stdin = stdin_pipe[1];
            close(stdin_pipe[0]);
        }
        
        if (pipe_flags & SPAWN_PIPE_STDOUT) {
            process.fd_stdout = stdout_pipe[0];
            close(stdout_pipe[1]);
        }
        
        if (pipe_flags & SPAWN_PIPE_STDERR) {
            process.fd_stderr = stderr_pipe[0];
            close(stderr_pipe[1]);
        }
        
        return process;
    }

    // child process

    if (pipe_flags & SPAWN_PIPE_STDIN) {
        close(STDIN_FILENO);
        dup2(stdin_pipe[0], STDIN_FILENO);
        close(stdin_pipe[1]);
    } else if (pipe_flags & SPAWN_REDIRECT_STDIN) {
        dup2(stdin_fd, STDIN_FILENO);
    }

    if (pipe_flags & SPAWN_PIPE_STDOUT) {
        close(STDOUT_FILENO);
        dup2(stdout_pipe[1], STDOUT_FILENO);
        close(stdout_pipe[0]);
    } else if (pipe_flags & SPAWN_REDIRECT_STDOUT) {
        dup2(stdout_fd, STDOUT_FILENO);
    }

    if (pipe_flags & SPAWN_PIPE_STDERR) {
        close(STDERR_FILENO);
        dup2(stderr_pipe[1], STDERR_FILENO);
        close(stderr_pipe[0]);
    } else if (pipe_flags & SPAWN_REDIRECT_STDERR) {
        dup2(stderr_fd, STDERR_FILENO);
    }

    if (usepath)
        execvp(file, argv);
    execv(file, argv);
    perror("Cannot start subprocess");
    exit(1);
}

child_process_t spawn_basic(const char *file, char *argv[], bool usepath)
{
    return spawn(file, argv, 0, usepath, 0, 0, 0);
}

child_process_t spawn_pipe(const char *file, char *argv[], int pipe_flags, bool usepath)
{
    return spawn(file, argv, pipe_flags, usepath, 0, 0, 0);
}

