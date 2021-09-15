#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <string.h>

int main(int ac,char** av){
	//av={"-8,-7, ... ,7,8","par=mtinv.par","gmt5","mtdegfree=5","use_snr","minsnr=3","shift", "ctol=0.85","maxshift=10"};
	//time counting
    struct timeval time1,time2;
    gettimeofday(&time1,0);
    char* ts0;
    int pid=0;
    int wait_pid;
    char ts0Arg[6];
	//separate "-8,-7, ... ,7,8"
    while((ts0 = strsep(&av[1],",")) != NULL){
        pid=fork();
        if(pid<0){
            fprintf(stderr, "fork error");
            fflush(stderr);
        }
        //child process
        else if(pid==0){
            sprintf(ts0Arg,"ts0=%s",ts0);
            execlp("/opt/mtinv.v3.0.5/bin/mtinv","mtinv",ts0Arg,av[2],av[3],av[4],av[5],av[6],av[7], av[8],av[9],(char *)0);
        }
    }
	//retrieve child process
    while( (wait_pid = waitpid( -1, NULL, 0 ))  )
    {
        if( errno == ECHILD ) {
            break;
        }
    }
    gettimeofday(&time2,0);
    double timeuse  = (1000000.0*(time2.tv_sec - time1.tv_sec) + time2.tv_usec - time1.tv_usec)/1000;
    fprintf(stderr,"multi process mtinv finished,total cost time:%fMS\n\n",timeuse);
}

