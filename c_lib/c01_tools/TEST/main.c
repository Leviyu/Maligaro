#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hongyulib.h>
#include <ESF.h>




int main()
{

char command[100];

double prem;
sprintf(command,"cat dd|awk 'NR==3 {print $2}'");
shell_pipe_double(command,&prem);

printf("prem is %lf \n", prem);



	return 0;
}
