#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hongyulib.h>




int main()
{

	int file_flag = file_exist("dd");
	printf("exist filag is %s \n", file_flag);

	return 0;
}
