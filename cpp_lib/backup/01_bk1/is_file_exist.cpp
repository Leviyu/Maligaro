#include<sys/stat.h>
#include<unistd.h>
#include "hongyulib.h"



bool is_file_exist(string name)
{
	ifstream f(name.c_str());
	return f.good();
}
