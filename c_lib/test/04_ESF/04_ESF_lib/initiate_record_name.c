#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<sacio.h>
#include<ESF.h>

int initiate_record_name(new_RECORD* my_record)
{
	my_record->EQ = (char*)malloc(sizeof(char)*20);
	my_record->COMP = (char*)malloc(sizeof(char)*10);
	my_record->PHASE = (char*)malloc(sizeof(char)*20);
	my_record->name = (char*)malloc(sizeof(char)*20);
	my_record->NET = (char*)malloc(sizeof(char)*10);

	return 0;
}
