#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<sacio.h>
#include<ESF.h>
int read_eventstation_list(new_RECORD* my_record, new_INPUT* my_input)
{
fprintf(my_input->out_logfile,"---> Read eventstation list info \n");
	FILE* in;
	in=fopen(my_input->event_file,"r");
	char buff[200];
	char buff3[200];
	char buff1[200];
	char sta_tmp[200];
	double polarity;
	int N = 200;
	int i,j,k;

	for(i=0;i<my_input->sta_num;i++)
	{
		initiate_record_name(&my_record[i]);
		fgets(buff,N,in);
		sscanf(buff,"%s %s %lf %s %lf %s %lf %s %lf %lf %lf %lf %lf %s %s %lf %s %s %s",
				my_record[i].name,my_record[i].NET,&my_record[i].DIST,buff1,&my_record[i].AZ,buff1,&my_record[i].BAZ,buff1,
				&my_record[i].sta_lat,&my_record[i].sta_lon,&my_record[i].eq_lat,&my_record[i].eq_lon,&my_record[i].eq_dep,
				buff1,buff1,&my_record[i].eq_mag,buff1,buff1,my_record[i].EQ);

		strcpy(my_record[i].PHASE, my_input->PHASE);
		
		//int DEBUG = 1;
		//if( DEBUG == 0 )
		//{
		 ////output for DUBUG
		//printf("%s %s %lf %s %lf %s %lf %s %lf %lf %lf %lf %lf %s %s %lf %s %s %s\n",
				//my_record[i].name,my_record[i].NET,my_record[i].DIST,buff1,my_record[i].AZ,buff1,my_record[i].BAZ,buff1,
				//my_record[i].sta_lat,my_record[i].sta_lon,my_record[i].eq_lat,my_record[i].eq_lon,my_record[i].eq_dep,
				//buff1,buff1,my_record[i].eq_mag,buff1,buff1,my_record[i].EQ);
		//}
	}

	fclose(in);


fprintf(my_input->out_logfile,"---> Read polarity  info \n");

	// read in polarity info
	FILE* in2;
	sprintf(my_input->polarity_file,"eventinfo.polarity.%s.%s",my_input->PHASE,my_input->COMP);
	in2=fopen(my_input->polarity_file,"r");
	int ista;
	// read in polarity file info
	for(i = 0; i< my_input->sta_num ; i++)
	{
		if(!in2) 
		{
			fprintf(my_input->out_logfile,"polarity file empty");
			break;
		}
		fgets(buff3,N,in2);
		sscanf(buff3, "%s %lf ", sta_tmp, &my_record[i].polarity);
		for(ista =0; ista<my_input->sta_num;ista ++)
		{
			if(strcmp(sta_tmp,my_record[ista].name)==0)
			{
			// set polar_flag according to the polarity info
			if( my_record[i].polarity > 0.15)
				my_record[i].polar_flag = 1;
			else if( my_record[i].polarity < -0.15)
				my_record[i].polar_flag = -1;
			else 
				my_record[i].polar_flag = 0;
			break;
			}
			else
				continue;
		}

	}
	//fclose(in2);

	// for SS we make the polarity flag all 0
	if(strcmp(my_input->PHASE, "SS") ==0 || 
			strcmp(my_input->PHASE, "SSSS") ==0||
			strcmp(my_input->PHASE, "SSSSm") ==0)
	{
		for(ista =0; ista<my_input->sta_num;ista ++)
			my_record[ista].polar_flag = 0;
	}

fprintf(my_input->out_logfile,"NOTE: polarity for SS S4 S4m are all set to NULL since we did a hilbert transform \n");

	return 0;
}
