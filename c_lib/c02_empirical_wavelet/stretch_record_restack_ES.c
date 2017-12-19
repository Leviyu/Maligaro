
#include "hongyulib.h"
#include "ESF.h"
/******************************************************************
 * stretch each record using stretch_facotr 
 * and then restack for new ES
 *
 *	Input:
 *
 *
 *	Output:
 *
 *
 *	DATE:				Keywords:
 *	Reference:
******************************************************************/

int stretch_record_restack_ES(new_RECORD* my_record, new_INPUT* my_input, double* ES)
{
	printf("---> stretch_record_restack_ES \n");
	int npts_phase = (int)( my_input->phase_len / my_input->delta  );
	double ES_tmp[npts_phase];
	int ista,i;
	int num_records_make_new_ES=0;
	double new_ES[npts_phase];
	double ccc_use = 1;

	double stretch_factor_record;
	char stretch_record_out_name[1000];
	double x_phase[npts_phase];

	// initiate ES
	for(i=0;i<npts_phase;i++)
		ES[i] = 0;
	


	for(ista=0;ista< my_input->sta_num; ista ++)
	{
		if(my_record[ista].quality == -1)
			continue;
		if(my_record[ista].traffic_phase_nearby == 1 )
			continue;
		if(fabs(my_record[ista].SNR) < my_input->SNR_CUT)
			continue;
		if(fabs(my_record[ista].ccc) < 0.8 )
			continue;



		// hongyu change
		// here we grab records from eventinfo.final
		// and we only use records that we picked already to construct the S E.W.
		if( strcmp( my_record[ista].PHASE, "S" )  == 0 
			&& my_input->ED_CHECKED_FLAG == 1 && 
			my_record[ista].checked_to_be_good == 1)
		{
			double kkk = 1;
		}
		else 
		{
			continue;
		}

		//printf(" %d th record weight is %lf \n", ista, my_record[ista].weight);
		if( my_record[ista].best_stretch_coefficient == 0)
			my_record[ista].best_stretch_coefficient = 1;
		stretch_factor_record = 1/my_record[ista].best_stretch_coefficient;
		//stretch phase win using stretch factor
		stretch_ES_function(my_record[ista].phase_win, npts_phase, stretch_factor_record, my_record[ista].stretched_phase_win);
		my_record[ista].stretched_phase_win_flag = 1;

		if(strstr(my_record[ista].name, "PPP") != NULL)
			continue;

		// output the stretched record
		//for(i=0;i<npts_phase;i++)
			//x_phase[i] = my_record[ista].phase_beg + i*my_input->delta;
		//sprintf(stretch_record_out_name,"%s.%s.%s.%s.phase_stretch", my_record[ista].EQ,my_record[ista].name,my_record[ista].PHASE,my_record[ista].COMP);
		//output_array2(stretch_record_out_name,x_phase,my_record[ista].stretched_phase_win, npts_phase, 0);

		if(my_record[ista].weight == 0)
			continue;
		for(i=0;i<npts_phase;i++)
		{
			if(my_record[ista].stretched_phase_win[i] != my_record[ista].stretched_phase_win[i] )
				my_record[ista].stretched_phase_win[i] = 0;

			//ES_tmp[i] += stretch_tmp[i]*my_record[ista].weight;
			ES[i] += my_record[ista].stretched_phase_win[i]*my_record[ista].weight;
		}
		num_records_make_new_ES++;
	}
	printf( "---> Restack stretch_record_restack_ES  E.W. Done %d records is included to remake the ES \n",num_records_make_new_ES);


	if(num_records_make_new_ES == 0)
	{
		printf(" if 0 record is needed to restack E.W. then we use everytrace with the same weight to stack for one \n");

		for(ista=0;ista< my_input->sta_num; ista ++)
		{
			//my_record[ista].weight = get_weight_from_SNR_CCC( my_record[ista].SNR, my_record[ista].ccc);
			//my_record[ista].weight = get_weight_from_SNR_CCC( my_record[ista].SNR, ccc_use);
			//if(my_record[ista].quality == -1)
				//continue;

			// if current record is not stretched yet
			if(my_record[ista].stretched_phase_win_flag != 1)
			{
				if(my_record[ista].best_stretch_coefficient == 0)
					my_record[ista].best_stretch_coefficient = 1;
				stretch_factor_record = 1/my_record[ista].best_stretch_coefficient;
				stretch_ES_function(my_record[ista].phase_win, npts_phase, stretch_factor_record, my_record[ista].stretched_phase_win);
			}

			for(i=0;i<npts_phase;i++)
			{
				if(my_record[ista].stretched_phase_win[i] != my_record[ista].stretched_phase_win[i] )
					my_record[ista].stretched_phase_win[i] = 0;

				//ES_tmp[i] += stretch_tmp[i]*my_record[ista].weight;
				ES[i] += my_record[ista].stretched_phase_win[i];
			}
			num_records_make_new_ES++;

		}
		printf(" ---> Restack E.W. Done %d records is included to remake the ES \n",num_records_make_new_ES);
	}

	//##output_array1("stretch_record_stacked_ES.tmp",ES, npts_phase);

	// we shrink the E.W. by 0.5 to make sure that we can match the skinny record
	//
	// int stretch_ES_function(double* current_ES, int npts_phase, double coeff,double*  tmp_ES)
	//double stretch_factor = 1.0;
	//stretch_ES_function(ES_tmp, npts_phase, stretch_factor,ES);

	return 0;
}
