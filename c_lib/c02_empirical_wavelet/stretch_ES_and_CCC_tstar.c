
#include "hongyulib.h"
#include "ESF.h"
// stretch ES to find the best waveform to match each record
int stretch_ES_and_CCC_tstar(new_RECORD* my_record, new_INPUT* my_input, double* current_ES)
{
	printf( "---> stretch_ES_and_CCC for tstar \n");
	int npts_phase = (int)(my_input->phase_len / my_input->delta);
	int ista,j;
	double stretched_ES[npts_phase];
	double best_ccc ;
	int best_time_shift;
	double best_coeff ;
	double coeff_min;
	double coeff_max;
	double coeff_delta;
	double coeff_array[50];
	int max_num = 4;
	int count;


	for(ista=0;ista<my_input->sta_num;ista++)
	{
		//if( my_record[ista].quality == -1)
			//continue;
		if(strcmp(my_input->stretch_flag,"tstar")==0 && my_record[ista].best_stretch_coefficient != 0 && my_record[ista].best_stretch_coefficient <1 )
		{
			coeff_min = 0.3;
			coeff_max = 1.5;
			coeff_array[0] = 0.2;
			coeff_array[1] = 0.1;
			coeff_array[2] = 0.05;
			coeff_array[3] = 0.02;
			coeff_array[4] = 0.01;
			//coeff_array[3] = 0.1;
			//coeff_array[4] = 0.0;
			//coeff_array[5] = 0.0;
			max_num = 5;
		}
		else
		{
			coeff_min = 0.1;
			coeff_max = 20.1;
			coeff_array[0] = 3;
			//##coeff_array[0] = 1.2;
			coeff_array[1] = 1.5;
			coeff_array[2] = 0.7;
			coeff_array[3] = 0.35;
			coeff_array[4] = 0.2;
			coeff_array[5] = 0.1;
			coeff_array[6] = 0.05;
			coeff_array[7] = 0.02;
			coeff_array[8] = 0.01;
			max_num = 9;
		}
		//printf(" Working on sta %d  %s \n ", ista, my_record[ista].name);	
		//output_array1("current_ES_tmp", current_ES,npts_phase);
		//output_array1("current_phase_win_tmp", my_record[ista].phase_win,npts_phase);

		stretch_ES_find_best_match_for_given_interval(&my_record[ista], current_ES, my_record[ista].phase_win, npts_phase, coeff_min, coeff_max, coeff_array[0], &best_ccc, &best_coeff,&best_time_shift,stretched_ES, my_input);

		//char kkk[500];
		////sprintf(kkk,"t_star.%s",my_record[ista].name);
		//output_array1(kkk,my_record[ista].phase_win,npts_phase);

		for(count = 1; count < max_num  ; count ++)
		{
			coeff_min = best_coeff - coeff_array[count-1]/2;
			coeff_max = best_coeff + coeff_array[count-1]/2;
			coeff_delta = coeff_array[count];
			stretch_ES_find_best_match_for_given_interval(&my_record[ista], current_ES, my_record[ista].phase_win, npts_phase, coeff_min, coeff_max, coeff_delta, &best_ccc, &best_coeff,&best_time_shift,stretched_ES, my_input);
			if(best_ccc > 0.96)
				break;
		}

		double shift_time = best_time_shift*my_input->delta;
		double shift_time_max = 1000;
		my_record[ista].phase_beg -= shift_time;
/*
			if( fabs(shift_time) > shift_time_max)
				shift_time = 0;
			if(my_record[ista].phase_beg < my_record[ista].long_beg) 
				my_record[ista].phase_beg += shift_time;

		int max_loc;
		double max_amp;
		amplitudeloc(current_ES, npts_phase, &max_loc, &max_amp,1);
		double EW_onset_relative_to_PREM = max_loc * my_input->delta + my_input->phase_beg - 5;
		


		my_record[ista].shift_time_recorder += shift_time;
		printf(" sta %s shift %lf total shift %lf EW_onset_relative_to_PREM %lf  \n", my_record[ista].name,
				shift_time, my_record[ista].shift_time_recorder, EW_onset_relative_to_PREM);
		double max_time = 30;
		double dt_shift = -1*my_record[ista].shift_time_recorder + EW_onset_relative_to_PREM;
		if( dt_shift > max_time || dt_shift < -1*max_time)
		{
			printf("%s sta shift to much \n", my_record[ista].name);
			my_record[ista].shift_time_recorder -= shift_time;
			my_record[ista].phase_beg += shift_time;
			my_record[ista].quality = -1;
		}


		*/

		// update phase win
		read_phase_window(&my_record[ista],my_input);

		my_record[ista].best_tstar_ccc = best_ccc;
		my_record[ista].best_tstar= best_coeff;

		for(j=0;j<npts_phase;j++)
			my_record[ista].stretched_ES_win[j] = stretched_ES[j];

		//printf("For ista %d stretching best ccc %lf best coeff %lf \n", ista, best_ccc, best_coeff);
	}



	return 0;
}

