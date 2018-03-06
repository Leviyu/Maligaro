
#include "ESF.h"
/******************************************************************
 * This is a c script to output the big info file
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
int output_final_pick_onset_info_for_record(new_RECORD* my_record, FILE* out)
{
	// output format
	// 1		2					3					4	
	// sta		time_phase_peak		npts_phase_peak		stretch_factor
	
	fprintf(out,"%7s %10.5lf %10d %5.2lf\n",my_record->name, my_record->time_phase_peak, my_record->npts_phase_peak, my_record->best_stretch_coefficient);

	return 0;
}
int output_final_event_info(new_RECORD* my_record, new_INPUT* my_input)
{
	printf("-> output final eventinfo \n");
	int output_final_event_info_for_record(new_RECORD* my_record,  FILE* out, new_INPUT* my_input);
	int ista;	
	int phase_flag;
	char phase[20];
	strcpy(phase,my_input->PHASE);
	phase_flag = 0;

	FILE* out;
	//FILE* out2;
	char out_name[200];
	char out_eventpick[200];
	sprintf(out_name,"eventinfo.%s.%s.%s",my_input->EQ, phase, my_input->COMP);
	out=fopen(out_name,"w");
	//sprintf(out_eventpick,"eventpick.%s.%s.%s",my_input->EQ, phase, my_input->COMP);
	//out2=fopen(out_eventpick,"w");
	

	for(ista = 0; ista< my_input->sta_num;ista++)
	{
		// output record info for each record
		//printf("%d / %d  %s \n ",ista, my_input->sta_num , my_record[ista].name);
		output_final_event_info_for_record(&my_record[ista], out, my_input);
		//output_final_pick_onset_info_for_record(&my_record[ista],out2);
	}
	fclose(out);
	//fclose(out2);

	return 0;
}

int output_final_event_info_for_record(new_RECORD* my_record, FILE* out, new_INPUT* my_input)
{

// output format
// 1 	2	3	4	5	6		7		8		9		10		11	12		13			14			   15   	16  			17  18 	 	19			   20 		21
// sta NET dist AZ BAZ STA_lat STA_lon  EQ_lat EQ_lon EQ_dep EQ_mag EQ_name polar_flag quality_flag   prem  phase_amplitude    ccc  SNR  dt_obs_PREM     phase_name best_stretch_ccc
// 22							23		24 					25					26 			27       28
// best_stretch_coefficient 	misfit 	COMP	time_phase_peak		npts_phase_peak		noise_beg 	noise_len

	
	double SNR = my_record->SNR;
	double SNR2 = my_record->SNR2;
	double best_tstar = my_record->best_tstar;
	double best_tstar_ccc = my_record->best_tstar_ccc;
	double best_ccc = my_record->best_stretch_ccc;
	double best_coefficient = my_record->best_stretch_coefficient;
	double ccc = my_record->ccc;

	if(ccc > 1 || ccc < -1 )
		ccc = 0;
	if(best_ccc > 1 || best_ccc < -1)
		ccc = 0;


	//printf("STA %s SNR %lf CCC %lf SNR_CUT %lf CCC_Cut %lf quality %d \n", my_record->name, my_record->SNR, my_record->best_stretch_ccc, my_input->SNR_CUT, my_input->CCC_CUT, my_record->quality);

	char sta[20];
	strcpy(sta,my_record->name);
	char NET[20];
	strcpy(NET,my_record->NET);
	double 	DIST = my_record->DIST;
	double AZ = my_record->AZ;
	double BAZ = my_record->BAZ;
	double STA_lat = my_record->sta_lat;
	double STA_lon = my_record->sta_lon;
	double EQ_lat = my_record->eq_lat;
	double EQ_lon = my_record->eq_lon;
	double EQ_dep = my_record->eq_dep;
	double EQ_mag = my_record->eq_mag;
double polarity = my_record->polarity;
	char EQ_name[20];
	strcpy(EQ_name,my_record->EQ);

	char phase_name[20];
	strcpy(phase_name,my_record->PHASE);
	int polar_flag = my_record->polar_flag;
	int quality_flag = my_record->quality;
	if (my_input->Reprocessing_Flag == 1)
		quality_flag = 1;
	double prem = my_record->prem;
	double amp = my_record->amplitude;
	double misfit = my_record->misfit;
	char COMP[20];
	strcpy(COMP,my_record->COMP);
	double time_phase_peak= my_record->time_phase_peak;
	int npts_phase_peak = my_record->npts_phase_peak;
	double noise_beg = my_record->noise_beg;
	double noise_len = my_record->noise_len;
	double phase_beg_time_relative_to_prem = my_record->phase_beg;
	double record_weight = my_record->weight;
	double misfit2 = my_record->misfit2;
	double misfit_pre = my_record->misfit_pre;
	double misfit_pre2T = my_record->misfit_pre2T;
	double misfit_pre3T = my_record->misfit_pre3T;
	double misfit_bak = my_record->misfit_bak;
	double misfit_bak2T = my_record->misfit_bak2T;
	double misfit_bak3T = my_record->misfit_bak3T;
	double record_gaussian_factor = my_record->record_gaussian_factor;
	double emp_gaussian_factor = my_record->emp_gaussian_factor;
	double gaussian_misfit = my_record->gaussian_misfit;
	int traffic_phase_nearby = my_record->traffic_phase_nearby;
	int polar_correct_flag = my_record->polar_correct_flag;

	double ONSET = my_record->ONSET;
	double ENDSET = my_record->ENDSET;
	double ccc3 = my_record->CCC3;
	double SNR3 = my_record->SNR3;
	double SNR4 = my_record->SNR4;
	int noise_too_short_flag = my_record->noise_too_short_flag;
	//printf("ccc2 \n");



	// for Lp Sdiff SS SSS data, we add delay time to the dt_residual
	//double dt_obs_prem = my_record->dt_obs_prem - my_input->bp_delay_time;
	double dt_obs_prem = my_record->dt_obs_prem;
	// for postprocessing
	// hardwire dt to be dt_picked shift
	if(my_input->Reprocessing_Flag == 1) 
		dt_obs_prem = my_record->dt_picked_shift;


	fprintf(out,"%6s %6s %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %13s %2d %2d %8.3lf %17.14lf %5.2lf %5.2lf %6.1lf %6s %5.2lf %5.2lf %5.2lf %5s %5.2lf %6d %5.2lf %5.2lf %8.3lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %6d %6d %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %3d\n",
			sta, NET, DIST, AZ, BAZ, STA_lat, STA_lon, EQ_lat, EQ_lon, EQ_dep, EQ_mag, EQ_name, polar_flag, quality_flag, prem, amp, ccc ,SNR, dt_obs_prem, phase_name,
			best_ccc, best_coefficient, misfit, COMP, time_phase_peak, npts_phase_peak, noise_beg, noise_len, phase_beg_time_relative_to_prem, record_weight,SNR2, misfit2,
			ONSET, ENDSET, best_tstar, best_tstar_ccc, ccc3, misfit_pre, misfit_bak, record_gaussian_factor, emp_gaussian_factor, gaussian_misfit, polarity, polar_correct_flag, 
			traffic_phase_nearby, misfit_pre2T, misfit_pre3T,misfit_bak2T, misfit_bak3T,SNR3,SNR4,noise_too_short_flag);
	return 0;
}



