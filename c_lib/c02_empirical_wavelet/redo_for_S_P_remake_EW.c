#include "ESF.h"



void redo_for_S_P_remake_EW(new_RECORD* my_record, new_INPUT* my_input)
{

	if( strcmp(my_input->PHASE, "S") == 0 || strcmp(my_input->PHASE, "P") == 0)
		printf("---> working on redo_for_S_P_remake_EW \n");
	else
		return ;

	int npts_phase = (int)(my_input->phase_len / my_input->delta);
	double EW[npts_phase];
	double stretched_EW[npts_phase];
	int ista;
	int count;

	
	// initiate
	for(count = 0; count < npts_phase ; count++ )
	{
		EW[count] = 0;
		stretched_EW[count] = 0;
	}

	for(ista = 0; ista < my_input->sta_num ; ista++)
	{
		if(my_record[ista].quality <= 0)
			continue;
		for(count = 0; count < npts_phase ; count++)
			EW[count] += my_record[ista].phase_win[count] * my_record[ista].weight;
		// 1. redefine the beyond_window_flag
		my_record[ista].beyong_window_flag = 1;
	}


	// 0. use all good records to stack new EW
	int loop_num;
	int loop_num_max = 3;
	for(loop_num = 2; loop_num <= loop_num_max ; loop_num ++)
	{
		printf(" loop %d / %d \n", loop_num , loop_num_max);

		my_input->iteration_flag = loop_num;
		empirical_source_for_each_record(my_record,my_input,EW,loop_num);
		get_ONSET_ENDSET_for_each_record_origional_phase(my_record,my_input,EW);
		restack_ES_for_phase_code_choice(my_record,my_input,EW);
	}

	store_ES_into_record(my_record,my_input,EW);
	find_best_match_gaussian_for_iterative_ES(my_record,my_input,EW);
	output_STD_of_first_ES(my_record,my_input, EW);
	output_current_ES_for_phase(my_input, EW);

	
	// 1. stretched records to EW
	strcpy(my_input->stretch_flag, "stretch");
	stretch_ES_and_CCC(my_record,my_input,EW);


	// 2. stretch records to restack
	stretch_record_restack_ES_based_on_code_choice(my_record,my_input,stretched_EW);
	output_STD_of_second_ES(my_record,my_input, stretched_EW);
	output_current_ES_for_phase_second(my_input, stretched_EW);


	//3. tstar
	strcpy(my_input->stretch_flag,"tstar");
	stretch_ES_and_CCC_tstar(my_record, my_input, stretched_EW);

	get_ONSET_ENDSET_for_each_record_stretched(my_record,my_input);
	define_goodness_of_record(my_record, my_input);
	output_STD_of_third_ES(my_record,my_input, stretched_EW);
	output_current_ES_for_phase_third(my_input, stretched_EW);


	define_stretch_EW_ONSET( my_record, my_input);
	output_ES_for_each_record(my_record, my_input);
	

	redefine_beyon_wind_flag(my_record,my_input,stretched_EW,stretched_EW);

}




