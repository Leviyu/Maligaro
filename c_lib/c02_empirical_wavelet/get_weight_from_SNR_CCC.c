
#include "ESF.h"
double get_weight_from_SNR_CCC(double SNR, double ccc_orig)
{
	double weight = 0;
	double ccc;
	double ccc_flag = 0;
	double SNR_flag = 0;
	if( fabs(ccc_orig) > 0.95 )
		ccc_flag = 1;
	else if( fabs(ccc_orig) > 0.90 ) 
		ccc_flag = 0.7;
	else if( fabs(ccc_orig) > 0.85 ) 
		ccc_flag = 0.5;
	else if( fabs(ccc_orig) > 0.8 ) 
		ccc_flag = 0.4;
	else if( fabs(ccc_orig) > 0.6 ) 
		ccc_flag = 0.2;
	else 
		ccc_flag = 0;



	if( SNR > 6 )
		SNR_flag = 1;
	else if ( SNR > 4)
		SNR_flag = 0.7;
	else if ( SNR > 3)
		SNR_flag = 0.5;
	else if ( SNR > 2.2)
		SNR_flag = 0.2;
	else
		SNR_flag = 0;

	weight = ccc_flag * SNR_flag;

//printf(" SNR %lf ccc %lf weight %lf \n", SNR,ccc_orig,weight);

	return weight;
}
