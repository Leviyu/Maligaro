#include "hongyulib.h"
#include "ESF.h"




int get_SNR3_and_4_for_record(double* phase_win,int phase_npts, double* noise_win,int noise_npts, double* SNR3, double* SNR4)
{
	int iphase, inoise;

	// 1. find the peak of noise/phase
	// int amplitudeloc(double* array, int len, int* max_amp_loc, double*
	//// amplitude, int flag)
	int npts_phase_max;
	int npts_noise_max;
	double noise_amp;
	double phase_amp;
	int flag = 1;
	amplitudeloc(phase_win, phase_npts, &npts_phase_max, &phase_amp, flag);
	amplitudeloc(noise_win, noise_npts, &npts_noise_max, &noise_amp, flag);

//printf("phase amp %lf noise amp %lf \n", )
	if(noise_amp == 0)
		noise_amp = 1;
	*SNR4 = fabs( phase_amp / noise_amp);

//printf("SNR4 %lf \n", *SNR4);

	
	// now we calculate the peak-to-trough amplitude
	int npts_phase_trough = 0;
	double phase_trough_amp = 0;
	for(iphase = npts_phase_max+1; iphase < phase_npts -1 ; iphase ++)
	{
		if( phase_win[iphase+1] - phase_win[iphase] > 0 )
		{
			npts_phase_trough = iphase;
			phase_trough_amp = fabs(phase_amp - phase_win[iphase]);
			break;
		}
	}

	int npts_noise_trough = 0;
	double noise_trough_amp = 0;
	for(inoise = npts_noise_max+1; inoise < noise_npts -1 ; inoise ++)
	{
		if( noise_win[inoise+1] - noise_win[inoise] > 0 )
		{
			npts_noise_trough = inoise;
			noise_trough_amp = fabs(noise_amp - noise_win[iphase]);
			break;
		}
	}
	
	if(noise_trough_amp == 0)
		noise_trough_amp = 1;
	*SNR3 = fabs(phase_trough_amp / noise_trough_amp);


	return 0;
}
