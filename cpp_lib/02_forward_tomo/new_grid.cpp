#include "forward_tomography.h"



new_grid::new_grid()
{

	cout << " new grid is declared! " << endl;
	// distribute space for record tag arrat
	int MAX = 10000;
	
	this->record_tag = new int[MAX];
	this->fix_BAZ_time = new double[MAX];
	this->fix_BAZ_slowness = new double[MAX];
	this->fix_BAZ_amp= new double[MAX];
	this->fix_slow_time = new double[MAX];
	this->fix_slow_BAZ = new double[MAX];
	this->npts_record_sum = 0;
	this->long_win = new double[MAX];



}
new_grid::~new_grid()
{

	
	delete[] this->record_tag;
	delete[] this->long_win;
}

void new_grid::initiate_grid()
{
	cout << " --> Build Virtual Station grid network" << endl;

	double ilat, ilon;
	double current_lat;

	// grid radius in km
	this->grid_radius = this->VS_LATITUDE_INC * 110;


	for(ilat = -89; ilat < 89; ilat+= this->VS_LATITUDE_INC)
	{
		current_lat = ilat;
		





	}




}


void new_grid::get_grid_dist()
{


	// get EQ lat lon
	
	double eq_lat = this->my_big_record->my_record[0].eq_lat;
	double eq_lon = this->my_big_record->my_record[0].eq_lon;
	double grid_lat , grid_lon;
	grid_lat = this->grid_lat;
	grid_lon = this->grid_lon;

	double grid_dist = dist_A_B(eq_lat ,eq_lon, grid_lat , grid_lon);
	grid_dist = grid_dist/111;
	this->grid_dist = grid_dist;






}


// This function first read in S_ES file 
// and then t*(S_ES) to find the best fit to virtual stack
// then find the best fit gaussian to t*(S_ES) and set ONSET on gaussian
int new_grid::find_stack_ONSET_time()
{

	// 1 . read in S_ES
	ifstream myfile;
	myfile.open(this->S_ES_file.c_str());
	int LINE = count_file_num(this->S_ES_file);
	int count;
	double X_TMP[LINE];
	double S_ES[LINE];

	cout << " S_ES is "<< this->S_ES_file<< endl;
	cout << " file count is "<< LINE << endl;

	for(count = 0; count < LINE; count++)
	{
		myfile >> X_TMP[count] >> S_ES[count];
	}


	myfile.close();


	// this->long is the stacked window, we cut the window to be the same size of S_ES
	// a find peak
	// int amplitudeloc(double* array, int len, int* max_amp_loc, double* amplitude, int flag)
	int max_loc;
	double amp;


	// we use a smaller window to find the phase max loc
	int long_zoom_beg_npts = (int) ( (fabs(this->LONG_BEG) - 15 ) / this->delta);
	double long_zoom_len = 50;
	int npts_long_zoom = long_zoom_len/this->delta;
	double long_zoom[npts_long_zoom];
	for(count = 0; count < npts_long_zoom ; count++ )
	{
		long_zoom[count] = this->long_win[ long_zoom_beg_npts + count];
	}



	//amplitudeloc( this->long_win, this->long_npts, & max_loc, &amp,1 );
	amplitudeloc( long_zoom  , npts_long_zoom, & max_loc, &amp,1 );
	if( max_loc > 50000 || max_loc < -50000 )
		max_loc = 0;
	if ( long_zoom_beg_npts  > 50000 || long_zoom_beg_npts < -50000 )
		long_zoom_beg_npts = 0;

	max_loc = max_loc + long_zoom_beg_npts;



	double phase_win[LINE];
	double phase_start_time = this->LONG_BEG + (max_loc - LINE/2) * this->delta;
	//double phase_start_time = 0;
	cout << "phase start time is "<< phase_start_time << endl;
	for(count = 0 ; count < LINE ; count++)
	{
		phase_win[count] = this->long_win[ max_loc - LINE/2 + count];
	}

	double xx[LINE];
	string phase_win_file = "phase_win."+to_string(this->ilat)+"."+to_string(this->ilon);
	for(count = 0 ; count < LINE; count++)
		xx[count] = phase_start_time+  count  *this->delta;

	output_array2(phase_win_file,xx,phase_win,LINE,0 );

	



	// 2. t* S_ES to fit virtual stack
	double coeff_min = 0.5;
	double coeff_max = 5;
	double coeff_delta = 0.5;
	double best_ccc;
	double best_coeff;
	int best_time_shift;
	double best_ES[LINE];
	double best_ES_gau[LINE];
	stretch_record_find_best_match_for_given_interval( S_ES,phase_win, LINE,  coeff_min, coeff_max, coeff_delta, &best_ccc, &best_coeff, & best_time_shift, best_ES );
	//stretch_record_find_best_match_for_given_interval( S_ES,this->long_win, LINE,  coeff_min, coeff_max, coeff_delta, &best_ccc, &best_coeff, & best_time_shift, best_ES );

	coeff_min = best_coeff - 0.3;
	coeff_max = best_coeff + 0.3;
	coeff_delta = 0.1;
	stretch_record_find_best_match_for_given_interval( S_ES,phase_win, LINE,  coeff_min, coeff_max, coeff_delta, &best_ccc, &best_coeff, & best_time_shift, best_ES );

	cout << " best ccc "<< best_ccc << " best coeff "<< best_coeff << " best time shift "<< best_time_shift << endl;
	
	// outout t* E_ES
	//best_time_shift = 0;
	string current_tstar_ES = "tstar_ES."+to_string(this->ilat)+"."+to_string(this->ilon);
	for(count = 0 ; count < LINE; count++)
		xx[count] = phase_start_time+ ( count - best_time_shift)  *this->delta;

	output_array2(current_tstar_ES,xx,best_ES,LINE,0 );

	int shift_time_tmp = best_time_shift;

	this->tstar_ccc = best_ccc;
	this->tstar_factor = best_coeff;

	// 3. find best fit gaussian 
	coeff_min = 0.1;
	coeff_max = 5;
	coeff_delta = 0.5;
	stretch_gaussian_find_best_match_for_given_interval( best_ES, LINE,  coeff_min, coeff_max, coeff_delta, &best_ccc, &best_coeff, & best_time_shift, best_ES_gau );
	coeff_min = best_coeff-0.3;
	coeff_max = best_coeff + 0.3;
	coeff_delta = 0.1;
	stretch_gaussian_find_best_match_for_given_interval( best_ES, LINE,  coeff_min, coeff_max, coeff_delta, &best_ccc, &best_coeff, & best_time_shift, best_ES_gau );
	
	cout << "gau best ccc "<< best_ccc << " best coeff "<< best_coeff << " best time shift "<< best_time_shift << endl;

	string current_gau = "gau_ES."+to_string(this->ilat)+"."+to_string(this->ilon);
	//best_time_shift = 0;
	for(count = 0 ; count < LINE; count++)
		xx[count] = phase_start_time+ ( count - best_time_shift - shift_time_tmp)   *this->delta;

	output_array2(current_gau,xx,best_ES_gau,LINE, 0 );

	this->gau_ccc = best_ccc;
	this->gau_factor = best_coeff;

	// store infor

	//find the gaussian function ONSET
	 //amplitudeloc( this->long_win, this->long_npts, & max_loc, &amp,1 );   	
	amplitudeloc( best_ES_gau ,LINE, &max_loc, &amp,1 );
	int ONSET;
	double ONSET_time;
	double threshold = 0.02;
	for(count = max_loc ; count > 0; count --)
	{
		if( best_ES_gau[count] < amp*threshold )
		{
			ONSET = count;
			ONSET_time = xx[count];
			this->virtual_stack_ONSET = ONSET_time;
			break;
		}



	}




}



void new_grid::find_records_within_range()
{
	
	int ista;



	for(ista = 0; ista < this->my_big_record->sta_num; ista ++)
	{
		double grid_lat = this->grid_lat;
		double grid_lon = this->grid_lon;
		double record_lat = this->my_big_record->my_record[ista].sta_lat;
		double record_lon = this->my_big_record->my_record[ista].sta_lon;
		double eq_lat = this->my_big_record->my_record[ista].eq_lat;
		double eq_lon = this->my_big_record->my_record[ista].eq_lon;




		// calculate the distance bewteen grid and record
		double distance;
		distance = dist_A_B( grid_lat, grid_lon, record_lat, record_lon);
		distance = distance/111;

		if(distance < this->grid_radius )
		{
			// add this record into grid
			int num = this->npts_record_sum;
			this->record_tag[ num ] = ista;
			this->npts_record_sum ++;
		}

	}


}

void new_grid::stack_records_from_one_EQ()
{

	int ista;
	int npts;
	//double weight_SUM = 0;


	// option #1 stacked records that are aligned to PHASE PREM time
	for(npts = 0 ; npts < this->long_npts ; npts ++)
	{
		this->long_win[npts] = 0;
		for(ista = 0 ; ista < this->npts_record_sum ; ista++)
		{
			int tag = this->record_tag[ista];
			// calculate the distance between station and grid center
			double weight;

			double dist;
			double sta_lat;
			double sta_lon;
			sta_lon = this->my_big_record->my_record[tag].sta_lon;
			sta_lat = this->my_big_record->my_record[tag].sta_lat;
			dist = dist_A_B( this->grid_lat, this->grid_lon,  sta_lat , sta_lon );
			dist = dist / 111;

			//gaussian_func(double a, double b, double c, double d, double x)
			weight = gaussian_func(1, 0, 4 , 0, dist);


			if(  this->my_big_record->my_record[tag].long_win[npts] !=  this->my_big_record->my_record[tag].long_win[npts] )
				continue;

			int current_record_polar = this->my_big_record->my_record[tag].polarity_flag;
			if(current_record_polar == 0)
				current_record_polar = 1;

			this->long_win[ npts] += this->my_big_record->my_record[tag].long_win[npts] * weight * current_record_polar;
		}

	}
	normalize_array( this->long_win, this->long_npts );



}


void new_grid::get_dt_ave_STD()
{

	double dt_SUM;
	double dt_ave;
	double dt_STD;
	int ista;
	dt_SUM = 0;
	for(ista = 0 ; ista < this->npts_record_sum ; ista++)
	{		
		int tag = this->record_tag[ista];
		double dt = this->my_big_record->my_record[tag].dt_obs_prem;
		dt_SUM += dt;
		cout << "+++ dt "<< dt << endl; 
	}


	dt_ave = dt_SUM / this->npts_record_sum;

	dt_STD = 0;
	for(ista = 0 ; ista < this->npts_record_sum ; ista++)
	{		
		int tag = this->record_tag[ista];
		double dt = this->my_big_record->my_record[tag].dt_obs_prem;

		dt_STD = (dt - dt_ave ) * (dt - dt_ave );
	}

	dt_STD = sqrt(dt_STD);

	this->dt_ave = dt_ave;
	this->dt_STD = dt_STD;

}

void new_grid::get_crust_correction_ave_STD()
{

	double dt_SUM;
	double dt_ave;
	double dt_STD;
	int ista;
	dt_SUM = 0;
	for(ista = 0 ; ista < this->npts_record_sum ; ista++)
	{		
		int tag = this->record_tag[ista];
		double dt = this->my_big_record->my_record[tag].crust_correction;
		dt_SUM += dt;
		//cout << "+++ dt "<< dt << endl; 
	}


	dt_ave = dt_SUM / this->npts_record_sum;

	dt_STD = 0;
	for(ista = 0 ; ista < this->npts_record_sum ; ista++)
	{		
		int tag = this->record_tag[ista];
		double dt = this->my_big_record->my_record[tag].crust_correction;

		dt_STD = (dt - dt_ave ) * (dt - dt_ave );
	}

	dt_STD = sqrt(dt_STD);

	this->crustal_correction_ave = dt_ave;
	this->crustal_correction_STD = dt_STD;

}

void new_grid::output_stacked_record()
{


	ofstream myfile;
	myfile.open(this->out_stacked_record_rel_PREM.c_str());

	double X[100000];

	int npts;

	normalize_array(this->long_win, this->long_npts);



	for(npts = 0; npts < this->long_npts; npts++)
	{
		X[npts] = this->LONG_BEG + npts * this->delta;
		myfile << X[npts] << "  " << this->long_win[npts] << endl;

	}


	myfile.close();
}


// 
// This code relocated the center of grid based on stations within range
// if stations within range is less or equal then 3, we skip
//
void new_grid::relocate_grid_center()
{
	
	if ( this->npts_record_sum <= this->station_num_threshold)
	{
		double skip;
		skip = 0;
	}
	else 
	{
		double new_lat;
		double new_lon;
		int count;
		int tag;
		new_lat = 0;
		for(count = 0; count < this->npts_record_sum ; count++)
		{
			// to avoid the 0-180 transition
			// when lon > 150 || < -150, we use 0-360 range and change it back to -180 ~180 range
			tag = this->record_tag[count];
			new_lat += this->my_big_record->my_record[tag].sta_lat;
		}
		new_lat = new_lat / this->npts_record_sum;


		// we use one record to decide whick lon range to use
		tag = 0;
		double lon_tmp;
		lon_tmp = this->my_big_record->my_record[tag].sta_lon;
		new_lon = 0;
		if(   lon_tmp > 150 || lon_tmp < -150  )
		{
			for(count = 0; count < this->npts_record_sum ; count++)
			{
				// to avoid the 0-180 transition
				// when lon > 150 || < -150, we use 0-360 range and change it back to -180 ~180 range
				tag = this->record_tag[count];
				new_lon += this->my_big_record->my_record[tag].sta_lon + 180;
			}
			new_lon = new_lon / this->npts_record_sum;
			new_lon = new_lon - 180;
		}
		else
		{
			for(count = 0; count < this->npts_record_sum ; count++)
			{
				// to avoid the 0-180 transition
				// when lon > 150 || < -150, we use 0-360 range and change it back to -180 ~180 range
				tag = this->record_tag[count];
				new_lon += this->my_big_record->my_record[tag].sta_lon;
			}
			new_lon = new_lon / this->npts_record_sum;

		}


		this->grid_lat = new_lat;
		this->grid_lon = new_lon;

	}





}


/******************************************************************
 * This is a c script to calculate the SNR 
 * 1. the average SNR for all stations winthin range
 * 2. the SNR of stacked record 
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
void new_grid::get_SNR_before_and_after_stack()
{

	double ave_SNR;
	double stack_SNR;
	ave_SNR = 0;
	stack_SNR = 0;
	int tag;

	int count;
	for(count = 0; count < this->npts_record_sum; count++     )
	{

		tag = this->record_tag[count];
		double tmp_SNR;
		tmp_SNR = this->my_big_record->my_record[tag].SNR;
		if(tmp_SNR != tmp_SNR || tmp_SNR == 0)
			cout << "--> EORROR tmp_SNR problem in get_SNR_before_and_after_stack "<< endl;
		ave_SNR += tmp_SNR;
	}

	ave_SNR = ave_SNR / this->npts_record_sum;
	this->ave_SNR = ave_SNR;


	// get SNR afeter 
	// make sure the long window exist
	if( this->long_win[0] != this->long_win[0])
		cout << "->> ERROR long win is not read in yet, SNR cant be calcualted" << endl;



	double noise_signal;
	double phase_signal;
	noise_signal = 0;
	phase_signal = 0;

	int npts_noise_beg;
	int npts_noise_len;
	int npts_phase_beg;
	int npts_phase_len;

	npts_noise_beg = (int) ( (this->noise_beg - this->LONG_BEG) / this->delta );
	npts_phase_beg = (int) ( (this->phase_beg - this->LONG_BEG ) / this->delta );
	npts_noise_len = (int) (this->noise_len  / this->delta);
	npts_phase_len = (int) (this->phase_len  / this->delta);

	cout << " npts_noise_beg is "<< npts_noise_beg << " npts_noise_len is "<< npts_noise_len << endl;
	for(count = 0; count < npts_noise_len ; count++)
	{
		noise_signal += fabs (  this->long_win[npts_noise_beg + count] );
	}

	for(count = 0; count < npts_phase_len ; count++)
	{
		phase_signal += fabs ( this->long_win[npts_phase_beg + count] );
	}

	if( noise_signal == 0)
	{
		cout << "ERROR when calculateing SNR, noise signal is zero!" << endl;
		this->stack_SNR = 1;
	}
	else
		this->stack_SNR= phase_signal / noise_signal;










}
