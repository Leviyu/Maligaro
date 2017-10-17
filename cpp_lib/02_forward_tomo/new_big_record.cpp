#include "forward_tomography.h"

big_new_record::big_new_record()
{

}


big_new_record::~big_new_record()
{
	cout << "big_new_record is deleted boys!" << endl;
	//delete[] this->my_record;
}



void big_new_record::big_record_read_cross_point_file(new_tomo* my_tomo)
{
	int count;
	count = 0;
	for(count = 0; count < this->sta_num; count++ )
	{

		//cout << "read cross file for "<< this->my_record[count].STA << endl;
		this->my_record[count].read_cross_point_info( my_tomo);
	}

}

void big_new_record::get_incident()
{
	int count;

	for(count = 0; count < this->sta_num ; count++)
		this->my_record[count].get_incident();


}


// Read INFILE
void big_new_record::read_INFILE()
{
	string infile = "./INFILE";
	ifstream myfile;
	myfile.open(infile.c_str());
	string tmp;

	while(getline(myfile,tmp))
	{
		istringstream ss(tmp);
		string sub1;
		string sub2;
		string sub3;
		string sub4;
		string flag;
		ss >> sub1 >> sub2 >> sub3 >> sub4;
		//cout << sub1 << sub2 << endl;
		flag = "<PHASE_LONG_WIN_LEN>";
		if(sub1.compare(flag) == 0)
			this->long_len = atof(sub2.c_str());

		flag = "<PHASE_LONG_WIN_BEG>";
		if(sub1.compare(flag) == 0)
			this->long_beg = atof(sub2.c_str());

		flag = "<PHASE_WIN_BEG>";
		if(sub1.compare(flag) == 0)
			this->phase_beg = atof(sub2.c_str());

		flag = "<PHASE_WIN_LEN>";
		if(sub1.compare(flag) == 0)
			this->phase_len = atof(sub2.c_str());

		flag = "<PHASE_NOISE_BEG>";
		if(sub1.compare(flag) == 0)
			this->noise_beg = atof(sub2.c_str());

		flag = "<PHASE_NOISE_LEN>";
		if(sub1.compare(flag) == 0)
			this->noise_len = atof(sub2.c_str());

		flag = "<DELTA>";
		if(sub1.compare(flag) == 0)
		{
			this->delta = atof(sub2.c_str());
			//cout << "delta "<< this->delta << endl;
		}

		flag = "<VS_LATITUDE_INC>";
		if(sub1.compare(flag) == 0)
			this->VS_LATITUDE_INC = atof(sub2.c_str());

		flag = "<VS_RADIUS>";
		if(sub1.compare(flag) == 0)
			this->VS_RADIUS_DEGREE = atof(sub2.c_str());
		flag = "<EXISTING_EVENTINFO>";
		if(sub1.compare(flag) == 0)
			this->EXISTING_EVENTINFO = sub2;
	}
}

void big_new_record::get_crustal_correction()
{
	int count;

	for(count = 0; count < this->sta_num ; count++)
		this->my_record[count].get_crustal_correction();


}

void big_new_record::read_in_polarity_file()
{
	ifstream myfile;
	string file = "eventinfo.polarity";
	myfile.open(file.c_str());

	int count;
	double polar;
	string tmp;
	for(count = 0; count < this->sta_num ; count++)
	{
		myfile >>  tmp >>  polar;
		if( polar < 0 )
			this->my_record[count].polarity_flag = -1;
		else 
			this->my_record[count].polarity_flag = 1;
		cout << " current polarity is "<< this->my_record[count].polarity_flag<<endl;
	}


}

void big_new_record::get_ellip_corr()
{
	int count;
	count = 0;
	ofstream out;
	string out_file = "ellip_cor.data";
	out.open(out_file.c_str());
	for(count = 0; count < this->sta_num; count++ )
	{
		//cout << "read cross file for "<< this->my_record[count].STA << endl;
		this->my_record[count].get_ellip_corr();
		out << this->my_record[count].STA << " " << this->my_record[count].ellip_corr << endl;
	}

	out.close();

}


void big_new_record::initiate_big_record()
{
	cout << " big_new_record is initiated ! " << endl;

}

void big_new_record::virtual_station_main()
{
	cout << "--> Running virtual Station stacking main " << endl;

	this->sta_num = count_file_num(this->EXISTING_EVENTINFO);
	cout << " --> Read in existing record num: " << this->sta_num << endl;
	this->existing_record = new new_record[this->sta_num];
	this->record_file = this->EXISTING_EVENTINFO;

	// read in eventinfo
	this->read_record_file(this->existing_record);

	// initiate virtual station
	this->virtual_station_grid_initiate();

	// catagorize existing eventinfo into virtual station

}


void big_new_record::catagorize_existing_eventinfo_to_VS()
{

	int ilat, ilon;
	for(ilat = 0; ilat < this->grid_lat_num ; ilat++)
		for(ilon = 0; ilon < this->grid_lon_num[ilat] ; ilon++)
		{
			double current_lat = this->my_grid[ilat][ilon].grid_lat;
			double current_lon = this->my_grid[ilat][ilon].grid_lon;

			// loop through stations





		}







}


void big_new_record::virtual_station_grid_initiate()
{
	new_grid** my_grid;
	this->my_grid = my_grid;
	// initiate virtual station grid
	this->grid_lat_num = (int)(180 / this->VS_LATITUDE_INC);
	this->grid_lon_num = (int*)malloc(sizeof(int)*this->grid_lat_num);
	this->my_grid = (new_grid**)malloc(sizeof(new_grid*)*this->grid_lat_num);


	cout << "grid lat num is " << this->grid_lat_num << endl;

	// calculate longitude grid num for each latitude
	int ilat,ilon;
	double current_lat = 0;
	double lat_inc_in_km = this->VS_LATITUDE_INC * 110;
	for( ilat = 0; ilat < this->grid_lat_num ; ilat++)
	{
		current_lat = -89 + ilat * this->VS_LATITUDE_INC;
		this->grid_lon_num[ilat] = floor( 2*3.1415926*6371*cos( current_lat * 3.1415926/180) / lat_inc_in_km );

		//cout << " current lat " << ilat << " lon num" << this->grid_lon_num[ilat] << endl;
		my_grid[ilat] = (new_grid*)malloc(sizeof(new_grid)*this->grid_lon_num[ilat]);
		
		for(ilon = 0; ilon < this->grid_lon_num[ilat] ; ilon++)
		{
			double lon_inc_in_degree = 360/this->grid_lon_num[ilat];
			double current_lon = ilon * lon_inc_in_degree - 180;
			//cout << " lat lon "<< current_lat << " "<< current_lon<< endl;
			this->my_grid[ilat][ilon].grid_lat = current_lat;
			this->my_grid[ilat][ilon].grid_lon = current_lon;
			this->my_grid[ilat][ilon].grid_radius = this->VS_RADIUS_DEGREE;

		}

	}


}



/*************************************************************
* This C function read in station information 
* for a certain phase
*	INPUT:
*	
*	Hongyu DATE: Aug 2016
*	Key words: read in station file
*************************************************************/
void big_new_record::read_record_file(new_record* my_record)
{

	int count;
	string tmp;
	ifstream myfile;
	myfile.open(this->record_file.c_str());
		string sub1;
		string flag;

	int line;
	line = 0;

	cout << "record file is "<< this->record_file << endl;
	cout <<"--> big_record  Read in "<< this->sta_num << " records" << endl;

	for(line = 0; line < this->sta_num ; line++)
	{
		getline(myfile,tmp);
		//cout << "read record "<< line+1 << endl;
		//cout << tmp << endl;
		istringstream ss(tmp);
		//cout << this->delta << " " << this->long_beg<< " "<< this->long_len << endl;

		// get delta value here 
		if(this->delta == this->delta && this->delta != 0)
			my_record[line].delta = this->delta;
		if(this->long_beg == this->long_beg && this->long_beg != 0)
			my_record[line].long_beg = this->long_beg;
		if(this->long_len == this->long_len && this->long_len != 0)
			my_record[line].long_len = this->long_len;



		my_record[line].record_file = this->record_file;
// get ista line
		for(count = 1 ; count <=45 ; count ++)
		{
			//cout << "ist " << count << endl;
			ss >> sub1;
			if(sub1.empty())
				continue;


			if( count == 1 )
			{
				my_record[line].STA = sub1;
			}
			else if( count == 2 )
				my_record[line].NET = sub1;
			else if( count == 3 )
				my_record[line].DIST = atof(sub1.c_str());
			else if( count == 4 )
				my_record[line].AZ = atof(sub1.c_str());
			else if( count == 5 )
				my_record[line].BAZ = atof(sub1.c_str());
			else if( count == 6 )
				my_record[line].sta_lat = atof(sub1.c_str());
			else if( count == 7 )
				my_record[line].sta_lon = atof(sub1.c_str());
			else if( count == 8 )
				my_record[line].eq_lat = atof(sub1.c_str());
			else if( count == 9 )
				my_record[line].eq_lon = atof(sub1.c_str());
			else if( count == 10 )
				my_record[line].eq_dep = atof(sub1.c_str());
			else if( count == 11)
				my_record[line].eq_mag = atof(sub1.c_str());
			else if( count == 12 )
				my_record[line].EQ = sub1;
			else if( count == 13 )
				my_record[line].polarity_flag = atoi(sub1.c_str());
			else if( count == 14 )
				my_record[line].quality_flag  = atoi(sub1.c_str());
			else if( count == 15 )
				my_record[line].PREM = atof(sub1.c_str());
			else if( count == 16 )
				my_record[line].phase_amplitude = atof(sub1.c_str());
			else if( count == 17 )
				my_record[line].CCC = atof(sub1.c_str());
			else if( count == 18 )
				my_record[line].SNR = atof(sub1.c_str());
			else if( count == 19 )
				my_record[line].dt_obs_prem  = atof(sub1.c_str());
			else if( count == 20 )
				my_record[line].PHASE = sub1;
			else if( count == 21 )
				my_record[line].stretch_ccc = atof(sub1.c_str());
			else if( count == 22 )
				my_record[line].stretch_coeff = atof(sub1.c_str());
			else if( count == 23 )
				my_record[line].misfit = atof(sub1.c_str());
			else if( count == 24 )
			{
				my_record[line].COMP = sub1;
				//cout << my_record[line].COMP<< endl;
			}
			else if( count == 25 )
				my_record[line].phase_peak_time_rel_PREM = atof(sub1.c_str());
			else if( count == 26 )
				my_record[line].npts_phase_peak_rel_start = atoi(sub1.c_str());
			else if( count == 27 )
				my_record[line].noise_beg = atof(sub1.c_str());
			else if( count == 28 )
				my_record[line].noise_len = atof(sub1.c_str());
			else if( count == 29 )
				my_record[line].phase_beg_rel_PREM = atof(sub1.c_str());
			else if( count == 30 )
				my_record[line].record_weight = atof(sub1.c_str());
			else if( count == 31 )
				my_record[line].SNR2 = atof(sub1.c_str());
			else if( count == 32 )
				my_record[line].misfit2 = atof(sub1.c_str());
			else if( count == 33 )
				my_record[line].ONSET = atof(sub1.c_str());
			else if( count == 34 )
				my_record[line].ENDSET = atof(sub1.c_str());
			else if( count == 35 )
				my_record[line].tstar_factor = atof(sub1.c_str());
			else if( count == 36 )
				my_record[line].tstar_ccc = atof(sub1.c_str());
			else if( count == 37 )
				my_record[line].ccc3 = atof(sub1.c_str());
			else if( count == 38 )
				my_record[line].misfit_pre = atof(sub1.c_str());
			else if( count == 39 )
				my_record[line].misfit_bak = atof(sub1.c_str());
			else if( count == 40 )
				my_record[line].record_gau_factor= atof(sub1.c_str());
			else if( count == 41 )
				my_record[line].EW_gau_factor = atof(sub1.c_str());
			else if( count == 42 )
				my_record[line].gau_misfit = atof(sub1.c_str());
			else if( count == 43 )
				my_record[line].polarity = atof(sub1.c_str());
			else if( count == 44 )
			{
				my_record[line].polarity_prediction= atof(sub1.c_str());
				//cout << my_record[line].polarity_prediction<< endl;
			}
			else if( count == 45 )
			{
				my_record[line].traffic_phase_nearby= atof(sub1.c_str());
				//cout << my_record[line].traffic_phase_nearby<< endl;
			}
		}



	}
	myfile.close();

	cout << "Big record read record file done" << endl;

}








/*************************************************************
* This C function read in eventStation
*	INPUT:
*	
*	Hongyu DATE: Aug 2016
*	Key words: read in station file
*************************************************************/
void big_new_record::read_eventStation()
{

	int count;
	string tmp;
	ifstream myfile;
	myfile.open(this->eventStation.c_str());
		string sub1;
		string flag;

	int line;
	line = 0;

	cout <<"--> big_record  Read in eventStation"<< this->sta_num << " records" << endl;

	for(line = 0; line < this->sta_num ; line++)
	{
		getline(myfile,tmp);
		//cout << "read record "<< line+1 << endl;
		//cout << tmp << endl;
		istringstream ss(tmp);

		// get delta value here 
		if(this->delta == this->delta && this->delta != 0)
			my_record[line].delta = this->delta;
		if(this->long_beg == this->long_beg && this->long_beg != 0)
			my_record[line].long_beg = this->long_beg;
		if(this->long_len == this->long_len && this->long_len != 0)
			my_record[line].long_len = this->long_len;



		this->my_record[line].record_file = this->record_file;
// get ista line
		for(count = 1 ; count < 36 ; count ++)
		{
			ss >> sub1;
			if(sub1.empty())
				continue;


			if( count == 1 )
				this->my_record[line].STA = sub1;
			else if( count == 2 )
				this->my_record[line].NET = sub1;
			else if( count == 3 )
				this->my_record[line].DIST = atof(sub1.c_str());
			else if( count == 5 )
				this->my_record[line].AZ = atof(sub1.c_str());
			else if( count == 7 )
				this->my_record[line].BAZ = atof(sub1.c_str());
			else if( count == 9 )
				this->my_record[line].sta_lat = atof(sub1.c_str());
			else if( count == 10 )
				this->my_record[line].sta_lon = atof(sub1.c_str());
			else if( count == 11 )
				this->my_record[line].eq_lat = atof(sub1.c_str());
			else if( count == 12 )
				this->my_record[line].eq_lon = atof(sub1.c_str());
			else if( count == 13 )
				this->my_record[line].eq_dep = atof(sub1.c_str());
			else if( count == 16 )
				this->my_record[line].eq_mag= atof(sub1.c_str());
			else if( count == 19 )
				this->my_record[line].EQ = sub1;
		}



	}
	myfile.close();

	cout << "Big record read record file done" << endl;

}











/*************************************************************
* This C function read in sac file for current station 
* for a certain phase
*	INPUT:
*	
*	Hongyu DATE: Sep 2016
*	Key words: read in station file
*************************************************************/
void big_new_record::read_sac_file()
{

	int ista;
	for(ista = 0; ista < this->sta_num;ista++)
		this->my_record[ista].read_sac_file();


}

/*************************************************************
* This C function calculate the SNR for each station
* for a certain phase
*	INPUT:
*	
*	Hongyu DATE: Sep 2016
*************************************************************/
void big_new_record::calculate_SNR()
{

	int ista;
	for(ista = 0; ista < this->sta_num;ista++)
	{

		this->my_record[ista].noise_beg = this->noise_beg;
		this->my_record[ista].noise_len = this->noise_len;
		this->my_record[ista].phase_beg = this->phase_beg;
		this->my_record[ista].phase_len = this->phase_len;


		this->my_record[ista].calculate_SNR();
	}


}
