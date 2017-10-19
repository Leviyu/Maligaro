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
		flag = "<VS_EXISTING_RECORD_NUM_THRESHOLD>";
		if(sub1.compare(flag) == 0)
			this->VS_EXISTING_RECORD_NUM_THRESHOLD = atof(sub2.c_str());
		flag = "<VS_RECORD_NUM_THRESHOLD>";
		if(sub1.compare(flag) == 0)
			this->VS_RECORD_NUM_THRESHOLD = atof(sub2.c_str());
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
	this->catagorize_existing_eventinfo_to_VS();


	// read eventStation file
	this->eventStation = "eventStation";
	int sta_num = count_file_num(this->eventStation);
	this->my_record = new new_record[sta_num];
	this->read_eventStation();
	this->catagorize_eventstation_to_VS();
	this->count_record_existance_for_grid_pair();




}

// ============================================
// loop through EQ_grid STA_grid PHASE pair and check
// 1. how many records exist for existing eventinfo
// 2. how many records exist for eventStation
void big_new_record::count_record_existance_for_grid_pair()
{
	cout << "---> count_record_existance_for_grid_pair " << endl;

	int ilat_eq = 0;
	int ilon_eq = 0;
	int ilat_sta = 0;
	int ilon_sta = 0;



	int flag = 0;
	//EQ grid loop
	for(ilat_eq = 0; ilat_eq<this->grid_lat_num; ilat_eq++)
		for(ilon_eq = 0; ilon_eq< this->grid_lon_num[ilat_eq]; ilon_eq++)
		{
			
			// if EQ num is 0
			if (this->my_grid[ilat_eq][ilon_eq].CU_EQ_NUM < 1)
				continue;
			cout << " grid CU EQ STA NUM " << this->my_grid[ilat_eq][ilon_eq].CU_EQ_NUM << endl;



			// sta grid loop
			for(ilat_sta = 0; ilat_sta<this->grid_lat_num; ilat_sta++)
				for(ilon_sta = 0; ilon_sta< this->grid_lon_num[ilat_sta]; ilon_sta++)
				{

					// if STA num < threshold
					// this is a restriction applied to single EQ processing
					if (this->my_grid[ilat_sta][ilon_sta].CU_STA_NUM < this->VS_RECORD_NUM_THRESHOLD )
						continue;


					// check eventinfo first
					// loop through current EQLIST STALIST
					//
					//int flag = this->loop_EX_EQ_STA_count_record(&this->my_grid[ilat_eq][ilon_eq], &this->my_grid[ilat_sta][ilon_sta]);
					flag +=1;
					cout << " current EX EQ STA count" << flag << endl;








				}


		}



}



int big_new_record::loop_EX_EQ_STA_count_record( new_grid* eq_grid, new_grid* sta_grid)
{

	if(eq_grid->EX_EQ_NUM < 1 || sta_grid->EX_STA_NUM < 1)
		return 0;

	int ieq = 0;
	int ista = 0;
	string eq_name = "";
	string sta_name = "";
	string phase_name = "";

	int return_value = 0;

	for(ieq = 0; ieq<eq_grid->EX_EQ_NUM ; ieq++)
		for(ista = 0; ista < sta_grid->EX_STA_NUM;ista++)
		{
			eq_name = eq_grid->EX_EQ_NAME[ieq];
			sta_name = sta_grid->EX_STA_NAME[ista];

			int irecord = 0;
			for(irecord = 0; irecord < this->sta_num; irecord++)
			{
				if( this->existing_record->EQ != eq_name  )
					continue;
				if( this->existing_record->STA != sta_name  )
					continue;
				if( this->PHASE != this->existing_record->PHASE )
					continue;

					return_value ++;



			}




		}

	
	return return_value;

}





// =========================================================
// Loop through each grid point
void big_new_record::catagorize_existing_eventinfo_to_VS()
{

	// step1, get unique EQ list
	string command;
	command = " cat eventinfo |awk '{print $12,$8,$9}'|sort|uniq > list.unique_EQ";
	exec(command);
	// step2 get unique sta list
	command = " cat eventinfo |grep -v PPP |awk '{print $1,$6,$7}'|sort|uniq > list.unique_STA";
	exec(command);


	int EQ_MAX = 1000;
	string EQ_LIST[EQ_MAX];
	double EQ_LAT[EQ_MAX];
	double EQ_LON[EQ_MAX];
	int EQ_NUM = 0;

	int STA_MAX = 10000;
	string STA_LIST[STA_MAX];
	double STA_LAT[STA_MAX];
	double STA_LON[STA_MAX];
	int STA_NUM = 0;

	ifstream myfile;
	string file_name = "list.unique_EQ";
	myfile.open(file_name);
	int file_count = 0;
	file_count = count_file_num(file_name);
	int count;
	for(count = 0; count < file_count ; count++)
	{
		myfile >> EQ_LIST[count] >> EQ_LAT[count] >> EQ_LON[count];
		EQ_NUM ++;
		//cout << EQ_LIST[count] << " " << EQ_LAT[count] << " " << EQ_LON[count] << endl;
	}
	myfile.close();

	//read in sta list
	file_name = "list.unique_STA";
	myfile.open(file_name);
	file_count = count_file_num(file_name);
	for(count = 0; count < file_count ; count++)
	{
		myfile >> STA_LIST[count] >> STA_LAT[count] >> STA_LON[count];
		//cout << STA_LIST[count] << endl;
		//if( STA_LIST[count].find( "PPP") != std::string::npos )
		//{
			//cout << "found PPP " << endl;
			//continue;
		//}
		STA_NUM ++;
		//cout << EQ_LIST[count] << " " << EQ_LAT[count] << " " << EQ_LON[count] << endl;
	}

	// step 3.
	// loop through EQLIST
	for(count = 0; count < EQ_NUM; count++)
	{
		//cout << "--> work on EQ NUM " << count  << " /" << EQ_NUM << endl;
		// loop through grid lat and lon and check grid that is within grid_radius in latitude
		int ilat = 0;
		int ilon = 0;
		for(ilat = 0; ilat < this->grid_lat_num ; ilat++)
			for(ilon = 0; ilon < this->grid_lon_num[ilat] ; ilon++)
			{
				double grid_lat = this->my_grid[ilat][ilon].grid_lat;
				if( fabs(grid_lat - EQ_LAT[count]) > this->VS_RADIUS_DEGREE)
					continue;

				double grid_lon = this->my_grid[ilat][ilon].grid_lon;
				// distance between EQ and grid
				double distance_km = dist_A_B( grid_lat, grid_lon, EQ_LAT[count], EQ_LON[count]);
				//cout << distance_km << endl;
				if( distance_km > 110*this->VS_RADIUS_DEGREE )
					continue;


				// add EQ info into current grid
				int current_index = this->my_grid[ilat][ilon].EX_EQ_NUM ;
				this->my_grid[ilat][ilon].EX_EQ_LAT[current_index] = EQ_LAT[count];
				this->my_grid[ilat][ilon].EX_EQ_LON[current_index] = EQ_LON[count];
				this->my_grid[ilat][ilon].EX_EQ_NAME[current_index] = EQ_LIST[count];
				this->my_grid[ilat][ilon].EX_EQ_NUM ++;
				//cout << " ilat lon EX EQ NUM" << ilat << " " << ilon << " " << this->my_grid[ilat][ilon].EX_EQ_NUM << endl;
			}

	}


	// step 4.
	// loop through STA_LIST
	for(count = 0; count < STA_NUM; count++)
	{
		//cout << "--> work on EQ NUM " << count  << " /" << STA_NUM << endl;
		// loop through grid lat and lon and check grid that is within grid_radius in latitude
		int ilat = 0;
		int ilon = 0;
		for(ilat = 0; ilat < this->grid_lat_num ; ilat++)
			for(ilon = 0; ilon < this->grid_lon_num[ilat] ; ilon++)
			{
				double grid_lat = this->my_grid[ilat][ilon].grid_lat;
				if( fabs(grid_lat - STA_LAT[count]) > this->VS_RADIUS_DEGREE)
					continue;

				double grid_lon = this->my_grid[ilat][ilon].grid_lon;
				// distance between EQ and grid
				double distance_km = dist_A_B( grid_lat, grid_lon, STA_LAT[count], STA_LON[count]);
				//cout << distance_km << endl;
				if( distance_km > 110*this->VS_RADIUS_DEGREE )
					continue;

				// add EQ info into current grid
				int current_index = this->my_grid[ilat][ilon].EX_STA_NUM ;
				this->my_grid[ilat][ilon].EX_STA_LAT[current_index] = STA_LAT[count];
				this->my_grid[ilat][ilon].EX_STA_LON[current_index] = STA_LON[count];
				this->my_grid[ilat][ilon].EX_STA_NAME[current_index] = STA_LIST[count];
				this->my_grid[ilat][ilon].EX_STA_NUM ++;
				//cout << " ilat lon EX STA NUM" << ilat << " " << ilon << " " << this->my_grid[ilat][ilon].EX_STA_NUM << endl;
			}

	}

}


void big_new_record::catagorize_eventstation_to_VS()
{
	cout << "-- > catagorize eventStation into VS " << endl;

	// step1, get unique EQ list
	string command;
	command = " cat eventStation |awk '{print $19,$11,$12}'|sort|uniq > list.unique_EQ";
	exec(command);
	// step2 get unique sta list
	command = " cat eventStation |grep -v PPP |awk '{print $1,$9,$10}'|sort|uniq > list.unique_STA";
	exec(command);


	int EQ_MAX = 1000;
	string EQ_LIST[EQ_MAX];
	double EQ_LAT[EQ_MAX];
	double EQ_LON[EQ_MAX];
	int EQ_NUM = 0;

	int STA_MAX = 10000;
	string STA_LIST[STA_MAX];
	double STA_LAT[STA_MAX];
	double STA_LON[STA_MAX];
	int STA_NUM = 0;

	ifstream myfile;
	string file_name = "list.unique_EQ";
	myfile.open(file_name);
	int file_count = 0;
	file_count = count_file_num(file_name);
	int count;
	for(count = 0; count < file_count ; count++)
	{
		myfile >> EQ_LIST[count] >> EQ_LAT[count] >> EQ_LON[count];
		EQ_NUM ++;
		//cout << EQ_LIST[count] << " " << EQ_LAT[count] << " " << EQ_LON[count] << endl;
	}
	myfile.close();

	//read in sta list
	file_name = "list.unique_STA";
	myfile.open(file_name);
	file_count = count_file_num(file_name);
	for(count = 0; count < file_count ; count++)
	{
		myfile >> STA_LIST[count] >> STA_LAT[count] >> STA_LON[count];
		//cout << STA_LIST[count] << endl;
		//if( STA_LIST[count].find( "PPP") != std::string::npos )
		//{
			//cout << "found PPP " << endl;
			//continue;
		//}
		STA_NUM ++;
		//cout << EQ_LIST[count] << " " << EQ_LAT[count] << " " << EQ_LON[count] << endl;
	}

	// step 3.
	// loop through EQLIST
	for(count = 0; count < EQ_NUM; count++)
	{
		// loop through grid lat and lon and check grid that is within grid_radius in latitude
		int ilat = 0;
		int ilon = 0;
		for(ilat = 0; ilat < this->grid_lat_num ; ilat++)
			for(ilon = 0; ilon < this->grid_lon_num[ilat] ; ilon++)
			{
				double grid_lat = this->my_grid[ilat][ilon].grid_lat;
				if( fabs(grid_lat - EQ_LAT[count]) > this->VS_RADIUS_DEGREE)
					continue;

				double grid_lon = this->my_grid[ilat][ilon].grid_lon;
				// distance between EQ and grid
				double distance_km = dist_A_B( grid_lat, grid_lon, EQ_LAT[count], EQ_LON[count]);
				//cout << distance_km << endl;
				if( distance_km > 110*this->VS_RADIUS_DEGREE )
					continue;

				// add EQ info into current grid
				int current_index = this->my_grid[ilat][ilon].CU_EQ_NUM ;
				this->my_grid[ilat][ilon].CU_EQ_LAT[current_index] = EQ_LAT[count];
				this->my_grid[ilat][ilon].CU_EQ_LON[current_index] = EQ_LON[count];
				this->my_grid[ilat][ilon].CU_EQ_NAME[current_index] = EQ_LIST[count];
				this->my_grid[ilat][ilon].CU_EQ_NUM ++;
			}

	}


	// step 4.
	// loop through STA_LIST
	for(count = 0; count < STA_NUM; count++)
	{
		// loop through grid lat and lon and check grid that is within grid_radius in latitude
		int ilat = 0;
		int ilon = 0;
		for(ilat = 0; ilat < this->grid_lat_num ; ilat++)
			for(ilon = 0; ilon < this->grid_lon_num[ilat] ; ilon++)
			{
				double grid_lat = this->my_grid[ilat][ilon].grid_lat;
				if( fabs(grid_lat - STA_LAT[count]) > this->VS_RADIUS_DEGREE)
					continue;

				double grid_lon = this->my_grid[ilat][ilon].grid_lon;
				// distance between EQ and grid
				double distance_km = dist_A_B( grid_lat, grid_lon, STA_LAT[count], STA_LON[count]);
				//cout << distance_km << endl;
				if( distance_km > 110*this->VS_RADIUS_DEGREE )
					continue;

				// add EQ info into current grid
				int current_index = this->my_grid[ilat][ilon].CU_STA_NUM ;
				this->my_grid[ilat][ilon].CU_STA_LAT[current_index] = STA_LAT[count];
				this->my_grid[ilat][ilon].CU_STA_LON[current_index] = STA_LON[count];
				this->my_grid[ilat][ilon].CU_STA_NAME[current_index] = STA_LIST[count];
				this->my_grid[ilat][ilon].CU_STA_NUM ++;
			}

	}

}
void big_new_record::virtual_station_grid_initiate()
{
	cout << "--> Initiate virtual station grid \n";
	new_grid** my_grid;
	// initiate virtual station grid
	this->grid_lat_num = (int)(180 / this->VS_LATITUDE_INC);
	this->grid_lon_num = (int*)malloc(sizeof(int)*this->grid_lat_num);
	my_grid = (new_grid**)malloc(sizeof(new_grid*)*this->grid_lat_num);


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
			my_grid[ilat][ilon].grid_lat = current_lat;
			my_grid[ilat][ilon].grid_lon = current_lon;
			my_grid[ilat][ilon].grid_radius = this->VS_RADIUS_DEGREE;

			// initiate EXisting EQ STA 
			my_grid[ilat][ilon].EX_EQ_NUM = 0;
			int MMM = 400;
			my_grid[ilat][ilon].EX_EQ_NAME = (string*)malloc(sizeof(string)*MMM);
			my_grid[ilat][ilon].EX_EQ_LAT = (double*)malloc(sizeof(double)*MMM);
			my_grid[ilat][ilon].EX_EQ_LON = (double*)malloc(sizeof(double)*MMM);
			my_grid[ilat][ilon].EX_STA_NUM = 0;
			MMM = 1000;
			my_grid[ilat][ilon].EX_STA_NAME = (string*)malloc(sizeof(string)*MMM);
			my_grid[ilat][ilon].EX_STA_LAT = (double*)malloc(sizeof(double)*MMM);
			my_grid[ilat][ilon].EX_STA_LON = (double*)malloc(sizeof(double)*MMM);

			// initiate current EQ STA
			my_grid[ilat][ilon].CU_EQ_NUM = 0;
			my_grid[ilat][ilon].CU_EQ_NAME = (string*)malloc(sizeof(string)*MMM);
			my_grid[ilat][ilon].CU_EQ_LAT = (double*)malloc(sizeof(double)*MMM);
			my_grid[ilat][ilon].CU_EQ_LON = (double*)malloc(sizeof(double)*MMM);
			my_grid[ilat][ilon].CU_STA_NUM = 0;
			MMM = 1000;
			my_grid[ilat][ilon].CU_STA_NAME = (string*)malloc(sizeof(string)*MMM);
			my_grid[ilat][ilon].CU_STA_LAT = (double*)malloc(sizeof(double)*MMM);
			my_grid[ilat][ilon].CU_STA_LON = (double*)malloc(sizeof(double)*MMM);


		}

	}
	this->my_grid = my_grid;


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
	int sta_num = count_file_num(this->eventStation);
	cout <<"--> big_record  Read in eventStation"<< sta_num << " records" << endl;

	for(line = 0; line < sta_num ; line++)
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



		//this->my_record[line].record_file = this->record_file;
// get ista line
		for(count = 1 ; count <= 19 ; count ++)
		{
			//cout << count << " / " << 19 << endl;
			ss >> sub1;
			if(sub1.empty())
				continue;

			//cout << sub1 << endl;

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
