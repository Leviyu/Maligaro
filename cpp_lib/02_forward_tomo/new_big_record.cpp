#include "forward_tomography.h"
#include <map>
big_new_record::big_new_record()
{

}


big_new_record::~big_new_record()
{
	cout << "big_new_record is deleted boys!" << endl;
	//delete[] this->my_record;
}

void big_new_record::virtual_station_main()
{
	cout << "=============> Running virtual Station stacking main " << endl;

	//this->sta_num = count_file_num(this->EXISTING_EVENTINFO);
	

	// lets define a virtual station array
	int vs_max = 10000;
	this->my_vs = new virtual_station[vs_max]; 
	this->my_vs_index = 0;
	


	// =======================================================
	// existing eventinfo block
	this->existing_sta_num = count_file_num(this->EXISTING_EVENTINFO);
	cout << " --> Read in existing record num: " << this->existing_sta_num << endl;
	this->existing_record = new new_record[this->existing_sta_num];
	this->record_file = this->EXISTING_EVENTINFO;
	// read in eventinfo
	this->sta_num = this->existing_sta_num;
	this->read_record_file(this->existing_record);
	// initiate virtual station
	this->virtual_station_grid_initiate();
	// catagorize existing eventinfo into virtual station
	//this->catagorize_records_into_VS_based_on_each_record( this->existing_record, this->existing_sta_num);
	this->catagorize_existing_eventinfo_to_VS();


	// =================================
	// possible eventStation block
	this->eventStation = "eventStation";
	this->sta_num = count_file_num(this->eventStation);
	this->my_record = new new_record[this->sta_num];
	this->read_eventStation();
	//this->catagorize_records_into_VS_based_on_each_record(this->my_record, this->sta_num);
	this->catagorize_eventstation_to_VS();
	

	// for each EQ_grid/STA_grid/PHASE, count
	// 1. good_record K1
	// 2. possible record K2
	this->count_record_existance_for_grid_pair();

	this->output_VS_info();


	//this->read_sac_file();
	// go ahead 
}

void big_new_record::output_VS_info()
{
	cout << "output_VS_info \n"<<endl;
	int ivs = 0;
	int ista = 0;
	int count;
	ofstream myfile;
	string outname = "out.station_stack.info."+this->PHASE;
	myfile.open(outname.c_str());

	for(ivs = 0; ivs < this->my_vs_index; ivs++)
	{	
		virtual_station my_vs = this->my_vs[ivs];
		int vs_index = ivs;
		int ilat_eq = my_vs.ilat_eq;
		int ilon_eq = my_vs.ilon_eq;
		double eq_lat = this->my_grid[ilat_eq][ilon_eq].grid_lat;
		double eq_lon = this->my_grid[ilat_eq][ilon_eq].grid_lon;
		double eq_radius = this->VS_EQ_RADIUS_DEGREE;
		
		int ilat_sta = my_vs.ilat_sta;
		int ilon_sta = my_vs.ilon_sta;
		double sta_lat = this->my_grid[ilat_sta][ilon_sta].grid_lat;
		double sta_lon = this->my_grid[ilat_sta][ilon_sta].grid_lon;
		double sta_radius = this->VS_STA_RADIUS_DEGREE;
		
		int record_num = my_vs.eventStation_index;
		string EQ_NAME = my_vs.EQ;

		// output eventStation list for current VS
		ofstream out;
		string eventout = "out.VS_eventStation_list."+this->PHASE+"."+to_string(ivs);
		out.open(eventout.c_str());
		for(count = 0; count < my_vs.eventStation_index;count++)
		{
			int sta_index = my_vs.eventStation_index_array[count];
			string long_file = "long_win."+my_record[sta_index].EQ+
			"."+my_record[sta_index].STA+"."+this->PHASE;
			double record_dist = my_record[sta_index].DIST;
			out<< long_file<<" "<< record_dist<<endl;
		}
		out.close();

		// cout << " sta "<< sta_lat << " "<< sta_lon<<endl; 
		myfile<< fixed
		<< setw(10) << setprecision(2) << ivs
		<< setw(15) << EQ_NAME
		<< setw(10) << setprecision(2) << eq_lat
		<< setw(10) << setprecision(2) << eq_lon
		<< setw(10) << setprecision(2) << eq_radius
		<< setw(10) << setprecision(2) << sta_lat
		<< setw(10) << setprecision(2) << sta_lon
		<< setw(10) << setprecision(2) << sta_radius
		<< setw(10) << setprecision(2) << record_num
		<< setw(10) << setprecision(2) << my_vs.grid_dist
		<<endl;
	}
	


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
		flag = "<SRCDIR>";
		if(sub1.compare(flag) == 0)
			this->SRCDIR = sub2;

		flag = "<DELTA>";
		if(sub1.compare(flag) == 0)
		{
			this->delta = atof(sub2.c_str());
			//cout << "delta "<< this->delta << endl;
		}

		flag = "<VS_LATITUDE_INC>";
		if(sub1.compare(flag) == 0)
			this->VS_LATITUDE_INC = atof(sub2.c_str());

		flag = "<VS_EQ_RADIUS>";
		if(sub1.compare(flag) == 0)
			this->VS_EQ_RADIUS_DEGREE = atof(sub2.c_str());
		flag = "<VS_STA_RADIUS>";
		if(sub1.compare(flag) == 0)
			this->VS_STA_RADIUS_DEGREE = atof(sub2.c_str());
		flag = "<EXISTING_EVENTINFO>";
		if(sub1.compare(flag) == 0)
			this->EXISTING_EVENTINFO = sub2;
		flag = "<VS_EXISTING_RECORD_NUM_THRESHOLD>";
		if(sub1.compare(flag) == 0)
			this->VS_EXISTING_RECORD_NUM_THRESHOLD = atof(sub2.c_str());
		flag = "<VS_RECORD_NUM_THRESHOLD>";
		if(sub1.compare(flag) == 0)
			this->VS_RECORD_NUM_THRESHOLD = atof(sub2.c_str());
		flag = "<eventinfo_max_threshold>";
		if(sub1.compare(flag) == 0)
			this->eventinfo_max_threshold = atoi(sub2.c_str());
		flag = "<eventStation_min_threshold>";
		if(sub1.compare(flag) == 0)
			this->eventStation_min_threshold = atoi(sub2.c_str());
		


	}
}


void big_new_record::get_PHASE_min_max_dist()
{

	// set default value
	this->phase_dist_min = 0;
	this->phase_dist_max = 180;

	map<string,double> dist_min;
	map<string,double> dist_max;


// ========= S ===================
	string myphase = "S";
	dist_min[myphase] = 30;
	dist_max[myphase] = 103;

	myphase = "Sdiff";
	dist_min[myphase] = 98;
	dist_max[myphase] = 180;

	myphase = "SS";
	dist_min[myphase] = 60;
	dist_max[myphase] = 180;

	myphase = "SSm";
	dist_min[myphase] = 154;
	dist_max[myphase] = 180;

	myphase = "S3";
	dist_min[myphase] = 90;
	dist_max[myphase] = 180;
	
	myphase = "S3m";
	dist_min[myphase] = 65;
	dist_max[myphase] = 180;

	myphase = "S4";
	dist_min[myphase] = 120;
	dist_max[myphase] = 180;

	myphase = "S4m";
	dist_min[myphase] = 0;
	dist_max[myphase] = 180;

	myphase = "S5";
	dist_min[myphase] = 130;
	dist_max[myphase] = 180;


	myphase = "S5m";
	dist_min[myphase] = 0;
	dist_max[myphase] = 180;


	myphase = "S6";
	dist_min[myphase] = 150;
	dist_max[myphase] = 180;


	myphase = "S6m";
	dist_min[myphase] = 0;
	dist_max[myphase] = 180;

// ============== ScS ================
	myphase = "ScS";
	dist_min[myphase] = 0;
	dist_max[myphase] = 83;

	myphase = "ScSScS";
	dist_min[myphase] = 0;
	dist_max[myphase] = 164;

	myphase = "ScS3";
	dist_min[myphase] = 0;
	dist_max[myphase] = 180;

	myphase = "ScS3m";
	dist_min[myphase] = 124;
	dist_max[myphase] = 180;

	myphase = "ScS4";
	dist_min[myphase] = 0;
	dist_max[myphase] = 180;

	myphase = "ScS4m";
	dist_min[myphase] = 32;
	dist_max[myphase] = 180;


	this->phase_dist_min = dist_min[this->PHASE];
	this->phase_dist_max = dist_max[this->PHASE];
	cout << " PHASE: " << this->PHASE << " dist_min_max " << this->phase_dist_min 
	<< " " << this->phase_dist_max << endl;

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
	//cout << " big_new_record is initiated ! " << endl;

}



// ============================================
// loop through EQ_grid STA_grid PHASE pair and check
// 1. how many records exist for existing eventinfo
// 2. how many records exist for eventStation
// to make it fast, we use SQL to implement the command
void big_new_record::count_record_existance_for_grid_pair()
{
	cout << "---> count_record_existance_for_grid_pair " << endl;


	// output record related 
	ofstream myfile;
	myfile.open("out.record_existance_for_EQ_STA_pair");
	
	/*=================================
	Part1 find existing record in eventinfo
	*/
	int ilat_eq = 0;
	int ilon_eq = 0;
	int ilat_sta = 0;
	int ilon_sta = 0;

	int vs_index = 0;


	int count;	
	for(ilat_eq = 0; ilat_eq<this->grid_lat_num; ilat_eq++)
		for(ilon_eq = 0; ilon_eq< this->grid_lon_num[ilat_eq]; ilon_eq++)
		{

			for(ilat_sta = 0; ilat_sta<this->grid_lat_num; ilat_sta++)
				for(ilon_sta = 0; ilon_sta< this->grid_lon_num[ilat_sta]; ilon_sta++)
				{
					//cout << "EQ_GRID_LAT_LON/STA_GRID_LAT_LON: "<< ilat_eq << " "<< ilon_eq<<" "<< ilat_sta << " "<<ilon_sta<<  endl;
					

					


					//int cu_eq_num = this->my_grid[ilat_eq][ilon_eq].CU_EQ_NUM;
					//int cu_sta_num = this->my_grid[ilat_sta][ilon_sta].CU_STA_NUM;
					int ex_eq_num = this->my_grid[ilat_eq][ilon_eq].EX_EQ_NUM;
					int ex_sta_num = this->my_grid[ilat_sta][ilon_sta].EX_STA_NUM;
					if(ex_eq_num == 0 || ex_sta_num == 0)
						continue;

					// cout << "EX eq sta num" << ex_eq_num << " " << ex_sta_num << endl;
					//myfile << "EQ_GRID_LAT_LON/STA_GRID_LAT_LON: "<< ilat_eq << " "<< ilon_eq<<" "<< ilat_sta << " "<<ilon_sta<<  endl;
					
				
					

					// IF eventinfo num > eventinfo_max_threshold, skip
					int exist_eventinfo = 0;
					count_eventinfo_existance_for_sing_grid_pair(ilat_eq,ilon_eq,ilat_sta,ilon_sta,&exist_eventinfo);
					// this->my_grid[i]
					// cout << "existed eventinfo "<< exist_eventinfo << endl;
					// cout << " eventinfo threshold "<< this->eventinfo_max_threshold<< endl;
					if(exist_eventinfo > this->eventinfo_max_threshold)
						continue;
					
					// IF eventStation < eventStation_min_threshold, skip

					int exist_eventStation = 0;
					count_eventStation_existance_for_sing_grid_pair(ilat_eq,ilon_eq,ilat_sta,ilon_sta,&exist_eventStation);
					
					if(exist_eventStation < this->eventStation_min_threshold)
						continue;
					// cout << "existed eventinfo "<< exist_eventinfo << endl;
					// cout << "existed eventStation "<< exist_eventStation << endl;
					// make virtual station for current EQ_grid/STA_grid
					// this->my_vs[this->my_vs_index].exist = "YES";

					make_virtual_station(ilat_eq,ilon_eq,ilat_sta,ilon_sta);



					// this->my_vs_index ++;

				}
		}

}


void big_new_record::make_virtual_station(int ilat_eq, int ilon_eq, int ilat_sta, int ilon_sta)
{
	cout << " --> Making virtualStation for " << ilat_eq << " "<< ilon_eq << " " << ilat_sta << " "<< ilon_sta
	<<  endl;  

	// scan each events in EQ_grid and for each event, find all possible EQ-STA-PHASE

	int ieq = 0;
	int ista = 0;
	string EQ_NAME;
	string STA_NAME;
	int eq_num = this->my_grid[ilat_eq][ilon_eq].EX_EQ_NUM;
	int sta_num = this->my_grid[ilat_sta][ilon_sta].EX_STA_NUM;
	for(ieq = 0; ieq < eq_num ; ieq++)
	{
		make_virtual_station_for_EQ(ilat_eq,ilon_eq,ilat_sta,ilon_sta,ieq);
	}



}

void big_new_record::make_virtual_station_for_EQ(int ilat_eq, int ilon_eq, int ilat_sta, int ilon_sta, int ieq_index)
{

	int ieq = 0;
	int ista = 0;
	string EQ_NAME;
	EQ_NAME = this->my_grid[ilat_eq][ilon_eq].EX_EQ_NAME[ieq_index];
	string STA_NAME;
	int eq_num = this->my_grid[ilat_eq][ilon_eq].EX_EQ_NUM;
	int sta_num = this->my_grid[ilat_sta][ilon_sta].EX_STA_NUM;

	int total_exist_eventStation_num = 0;

	// initiate virtual staion
	// if current vs does not have enough stations, then it is overwrite by the next 
	int vs_index = this->my_vs_index;
	this->my_vs[vs_index].initiate();
	this->my_vs[vs_index].my_big_record = this;
	this->my_vs[vs_index].ilat_eq = ilat_eq;
	this->my_vs[vs_index].ilon_eq = ilon_eq;
	this->my_vs[vs_index].ilat_sta = ilat_sta;
	this->my_vs[vs_index].ilon_sta = ilat_sta;
	this->my_vs[vs_index].EQ = EQ_NAME;
	this->my_vs[vs_index].eventStation_index = 0;
	this->my_vs[vs_index].long_npts = (int)(this->long_len / this->delta);
	this->my_vs[vs_index].LONG_BEG = this->long_beg;
	this->my_vs[vs_index].delta = this->delta;



	int max = 1000;
	int stack_record_array[max];
	int istack = 0;

	for(ista = 0; ista < sta_num; ista++)
	{
		STA_NAME = this->my_grid[ilat_sta][ilon_sta].EX_STA_NAME[ista];
		int ista_eventStation = this->find_EQ_STA_PHASE_number_in_eventStation(EQ_NAME,STA_NAME);
		if(ista_eventStation == 0) 
			continue;
		stack_record_array[istack] = ista_eventStation;
		istack++;
		total_exist_eventStation_num ++;

		int iindex = this->my_vs[vs_index].eventStation_index;
		this->my_vs[vs_index].eventStation_index_array[iindex] = ista_eventStation;
		this->my_vs[vs_index].eventStation_index ++;
	}
	
	// if there is not enough records to stack for current event, then skip
	if(total_exist_eventStation_num < this->eventStation_min_threshold)
		return;


	// store EQ/STA grid info into virtual station



	// 1.0 download sac file for EQ-STA-PHASE
	this->individual_VS_processing();

	

	

	// 5. input S_ES file and find virtual station ONSET

	// 6. make code choice for virtual stations

	// 7. 
	cout << "Current VS index " << this->my_vs_index << endl;
	this->my_vs_index++;


}

void big_new_record::individual_VS_processing()
{
	int vs_index = this->my_vs_index;
	int count;
	this->my_vs[vs_index].npts_record_sum = this->my_vs[vs_index].eventStation_index;

	int tag_index = 0;
	for(count = 0; count < this->my_vs[vs_index].eventStation_index; count ++)
	{
		int station_index = this->my_vs[vs_index].eventStation_index_array[count];
		// cout << " eventStation array num "<< stack_record_array[count] << endl;
		
		cout << "my vs index is "<< vs_index << " eventStation index is  "<< station_index << " eventstation name is "<< this->my_record[station_index].STA << endl;
		// cout << "record lat lon is "<< this->my_record[station_index].sta_lat<< endl;
		this->my_record[station_index].download_sac_file();
		this->my_record[station_index].PHASE = this->PHASE;
		this->my_record[station_index].delta = 0.1;
		this->my_record[station_index].noise_beg = this->noise_beg;
		this->my_record[station_index].noise_len = this->noise_len;
		
		this->my_record[station_index].phase_len = this->phase_len;
		this->my_record[station_index].phase_beg = this->phase_beg;
		this->my_record[station_index].long_len = this->long_len;
		this->my_record[station_index].long_beg = this->long_beg;
		// cout << "long beg is "<< this->long_beg<<endl;

		// this->my_record[count].
		// this->my_record[count].
		// 1.1 make and read in polarity file

		// 1.2 read in sac file
		this->my_record[station_index].read_sac_file();
		

		// 1.3 calculate SNR before stack
		this->my_record[station_index].calculate_SNR();

		this->my_vs[vs_index].record_tag[tag_index] = station_index;
		tag_index ++;

	}


	this->my_vs[vs_index].stack_records_from_one_EQ();
	this->my_vs[vs_index].output_stacked_record();

	// 4. get SNR  after stack


}

void big_new_record::get_unique_EQ_latlon(string EQ_NAME, double* lat, double* lon)
{
	int count = 0;
	for(count = 0; count < this->unique_EQ_num; count++)
	{
		if(EQ_NAME.compare(this->EQ_LIST[count]) == 0)
		{
			*lat = this->EQ_LAT[count];
			*lon = this->EQ_LON[count];
			break;
		}

	}
}
void big_new_record::get_unique_STA_latlon(string STA_NAME, double* lat, double* lon)
{
	int count = 0;
	for(count = 0; count < this->unique_STA_num; count++)
	{
		if(STA_NAME.compare(this->STA_LIST[count]) == 0)
		{
			*lat = this->STA_LAT[count];
			*lon = this->STA_LON[count];
			break;
		}
	}
}


void big_new_record::count_eventStation_existance_for_sing_grid_pair(int ilat_eq, int ilon_eq, int ilat_sta, int ilon_sta, int* exist_eventStation)
{
	//cout << "count_record_existance_for_sing_grid_pair !\n"<<endl;

	int ieq = 0;
	int ista = 0;

	// 1. for all EQs in EQ_grid and STAs in STA_grid, check possible eventStation records
	int ex_eq_num = this->my_grid[ilat_eq][ilon_eq].EX_EQ_NUM;
	int ex_sta_num = this->my_grid[ilat_sta][ilon_sta].EX_STA_NUM;
	// if eq_num is zero or sta num is 0, then skip
	if( ex_eq_num == 0 || ex_sta_num == 0)
	{
		*exist_eventStation = 0;
	}
	else
	{
		for(ieq = 0; ieq<ex_eq_num; ieq ++)
			for(ista = 0; ista < ex_sta_num;ista++)
			{
				// calculate distance between EQ and STA
				string EQ_NAME = this->my_grid[ilat_eq][ilon_eq].EX_EQ_NAME[ieq];
				string STA_NAME = this->my_grid[ilat_sta][ilon_sta].EX_STA_NAME[ista];
				int num = find_EQ_STA_PHASE_number_in_eventStation(EQ_NAME,STA_NAME);
				// int array_index = this->my_vs[this->my_vs_index].EQ_index;
				// this->my_vs[this->my_vs_index].EQ_NAME_array[array_index] = EQ_NAME;
				// array_index = this->my_vs[this->my_vs_index].eventStation_index;
				// this->my_vs[this->my_vs_index].eventStation_index_array[array_index] = num;
				if (num != 0)
					*exist_eventStation += 1;
			}
	}
}
void big_new_record::count_eventinfo_existance_for_sing_grid_pair(int ilat_eq, int ilon_eq, int ilat_sta, int ilon_sta, int* exist_eventinfo)
{
	int ieq = 0;
	int ista = 0;
	// 2. for all EQs in EQ_grid and STAs in STA_grid, check existing eventinfo records
	int cu_eq_num = this->my_grid[ilat_eq][ilon_eq].CU_EQ_NUM;
	int cu_sta_num = this->my_grid[ilat_sta][ilon_sta].CU_STA_NUM;
	//cout << " existing eq num sta num "<< ex_eq_num << " " << ex_sta_num << endl;
	if(cu_eq_num == 0 || cu_sta_num == 0)
	{
		*exist_eventinfo = 0;
	}
	else 
	{
		*exist_eventinfo = 0;
		for(ieq = 0; ieq<cu_eq_num; ieq ++)
			for(ista = 0; ista < cu_sta_num;ista++)
			{
				string EQ_NAME = this->my_grid[ilat_eq][ilon_eq].CU_EQ_NAME[ieq];
				string STA_NAME = this->my_grid[ilat_sta][ilon_sta].CU_STA_NAME[ista];
				// find how many if EQ-STA-PHASE exist in eventinfo
				int num = this->find_EQ_STA_PHASE_number_in_eventinfo(EQ_NAME,STA_NAME,this->PHASE);
				if(num != 0)
					*exist_eventinfo += 1 ;
			}

		

	}
	//cout << "eventinfo num is "<< *exist_eventinfo<< endl;
}

int big_new_record::find_EQ_STA_PHASE_number_in_eventinfo(string EQ_NAME, string STA, string PHASE)
{
	int return_num = 0;
	int count;
	for(count = 0; count < this->existing_sta_num; count ++)
	{

		if(EQ_NAME.compare(this->existing_record[count].EQ) != 0)
			continue;
		if(STA.compare(this->existing_record[count].STA) != 0)
			continue;
		if(PHASE.compare(this->existing_record[count].PHASE) != 0)
			continue;
		return_num = count;
		break;

	}
	return return_num;

}
int big_new_record::find_EQ_STA_PHASE_number_in_eventStation(string EQ_NAME, string STA)
{
	int return_num = 0;
	int count;
	for(count = 0; count < this->sta_num; count ++)
	{
	

		if(EQ_NAME.compare(this->my_record[count].EQ) != 0)
			continue;

		if(STA.compare(this->my_record[count].STA) != 0)
			continue;
		 // cout << this->my_record[count].EQ << " sta "<< this->my_record[count].STA << endl;
		double dist = this->my_record[count].DIST;
		//cout << " dist is " << dist << endl;
		// cout << "distminmax " << this->phase_dist_min << " " << this->phase_dist_max << endl;
		if ( dist < this->phase_dist_min || dist > this->phase_dist_max)
			continue;

		return_num = count;
		break;

	}
	return return_num;

}


int big_new_record::find_EQ_STA_bin_record_existance_num(new_record* my_record, int sta_num, int ilat_eq, int ilon_eq, int ilat_sta, int ilon_sta)
{

	int exist_record_num =0;
	int irecord = 0;
	for(irecord =0; irecord < sta_num; irecord++)
	{
		//cout << "--> Work on irecord :"<< irecord<< endl;
		int ieq = 0;
		int ista = 0;
		for(ieq = 0; ieq < my_record[irecord].EQ_BIN_NUM;ieq++)
			for(ista = 0; ista < my_record[irecord].STA_BIN_NUM;ista++)
			{
				//cout << " on "
				if(ilat_eq == my_record[irecord].EQ_BIN_LAT[ieq] 
						&& ilon_eq == my_record[irecord].EQ_BIN_LON[ieq]
						&& ilon_sta == my_record[irecord].STA_BIN_LON[ista]
						&& ilat_sta == my_record[irecord].STA_BIN_LAT[ista]
						)
					exist_record_num ++;
				break;

			}
	}

	return existing_sta_num;

}



void big_new_record::for_each_EQ_check_num_records_for_stacking(new_record* my_record, int sta_num, int ilat_eq, int ilon_eq, int ilat_sta, int ilon_sta  )
{

	int ieq = 0;
	for(ieq = 0; ieq < this->my_grid[ilat_eq][ilon_eq].CU_EQ_NUM;ieq++)
	{
		int current_EQ_record_num = 0;
		string current_eq = this->my_grid[ilat_eq][ilon_eq].CU_EQ_NAME[ieq];
		int irecord;
		for(irecord = 0; irecord < sta_num ; irecord ++)
		{
			string record_EQ = my_record[irecord].EQ;
			if( record_EQ != current_eq )
				continue;
			int ista;
			for(ista = 0; ista < my_record[irecord].STA_BIN_NUM;ista++)
			{
				if( ilon_sta == my_record[irecord].STA_BIN_LON[ista]
						&& ilat_sta == my_record[irecord].STA_BIN_LAT[ista]
						)
					current_EQ_record_num ++;

			}

		}
		if(current_EQ_record_num > this->VS_RECORD_NUM_THRESHOLD)
		{
			// stack for current EQ and STA_BIN


		}

	}

}






// for all possible EQ -STA pair, count the existance number in eventinfo

int big_new_record::loop_EX_EQ_STA_count_record( virtual_station* eq_grid, virtual_station* sta_grid)
{

	if( eq_grid->eq_skip_flag ==1  || sta_grid->sta_skip_flag == 1)
		return 0;
cout << "check for " << endl;




	int ieq = 0;
	int ista = 0;
	string eq_name = "";
	string sta_name = "";
	string phase_name = "";

	int return_value = 0;
	ifstream myfile;


/*
			string command = " cat eventinfo.EQ.STA.PHASE |grep -w " + eq_name + " |grep -w  " + sta_name + " |wc -l > tmp.exist_num";
			exec(command);
			myfile.open("tmp.exist_num");
			string input;
			myfile >> input;

			double exist_num = atof(input.c_str());
			if( exist_num != 0)
				return_value ++;
			//cout << "exist num is " << exist_num << endl;
*/


			int irecord = 0;
			for(irecord = 0; irecord < this->existing_sta_num; irecord++)
			{

				for(ieq = 0; ieq<eq_grid->EX_EQ_NUM ; ieq++)
					for(ista = 0; ista < sta_grid->EX_STA_NUM;ista++)
					{
						eq_name = eq_grid->EX_EQ_NAME[ieq];
						sta_name = sta_grid->EX_STA_NAME[ista];


							//if( this->existing_record[irecord].EQ != eq_name  )
								//continue;
							//if( this->existing_record[irecord].STA != sta_name  )
								//continue;
							//if( this->PHASE != this->existing_record[irecord].PHASE )
								//continue;
			//
								//return_value ++;
								//break;
					}
			}






	
	return return_value;

}


void big_new_record::catagorize_records_into_VS_based_on_each_record( new_record* my_record, int sta_num)
{

	int ilat = 0;
	int ilon = 0;
	int ista = 0;
	double eq_lat = 0;
	double eq_lon = 0;
	double sta_lat = 0;
	double sta_lon = 0;

	// initiate record EQ_STA bin
	int MAX = 10;

	for(ista = 0 ; ista < sta_num ; ista++)
	{
		my_record[ista].EQ_BIN_LAT = (int*)malloc(sizeof(int)*MAX);
		my_record[ista].EQ_BIN_LON = (int*)malloc(sizeof(int)*MAX);
		my_record[ista].STA_BIN_LAT = (int*)malloc(sizeof(int)*MAX);
		my_record[ista].STA_BIN_LON = (int*)malloc(sizeof(int)*MAX);
		my_record[ista].EQ_BIN_NUM = 0;
		my_record[ista].STA_BIN_NUM = 0;

		eq_lat = my_record[ista].eq_lat;
		eq_lon = my_record[ista].eq_lon;
		sta_lat = my_record[ista].sta_lat;
		sta_lon = my_record[ista].sta_lon;

		// find and store EQ info
		for(ilat = 0; ilat < this->grid_lat_num ; ilat++)
			for(ilon = 0; ilon < this->grid_lon_num[ilat] ; ilon++)
			{
				double grid_lat = this->my_grid[ilat][ilon].grid_lat;
				if( fabs(grid_lat - eq_lat) > this->VS_EQ_RADIUS_DEGREE)
					continue;

				double grid_lon = this->my_grid[ilat][ilon].grid_lon;
				// distance between EQ and grid
				double distance_km = dist_A_B( grid_lat, grid_lon, eq_lat, eq_lon);
				//cout << distance_km << endl;
				if( distance_km > 110*this->VS_EQ_RADIUS_DEGREE )
					continue;
					my_record[ista].EQ_BIN_LAT[my_record[ista].EQ_BIN_NUM] = ilat;
					my_record[ista].EQ_BIN_LON[my_record[ista].EQ_BIN_NUM] = ilon;
					my_record[ista].EQ_BIN_NUM ++;
			}
		// find and store STA info
		for(ilat = 0; ilat < this->grid_lat_num ; ilat++)
			for(ilon = 0; ilon < this->grid_lon_num[ilat] ; ilon++)
			{
				double grid_lat = this->my_grid[ilat][ilon].grid_lat;
				if( fabs(grid_lat - sta_lat) > this->VS_EQ_RADIUS_DEGREE)
					continue;

				double grid_lon = this->my_grid[ilat][ilon].grid_lon;
				// distance between EQ and grid
				double distance_km = dist_A_B( grid_lat, grid_lon, sta_lat, sta_lon);
				//cout << distance_km << endl;
				if( distance_km > 110*this->VS_STA_RADIUS_DEGREE )
					continue;
					my_record[ista].STA_BIN_LAT[my_record[ista].STA_BIN_NUM] = ilat;
					my_record[ista].STA_BIN_LON[my_record[ista].STA_BIN_NUM] = ilon;
					my_record[ista].STA_BIN_NUM ++;
			}
	}

	// output record related 
	ofstream myfile;
	myfile.open("output.record_EQ_STA_info");
	for(ista = 0; ista < sta_num;ista++)
	{
		myfile<<"Record_index: " << ista << "/" << sta_num << " EQ_NUM " << my_record[ista].EQ_BIN_NUM 
			<< " STA_NUM " << my_record[ista].STA_BIN_NUM << endl;
	}

	myfile.close();



}



// =========================================================
// Loop through each grid point
void big_new_record::catagorize_existing_eventinfo_to_VS()
{

	// step1, get unique EQ list
	string command;
	command = " cat eventinfo |awk '{print $12,$8,$9}'|sort|uniq |sort -u -k1,1>> list.unique_EQ";
	exec(command);
	// step2 get unique sta list
	command = " cat eventinfo |grep -v PPP |awk '{print $1,$6,$7}'|sort|uniq > list.unique_STA";
	exec(command);




	int EQ_NUM = 0;
	int STA_NUM = 0;

	ifstream myfile;
	string file_name = "list.unique_EQ";
	myfile.open(file_name);
	int file_count = 0;
	file_count = count_file_num(file_name);
	int count;
	for(count = 0; count < file_count ; count++)
	{
		myfile >> this->EQ_LIST[count] >> this->EQ_LAT[count] >> this->EQ_LON[count];
		this->unique_EQ_num++;
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
		myfile >> this->STA_LIST[count] >> this->STA_LAT[count] >> this->STA_LON[count];
		//cout << STA_LIST[count] << endl;
		//if( STA_LIST[count].find( "PPP") != std::string::npos )
		//{
			//cout << "found PPP " << endl;
			//continue;
		//}
		STA_NUM ++;
		this->unique_STA_num++;
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
				if( fabs(grid_lat - EQ_LAT[count]) > this->VS_EQ_RADIUS_DEGREE)
					continue;

				double grid_lon = this->my_grid[ilat][ilon].grid_lon;
				// distance between EQ and grid
				double distance_km = dist_A_B( grid_lat, grid_lon, EQ_LAT[count], EQ_LON[count]);
				//cout << distance_km << endl;
				if( distance_km > 110*this->VS_EQ_RADIUS_DEGREE )
					continue;


				// add EQ info into current grid
				int current_index = this->my_grid[ilat][ilon].EX_EQ_NUM ;
				this->my_grid[ilat][ilon].EX_EQ_LAT[current_index] = EQ_LAT[count];
				this->my_grid[ilat][ilon].EX_EQ_LON[current_index] = EQ_LON[count];
				this->my_grid[ilat][ilon].EX_EQ_NAME[current_index] = EQ_LIST[count];
				this->my_grid[ilat][ilon].EX_EQ_NUM ++;
				// cout << " ilat lon EX EQ NUM" << ilat << " " << ilon << " " << this->my_grid[ilat][ilon].EX_EQ_NUM << endl;
				// cout << " ilat lon " << this->my_grid[ilat][ilon].grid_lat << " " << this->my_grid[ilat][ilon].grid_lon <<  endl;
				// cout << " eventinf EQ lat lon " << EQ_LAT[count]<< " " << EQ_LON[count]<<  endl;
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
				if( fabs(grid_lat - STA_LAT[count]) > this->VS_STA_RADIUS_DEGREE)
					continue;

				double grid_lon = this->my_grid[ilat][ilon].grid_lon;
				// distance between EQ and grid
				double distance_km = dist_A_B( grid_lat, grid_lon, STA_LAT[count], STA_LON[count]);
				//cout << distance_km << endl;
				if( distance_km > 110*this->VS_STA_RADIUS_DEGREE )
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

	// output grid info for existing eventinfo
	ofstream myfile2;
	string filename = "out.existing_eventinfo.EQ_grid.STA_grid";
	myfile2.open(filename.c_str());
	int ilat = 0;
	int ilon = 0;
	for(ilat = 0; ilat < this->grid_lat_num ; ilat++)
		for(ilon = 0; ilon < this->grid_lon_num[ilat] ; ilon++)
		{
			int ex_EQ_num = this->my_grid[ilat][ilon].EX_EQ_NUM;
			int ex_sta_num = this->my_grid[ilat][ilon].EX_STA_NUM;
			if(ex_EQ_num > 0 || ex_sta_num > 0)
			{
				myfile2<< "GRID_LAT_LON: "<< ilat << " " << ilon << " lat_lon" << 
				this->my_grid[ilat][ilon].grid_lat<< " " << this->my_grid[ilat][ilon].grid_lon <<endl;
				myfile2<< "GRID_EQ_NUM: " << this->my_grid[ilat][ilon].EX_EQ_NUM<< endl;
				myfile2<< "GRID_STA_NUM: " << this->my_grid[ilat][ilon].EX_STA_NUM<< endl;
			}
		

		}
	myfile2.close();

}


void big_new_record::catagorize_eventstation_to_VS()
{
	cout << "-- > catagorize eventStation into VS " << endl;

	// step1, get unique EQ list
	string command;
	command = " cat eventStation |awk '{print $19,$11,$12}'|sort|uniq |sort -u -k1,1> list.unique_EQ";
	exec(command);
	// step2 get unique sta list
	command = " cat eventStation |grep -v PPP |awk '{print $1,$9,$10}'|sort|uniq > list.unique_STA";
	exec(command);


	int EQ_MAX = 2000;
	string EQ_LIST[EQ_MAX];
	double EQ_LAT[EQ_MAX];
	double EQ_LON[EQ_MAX];
	int EQ_NUM = 0;

	int STA_MAX = 20000;
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
				if( fabs(grid_lat - EQ_LAT[count]) > this->VS_EQ_RADIUS_DEGREE)
					continue;

				double grid_lon = this->my_grid[ilat][ilon].grid_lon;
				// distance between EQ and grid
				double distance_km = dist_A_B( grid_lat, grid_lon, EQ_LAT[count], EQ_LON[count]);
				//cout << distance_km << endl;
				if( distance_km > 110*this->VS_EQ_RADIUS_DEGREE )
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
				if( fabs(grid_lat - STA_LAT[count]) > this->VS_STA_RADIUS_DEGREE)
					continue;

				double grid_lon = this->my_grid[ilat][ilon].grid_lon;
				// distance between EQ and grid
				double distance_km = dist_A_B( grid_lat, grid_lon, STA_LAT[count], STA_LON[count]);
				//cout << distance_km << endl;
				if( distance_km > 110*this->VS_STA_RADIUS_DEGREE )
					continue;

				// add EQ info into current grid
				int current_index = this->my_grid[ilat][ilon].CU_STA_NUM ;
				this->my_grid[ilat][ilon].CU_STA_LAT[current_index] = STA_LAT[count];
				this->my_grid[ilat][ilon].CU_STA_LON[current_index] = STA_LON[count];
				this->my_grid[ilat][ilon].CU_STA_NAME[current_index] = STA_LIST[count];
				this->my_grid[ilat][ilon].CU_STA_NUM ++;
			}

	}

	// output grid info for existing eventinfo
	ofstream myfile2;
	string filename = "out.evetStation.EQ_grid.STA_grid";
	myfile2.open(filename.c_str());
	int ilat = 0;
	int ilon = 0;
	for(ilat = 0; ilat < this->grid_lat_num ; ilat++)
		for(ilon = 0; ilon < this->grid_lon_num[ilat] ; ilon++)
		{
			int cu_EQ_num = this->my_grid[ilat][ilon].CU_EQ_NUM;
			int cu_sta_num = this->my_grid[ilat][ilon].CU_STA_NUM;
			if(cu_EQ_num > 0 || cu_sta_num > 0)
			{
				myfile2<< "GRID_LAT_LON: "<< ilat << " " << ilon << " lat_lon" << 
				this->my_grid[ilat][ilon].grid_lat<< " " << this->my_grid[ilat][ilon].grid_lon <<endl;
				myfile2<< "GRID_EQ_NUM: " << this->my_grid[ilat][ilon].CU_EQ_NUM<< endl;
				myfile2<< "GRID_STA_NUM: " << this->my_grid[ilat][ilon].CU_STA_NUM<< endl;
			}


		}
	myfile2.close();



}
void big_new_record::virtual_station_grid_initiate()
{
	cout << "--> Initiate virtual station grid \n";
	virtual_station** my_grid;
	// initiate virtual station grid
	this->grid_lat_num = (int)(180 / this->VS_LATITUDE_INC);
	this->grid_lon_num = (int*)malloc(sizeof(int)*this->grid_lat_num);
	my_grid = (virtual_station**)malloc(sizeof(virtual_station*)*this->grid_lat_num);


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
		my_grid[ilat] = (virtual_station*)malloc(sizeof(virtual_station)*this->grid_lon_num[ilat]);
		
		for(ilon = 0; ilon < this->grid_lon_num[ilat] ; ilon++)
		{
			double lon_inc_in_degree = 360/this->grid_lon_num[ilat];
			double current_lon = ilon * lon_inc_in_degree - 180;
			//cout << " lat lon "<< current_lat << " "<< current_lon<< endl;
			my_grid[ilat][ilon].grid_lat = current_lat;
			my_grid[ilat][ilon].grid_lon = current_lon;
			//my_grid[ilat][ilon].grid_radius = this->VS_RADIUS_DEGREE;

			// initiate EXisting EQ STA 
			my_grid[ilat][ilon].EX_EQ_NUM = 0;
			int MMM = 1000;
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


	// also initiate EQ_LIST
	int max_EQ = 2000;
	int max_STA = 20000;
	this->unique_EQ_num = 0;
	this->EQ_LIST = (string*)malloc(sizeof(string)*max_EQ);
	this->EQ_LAT = (double*)malloc(sizeof(double)*max_EQ);	
	this->EQ_LON = (double*)malloc(sizeof(double)*max_EQ);

	this->unique_STA_num = 0;
	this->STA_LIST = (string*)malloc(sizeof(string)*max_STA);
	this->STA_LAT = (double*)malloc(sizeof(double)*max_STA);	
	this->STA_LON = (double*)malloc(sizeof(double)*max_STA);


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
	cout << "--> Read in eventStation "<< sta_num << " stations is read in " << endl;


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
	{
		cout << " -> Read sac file " << ista << "/" << this->sta_num << endl; 
		this->my_record[ista].PHASE = this->PHASE;
		this->my_record[ista].download_sac_file();
		this->my_record[ista].read_sac_file();
	}


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
