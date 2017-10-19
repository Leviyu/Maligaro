#include "../02_forward_tomo/forward_tomography.h"
#include "../01_cpp_lib/hongyulibcpp.h"

using namespace std;
int main()
{
	string PHASE;
	cin >> PHASE ;
	cout << "============= BEGIN" << endl;
	string infile="INFILE";

	big_new_record my_big_record;
	my_big_record.PHASE = PHASE;
	my_big_record.read_INFILE();
	my_big_record.initiate_big_record();
	my_big_record.virtual_station_main();
/*
	cin >> eventinfo_file >> PHASE>> S_ES_file;

	my_big_record.sta_num = count_file_num(eventinfo_file);


	my_big_record.eventStation= eventinfo_file;
	my_big_record.read_eventStation();

	// set phase 
	for(count =0; count < my_big_record.sta_num ; count++ )
		my_big_record.my_record[count].PHASE = PHASE;


	my_big_record.read_sac_file();
	my_big_record.calculate_SNR();
	my_big_record.read_in_polarity_file();


	cout << " 1. define grid file "<< endl;
	// 1. define grid file
	int grid_lat_min = -80;
	int grid_lat_max = 80;
	int grid_lat_delta = 8;
	int grid_lon_min = -180;
	int grid_lon_max = 180;
	int grid_lon_delta = 5;	// degree
	int station_num_threshold = 5;


	double grid_radius = 5;

	int npts_grid;
	int npts_grid_lat =  (grid_lat_max - grid_lat_min)/ grid_lat_delta +1;
	int npts_grid_lon = (grid_lon_max - grid_lon_min) / grid_lon_delta  +1;
	npts_grid = npts_grid_lat * npts_grid_lon ;

	new_grid my_grid[npts_grid_lat][npts_grid_lon];


	
	int ilat, ilon;
	for(ilat = 0; ilat < npts_grid_lat ; ilat++)
		for( ilon = 0; ilon < npts_grid_lon ; ilon++)
		{
			//cout << "initiate grid " << ilat << " " << ilon << endl;
			my_grid[ilat][ilon].initiate_grid();
			my_grid[ilat][ilon].grid_lat = grid_lat_min + ilat * grid_lat_delta;
			my_grid[ilat][ilon].grid_lon = grid_lon_min + ilon * grid_lon_delta;
			my_grid[ilat][ilon].grid_radius = grid_radius;
			my_grid[ilat][ilon].station_num_threshold = station_num_threshold;
			my_grid[ilat][ilon].long_npts = (int) (phase_long_win_len / delta);
			my_grid[ilat][ilon].my_big_record = &my_big_record;
			my_grid[ilat][ilon].LONG_BEG = phase_long_win_beg;
			my_grid[ilat][ilon].LONG_LEN = phase_long_win_len;
			my_grid[ilat][ilon].noise_beg = phase_noise_win_beg;
			my_grid[ilat][ilon].noise_len = phase_noise_win_len;
			my_grid[ilat][ilon].phase_beg = phase_win_beg;
			my_grid[ilat][ilon].phase_len = phase_win_len;

			my_grid[ilat][ilon].delta = delta;

			//cout << " work on ilat ilon" << ilat << " " << ilon << " grid radius  "<< my_grid[ilat][ilon].grid_radius << endl;
		}


	// 2. go through stations to find stations within range for each grid
	cout << " STEP 2 "<< endl;
	for(ilat = 0; ilat < npts_grid_lat ; ilat++)
		for( ilon = 0; ilon < npts_grid_lon ; ilon++)
		{
			my_grid[ilat][ilon].find_records_within_range();
			// 
			// we find the relocated grid center based on stations within range
			my_grid[ilat][ilon].relocate_grid_center();
			// find station within range again	
			//my_grid[ilat][ilon].find_records_within_range();
			//my_grid[ilat][ilon].relocate_grid_center();
			my_grid[ilat][ilon].get_grid_dist();
		}

	


	// 2. Iteratively go through all grid points and decide what records are in range
cout << "STEP 3" << endl;	
	ofstream myfile;
	string station_stack_file = "out.station_stack.info." + PHASE;
	myfile.open(station_stack_file.c_str());
	myfile<< " ilat ilon lat lon radius  num_records grid_distance ave_SNR stack_SNR"  << endl;

	for(ilat = 0; ilat < npts_grid_lat ; ilat++)
		for( ilon = 0; ilon < npts_grid_lon ; ilon++)
		{
			if(  my_grid[ilat][ilon].npts_record_sum <= station_num_threshold )
				continue;

			//// skip grid that has mid point in the north hemisphere
			//double current_eq_lat = my_big_record.my_record[2].eq_lat;
			//double current_eq_lon = my_big_record.my_record[2].eq_lon;
			//double current_grid_lat = my_grid[ilat][ilon].grid_lat;
			//double current_grid_lon = my_grid[ilat][ilon].grid_lon;
			//double current_DIST = my_grid[ilat][ilon].grid_dist;
			//double mid_lat = 0;
			//double mid_lon = 0;
			//double mid_az = 0;
			//mid_az = az_A_B(current_eq_lat, current_eq_lon, current_grid_lat, current_grid_lon);
			//point_AZ_dist_point( current_eq_lat, current_eq_lon, mid_az, current_DIST/2, &mid_lat, & mid_lon );
			//if( mid_lat > 0 )
				//continue;


			my_grid[ilat][ilon].ilat = ilat;
			my_grid[ilat][ilon].ilon = ilon;


			//cout << " --> stacking for ilat ilon"<< ilat << " "<< ilon << endl;
			// 3. Stack records within range for given window
			my_grid[ilat][ilon].stack_records_from_one_EQ();

			// 4. Output
			my_grid[ilat][ilon].out_stacked_record_rel_PREM = "out.station_stack."+PHASE+"." +to_string(ilat)+"."+to_string(ilon);

			// 5. get SNR
			my_grid[ilat][ilon].get_SNR_before_and_after_stack();

			//cout << "out stacked name is " <<  my_grid[ilat][ilon].out_stacked_record_rel_PREM << endl;
			my_grid[ilat][ilon].output_stacked_record();

			// find the ONSET time of virtual stack by t*(S.E.W)
			my_grid[ilat][ilon].S_ES_file = S_ES_file;
			my_grid[ilat][ilon].find_stack_ONSET_time();
			


			// lets define a record quality standard
			double ccc_threshold = 0.9;
			double SNR_threshold = 1.5;
			double tstar_factor_threshold = 1.5;

			if( my_grid[ilat][ilon].stack_SNR > SNR_threshold 
					&& my_grid[ilat][ilon].tstar_ccc > ccc_threshold 
					&& my_grid[ilat][ilon].tstar_factor > tstar_factor_threshold )
				my_grid[ilat][ilon].quality = 1;
			else 
				my_grid[ilat][ilon].quality = 0;
					



			myfile << " "
				<< fixed << setw(10) << ilat
				<< fixed << setw(10) << ilon
				<< fixed << setw(10) << setprecision(3)<< my_grid[ilat][ilon].grid_lat
				<< fixed << setw(10) << setprecision(3)<< my_grid[ilat][ilon].grid_lon
				<< fixed << setw(10) << setprecision(3)<< my_grid[ilat][ilon].grid_radius
				<< fixed << setw(10) << my_grid[ilat][ilon].npts_record_sum
				<< fixed << setw(10) << setprecision(3)<< my_grid[ilat][ilon].grid_dist
				<< fixed << setw(10) << setprecision(3)<< my_grid[ilat][ilon].ave_SNR
				<< fixed << setw(10) << setprecision(3)<< my_grid[ilat][ilon].stack_SNR
				<< fixed << setw(10) << setprecision(3)<< my_grid[ilat][ilon].virtual_stack_ONSET
				<< fixed << setw(10) << setprecision(3)<< my_grid[ilat][ilon].tstar_ccc
				<< fixed << setw(10) << setprecision(3)<< my_grid[ilat][ilon].tstar_factor
				<< fixed << setw(10) << setprecision(3)<< my_grid[ilat][ilon].gau_ccc
				<< fixed << setw(10) << setprecision(3)<< my_grid[ilat][ilon].gau_factor
				<< fixed << setw(10) << my_grid[ilat][ilon].quality
				<< endl;


				ofstream out;
				string grid_station_list_file = "out.grid_station_list."+ PHASE + "."+to_string(ilat)+"."+to_string(ilon);
				out.open(grid_station_list_file.c_str());
				for(count = 0; count < my_grid[ilat][ilon].npts_record_sum; count++)
				{
					int tag = my_grid[ilat][ilon].record_tag[count];
					string sta = my_big_record.my_record[tag].STA;
					out << sta << endl;
				}

				out.close();
				


			
		}
	




	// read in eventinfo
	//string eventinfo_file;
	//cin >> eventinfo_file;

	//int num_file;
	//num_file = count_file_num(eventinfo_file);
	//cout << " num of record is " << num_file << endl;

	//big_new_record my_big_record;
	//my_big_record.sta_num = num_file;



	//


	////// read record file
	//my_big_record.initiate_big_record();
	//my_big_record.record_file = eventinfo_file;
	//my_big_record.read_record_file();


	//// read in sac file 
	//my_big_record.read_sac_file();
	

*/





	return 0;
}
