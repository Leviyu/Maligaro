#ifndef __FORWARD_TOMOGRAPHY_H
#define __FORWARD_TOMOGRAPHY_H
// forward tomography data structure declaration here
#include<cstdlib>
#include<string>
#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
#include<iomanip>
#include "../01_cpp_lib/hongyulibcpp.h"
#include<cmath>
#include<complex>
//extern "C"{
//#include "ESF.h"
//}
using namespace std;


class big_new_record;
class new_record;
class new_tomo;
class CRUST_MODEL;
class new_grid;
class virtual_station;


// strucutre of cell
class new_cell
{
	public:
	double dvs;
	double STD;							// standard deviation 
	int idep, ilon, ilat;


	// find boundary


	// sumation of number of ray going through cell
	int 	sum_num_in_cell;
	double 	sum_dvs_in_cell;
	double 	delta_dvs;					// delta of updated_tomo - before_updated_tomo
	double 	gradient;					// the 3D gradient structure model
	double 	BAZ;
	int		num_BAZ;
	int 	vtk_index;
	int*	vtk_index_array;

	new_cell();
	~new_cell();


};

// Model Structure
class new_tomo
{

	public:
		string INFILE;
		string MODEL_INFILE;
		string MODEL_model;
		string MODEL;
		string MODEL_DIR;
		string dv_type;

		string Iteration_Scheme;
		int Iteration_M;
		int Iteration_N;
		int Iteration_MN;



		double weight_ratio_path_length_RMS_tomo;
		int		num_dep,num_lat,num_lon;
		int* num_lon2;
		double dep_min, dep_max, dep_delta;
		double lat_min, lat_max, lat_delta;
		double lon_min, lon_max, lon_delta;


		int dep_profile_num;
		double* depth_profile;

		// plotting parameter
		int plot_profile_num;
		double* plot_profile_dep;

		// iteration info
		int current_iteration;
		int current_layer_num;
		int tomo_iterations;
		int tomo_layer_number;
		double* tomo_layer_array;

		double current_layer_dep_min;		// depth of current layer min value
		double current_layer_dep_max;		// depth of current layer max value

		int CMB_lat_cell_number;

		string taup_path_dir;
		string cross_point_dir;

		double* dep;
		double* lat;
		double* lon;
		double** lon2;

		new_cell*** my_cell;



		double 	  layer_thickness;			// the layer_thickness to look for when averaging for dvs

		string tomo_model_file;				// netcdf tomography file

		void initiate_tomo();				// allocate space for model
		int read_INFILE();
		int read_MODEL_INFILE();
		int read_tomography();


		string output_tomo_file;
		int output_tomography_info();
		int output_tomography_info2();
		int output_tomography_info3();
		int output_starting_tomography();


		int convert_to_new_tomo(new_tomo* old_tomo);
		int define_vertical_profile();
		int define_horizontal_profile();


		double find_dvs_with_dep_lat_lon(double dep, double lat, double lon);
		int    find_index_lat(double lat);
		int    find_index_lon(double lon, int index_lat);
		int    find_index_dep(double dep,int* index_min, int* index_max);
		int    find_index_dep2(double dep);


		string flag_is_new_tomo;
		bool	is_new_tomo;

		int getDep(int idep);

		int print();

		int RMS_layer_num;
		double* RMS_dep;
		double* RMS_weight;
		double get_delta_lat_lon_distance_with_depth(int index_dep);


		void construct_RMS_profile();


		void pointA_pointB_find_cross_pointC(double dep1, double lat1, double lon1,
				double dep2, double lat2, double lon2,
				int index_dep1,int index_lat1, int index_lon1,
				int index_dep2,int index_lat2, int index_lon2,
				double* mid_dep, double* mid_lat, double* mid_lon,
				double* delta_angle_d1_d2, double* radius);

		void forward_tomography_func(big_new_record* my_big_record);
		void forward_tomo_for_station(big_new_record* my_big_record);
		int forward_tomo_for_one_station(new_record* my_record);
		int distribute_dt_residual_and_convert_to_dvs(new_record* record);
		void update_tomo_for_current_iteration();
		void output_vtk_format_model();
		void output_time_info(big_new_record* my_big_record);
		void output_record_path_crosssection(big_new_record* my_big_record);
		void output_record_path_line(big_new_record* my_big_record);
		void output_record_path_line_in_one_file(big_new_record* my_big_record);
		void output_coverage_hit_count();
		void output_delta_tomography(new_tomo* my_tomo);
		void output_gradient_tomography();
		void output_vertical_gradient_tomograph();
		void output_STD_tomography();
		void output_smoothed_tomography();
		void S1_update(new_record* my_record);
		//void S2_update(big_new_record* my_big_record);
		//void S3_update(big_new_record* my_big_record);
		//void S4_update(big_new_record* my_big_record);
		//void S5_update(big_new_record* my_big_record);
		void calculate_turning_depth(new_record* my_record);
		void output_cross_section(big_new_record* my_big_record);
		void output_cross_section_for_one_record(new_record* my_record);
		//void output_cross_section_hit_count(big_new_record* my_big_record);
		//void output_cross_section_for_one_record_hit_count(new_record* my_record);
		double find_hit_with_dep_lat_lon(double dep, double lat , double lon);

		int CONVERT_write_to_netcdf(char* netcdf_file);
		int construct_1D_reference();
		void output_RMS_profile();
		void output_tomo_for_media();

		new_tomo();
		~new_tomo();
};



// Record Structure
class new_record
{
	public:
		string EQ;
		string STA;
		string NET;
		string PHASE;
		string COMP;


		double eq_lat;
		double eq_lon;
		double eq_dep;
		double sta_lat;
		double sta_lon;
		double eq_mag;
		double DIST;
		double AZ;
		double BAZ;
		int polarity_flag;	
		int  quality_flag;			// 1 is good -1 is bad 0 is not known
		double CCC;
		double SNR;
		double stretch_ccc;
		double stretch_coeff;
		double misfit;
		double phase_amplitude;

		double dt_obs_prem;				// travel time anomaly relative to PREM from Empirical Waveform
		double PREM;
		double code_PREM;
		double PREM_tomo;
		double dt_tomo_PREM;
		//double dt_residual_before;		// dt(T_obs - T_tomo)
		double dt_residual_after;		// dt(T_obs - T_new_tomo)
		double dt_residual_for_current_iteration;	// dt(T_obs - T_tomo)

		double turning_depth;
		double dt_tomo_correction;
		string record_file;


		double phase_peak_time_rel_PREM;
		int 	npts_phase_peak_rel_start;
		double 	noise_beg;
		double 	noise_len;
		double phase_beg;
		double phase_len;
		double  long_beg;
		double  long_len;
		double* long_win;
		string read_sac_flag;		
		double delta;
		double phase_beg_rel_PREM;		// phase begin time relative to PREM
		double record_weight;
		double SNR2;
		double misfit2;
		double ONSET;
		double ENDSET;
		double ellip_corr;			// ellipticity correction for current record
		double incident;			//incident angle

		double current_iteration_coeff;

		// store info for taup_path
		double* angle;
		double* radius;
		int taup_path_max_num;

		// here is the array to store the cross-points for each record
		int CP_num;
		double* CP_lat;
		double* CP_lon;
		double* CP_dep;
		int* CP_ilat;
		int* CP_ilon;
		int* CP_idep;
		double*	CP_dl;
		double*	CP_v_PREM;			// PREM velocity 
		double* CP_weight;			// weight from CCC and SNR
		double* current_CP_weight;
		void initiate_CP();
		void free_CP();

		double crust_correction;

		// extra info read in
		// taup_path and cross-point path
		string taup_path_file;
		string cross_point_file;

		void read_cross_point_info( new_tomo* my_tomo);
		void read_taup_path_info(string taup_path_dir);
		void find_cross_points(new_tomo* my_tomo);

		void calculate_tomo_correction(new_tomo* my_tomo);
		void get_ellip_corr();
		void read_sac_file_relative_to_PREM();
		void read_sac_file();
		void calculate_SNR();
		void get_crustal_correction();
		void get_incident();
		string sac_file;

		CRUST_MODEL* my_crust;




		new_record();
		~new_record();

};

// declare big_new_record
class big_new_record
{
	public:
		big_new_record();
		~big_new_record();

		double VS_LATITUDE_INC;


		int sta_num;
		string record_file;
		string eventStation;
		new_record* my_record;
		double delta;
		double long_len;
		double long_beg;
		double noise_beg;
		double noise_len;
		double phase_beg;
		double phase_len;


		string timeinfo_outfile;
		
		void read_record_file();
		void initiate_big_record();
		void read_INFILE();


		void big_record_read_cross_point_file(new_tomo* my_tomo);
		void get_ellip_corr();
		void get_incident();
		void read_sac_file_and_store_relative_to_phase();
		void read_eventStation();
		void get_crustal_correction();

		void read_sac_file();
		void calculate_SNR();
		void read_in_polarity_file();
		new_grid* my_grid;


		// grid related parameter
		int grid_lat_num;
		int* grid_lon_num;


};


// declare virtual_station
class virtual_station
{

};



class new_grid
{
	public:
		// grid basic information
		double grid_lat;
		double grid_lon;
		int ilat;
		int ilon;
		double grid_radius;
		double grid_height;
		double grid_time;					// centra grid travel time for one EQ, as a reference time

		double VS_LATITUDE_INC;

		// for records within grid range, store record tag info
		int npts_record_sum;
		int* record_tag;				// store the record line number in eventinfo to make sure that we can find each record

		// for given EQ, when we do stacking, we need az baz slowness info
		double grid_dist;					// distance from EQ to grid center
		double AZ;							//az from EQ to grid center
		double BAZ;
		double incident_angle;
		double PREM;						//PREM time from EQ to grid center
		// the stacking window
		double delta;
		double WIN_BEG;
		double WIN_LEN;
		double LONG_BEG;
		double LONG_LEN;
		double noise_beg;
		double noise_len;
		double phase_beg;
		double phase_len;
		int     long_npts;
		double* long_win;			// the stacked record for current grid
		int station_num_threshold;			//  if station number in range is not greater then this value, grid is skipped
		double ave_SNR;						// average SNR of records within range
		double stack_SNR;
		string S_ES_file;
		int find_stack_ONSET_time();
		double virtual_stack_ONSET;
		double tstar_ccc;
		double tstar_factor;
		double gau_ccc;
		double gau_factor;

		int quality;

		// when we do stacking, arrays to store time-slowness-amplitude info
		
		double dt_ave;
		double dt_STD;
		double crustal_correction_ave;
		double crustal_correction_STD;



		double* fix_BAZ_time;
		double* fix_BAZ_slowness;
		double* fix_BAZ_amp;

		double* fix_slow_time;
		double* fix_slow_BAZ;
		double* fix_slow_amp;

		double fix_BAZ_delay_time;
		double fix_slow_delay_time;

		big_new_record* my_big_record;
		void find_records_within_range();
		void stack_records_from_one_EQ();
		void output_stacked_record();
		void get_dt_ave_STD();
		void get_crust_correction_ave_STD();

		void relocate_grid_center();
		void get_SNR_before_and_after_stack();
		void get_grid_dist();


		new_grid();
		~new_grid();

		void initiate_grid();

		// station stacking usage
		string 	out_stacked_record_rel_PREM;
		string 	out_stacked_record_raw;
		double* record_win_rel_PREM;
		double* record_win_rel_PREM_shifted;
		double* record_win_raw;
		double* record_win_raw_shifted;

};






#endif
