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
	  this->my_record = new new_record[this->sta_num]; 

}



/*************************************************************
* This C function read in station information 
* for a certain phase
*	INPUT:
*	
*	Hongyu DATE: Aug 2016
*	Key words: read in station file
*************************************************************/
void big_new_record::read_record_file()
{

	int count;
	string tmp;
	ifstream myfile;
	myfile.open(this->record_file.c_str());
		string sub1;
		string flag;

	int line;
	line = 0;

	cout <<"--> big_record  Read in "<< this->sta_num << " records" << endl;

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
			{
				this->my_record[line].STA = sub1;
			}
			else if( count == 2 )
				this->my_record[line].NET = sub1;
			else if( count == 3 )
				this->my_record[line].DIST = atof(sub1.c_str());
			else if( count == 4 )
				this->my_record[line].AZ = atof(sub1.c_str());
			else if( count == 5 )
				this->my_record[line].BAZ = atof(sub1.c_str());
			else if( count == 6 )
				this->my_record[line].sta_lat = atof(sub1.c_str());
			else if( count == 7 )
				this->my_record[line].sta_lon = atof(sub1.c_str());
			else if( count == 8 )
				this->my_record[line].eq_lat = atof(sub1.c_str());
			else if( count == 9 )
				this->my_record[line].eq_lon = atof(sub1.c_str());
			else if( count == 10 )
				this->my_record[line].eq_dep = atof(sub1.c_str());
			else if( count == 11)
				this->my_record[line].eq_mag = atof(sub1.c_str());
			else if( count == 12 )
				this->my_record[line].EQ = sub1;
			else if( count == 13 )
				this->my_record[line].polarity_flag = atoi(sub1.c_str());
			else if( count == 14 )
				this->my_record[line].quality_flag  = atoi(sub1.c_str());
			else if( count == 15 )
				this->my_record[line].PREM = atof(sub1.c_str());
			else if( count == 16 )
				this->my_record[line].phase_amplitude = atof(sub1.c_str());
			else if( count == 17 )
				this->my_record[line].CCC = atof(sub1.c_str());
			else if( count == 18 )
				this->my_record[line].SNR = atof(sub1.c_str());
			else if( count == 19 )
				this->my_record[line].dt_obs_prem  = atof(sub1.c_str());
			else if( count == 20 )
				this->my_record[line].PHASE = sub1;
			else if( count == 21 )
				this->my_record[line].stretch_ccc = atof(sub1.c_str());
			else if( count == 22 )
				this->my_record[line].stretch_coeff = atof(sub1.c_str());
			else if( count == 23 )
				this->my_record[line].misfit = atof(sub1.c_str());
			else if( count == 24 )
				this->my_record[line].COMP = sub1;
			else if( count == 25 )
				this->my_record[line].phase_peak_time_rel_PREM = atof(sub1.c_str());
			else if( count == 26 )
				this->my_record[line].npts_phase_peak_rel_start = atoi(sub1.c_str());
			else if( count == 27 )
				this->my_record[line].noise_beg = atof(sub1.c_str());
			else if( count == 28 )
				this->my_record[line].noise_len = atof(sub1.c_str());
			else if( count == 29 )
				this->my_record[line].phase_beg_rel_PREM = atof(sub1.c_str());
			else if( count == 30 )
				this->my_record[line].record_weight = atof(sub1.c_str());
			else if( count == 31 )
				this->my_record[line].SNR2 = atof(sub1.c_str());
			else if( count == 32 )
				this->my_record[line].misfit2 = atof(sub1.c_str());
			else if( count == 33 )
				this->my_record[line].ONSET = atof(sub1.c_str());
			else if( count == 34 )
				this->my_record[line].ENDSET = atof(sub1.c_str());
			//else if( count == 35 )
				//this->my_record[line].STA = sub1;
			//else if( count == 36 )
				//this->my_record[line].STA = sub1;
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