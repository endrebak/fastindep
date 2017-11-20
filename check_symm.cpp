// This is an auxilliary program for checking the symmetry of the data matrix used 
// in FastIndep.   
// Author: Joseph Abraham November 2013

#include <istream>
#include <ostream>
#include <sstream>
#include <iostream>
#include <fstream>
#include <ostream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>
#include <limits>
#include <cstdlib>
#include <time.h>
#include <cerrno>
#include <algorithm>
#include <functional>

using namespace std;


int main(int argc, char* argv[])
{
   if (argc < 2) {
    cout << " Must specify datafile " << endl;
    cout << " program aborting " << endl;
    return 0;
   }

	clock_t t_begin,t_end; 
	t_begin = clock();


   unsigned j,linenum;
	
   stringstream strstr;
   string this_ent,s;
	

   string homedir = "./";
   string datafile = homedir + argv[1]; // original data file
   ifstream gg(datafile.c_str());
    if (!gg) 
    {
     cout << " File " << datafile << " not found ... program aborting " << endl;
     return 0;
    }
	

	linenum = 1;
		
	
	cout << "Reading in data from file " << datafile << endl;
	
	linenum = 0;
	vector<string> names;
	names.clear();
	vector<unsigned> sep_count,tempseps;
	vector<unsigned>::iterator first_pos;
	vector<vector<double> > all_dat;
	vector<double> this_line;
	this_line.clear();
	all_dat.clear();
	unsigned numents = 0;

    while(getline(gg,s)){
	 linenum++; // linenum is a counter for lines
	  strstr.clear();
	  strstr.str(s);
	   if(linenum == 1){ // first line has names only 
		   while(strstr >> this_ent){
			   names.push_back(this_ent);
		   } // Initialize the names 
		  numents = names.size();
		  this_line.resize(numents,numeric_limits<double>::quiet_NaN());
		   all_dat.resize(numents);
		   for(j =0; j != numents; j++){
			   all_dat[j] = this_line;
		   }
        }
	    if(linenum == 1) continue; // have read in all the names of elements		
		j = 0; // index as we move across each line 
		this_line.resize(numents,numeric_limits<double>::quiet_NaN());
		while(strstr >> this_ent){
			if(j == 0){ // name check 
				if(this_ent != names[linenum-2]){
					cout << "Discrepancy in names in line " << linenum << " of datafile :: " << 
					this_ent << " " << names[linenum-2];
					cout << endl << "Program Terminating " << endl;
					exit(-1);
				}
			}else{
				 if(j != (linenum -1)){ // skip diagonal element 
				   this_line[j-1] = atof(this_ent.c_str());
					 if(j < (linenum -1) ){ // left of diagonal elements check all others 
					  if(this_line[j-1] != all_dat[j-1][linenum -2]){
						cout << "Lack of symmetry in data " << endl;
						cout << "Row " << (linenum -1) << " element " << j << " = " << this_line[j-1] ;
						cout << ";  Row " << j << "  " << " element " << (linenum -1) <<  " = "
						  << all_dat[j-1][linenum -2] << endl;
						cout << "Program Terminating " << endl;
						exit(-1);
					   }
			     }
			 }
		  }
		 j++;
		} // end of line 
		if(j != (numents+1) ){
		    cout << " Expecting " << (numents+1) << " entries found " << j << "  ";
			cout << " in line " << linenum << " program terminating " << endl;
			throw exception();
		}
		all_dat[linenum -2] = this_line;
	}

	gg.close();

   cout << "Finished processing data " << endl; 
	
	
   cout << " Data read in from file " << datafile << " with " << (numents +1) << " lines " << endl;
	cout << " Data symmetric " << endl;
   
	
	t_end = clock();
	cout << " Number of Processor Clock Ticks is " << double(t_end -t_begin) << endl;
	cout << " Number of seconds is " << double(t_end -t_begin)/CLOCKS_PER_SEC << endl;
	
    cout << " Program terminating " << endl;
	
	return 0;
} 

	

