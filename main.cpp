// This is the main program for FastIndep, a C++ program for 
// finding maximal independent sets given an input matrix and a threshold. 
// This program may be freely distributed or modified, provided these 
// comment lines are retained. In particular, there are no restrictions
// on non-commercial use. If this code (either in the original form or 
// modified) is used in academic publications, a citation would be appreciated.  
// Please read the accompanying ReadMe.txt file for instructions on how to compile 
// and run the program.  
// Author: Joseph Abraham May 2013

#include "./headers.h"
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

bool setcomp(vector<unsigned> vec1, vector<unsigned> vec2)
{
  if (vec1.size() < vec2.size()) 
     {
      return true;
     }
   else 
     {
      return false;
       } 
}

int main(int argc, char* argv[])
{
   if (argc < 3) {
    cout << " Must specify datafile and parameter input file " << endl;
    cout << " program aborting " << endl;
    return 0;
   }

	clock_t t_begin,t_end; 
	t_begin = clock();


   unsigned i,j,linenum;
	
   stringstream strstr;
   string this_ent,s;
	
   double lowLDlimit = numeric_limits<double>::quiet_NaN();
   bool seed_specified = false;
   int num_sweeps = 0;
   unsigned long longseed = 0;
	

   string homedir = "./";
   string datafile = homedir + argv[1]; // original data file
   ifstream gg(datafile.c_str());
    if (!gg) 
    {
     cout << " File " << datafile << " not found ... program aborting " << endl;
     return 0;
    }
	

	linenum = 1;
	
	string param_input_file = homedir + argv[2];
    cout << " Reading in data from file " << param_input_file << endl;
	ifstream param_input_file_stream(param_input_file.c_str());
	if (!param_input_file_stream) 
    {
		cout << " File " << param_input_file << " not found ... program aborting " << endl;
		return 0;
    }
    while(getline(param_input_file_stream,s)){
		strstr.clear();
		strstr.str(s);
		j = 0;
		while(strstr >> this_ent){
			j++;
			if(j > 1) continue;
			if(linenum == 1){
				lowLDlimit = (double) atof(this_ent.c_str());
				if(lowLDlimit <= 0) 
				{ 
					cout << lowLDlimit << " not acceptable cut-off " << endl;
					cout << " Check parameter.txt file. Program aborting " << endl; 
					throw exception();
				} 
			} else if(linenum == 2){
				num_sweeps = (unsigned) atof(this_ent.c_str());
				if(num_sweeps <= 0) 
				{
					cout <<  num_sweeps << " Defective value for number of sweeps " << endl;
					cout << " Check parameter.txt file. Program aborting " << endl; 
					throw exception();
				} 
			} else if(linenum == 3){
				longseed = atol(this_ent.c_str()); 
				if(longseed > 0 )
				{
					seed_specified = true;
				}else{ 
					cout << " Invalid choice of seed for random number generator " << endl;
					cout << " Will use program generated seed instead " << endl;
					seed_specified = false;
				}
			} else {continue;}
			linenum++;
		}
	}
			 

// Define data object 
	
	Data OurData;
	OurData.NameIdMap.clear();
	OurData.SNP_vector.clear();
	OurData.coveredVec.clear();
	OurData.allIds.clear();
	OurData.singletons.clear();	
	
	vector<unsigned> singletons;
	singletons.clear();
	vector<unsigned> all_connectors;
	all_connectors.clear();
	set<unsigned>::iterator unsigit;
		
	// Data object defined
	// Transfer File information to Data object
	// First line of datafile is names only 
	// Subsequent lines contain association between elements 
	
	//	string outprefix = "indep_set_";
	//	string outputfile = outprefix + argv[1];
	string outputfile = "outfile.txt";
	ofstream outfile(outputfile.c_str());
	
	
	cout << " Upper limit of association between separated elements is " << lowLDlimit << endl;
    cout << " " << num_sweeps << " sets of unassociated elements requested" << endl << endl;

	
	cout << "Reading in data from file " << datafile << endl;
	
	linenum = 0;
	vector<string> names;
	names.clear();
	vector<unsigned> sep_count,tempseps;
	vector<unsigned>::iterator first_pos;
	double currentval;
	vector<unsigned> search_elem;
	search_elem.push_back(0); // important later 
	unsigned newsize;
	unsigned numents = 0;
	unsigned largest_name_size = 0;
	set<unsigned> tempSet;
	
	unsigned num_edges = 0;

    while(getline(gg,s)){
	 linenum++;
	  strstr.clear();
	  strstr.str(s);
	   if(linenum == 1){ // first line has names only 
		   while(strstr >> this_ent){
			   names.push_back(this_ent);
		   } // Inditialize the data object 
		  numents = names.size();
		  OurData.SNP_vector.resize(numents);
		   OurData.allIds.resize(numents);
		   for (j = 0; j != numents; j++){
			SNP* newSNP = new SNP;
			 newSNP->SNPId = (j +1);
			  newSNP->SNPname = names[j];
			   if(names[j].size() > largest_name_size) largest_name_size = names[j].size();
			   newSNP->seps.clear(); 
				OurData.SNP_vector[j] = newSNP;
				 OurData.NameIdMap[newSNP->SNPname]= (j+1);
				  OurData.allIds[j] = (j+1);
			}
		  names.clear();
		  sep_count.resize(numents,0); // zero is not a valid node id 
        }
	    if(linenum == 1) continue; // have read in all the names of elements		
		j = 0;
		while(strstr >> this_ent){
			if(j == 0){ // name check 
				if(this_ent != (OurData.SNP_vector[linenum-2])->SNPname){
					cout << "Discrepancy in names entry in line " << linenum << "of datafile :: " << 
					this_ent << " " << OurData.SNP_vector[linenum-2]->SNPname;
					cout << " Program Terminating " << endl;
					exit(-1);
				}
			}else{
				 if(j != (linenum -1)){ // skip diagonal element 
					 currentval = atof(this_ent.c_str());
					 if(currentval <= lowLDlimit){
						 sep_count[j-1] = j;
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
		sort(sep_count.begin(), sep_count.end(), greater<unsigned>() ); // descending sort 	
		first_pos = find_first_of(sep_count.begin(), sep_count.end(), 
					 search_elem.begin(), search_elem.end() ); // first occurrence of zero 		
		newsize = unsigned(first_pos - sep_count.begin() );	// number of non zero entries 
		if(newsize > 0){
			tempseps.resize(newsize,0);
			copy(sep_count.begin(), first_pos, tempseps.begin() );
			sort(tempseps.begin(), tempseps.end() );
			OurData.SNP_vector[linenum-2]->seps = tempseps;
			if(newsize == (j -2) ) singletons.push_back(linenum -1); // (j-2) b/c skip name & diag
		}else{ // all entries are zero no distantly related nodes 
			all_connectors.push_back(linenum -1);
		}
		num_edges += (numents -1 -newsize);
		sep_count.clear();
		sep_count.resize(numents,0);
	}

	gg.close();


// end of data entry
   cout << "Finished processing data " << endl; 
	
   outfile << " Data read in from file " << datafile << endl;
   outfile << " Total number of elements in dataset = "<< numents << endl;
   outfile << " "<< num_sweeps << " sets of elements with pairwise association <= ";
   outfile << lowLDlimit << " requested" << endl ;

	
    cout << " Number of singletons is " << singletons.size() << endl;
	cout << " Number related to all " << all_connectors.size() << endl;
	
	outfile << " Number of singletons is " << singletons.size() << endl;
	outfile << " Number related to all " << all_connectors.size() << endl;
	
	double connectivity = double(num_edges)/(double(numents*(numents -1.0)));
	
	outfile << " Connectivity = " << connectivity << endl;

	
	if(singletons.size() == OurData.allIds.size() ){
		cout << " Data contains only singletons " << endl;
		cout << " Program exiting " << endl;
		outfile << " Data contains only singletons " << endl;
		outfile << " Program exiting " << endl;
		outfile << "-----------------------------------------------------------------" << endl;
		return 0;
	}
	
	if(all_connectors.size() == OurData.allIds.size() ){
		cout << " Data contains one large clique " << endl;
		cout << " Program exiting " << endl;
		outfile << " Data contains one large clique " << endl;
		outfile << " Program exiting " << endl;
		outfile << "-----------------------------------------------------------------" << endl;
		return 0;
	}
	
	if(singletons.size() > 0){
		sort(singletons.begin(), singletons.end() );
		OurData.singletons.resize(singletons.size() );
		copy(singletons.begin(), singletons.end(), OurData.singletons.begin() );
	}
	
	if(all_connectors.size() > 0){
		sort(all_connectors.begin(), all_connectors.end() );
		OurData.all_connectors.resize(all_connectors.size() );
		copy(all_connectors.begin(), all_connectors.end(), OurData.all_connectors.begin() );
	}
	
	
   sort(OurData.allIds.begin(), OurData.allIds.end() );
	tempSet.clear();
	
	unsigned datasize = OurData.SNP_vector.size();
	
	
    MTRand main_rand;
	
    if (!seed_specified) 
    {
		longseed = main_rand.randInt();
		main_rand.seed(longseed);
    } else 
	{ 
		main_rand.seed(longseed);
	}
	
	if (!seed_specified){
		cout << "Program generated seed for Mersenne Twister is " << longseed << endl;
	}else{
		cout << "Using supplied seed  "<< longseed << " for Mersenne Twister is " << endl;
	} 
	
	if (!seed_specified){
		outfile << " Program generated seed for Mersenne Twister is " << longseed << endl;
	}else{
		outfile << " Using user supplied seed  "<< longseed << " for Mersenne Twister" << endl;
	} 
	
	OurData.longseed = longseed;
	OurData.initialize_rand_gen();
		
    cout << " Calculating "<< num_sweeps << " sets of unassociated elements " << endl ;
    vector<vector<unsigned> > vecSNPIds;
    vecSNPIds.clear();
	map<unsigned,unsigned> size_dist; // distribution of output set sizes
	size_dist.clear();
	map<unsigned,unsigned> Id_freq; // distribution of output elements 
	Id_freq.clear();
	map<unsigned,unsigned>::iterator mapit; 
	bool found_already = false;
	bool greedy_found_already = false;
	
    if (num_sweeps > 1)
    {
		vector<unsigned> ReturnVec;
		ReturnVec.clear();
		for(j = 0; j != (unsigned) (num_sweeps-1); j++){
			OurData.data_rand = true;
			found_already = false;
			OurData.FindSeparatedVec(ReturnVec);
			if (ReturnVec.size() <= 0) cout << " cannot be ! " << endl;
			cout << " Found " << (j+1) << " separated sets " << endl;
			sort(ReturnVec.begin(), ReturnVec.end() );
			if(j > 0){ // Have we found this one already ?? 
				for(i = 0; i != vecSNPIds.size(); i++){
					if(ReturnVec.size() != vecSNPIds[i].size() ) continue;
					found_already = equal(ReturnVec.begin(),ReturnVec.end(),vecSNPIds[i].begin());
					if(found_already) break;
				}
			}
			if(found_already) continue; // Found this one already 
			vecSNPIds.push_back(ReturnVec);
			if(OurData.coveredVec.size() != OurData.SNP_vector.size() ) 
			{ 
				cout << " Size of Independent Set is " << ReturnVec.size() << endl;
				cout << "Covered " << OurData.coveredVec.size() << endl;
				cout << "Full " << OurData.SNP_vector.size() << endl;
				set<unsigned> disc_set;
				disc_set.clear();
				set_difference(OurData.allIds.begin(),OurData.allIds.end(),\
							   OurData.coveredVec.begin(), OurData.coveredVec.end(),\
							   inserter(disc_set,disc_set.begin() ) );
				cout << " Missing SNPs : " << endl;
				for(unsigit = disc_set.begin(); unsigit != disc_set.end(); unsigit++){
					cout << OurData.SNP_vector[(*unsigit)-1]->SNPname << endl;
				}
				cout << "Inconsistent set of SNPs ... program aborting" << endl;
				throw exception();
			}
			if(singletons.size() > 0){
				tempSet.clear();			
				sort(ReturnVec.begin(), ReturnVec.end() );
				set_intersection(ReturnVec.begin(), ReturnVec.end(), singletons.begin(), \
								 singletons.end(), inserter(tempSet,tempSet.begin() ) );
				if( (ReturnVec.size() < singletons.size())||(tempSet.size() < singletons.size()) ){
					set<unsigned> diff_set; diff_set.clear();
					cout << "Size of Independent Set found " << ReturnVec.size() << endl;
					cout << " Number of singletons " << singletons.size() << endl;
					cout << " Number of Common Elements " << tempSet.size() << endl;
					cout << " Following singletons missing from Return Set: " << endl; 
					set_difference(singletons.begin(), singletons.end(),tempSet.begin(),\
								   tempSet.end(),inserter(diff_set,diff_set.begin() ) );
					for(unsigit = diff_set.begin(); unsigit != diff_set.end(); unsigit++){
						cout << OurData.SNP_vector[(*unsigit)-1]->SNPname << endl;
					}
					cout << "Inconsistent set of SNPs ... program aborting" << endl;
					throw exception();				
				}
			}	
			mapit=size_dist.find(ReturnVec.size());
			if(mapit == size_dist.end() ){
				size_dist[ReturnVec.size()] = 1;
			}else{
				size_dist[ReturnVec.size()]++;
			}
			for(unsigned jj = 0; jj != ReturnVec.size(); jj++){
				mapit = Id_freq.find(ReturnVec[jj]);
				if(mapit == Id_freq.end() ){
					Id_freq[ReturnVec[jj]] =1;
				}else{
					Id_freq[ReturnVec[jj]]++;
				}
			}
			ReturnVec.clear();
		}
		stable_sort(vecSNPIds.begin(), vecSNPIds.end(), setcomp);
	} 
	
	// now the greedy set 
	OurData.data_rand = false;
	vector<unsigned> GreedyVec;
	GreedyVec.clear();
	cout << " Finding Greedy Set " << endl;
	OurData.FindSeparatedVec(GreedyVec);
	cout << " Found Greedy Set " << endl;
	if (OurData.coveredVec.size() != OurData.SNP_vector.size() ) 
	{ 
		cout << " Size of Independent Set is " << GreedyVec.size() << endl;
		cout << "Covered " << OurData.coveredVec.size() << endl;
		cout << "Full " << OurData.SNP_vector.size() << endl;
		cout.flush();
		set<unsigned> disc_set;
		disc_set.clear();
		set_difference(OurData.allIds.begin(),OurData.allIds.end(),\
					   OurData.coveredVec.begin(), OurData.coveredVec.end(),\
					   inserter(disc_set,disc_set.begin() ) );
		cout << " Missing SNPs : " << endl;
		for(unsigit = disc_set.begin(); unsigit != disc_set.end(); unsigit++){
			cout << (*unsigit) << "  " << OurData.SNP_vector[(*unsigit)-1]->seps.size() << endl;
		}
		cout << "Inconsistent set of SNPs ... program aborting" << endl;
		throw exception();
	}		
	if(singletons.size() > 0){
		tempSet.clear();	
		sort(GreedyVec.begin(), GreedyVec.end() );
		set_intersection(GreedyVec.begin(), GreedyVec.end(), singletons.begin(), \
						 singletons.end(), inserter(tempSet,tempSet.begin() ) );
		if( (GreedyVec.size() < singletons.size())||(tempSet.size() < singletons.size()) ){
			set<unsigned> diff_set; diff_set.clear();
			cout << "Size of Independent Set found " << GreedyVec.size() << endl;
			cout << " Number of singletons " << singletons.size() << endl;
			cout << " Number of Common Elements " << tempSet.size() << endl;
			cout << " Following singletons missing from Return Set: " << endl; 
			set_difference(singletons.begin(), singletons.end(),tempSet.begin(),\
						   tempSet.end(),inserter(diff_set,diff_set.begin() ) );
			for(unsigit = diff_set.begin(); unsigit != diff_set.end(); unsigit++){
				cout << OurData.SNP_vector[(*unsigit)-1]->SNPname << "  ";
				cout << OurData.SNP_vector[(*unsigit)-1]->SNPId << endl;
			}
			cout << "Inconsistent set of SNPs ... program aborting" << endl;
			throw exception();				
		}
	}
	if(num_sweeps > 1){
		greedy_found_already = false;
                sort(GreedyVec.begin(), GreedyVec.end() );
		for(i = 0; i != vecSNPIds.size(); i++){
			if(GreedyVec.size() != vecSNPIds[i].size() ) continue;
                        greedy_found_already = equal(GreedyVec.begin(),GreedyVec.end(),vecSNPIds[i].begin());
			if(greedy_found_already) break;
		}
	}
	cout << " Greedy Set Consistent " << endl;
	
        if(!greedy_found_already){
	  mapit=size_dist.find(GreedyVec.size());
	  if(mapit == size_dist.end() ){
		size_dist[GreedyVec.size()] = 1;
	  }else{
		size_dist[GreedyVec.size()]++;
	  }
        }

	for(j = 0; j != GreedyVec.size(); j++){
		mapit = Id_freq.find(GreedyVec[j]);
		if(mapit == Id_freq.end() ){
			Id_freq[GreedyVec[j]] =1;
		}else{
			Id_freq[GreedyVec[j]]++;
		}
	}
	
	if(num_sweeps == 1){
        outfile << " Data read in from file " << datafile << endl;
        outfile << " Total number of elements in dataset = "<< datasize << endl;
        if (!seed_specified){
			outfile << " Program generated seed for Mersenne Twister is " << longseed << endl;
		}else{
			outfile << " Using user supplied seed  "<< longseed << " for Mersenne Twister" << endl;
		}
		outfile << " "<< num_sweeps << " sets of elements with pairwise association <= " <<\
		lowLDlimit << " requested" << endl ;
		outfile << " Number of singletons is " << singletons.size() << endl;
		outfile << " Output from the deterministic heuristic " << endl << endl;
		outfile << " size = " << setw(8) << GreedyVec.size() << "   ";
		outfile << endl;
		OurData.cleanup();
		
		t_end = clock();
		cout << " Program terminating:";
		t_end = clock();
		cout << " Number of Processor Clock Ticks is " << double(t_end -t_begin) << endl;
		cout << " Number of seconds is " << double(t_end -t_begin)/CLOCKS_PER_SEC << endl;
		outfile << " Number of Processor Clock Ticks is " << double(t_end -t_begin) << endl;
		outfile << " Number of seconds is " << double(t_end -t_begin)/CLOCKS_PER_SEC << endl;		
		return 0;
	}
	
	vecSNPIds.push_back(GreedyVec);
	vector<unsigned> this_set;
	

	vector<unsigned> currentVec;

	unsigned largest_setsize;
	unsigned distinct_sweeps = vecSNPIds.size() -1; 
	// How many distinct return sets found from randomized algorithm 
	
	//   greedy set may be largest or not 
	//   if largest, then last element is largest, if not then second last is
	
	if (vecSNPIds[distinct_sweeps].size() > vecSNPIds[distinct_sweeps -1].size()){
		largest_setsize = vecSNPIds[distinct_sweeps].size();
	}else{
		largest_setsize = vecSNPIds[distinct_sweeps-1].size();
	}
	
	unsigned num_unique =0;
		
	if(greedy_found_already){
		num_unique = distinct_sweeps; // greedy set found during random sweeps
	}else{
		num_unique = distinct_sweeps +1; // greedy set not discovered earlier 
	}
	
	unsigned multip = 0;
	
	for(mapit = size_dist.begin();mapit != size_dist.end();mapit++){
		if(mapit->first == largest_setsize) multip = mapit->second;
	}
	
	vector<vector<string> >transfer;
	vector<string> tempvec;
	transfer.clear();
	transfer.resize(vecSNPIds.size() ); 
	vector<vector<unsigned> >::reverse_iterator revit;
	
	// reverse the order so that the greedy set appears first
	unsigned k = 0;
    for (revit = vecSNPIds.rbegin();revit != vecSNPIds.rend(); revit++) {
		tempvec.clear();
		this_set = *revit;
		tempvec.resize(this_set.size() );
		j = 0; 
        for (vector<unsigned>::iterator it = this_set.begin(); it != this_set.end(); it++) { 
            tempvec[j] = OurData.SNP_vector[(*it)-1]->SNPname;
			j++;
		} 
		transfer[k] = tempvec;
		k++;
	}
	

	outfile << " From " << num_sweeps << " attempts ";
	outfile << num_unique << " distinct sets found " << endl;
	if(greedy_found_already) outfile << " Greedy Set also found during random sweeps " << endl;
    outfile << " The sets of elements are arranged row by row with the deterministic output first " << endl
	<< " subsequent rows are ordered in decreasing size" << endl << endl;
    outfile << " "\
	<< "Number of elements from deterministic heuristic is " << vecSNPIds[distinct_sweeps].size() << endl;
    outfile << " "\
	<< "Largest Number of elements from random heuristic is " << vecSNPIds[distinct_sweeps -1].size() << endl;
    outfile << " "<< "Smallest Number of elements from random heuristic is " << vecSNPIds[0].size() << endl;
	outfile << " "<< "Largest Maximal independent set has size " << largest_setsize;
	outfile << " and multiplicity " << multip << endl;
	outfile << " "<< "Size distribution of independent sets is: " << endl;
	for(mapit = size_dist.begin();mapit != size_dist.end();mapit++){
		outfile << " Map Size " << mapit->first << " occurs " << mapit->second << " times " << endl;
	}
	std::pair<string,unsigned> temp_pair;
	unsigned sum1 = 0;
	
	
	for(mapit = Id_freq.begin();mapit != Id_freq.end();mapit++){
		temp_pair.first = OurData.SNP_vector[mapit->first -1]->SNPname ;
		temp_pair.second = mapit->second;
		sum1 += mapit->second;
	}
		
	outfile << endl;
	outfile << " Total number of distinct elements covered is " << Id_freq.size() << endl;
	outfile << " Total number of elements covered is " << sum1 << endl;
	outfile << " Printing Independent Sets first set is the Greedy Set" << endl;
	
	
	string offsetstring = "";
	unsigned sum2 = 0;
	
	for(i = 0; i != largest_name_size; i++){
		string shift_string = " ";
	    offsetstring = offsetstring + shift_string;
	}
	
	unsigned fieldwidth = 0;
	
	for(i = 0; i != transfer.size(); i++){
		tempvec = transfer[i];
		fieldwidth = unsigned(ceil(1.0 + log10((double(tempvec.size())))) );
	        outfile << endl << " Set Size = " << setw(fieldwidth) << tempvec.size();
		outfile << "            ";
		for(j = 0; j != tempvec.size(); j++){
	                 outfile << setw(largest_name_size) << tempvec[j] << "      ";
			sum2++;
		}
		outfile << endl;
	}
	
	
	if(sum1 != sum2){
		cout << " Discrepancy in elements found and printed : ";
		cout << sum1 << " found by searching data " << sum2 << " printed out " << endl;
		cout << " Program terminating " << endl;
		outfile << " Discrepancy in elemnts found and printed : ";
		outfile << sum1 << " found by searching data " << sum2 << " printed out " << endl;
		outfile << " Program terminating " << endl;
		throw exception();
	}
	

	
    for (i = 0; i != vecSNPIds.size(); i++) {
		this_set = vecSNPIds[i];
		this_set.clear();
    }
	
    vecSNPIds.clear();
    OurData.cleanup();
	
	t_end = clock();
	cout << " Number of Processor Clock Ticks is " << double(t_end -t_begin) << endl;
	cout << " Number of seconds is " << double(t_end -t_begin)/CLOCKS_PER_SEC << endl;
	outfile << endl << " Number of Processor Clock Ticks is " << double(t_end -t_begin) << endl;
	outfile << " Number of seconds is " << double(t_end -t_begin)/CLOCKS_PER_SEC << endl;
	
	
    outfile << "-----------------------------------------------------------------" << endl;	
	
    cout << " Program terminating " << endl;
	
	return 0;
} 

	

