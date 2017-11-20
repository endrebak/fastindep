// This file contains graph theoretical subroutines for FastIndep, a C++ program 
// for finding maximal independent sets given an input matrix and a threshold. 
// Methods with "rand" in the title are randomized greedy algorithms
// This program may be freely distributed or modified, provided these 
// comment lines are retained. In particular, there are no restrictions
// on non-commercial use. If this code (either in the original form or 
// modified) is used in academic publications, a citation would be appreciated.  
// Please read the accompanying ReadMe.txt file for instructions on how to compile 
// and run the program.  
// Author: Joseph Abraham May 2013

#include <iostream>
#include <fstream>
#include <ostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <limits>
#include "./headers.h"
	
	using namespace std;
	
	void Data::FindSeparatedVec(vector<unsigned>& ReturnVec)
{
		vector<unsigned> interVec;
		vector<unsigned>::iterator endit,vecit;
		vector<unsigned> candVec,tempVec,next_vec;
		unsigned i;
		unsigned nextSNP = 0;
		unsigned startSNP = 0;
		
		// initialization
		candVec.clear();
		ReturnVec.clear();
		tempVec.clear();
		coveredVec.clear();
		interVec.clear();
		
		// Finding a good starting SNP
		
		if(data_rand){
			startSNP = FindstartSNP_rand();
		} else{
			startSNP = FindstartSNP();
		}
		
		if (startSNP == 0){
			std::cout << startSNP << " Bad Return " << endl;
			std::cout.flush();
			throw exception();
		}
		
		tempVec.push_back(startSNP);
		
		candVec = SNP_vector[startSNP -1]->seps;
		// experimental code, remove singletons at start 	
		
		if(singletons.size() > 0){
			interVec.resize(candVec.size() );
			endit = set_difference(candVec.begin(), candVec.end(), singletons.begin(), \
								   singletons.end(), interVec.begin()  );
			candVec.resize(unsigned(endit -interVec.begin()));
			copy(interVec.begin(), endit, candVec.begin() );
			sort(candVec.begin(), candVec.end() );
		}
		
		AugmentcoveredVec(startSNP);
		
		if(singletons.size() > 0){
			copy(singletons.begin(), singletons.end(), back_inserter(coveredVec) );
			sort(coveredVec.begin(), coveredVec.end() );
			copy(singletons.begin(), singletons.end(), back_inserter(tempVec) );
		}
		
		
		if(candVec.size() == 0){
			if(singletons.size() == 0){
				std::cout << "Should not happen ... Check entries in parameters.txt !" << endl;
				std::cout << "Started from " << startSNP << endl;
				std::cout.flush();
				throw exception();
			}else{
				ReturnVec.resize(tempVec.size() );
				copy(tempVec.begin(), tempVec.end(),ReturnVec.begin() );
				return; 
			}
		}
		
		//     Main cycle of program
		do {
			
			if(data_rand){
				nextSNP = FindnextSNP_rand(candVec);
			}else{
				nextSNP = FindnextSNP(candVec);
			}
			
			if(nextSNP == 0){
				ReturnVec.resize(tempVec.size() );
				copy(tempVec.begin(),tempVec.end(),ReturnVec.begin() );
				return;
			}
			
			
			next_vec = SNP_vector[nextSNP -1]->seps; 
			
			sort(next_vec.begin(), next_vec.end() );
			sort(tempVec.begin(), tempVec.end() );
			interVec.resize(tempVec.size() );
			
			endit =  set_intersection(next_vec.begin(), next_vec.end(), 
									  tempVec.begin(), tempVec.end(), interVec.begin() );
			
			if(unsigned(endit -interVec.begin() )  != tempVec.size() ){		
				for(i = 0;i != tempVec.size();i++){
					vecit = find(next_vec.begin(), next_vec.end(), tempVec[i]);
					if(vecit == next_vec.end() ) {
						std::cout << " Inconsistent SNP Added " << endl;
						std::cout << SNP_vector[nextSNP-1]->SNPname \
						<< " and " << SNP_vector[tempVec[i]-1]->SNPname << " correlated " << endl;
						std::cout.flush();
						std::cout << " Program Terminating " << endl;
						throw exception();		
					}
				}
			}
			
			AugmentcoveredVec(nextSNP);
			tempVec.push_back(nextSNP);	
			
			candVec.erase(find(candVec.begin(),candVec.end(), nextSNP));
			
			sort(candVec.begin(), candVec.end() );
			sort(next_vec.begin(), next_vec.end() );
			
			interVec.resize(candVec.size() );
			
			endit = set_intersection(candVec.begin(),candVec.end(), next_vec.begin(), next_vec.end(),
									 interVec.begin() );
			
			
			if(unsigned(endit -interVec.begin())== 0){
				ReturnVec.resize(tempVec.size() );
				copy(tempVec.begin(), tempVec.end(),ReturnVec.begin() );
				return;
			}
			
			candVec.resize(unsigned(endit -interVec.begin()));
			copy(interVec.begin(), endit, candVec.begin() );
			
		}while(tempVec.size() < SNP_vector.size()); // dummy condition
		
		return;
} 
	
	unsigned Data::FindstartSNP()
{	
		// Greedy Heuristic to find starting SNP
		
		unsigned nextSNP = 0;
		unsigned i;
		vector<unsigned>::iterator vecit;
		bool is_connector,is_singleton;
		double num_choices,largest_so_far;
		largest_so_far = 0;
		num_choices = 0;
		
		for(i=0;i != SNP_vector.size();i++){
			is_connector = false;
			if(all_connectors.size() > 0) {
				vecit = find(all_connectors.begin(), all_connectors.end(), (i+1) );
				if(vecit != all_connectors.end() ) is_connector = true;
			}
			if(is_connector) continue;// do not include connectors among possible starts
			is_singleton = false;
			if(singletons.size() > 0){
				vecit = find(singletons.begin(), singletons.end(), (i+1) );
				if(vecit != singletons.end() ) is_singleton = true;
			}
			if(is_singleton) continue; // do not include singletons among possible starts 
			num_choices = SNP_vector[i]->seps.size(); 
			if(num_choices > largest_so_far){
				largest_so_far = num_choices;
				nextSNP = (i+1);
			} 
		}
		return nextSNP;
}
	
	
	unsigned Data::FindnextSNP(const vector<unsigned>& CandVec)
{	
		// Greedy Heuristic to find next SNP
		
		unsigned nextSNP = 0;
		unsigned i;
		double num_choices,largest_so_far;
		largest_so_far = 0;
		num_choices = 0;
		
		for(i=0;i != CandVec.size();i++){
			num_choices = SNP_vector[CandVec[i] -1]->seps.size();
			if(num_choices > largest_so_far){
				largest_so_far = num_choices;
				nextSNP = CandVec[i];
			} 
		}
		
		if(nextSNP == 0) {
			std::cout << " Possible Bug " << endl;
			std::cout.flush();
			throw exception();
		}
		
		return nextSNP;
}
	
	
	unsigned Data::FindnextSNP_rand(const vector<unsigned>& candVec)
{
		unsigned nextSNP = 0;
		unsigned i;
		double num_choices = 0;
		double total_num = 0;
	    double incr;
		vector<std::pair<unsigned,double> > neighvec;
		neighvec.resize(candVec.size() );
		
		for(i=0; i != candVec.size(); i++){
			num_choices = SNP_vector[candVec[i] -1]->seps.size(); 
			incr = (double)num_choices;
			std::pair<unsigned,double> thispair(candVec[i],incr);
			neighvec[i] =thispair; 
			total_num += incr;
		}
		
		
		if (total_num == 0) {
			{ 
				std::cout << " Possible Bug " << endl; 
				std::cout.flush();
				throw exception();
			}
		}
		
		for(i =0 ; i != neighvec.size(); i++){
			neighvec[i].second = (neighvec[i].second/total_num);
		}
		
		double previous = 0;
		double u = ran_obj.rand();
		
		for(i = 0; i != neighvec.size(); i++){
			if(u < (neighvec[i].second + previous)) {nextSNP = neighvec[i].first;}
			else{previous += neighvec[i].second;}     
			if (nextSNP != 0) break;
		}
		
		return nextSNP;
}
	
	unsigned Data::FindstartSNP_rand()
{
		unsigned nextSNP = 0;
		unsigned i;
		double num_choices = 0;
		double total_num = 0;
		vector<std::pair<unsigned,double> > neighvec;
		vector<unsigned>::iterator vecit;
		bool is_connector;
		bool is_singleton;
		
		for(i=0; i != SNP_vector.size(); i++){
			double incr = numeric_limits<double>::quiet_NaN();
			is_connector = false;
			if(all_connectors.size() > 0){
				vecit = find(all_connectors.begin(), all_connectors.end(), (i+1) );
				if(vecit != all_connectors.end() ) is_connector = true;
			}
			if(is_connector) continue; // do not include connectors among possible starts
			is_singleton = false;
			if(singletons.size() > 0){
				vecit = find(singletons.begin(), singletons.end(), (i+1) );
				if(vecit != singletons.end() ) is_singleton = true;
			}
			if(is_singleton) continue; // do not include singletons among possible starts 
			num_choices = SNP_vector[i]->seps.size();
			if(num_choices <= 0){
				std::cout << " Entries for element " << (i+1) 
				<<	" dubious ... program terminating " << endl;
				std::cout.flush();
				throw exception();
			}
			incr = pow((num_choices*1.0),1.0);
			std::pair<unsigned,double> thispair((i+1),incr);
			neighvec.push_back(thispair); 
			total_num += incr;
		}
		
		
		if(total_num == 0) {
			std::cout << " Total sum of weights vanishes ... error ... program terminating " << endl;
			std::cout.flush();
			throw exception();
		}
		
		for(i =0 ; i != neighvec.size(); i++){
			neighvec[i].second = (neighvec[i].second/total_num);
		}
		
		double previous = 0;
		double u = ran_obj.rand();
		
		for(i = 0; i != neighvec.size(); i++){
			if(u < (neighvec[i].second + previous)) {nextSNP = neighvec[i].first;}
			else{previous += neighvec[i].second;}     
			if (nextSNP != 0) break;
		}
		
		if(nextSNP == 0){
			std::cout << " Possible Problems with biased coin toss ... program terminating " << endl;
			std::cout.flush();
			throw exception();
		}
		
		return nextSNP;
}
	
	// Used to check for dominating property 
	void Data::AugmentcoveredVec(const unsigned nextSNP){
		
		vector<unsigned> diffVec;
		diffVec.resize(allIds.size(),0);
		//		set<unsigned> newSet;
		//		newSet.clear();
		vector<unsigned> newVec;
		newVec.clear();
		vector<unsigned>::iterator endit,endit1;
		
		vector<unsigned> next_vec = SNP_vector[nextSNP -1]->seps;
		
		endit = set_difference(allIds.begin(),allIds.end(),\
							   next_vec.begin(),next_vec.end(), diffVec.begin() );
		
		if(unsigned(endit -diffVec.begin() ) != unsigned(allIds.size() -next_vec.size()) ){
			std::cout << " Size mismatch with cover of " << nextSNP 
			<< "  " << next_vec.size() << "  " << unsigned(endit -diffVec.begin() ) << endl;
			std::cout.flush();
			throw exception();
		}
		
		
		if(coveredVec.size() == 0){			
			coveredVec.resize(unsigned(endit -diffVec.begin()) );
			copy(diffVec.begin(),endit, coveredVec.begin() );
		}else{	
			newVec.resize(coveredVec.size() + unsigned(endit -diffVec.begin() ) );
			endit1 = set_union(coveredVec.begin(), coveredVec.end(), diffVec.begin(), \
							   endit,newVec.begin() );
			coveredVec.resize(unsigned(endit1 -newVec.begin()) );
			copy(newVec.begin(),endit1,coveredVec.begin() );
		}
		sort(coveredVec.begin(), coveredVec.end());
	    return;
		
	}
	
	void Data::cleanup()
{
		allIds.clear();
		singletons.clear();
		all_connectors.clear();
		coveredVec.clear();
		
		for(unsigned i = 0; i != SNP_vector.size();i++){
			SNP_vector[i]->seps.clear();
			SNP_vector[i]->SNPname.clear();
			delete SNP_vector[i];
		}
		
		SNP_vector.clear();
		NameIdMap.clear(); 
		
		return;
}
	
	void Data::initialize_rand_gen(){
		ran_obj = MTRand();
		ran_obj.seed(longseed);
}



