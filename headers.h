// This is the header file for FastIndep, a C++ program for 
// finding maximal independent sets given an input matrix and a threshold. 
// This program may be freely distributed or modified, provided these 
// comment lines are retained. In particular, there are no restrictions
// on non-commercial use. If this code (either in the original form or 
// modified) is used in academic publications, a citation would be appreciated.  
// Please read the accompanying ReadMe.txt file for instructions on how to compile 
// and run the program.  
// Author: Joseph Abraham May 2013


#include <map>
#include <set>
#include <vector>
#include <algorithm> 
#include <string>
#include "./MersenneTwister.h"
using namespace std;

class SNP{
public:
      SNP(){};
      string SNPname; // name from user supplied file
      unsigned SNPId; // Label for each node
      vector<unsigned> seps; // Ids of SNPs in low LD with SNP
    };

class Data{
public:
      Data(){};
      vector<SNP*> SNP_vector;
	  vector<unsigned> allIds;
	  vector<unsigned> singletons;  // singletons 
	  vector<unsigned> all_connectors; // connected to all 
	  bool data_rand;
      map<string, unsigned> NameIdMap; // could be useful later
	  vector<unsigned> coveredVec;
      MTRand ran_obj; // needed to run Mersenne twister
      unsigned long longseed;
      void AugmentcoveredVec(const unsigned nextSNP);
      void FindSeparatedVec(vector<unsigned>& ReturnVec); // Finds an independent set 
      unsigned FindnextSNP(const vector<unsigned>& candSet); 
      // Used for the finding unassociated nodes
      unsigned FindnextSNP_rand(const vector<unsigned>& candSet); // see above with randomization
      unsigned FindstartSNP_rand();
      unsigned FindstartSNP(); // greedy heuristic for finding starting SNPent
      void cleanup(); // invoked at end of main program to delete pointers etc.
      void initialize_rand_gen();
    };

