/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef NONABELIANGROUP_HEADER
#define NONABELIANGROUP_HEADER

#include <vector>
#include <string>
#include "ObjectMatrix.h"
#include <iostream>
#include <fstream>

//using namespace std;

namespace SpinAdapted {
  class nonAbelianGroup {
  protected:
    int numIrreps;
    vector<string> irrepNames;
    vector<int> irrepSizes;
    ObjectMatrix< vector<int> > multiplicationTable;
    ObjectMatrix3D< double> CGcoeffs;
    
  public :
    nonAbelianGroup() {};
    virtual ~nonAbelianGroup() {};
    string getIrrepName(int i) const {return irrepNames[i];}
    int getIrrepSize(int i) const {return irrepSizes[i];}
    int getNumIrreps() const {return numIrreps;} 
    vector<int> getProduct(int i, int j) const {return multiplicationTable(i,j);}
    int getCombinedIndex(int A, int rowa) const
    {
      int index = 0;
      for (int i=0; i<A; i++) index += irrepSizes[i];
      return index+rowa;
    }
    double getCG(int A, int B, int C, int i, int j, int k) const {
      int Ai = getCombinedIndex(A,i), Bj = getCombinedIndex(B,j), Ck = getCombinedIndex(C,k);
      return CGcoeffs(Ai, Bj, Ck);
    }
    virtual string getGroupName() const {return "nonabelian";}
  };

  class C3v : public nonAbelianGroup {
  public:
    C3v() 
    {
      numIrreps = 3;
      irrepNames.push_back("A1"); irrepNames.push_back("A2"); irrepNames.push_back("E");
      irrepSizes.push_back(1); irrepSizes.push_back(1); irrepSizes.push_back(2);

      //read_multiplication Table
      multiplicationTable.ReSize(3,3);
      vector<int> multentry(1,0);
      multiplicationTable(0,0) = multentry;
      multiplicationTable(1,0) = vector<int>(1,1); multiplicationTable(0,1) = multiplicationTable(1,0);
      multiplicationTable(2,0) = vector<int>(1,2); multiplicationTable(0,2) = multiplicationTable(2,0);

      multiplicationTable(1,1) = vector<int>(1,0);
      multiplicationTable(2,1) = vector<int>(1,2); multiplicationTable(1,2) = multiplicationTable(2,1);

      int multarray[] = {2,1,0};
      multiplicationTable(2,2) = vector<int>(multarray, multarray + sizeof(multarray)/sizeof(int));

      readCGfromDisk();

    }

    string getGroupName() {return "C3v";}

    void readCGfromDisk() {
      int largestIrrep = numIrreps-1;
      int largestRow = irrepSizes[largestIrrep]-1;
      int loopOver = getCombinedIndex(largestIrrep, largestRow) + 1;
      
      CGcoeffs.ReSize(loopOver, loopOver, loopOver);
      std::ifstream instream("c3v.cgcoeffs", std::ifstream::in);

      int sizeofCG = 0;
      for (int i=0; i<loopOver; i++)
	for (int j=0; j<loopOver; j++)
	  for (int k=0; k<loopOver; k++) {
	    instream >> CGcoeffs(i,j,k);
	    //cout << i<<"  "<<j<<"  "<<k<<"  "<<CGcoeffs(i,j,k)<<endl;
	  }


    }
  };


  class C5v : public nonAbelianGroup {
  public:
    C5v() 
    {
      numIrreps = 4;
      irrepNames.push_back("A1"); irrepNames.push_back("A2"); irrepNames.push_back("E1"), irrepNames.push_back("E2");
      irrepSizes.push_back(1); irrepSizes.push_back(1); irrepSizes.push_back(2); irrepSizes.push_back(2);

      int multarray1[] = {0,1,2};//{2,1,0};
      int multarray2[] = {0,1,3};//{3,1,0};
      int multarray3[] = {2,3};

      //read_multiplication Table
      multiplicationTable.ReSize(4,4);
      multiplicationTable(0,0) = vector<int>(1,0);
      multiplicationTable(1,0) = vector<int>(1,1); 
      multiplicationTable(2,0) = vector<int>(1,2); 
      multiplicationTable(3,0) = vector<int>(1,3); 

      multiplicationTable(0,1) = vector<int>(1,1);
      multiplicationTable(1,1) = vector<int>(1,0); 
      multiplicationTable(2,1) = vector<int>(1,2); 
      multiplicationTable(3,1) = vector<int>(1,3); 

      multiplicationTable(0,2) = vector<int>(1,2);
      multiplicationTable(1,2) = vector<int>(1,2); 
      multiplicationTable(2,2) = vector<int>(multarray2, multarray2 + sizeof(multarray2)/sizeof(int));
      multiplicationTable(3,2) = vector<int>(multarray3, multarray3 + sizeof(multarray3)/sizeof(int));

      multiplicationTable(0,3) = vector<int>(1,3);
      multiplicationTable(1,3) = vector<int>(1,3); 
      multiplicationTable(2,3) = vector<int>(multarray3, multarray3 + sizeof(multarray3)/sizeof(int));
      multiplicationTable(3,3) = vector<int>(multarray1, multarray1 + sizeof(multarray1)/sizeof(int));

      readCGfromDisk();

    }

    string getGroupName() {return "C5v";}

    void readCGfromDisk() {
      int largestIrrep = numIrreps-1;
      int largestRow = irrepSizes[largestIrrep]-1;
      int loopOver = getCombinedIndex(largestIrrep, largestRow) + 1;
      
      CGcoeffs.ReSize(loopOver, loopOver, loopOver);
      std::ifstream instream("c5v.cgcoeffs", std::ifstream::in);

      int sizeofCG = 0;
      for (int i=0; i<loopOver; i++)
	for (int j=0; j<loopOver; j++)
	  for (int k=0; k<loopOver; k++) {
	    instream >> CGcoeffs(i,j,k);
	    //cout << i<<"  "<<j<<"  "<<k<<"  "<<CGcoeffs(i,j,k)<<endl;
	  }


    }
  };

}

#endif
