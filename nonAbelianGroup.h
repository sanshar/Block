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
#include "pario.h"

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
	    //pout << i<<"  "<<j<<"  "<<k<<"  "<<CGcoeffs(i,j,k)<<endl;
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
	    //pout << i<<"  "<<j<<"  "<<k<<"  "<<CGcoeffs(i,j,k)<<endl;
	  }


    }
  };


  class D5h: public nonAbelianGroup {
  public:
    D5h() 
    {
      numIrreps = 8;
      irrepNames.push_back("A1'"); irrepNames.push_back("A2'"); irrepNames.push_back("E1'"), irrepNames.push_back("E2'");
      irrepNames.push_back("A1\""); irrepNames.push_back("A2\""); irrepNames.push_back("E1\""), irrepNames.push_back("E2\"");
      irrepSizes.push_back(1); irrepSizes.push_back(1); irrepSizes.push_back(2); irrepSizes.push_back(2);
      irrepSizes.push_back(1); irrepSizes.push_back(1); irrepSizes.push_back(2); irrepSizes.push_back(2);

      int multarray1[] = {0,1,2};//{2,1,0};
      int multarray2[] = {0,1,3};//{3,1,0};
      int multarray3[] = {2,3};
      int multarray1b[] = {4,5,6};//{2,1,0};
      int multarray2b[] = {4,5,7};//{3,1,0};
      int multarray3b[] = {6,7};

      //read_multiplication Table
      multiplicationTable.ReSize(8,8);
      multiplicationTable(0,0) = vector<int>(1,0);
      multiplicationTable(1,0) = vector<int>(1,1); 
      multiplicationTable(2,0) = vector<int>(1,2); 
      multiplicationTable(3,0) = vector<int>(1,3); 
      multiplicationTable(4,0) = vector<int>(1,4); 
      multiplicationTable(5,0) = vector<int>(1,5); 
      multiplicationTable(6,0) = vector<int>(1,6); 
      multiplicationTable(7,0) = vector<int>(1,7); 

      multiplicationTable(0,1) = vector<int>(1,1);
      multiplicationTable(1,1) = vector<int>(1,0); 
      multiplicationTable(2,1) = vector<int>(1,2); 
      multiplicationTable(3,1) = vector<int>(1,3); 
      multiplicationTable(4,1) = vector<int>(1,5); 
      multiplicationTable(5,1) = vector<int>(1,4); 
      multiplicationTable(6,1) = vector<int>(1,6); 
      multiplicationTable(7,1) = vector<int>(1,7); 

      multiplicationTable(0,2) = vector<int>(1,2);
      multiplicationTable(1,2) = vector<int>(1,2); 
      multiplicationTable(2,2) = vector<int>(multarray2, multarray2 + sizeof(multarray2)/sizeof(int));
      multiplicationTable(3,2) = vector<int>(multarray3, multarray3 + sizeof(multarray3)/sizeof(int));
      multiplicationTable(4,2) = vector<int>(1,6);
      multiplicationTable(5,2) = vector<int>(1,6); 
      multiplicationTable(6,2) = vector<int>(multarray2b, multarray2b + sizeof(multarray2b)/sizeof(int));
      multiplicationTable(7,2) = vector<int>(multarray3b, multarray3b + sizeof(multarray3b)/sizeof(int));

      multiplicationTable(0,3) = vector<int>(1,3);
      multiplicationTable(1,3) = vector<int>(1,3); 
      multiplicationTable(2,3) = vector<int>(multarray3, multarray3 + sizeof(multarray3)/sizeof(int));
      multiplicationTable(3,3) = vector<int>(multarray1, multarray1 + sizeof(multarray1)/sizeof(int));
      multiplicationTable(4,3) = vector<int>(1,7);
      multiplicationTable(5,3) = vector<int>(1,7); 
      multiplicationTable(6,3) = vector<int>(multarray3b, multarray3b + sizeof(multarray3b)/sizeof(int));
      multiplicationTable(7,3) = vector<int>(multarray1b, multarray1b + sizeof(multarray1b)/sizeof(int));

      multiplicationTable(0,4) = vector<int>(1,4);
      multiplicationTable(1,4) = vector<int>(1,5); 
      multiplicationTable(2,4) = vector<int>(1,6); 
      multiplicationTable(3,4) = vector<int>(1,7); 
      multiplicationTable(4,4) = vector<int>(1,0); 
      multiplicationTable(5,4) = vector<int>(1,1); 
      multiplicationTable(6,4) = vector<int>(1,2); 
      multiplicationTable(7,4) = vector<int>(1,3); 

      multiplicationTable(0,5) = vector<int>(1,5);
      multiplicationTable(1,5) = vector<int>(1,4); 
      multiplicationTable(2,5) = vector<int>(1,6); 
      multiplicationTable(3,5) = vector<int>(1,7); 
      multiplicationTable(4,5) = vector<int>(1,1); 
      multiplicationTable(5,5) = vector<int>(1,0); 
      multiplicationTable(6,5) = vector<int>(1,2); 
      multiplicationTable(7,5) = vector<int>(1,3); 

      multiplicationTable(0,6) = vector<int>(1,6);
      multiplicationTable(1,6) = vector<int>(1,6); 
      multiplicationTable(2,6) = vector<int>(multarray2b, multarray2b + sizeof(multarray2b)/sizeof(int));
      multiplicationTable(3,6) = vector<int>(multarray3b, multarray3b + sizeof(multarray3b)/sizeof(int));
      multiplicationTable(4,6) = vector<int>(1,2);
      multiplicationTable(5,6) = vector<int>(1,2); 
      multiplicationTable(6,6) = vector<int>(multarray2, multarray2 + sizeof(multarray2)/sizeof(int));
      multiplicationTable(7,6) = vector<int>(multarray3, multarray3 + sizeof(multarray3)/sizeof(int));

      multiplicationTable(0,7) = vector<int>(1,7);
      multiplicationTable(1,7) = vector<int>(1,7); 
      multiplicationTable(2,7) = vector<int>(multarray3b, multarray3b + sizeof(multarray3b)/sizeof(int));
      multiplicationTable(3,7) = vector<int>(multarray1b, multarray1b + sizeof(multarray1b)/sizeof(int));
      multiplicationTable(4,7) = vector<int>(1,3);
      multiplicationTable(5,7) = vector<int>(1,3); 
      multiplicationTable(6,7) = vector<int>(multarray3, multarray3 + sizeof(multarray3)/sizeof(int));
      multiplicationTable(7,7) = vector<int>(multarray1, multarray1 + sizeof(multarray1)/sizeof(int));

      readCGfromDisk();

    }

    string getGroupName() {return "D5h";}

    void readCGfromDisk() {
      int largestIrrep = numIrreps-1;
      int largestRow = irrepSizes[largestIrrep]-1;
      int loopOver = getCombinedIndex(largestIrrep, largestRow) + 1;
      
      CGcoeffs.ReSize(loopOver, loopOver, loopOver);
      std::ifstream instream("d5h.cgcoeffs", std::ifstream::in);

      int sizeofCG = 0;
      for (int i=0; i<loopOver; i++)
	for (int j=0; j<loopOver; j++)
	  for (int k=0; k<loopOver; k++) {
	    instream >> CGcoeffs(i,j,k);
	    //pout << i<<"  "<<j<<"  "<<k<<"  "<<CGcoeffs(i,j,k)<<endl;
	  }


    }
  };


  class D4h: public nonAbelianGroup {
  public:
    D4h() 
    {
      enum{A1g, A2g, B1g, B2g, Eg, A1u, A2u, B1u, B2u, Eu};
      numIrreps = 10;
      irrepNames.push_back("A1g"); irrepNames.push_back("A2g"); irrepNames.push_back("B1g"), irrepNames.push_back("B2g"), irrepNames.push_back("Eg");
      irrepNames.push_back("A1u"); irrepNames.push_back("A2u"); irrepNames.push_back("B1u"), irrepNames.push_back("B2u"), irrepNames.push_back("Eu");
      irrepSizes.push_back(1); irrepSizes.push_back(1); irrepSizes.push_back(1); irrepSizes.push_back(1); irrepSizes.push_back(2);
      irrepSizes.push_back(1); irrepSizes.push_back(1); irrepSizes.push_back(1); irrepSizes.push_back(1); irrepSizes.push_back(2);


      int multarray1[] = {0,1,2,3};//{2,1,0};
      int multarray2[] = {5,6,7,8};//{3,1,0};

      //read_multiplication Table
      multiplicationTable.ReSize(10,10);
      multiplicationTable(0,0) = vector<int>(1,0);
      multiplicationTable(1,0) = vector<int>(1,1); 
      multiplicationTable(2,0) = vector<int>(1,2); 
      multiplicationTable(3,0) = vector<int>(1,3); 
      multiplicationTable(4,0) = vector<int>(1,4); 
      multiplicationTable(5,0) = vector<int>(1,5); 
      multiplicationTable(6,0) = vector<int>(1,6); 
      multiplicationTable(7,0) = vector<int>(1,7); 
      multiplicationTable(8,0) = vector<int>(1,8); 
      multiplicationTable(9,0) = vector<int>(1,9); 

      multiplicationTable(0,1) = vector<int>(1,A2g);
      multiplicationTable(1,1) = vector<int>(1,A1g); 
      multiplicationTable(2,1) = vector<int>(1,B2g); 
      multiplicationTable(3,1) = vector<int>(1,B1g); 
      multiplicationTable(4,1) = vector<int>(1,Eg); 
      multiplicationTable(5,1) = vector<int>(1,A2u); 
      multiplicationTable(6,1) = vector<int>(1,A1u); 
      multiplicationTable(7,1) = vector<int>(1,B2u); 
      multiplicationTable(8,1) = vector<int>(1,B1u); 
      multiplicationTable(9,1) = vector<int>(1,Eu); 

      multiplicationTable(0,2) = vector<int>(1,B1g);
      multiplicationTable(1,2) = vector<int>(1,B2g); 
      multiplicationTable(2,2) = vector<int>(1,A1g); 
      multiplicationTable(3,2) = vector<int>(1,A2g); 
      multiplicationTable(4,2) = vector<int>(1,Eg); 
      multiplicationTable(5,2) = vector<int>(1,B1u); 
      multiplicationTable(6,2) = vector<int>(1,B2u); 
      multiplicationTable(7,2) = vector<int>(1,A1u); 
      multiplicationTable(8,2) = vector<int>(1,A2u); 
      multiplicationTable(9,2) = vector<int>(1,Eu); 

      multiplicationTable(0,3) = vector<int>(1,B2g);
      multiplicationTable(1,3) = vector<int>(1,B1g); 
      multiplicationTable(2,3) = vector<int>(1,A2g); 
      multiplicationTable(3,3) = vector<int>(1,A1g); 
      multiplicationTable(4,3) = vector<int>(1,Eg); 
      multiplicationTable(5,3) = vector<int>(1,B2u); 
      multiplicationTable(6,3) = vector<int>(1,B1u); 
      multiplicationTable(7,3) = vector<int>(1,A2u); 
      multiplicationTable(8,3) = vector<int>(1,A1u); 
      multiplicationTable(9,3) = vector<int>(1,Eu); 

      multiplicationTable(0,4) = vector<int>(1,Eg);
      multiplicationTable(1,4) = vector<int>(1,Eg); 
      multiplicationTable(2,4) = vector<int>(1,Eg); 
      multiplicationTable(3,4) = vector<int>(1,Eg); 
      multiplicationTable(4,4) = vector<int>(multarray1, multarray1+sizeof(multarray1)/sizeof(int)); 
      multiplicationTable(5,4) = vector<int>(1,Eu); 
      multiplicationTable(6,4) = vector<int>(1,Eu); 
      multiplicationTable(7,4) = vector<int>(1,Eu); 
      multiplicationTable(8,4) = vector<int>(1,Eu); 
      multiplicationTable(9,4) = vector<int>(multarray2, multarray2+sizeof(multarray2)/sizeof(int));  

      multiplicationTable(0,5) = vector<int>(1,A1u);
      multiplicationTable(1,5) = vector<int>(1,A2u); 
      multiplicationTable(2,5) = vector<int>(1,B1u); 
      multiplicationTable(3,5) = vector<int>(1,B2u); 
      multiplicationTable(4,5) = vector<int>(1,Eu); 
      multiplicationTable(5,5) = vector<int>(1,A1g); 
      multiplicationTable(6,5) = vector<int>(1,A2g); 
      multiplicationTable(7,5) = vector<int>(1,B1g); 
      multiplicationTable(8,5) = vector<int>(1,B2g); 
      multiplicationTable(9,5) = vector<int>(1,Eg); 

      multiplicationTable(0,6) = vector<int>(1,A2u);
      multiplicationTable(1,6) = vector<int>(1,A1u); 
      multiplicationTable(2,6) = vector<int>(1,B2u); 
      multiplicationTable(3,6) = vector<int>(1,B1u); 
      multiplicationTable(4,6) = vector<int>(1,Eu); 
      multiplicationTable(5,6) = vector<int>(1,A2g); 
      multiplicationTable(6,6) = vector<int>(1,A1g); 
      multiplicationTable(7,6) = vector<int>(1,B2g); 
      multiplicationTable(8,6) = vector<int>(1,B1g); 
      multiplicationTable(9,6) = vector<int>(1,Eg); 

      multiplicationTable(0,7) = vector<int>(1,B1u);
      multiplicationTable(1,7) = vector<int>(1,B2u); 
      multiplicationTable(2,7) = vector<int>(1,A1u); 
      multiplicationTable(3,7) = vector<int>(1,A2u); 
      multiplicationTable(4,7) = vector<int>(1,Eu); 
      multiplicationTable(5,7) = vector<int>(1,B1g); 
      multiplicationTable(6,7) = vector<int>(1,B2g); 
      multiplicationTable(7,7) = vector<int>(1,A1g); 
      multiplicationTable(8,7) = vector<int>(1,A2g); 
      multiplicationTable(9,7) = vector<int>(1,Eg); 

      multiplicationTable(0,8) = vector<int>(1,B2u);
      multiplicationTable(1,8) = vector<int>(1,B1u); 
      multiplicationTable(2,8) = vector<int>(1,A2u); 
      multiplicationTable(3,8) = vector<int>(1,A1u); 
      multiplicationTable(4,8) = vector<int>(1,Eu); 
      multiplicationTable(5,8) = vector<int>(1,B2g); 
      multiplicationTable(6,8) = vector<int>(1,B1g); 
      multiplicationTable(7,8) = vector<int>(1,A2g); 
      multiplicationTable(8,8) = vector<int>(1,A1g); 
      multiplicationTable(9,8) = vector<int>(1,Eg); 

      multiplicationTable(0,9) = vector<int>(1,Eu);
      multiplicationTable(1,9) = vector<int>(1,Eu); 
      multiplicationTable(2,9) = vector<int>(1,Eu); 
      multiplicationTable(3,9) = vector<int>(1,Eu); 
      multiplicationTable(4,9) = vector<int>(multarray2, multarray2+sizeof(multarray1)/sizeof(int)); 
      multiplicationTable(5,9) = vector<int>(1,Eg); 
      multiplicationTable(6,9) = vector<int>(1,Eg); 
      multiplicationTable(7,9) = vector<int>(1,Eg); 
      multiplicationTable(8,9) = vector<int>(1,Eg); 
      multiplicationTable(9,9) = vector<int>(multarray1, multarray1+sizeof(multarray2)/sizeof(int));  


      readCGfromDisk();

    }

    string getGroupName() {return "D4h";}

    void readCGfromDisk() {
      int largestIrrep = numIrreps-1;
      int largestRow = irrepSizes[largestIrrep]-1;
      int loopOver = getCombinedIndex(largestIrrep, largestRow) + 1;
      
      CGcoeffs.ReSize(loopOver, loopOver, loopOver);
      std::ifstream instream("d4h.cgcoeffs", std::ifstream::in);

      int sizeofCG = 0;
      for (int i=0; i<loopOver; i++)
	for (int j=0; j<loopOver; j++)
	  for (int k=0; k<loopOver; k++) {
	    instream >> CGcoeffs(i,j,k);
	    //pout << i<<"  "<<j<<"  "<<k<<"  "<<CGcoeffs(i,j,k)<<endl;
	  }


    }
  };

}

#endif
