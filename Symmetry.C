#include "Symmetry.h"

// some explicit instantiations to prevent IBM linker from complaining
#include <string>
#include <newmat.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include "global.h"
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
template class multiarray<int>;
template class array_2d<int>;
array_2d<int> groupTable;

namespace SpinAdapted {

void Symmetry::InitialiseTable(string psym)
{
  if (sym == "c1")
    {
      groupTable.resize(1, 1);
      groupTable(0, 0) = 0;
    }
  else if (sym == "ci")
    {
      /*
	inversion symmetry table
	0 1 
	1 0
      */
      groupTable.resize(2, 2);
      groupTable(0, 0) = 0;
      groupTable(0, 1) = 1;
      groupTable(1, 0) = 1;
      groupTable(1, 1) = 0;
    }
  else if (sym == "c2v" || sym == "c2h")
    {
      /*
	c2v symmetry table
	0 1 2 3
	1 0 3 2
	2 3 0 1
	3 2 1 0
      */
      groupTable.resize(4, 4);
      groupTable(0, 0) = 0;
      groupTable(0, 1) = 1;
      groupTable(0, 2) = 2;
      groupTable(0, 3) = 3;

      groupTable(1, 0) = 1;
      groupTable(1, 1) = 0;
      groupTable(1, 2) = 3;
      groupTable(1, 3) = 2;

      groupTable(2, 0) = 2;
      groupTable(2, 1) = 3;
      groupTable(2, 2) = 0;
      groupTable(2, 3) = 1;

      groupTable(3, 0) = 3;
      groupTable(3, 1) = 2;
      groupTable(3, 2) = 1;
      groupTable(3, 3) = 0;

    }
  else if (sym == "d2h")
    {   
      /*
      data multt /1,2,3,4,5,6,7,8,
     1            2,1,4,3,6,5,8,7,
     2            3,4,1,2,7,8,5,6,
     3            4,3,2,1,8,7,6,5,
     4            5,6,7,8,1,2,3,4,
     5            6,5,8,7,2,1,4,3,
     6            7,8,5,6,3,4,1,2,
     7            8,7,6,5,4,3,2,1/
      */
      Matrix d2h(8, 8);

      d2h << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7
	  << 1 << 0 << 3 << 2 << 5 << 4 << 7 << 6
	  << 2 << 3 << 0 << 1 << 6 << 7 << 4 << 5
	  << 3 << 2 << 1 << 0 << 7 << 6 << 5 << 4
	  << 4 << 5 << 6 << 7 << 0 << 1 << 2 << 3
	  << 5 << 4 << 7 << 6 << 1 << 0 << 3 << 2
	  << 6 << 7 << 4 << 5 << 2 << 3 << 0 << 1 
	  << 7 << 6 << 5 << 4 << 3 << 2 << 1 << 0;

      groupTable.resize(8, 8);
      for (int i = 0; i < 8; ++i)
	for (int j = 0; j < 8; ++j)
	  groupTable(i, j) = static_cast<int>(d2h.element(i, j));
    }
  else if (sym == "dinfh") {
    groupTable.resize(4, 4);
    groupTable(0, 0) = 0;
    groupTable(0, 1) = 1;
    groupTable(0, 2) = 2;
    groupTable(0, 3) = 3;
    
    groupTable(1, 0) = 1;
    groupTable(1, 1) = 0;
    groupTable(1, 2) = 3;
    groupTable(1, 3) = 2;
    
    groupTable(2, 0) = 2;
    groupTable(2, 1) = 3;
    groupTable(2, 2) = 0;
    groupTable(2, 3) = 1;
    
    groupTable(3, 0) = 3;
    groupTable(3, 1) = 2;
    groupTable(3, 2) = 1;
    groupTable(3, 3) = 0;
  }
  else
    abort();

}

int Symmetry::sizeofIrrep(int irrep)
{
  if (sym == "dinfh")
    return irrep > 3 ? 2 : 1;
  else
    return 1;
}

std::vector<int> Symmetry::add(int irrepl, int irrepr)
{
  if (sym == "dinfh") {
    std::vector<int> vec;
    int goru = ((irrepl%2==0 && irrepr%2==0) || (irrepl%2==1 && irrepr%2==1)) ? 0 : 1;
    
    if (irrepl < 4 && irrepr< 4) {
      vec.push_back( groupTable(irrepl, irrepr));
      return vec;
    }
    else if (irrepl <4 && irrepr >= 4) {
      vec.push_back(2*abs(irrepr/2) + goru);
      return vec;
    } 
    else if (irrepl >= 4 && irrepr <4) {
      vec.push_back(2*abs(irrepl/2) + goru);
      return vec;
    }
    else {
      int irrep1 = 2*abs(irrepl/2 - irrepr/2) + goru;
      int irrep2 = 2*abs(irrepl/2 + irrepr/2) + goru - 2;
      
      if (irrep1 >=2) irrep1 += 2;
      vec.push_back(irrep1);
      if (irrep1 < 3)
	vec.push_back(irrep1+2);
      vec.push_back(irrep2);
      
      return vec;
    }
  }
  else if (sym == "c1") {
    std::vector<int> vec;
    vec.push_back(0);
    return vec;
  }
  else {
    std::vector<int> vec;
    vec.push_back( groupTable(irrepl, irrepr));
    return vec;
  }
  
}



double Symmetry::spatial_sixj(int j1, int j2, int j3, int j5, int j4, int j7) {
    double out = spatial_ninej(j1, j2, j3, j4, j5, j3, j7, j7, 0);
    return out;
}

double Symmetry::spatial_ninej(int j1, int j2, int j12, int j3, int j4, int j34, int j13, int j24, int j) {
  
    // all the numbers are irreps
    int m = 0; //since 9-j does not depend on m, we use m=j/2 to calculate the coefficient
    
    double out = 0.0;

    int m1step = sizeofIrrep(j1);
    int m2step = sizeofIrrep(j2);
    int m3step = sizeofIrrep(j3);
    int m4step = sizeofIrrep(j4);
    int m12step = sizeofIrrep(j12);
    int m13step = sizeofIrrep(j13);
    int m34step = sizeofIrrep(j34);
    int m24step = sizeofIrrep(j24);
    
    for (int m1=0;m1<m1step;m1++)
      for (int m2=0;m2<m2step;m2++)
	for (int m3=0;m3<m3step;m3++)
	  for (int m4=0;m4<m4step;m4++)
	    {
	      double first = 0.0, second = 0.0;
	      for (int m12=0;m12<m12step;m12++)
		for (int m34=0;m34<m34step;m34++)
		  {
		    first +=  spatial_cg(j1, j2, j12, m1, m2, m12) *
		      spatial_cg(j3, j4, j34, m3, m4, m34) *
		      spatial_cg(j12, j34, j, m12, m34, m) ;
		  } 
	      
	      for (int m13=0;m13<m13step;m13++)
		for (int m24=0;m24<m24step;m24++)
		  {
		    second +=  spatial_cg(j1, j3, j13, m1, m3, m13) *
		      spatial_cg(j2, j4, j24, m2, m4, m24) *
		      spatial_cg(j13, j24, j, m13, m24, m) ;
		  } 
	      out += first*second;
	      //cout << m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<first<<"  "<<second<<" "<<out<<endl;
	    }
    return out;
}

double Symmetry::spatial_cg(int a, int b, int c, int rowa, int rowb, int rowc) {
  if (sym == "dinfh") {
    if (a<4 && rowa != 0)
      { cerr<<"a= "<<a<<" and row = "<<rowa<<endl; exit(0);} 
    if (b<4 && rowb != 0)
      { cerr<<"b= "<<b<<" and row = "<<rowb<<endl; exit(0);} 
    if (c<4 && rowc != 0)
      { cerr<<"c= "<<c<<" and row = "<<rowc<<endl; exit(0);} 
    int la, lb, lc;
    la = (2*rowa-1) * (max(0,a-2))/2; 
    lb = (2*rowb-1) * (max(0,b-2))/2; 
    lc = (2*rowc-1) * (max(0,c-2))/2; 
    //a, b, c, are irreps
    // la, lb, lc are z- angular momentums e.g. la = a/2, -a/2
    if (a<4 && b<4 && c<4) {
      if (c == groupTable(a,b))
	return 1.0;
      else
	return 0.0;
    }
    double out = 1.0;
    if (lc != la+lb)
      return 0.0;
    if ((a+b+c)%2 != 0 )
      return 0.0;
    if(lc == la+lb && lc == 0 && la != 0) {
      if (c == 0 || c==1 || (c==2&&la>0) || (c==3&&la>0))
	out = 1.0/sqrt(2.0);
      else
	out = -1.0/sqrt(2.0);
      
    }
    if (lc <0) {
      if (a == 2 || a==3 || b==2 || b==3)
	return -1.0;
    }
    
    return out;
  }
  else if (sym == "c1") {
    return 1.0;
  }
  else {
    if (c == groupTable(a,b))
      return 1.0;
    else
      return 0.0;
  }
}

  
}
