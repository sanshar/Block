/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "nonAbelianGroup.h"
#include "Symmetry.h"
#include <boost/lexical_cast.hpp>
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
SpinAdapted::nonAbelianGroup nonAbelianGrp;

namespace SpinAdapted {

void Symmetry::InitialiseTable(string psym)
{
  if (sym == "c1")
    {
      groupTable.resize(1, 1);
      groupTable(0, 0) = 0;
    }
  else if (sym == "ci" || sym == "cs" || sym == "c2")
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
  else if (sym == "c2v" || sym == "c2h" || sym == "d2")
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
    NonabelianSym = true;
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

  else if (sym == "c3v") {
    NonabelianSym = true;
    nonAbelianGrp = C3v();
  }
  else if (sym == "c5v") {
    NonabelianSym = true;
    nonAbelianGrp = C5v();
  }
  else if (sym == "d5h") {
    NonabelianSym = true;
    nonAbelianGrp = D5h();
  }
  else if (sym == "d4h") {
    NonabelianSym = true;
    nonAbelianGrp = D4h();
  }

  else if (sym == "trans") {
    //do nothing;
  }
  else if (sym == "lzsym") {
    //do nothing;
  }
  else if (sym == "dinfh_abelian") {
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
  else {
    pout << "Symmetry of the molecule has to be one of c1, ci, cs, c2, c2h, c2v, d2, d2h or dinfh"<<endl;
    pout << "Symmetry provided in the input file "<<sym<<endl;
    abort();
  }

}

bool Symmetry::irrepAllowed(int irrep)
{
  if ( (sym == "dinfh"|| sym == "dinfh_abelian")) {
    if ((irrep<0 && irrep >-4) || irrep == 2 || irrep == 3) {
      pout << "Orbital cannot have an irreducible representation of "<<irrep+1<<"  with dinfh symmetry"<<endl;
      abort();
    }
    else
      return true;
  }
  if (sym == "d2h" && (irrep<0 || irrep >= 8)) {
    pout << "Orbital cannot have an irreducible representation of "<<irrep+1<<"  with "<<sym<<" symmetry"<<endl;
    abort();
  }
  if ((sym == "c2v" || sym == "c2h" || sym == "d2") && (irrep<0 || irrep >= 4)) {
    pout << "Orbital cannot have an irreducible representation of "<<irrep+1<<"  with "<<sym<<" symmetry"<<endl;
    abort();
  }
  if ((sym == "ci" || sym == "c2" || sym == "cs" ) && (irrep <0 || irrep >=2)) {
    pout << "Orbital cannot have an irreducible representation of "<<irrep+1<<"  with "<<sym<<" symmetry"<<endl;
    abort();
  }
  if (sym == "c1" && irrep != 0) {
    pout << "Orbital cannot have an irreducible representation of "<<irrep+1<<"  with "<<sym<<" symmetry"<<endl;
    abort();
  }
  if ( (NonabelianSym) && (irrep < 0 || irrep >= nonAbelianGrp.getNumIrreps())) {
    pout << "Orbital cannot have an irreducible representation of "<<irrep+1<<"  with "<<sym<<" symmetry"<<endl;
    abort();
  }
  if (sym == "trans") {
    std::vector<int> irreps = decompress(irrep);
    if (irreps[0] >= NPROP[0] || irreps[0]< 0 ||
	irreps[1] >= NPROP[1] || irreps[1]< 0 ||
	irreps[2] >= NPROP[2] || irreps[2]< 0 ) {
      
      pout << "decompressing the irrep "<<irrep<<" leads to k points "<<irreps[0]<<"  "<<irreps[1]<<"  "<<irreps[2]<<endl;
      abort();
    }
  }
  if (sym == "lzsym") {
    return true;
  }
  return true;
}

std::vector<int> Symmetry::decompress(int pirrep) 
{
  //this is used to decompress the irrep to 3 k points
  std::vector<int> out(3,0);
  int irrep = abs(pirrep);
  out[2] = irrep/PROPBITLEN/PROPBITLEN;
  out[1] = (irrep - out[2]*PROPBITLEN*PROPBITLEN)/PROPBITLEN;
  out[0] = irrep - out[2]*PROPBITLEN*PROPBITLEN - out[1]*PROPBITLEN;
  if (irrep == -pirrep) {
    out[0] *=-1;
    out[1] *=-1;
    out[2] *=-1;
  }
  return out;
}

int Symmetry::compress(std::vector<int>& irreps) 
{
  return irreps[0] + irreps[2]*PROPBITLEN*PROPBITLEN + irreps[1]*PROPBITLEN;
}


string Symmetry::stringOfIrrep(int irrep) 
{
   string symbol;
  if (sym == "d2h") {
    switch(irrep)
      {
      case(0): 
	symbol = "Ag"; break;
      case(1):
	symbol = "B3u"; break;
      case(2):
	symbol = "B2u"; break;
      case(3):
	symbol = "B1g"; break;
      case(4):
	symbol = "B1u"; break;
      case(5):
	symbol = "B2g"; break;
      case(6):
	symbol = "B3g"; break;
      case(7):
	symbol = "Au"; break;
      }
  }
  else if (sym == "c2v") {
    switch(irrep)
      {
      case(0): 
	symbol = "A1"; break;
      case(1):
	symbol = "B1"; break;
      case(2):
	symbol = "B2"; break;
      case(3):
	symbol = "A2"; break;
      }
    }
  else if (sym == "c2h") {
    switch(irrep)
      {
      case(0): 
	symbol = "Ag"; break;
      case(1):
	symbol = "Au"; break;
      case(2):
	symbol = "Bu"; break;
      case(3):
	symbol = "Bg"; break;
      }
  }
  else if (sym == "d2") {
    switch(irrep)
      {
      case(0): 
	symbol = "A"; break;
      case(1):
	symbol = "B3";//B3
   break;
      case(2):
	symbol = "B2";//B1
   break;
      case(3):
	symbol = "B1";//B2
   break;
      }
    }
  else if (sym == "cs") 
    symbol = (irrep == 0) ? "A'" : "A''";
  else if (sym == "c2") 
    symbol = (irrep == 0) ? "A" : "B";
  else if (sym == "ci") 
    symbol = (irrep == 0) ? "Ag" : "Au";
  else if (sym == "dinfh") {
    string output = "";
    char goru = irrep%2 == 0 ? 'g' : 'u';
    int lz = max(0,(abs(irrep)-2)/2);
    lz *= irrep<0 ? -1 : 1; 
    output+= boost::lexical_cast<string>(lz);
    output+=goru;
    if (irrep <2) output+= '+';
    else if (irrep >=2 && irrep <4 ) output+= '-';
    symbol = output;
  }
  else if (NonabelianSym) {
    return nonAbelianGrp.getIrrepName(irrep);
  }
  else if(sym == "trans") {
    std::vector<int> irreps = decompress(irrep);
    string output = "";
    output+=boost::lexical_cast<string>(irreps[2]);
    output+=boost::lexical_cast<string>(irreps[1]);
    output+=boost::lexical_cast<string>(irreps[0]);
    symbol = output;
  }
  else if (sym == "lzsym") {
    symbol = boost::lexical_cast<string>(irrep);
  }
  else if (sym == "dinfh_abelian") {
    string output = "";
    char goru = irrep%2 == 0 ? 'g' : 'u';
    int lz = max(0,(abs(irrep)-2)/2);
    lz *= irrep<0 ? -1 : 1; 
    output+= boost::lexical_cast<string>(lz);
    output+=goru;
    symbol = output;
  }
  else 
    symbol = "A";
  return symbol;
}

int Symmetry::sizeofIrrep(int irrep)
{
  if (sym == "dinfh")
    return irrep > 3 ? 2 : 1;
  else if (NonabelianSym)
    return nonAbelianGrp.getIrrepSize(irrep);
  else
    return 1;

}

int Symmetry::negativeof(int irrep)
{
  if (sym == "trans") {
    std::vector<int> lirrep = decompress(irrep);
    for (int i=0; i<lirrep.size(); i++) {
      lirrep[i] = lirrep[i] == 0 ? lirrep[i] : (NPROP[i] - lirrep[i]);//%NPROP[i];
      if(lirrep[i] >= NPROP[i] || lirrep[i] < 0) {
	pout << "cannot find the negative of "<<i<<" component of lirrep "<<lirrep[i]<<endl;
	pout << "it is not in the first bruillion zone"<<endl;
	exit(0);
      }
    }
    int outirrep = compress(lirrep);
    return outirrep;
  }
  else if (sym == "lzsym") {
    return -irrep;
  }
  else if (sym == "dinfh_abelian") {
    if (irrep >= 0 && irrep < 4)
      return irrep;
    else
      return -irrep;
  }
  else
    return irrep;
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
  else if (NonabelianSym) {
    return nonAbelianGrp.getProduct(irrepl, irrepr);
  }
  else if (sym == "dinfh_abelian") {
    std::vector<int> vec;
    int goru = ((abs(irrepl)%2==0 && abs(irrepr)%2==0) || (abs(irrepl)%2==1 && abs(irrepr)%2==1)) ? 0 : 1;
    
    if (abs(irrepl) < 4 && abs(irrepr)< 4) {
      vec.push_back( groupTable(irrepl, irrepr));
      return vec;
    }
    else if (abs(irrepl) <4 && abs(irrepr) >= 4) {
      int irrepout = 2*(abs(irrepr)/2) + goru;
      if (irrepr <0)
	vec.push_back(-irrepout);
      else 
	vec.push_back(irrepout);
      return vec;
    } 
    else if (abs(irrepl) >= 4 && abs(irrepr) <4) {
      int irrepout = 2*(abs(irrepl)/2) + goru;
      if (irrepl <0)
	vec.push_back(-irrepout);
      else 
	vec.push_back(irrepout);
      return vec;
    }
    else {
      int irrep1 = 2*abs(abs(irrepl)/2 - abs(irrepr)/2) + goru;
      int irrep2 = 2*abs(abs(irrepl)/2 + abs(irrepr)/2) + goru - 2;

      int irrep3 = 0;
      if (irrep1 >=2) irrep1 += 2;

      if (irrepl*irrepr > 0) {
	if (irrepl <0)
	  vec.push_back(-irrep2);
	else
	  vec.push_back(irrep2);
      }
      else {
	if (irrepl + irrepr > 0) {
	  if (irrep1 >=2 ) vec.push_back(irrep1);
	  else vec.push_back(abs(irrep1));
	}
	else {
	  if (irrep1 >=2 ) vec.push_back(-irrep1);
	  else vec.push_back(abs(irrep1));
	}
	
      }
      
      return vec;

    }
  }
  else if (sym == "c1") {
    std::vector<int> vec;
    vec.push_back(0);
    return vec;
  }
  else if(sym == "lzsym") {
    std::vector<int> vec;
    vec.push_back(irrepl+irrepr);
    return vec;
  }
  else if (sym == "trans") {
    std::vector<int> vec;
    std::vector<int> lirrep = decompress(irrepl);
    std::vector<int> rirrep = decompress(irrepr);
    std::vector<int> out(3,0);
    out[0] = (lirrep[0]+rirrep[0])%NPROP[0];
    out[1] = (lirrep[1]+rirrep[1])%NPROP[1];
    out[2] = (lirrep[2]+rirrep[2])%NPROP[2];
    int outirrep = compress(out);
    vec.push_back(outirrep);
    return vec;
  }
  else {
    std::vector<int> vec;
    vec.push_back( groupTable(irrepl, irrepr));
    return vec;
  }
}



double Symmetry::spatial_sixj(int j1, int j2, int j3, int j5, int j4, int j7) {
  if (! (NonabelianSym) ) {
    if (j3 != add(j1,j2)[0]) return 0.0;
    if (j7 != add(j2,j5)[0]) return 0.0;
    if (j4 != add(j3,j5)[0]) return 0.0;
    return 1.0;
  }
  else {
    //double out = spatial_ninej(j1, j2, j3, j4, j5, j3, j7, j7, 0);
    double out = spatial_ninej(j1, j2, j3, 0, j5, j5, j1, j7, j4);
    return out;
  }
}

double Symmetry::spatial_ninej(int j1, int j2, int j12, int j3, int j4, int j34, int j13, int j24, int j) {
  
  if (!(NonabelianSym)) {
    return 1.0;
  }

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
	      //pout << m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<first<<"  "<<second<<" "<<out<<endl;
	    }
    return out;
}

double Symmetry::spatial_cg(int a, int b, int c, int rowa, int rowb, int rowc) {

  if (sym == "dinfh") {
    if (a<4 && rowa != 0)
      { pout<<"a= "<<a<<" and row = "<<rowa<<endl; exit(0);} 
    if (b<4 && rowb != 0)
      { pout<<"b= "<<b<<" and row = "<<rowb<<endl; exit(0);} 
    if (c<4 && rowc != 0)
      { pout<<"c= "<<c<<" and row = "<<rowc<<endl; exit(0);} 
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
  else if (NonabelianSym)
    return nonAbelianGrp.getCG(a, b, c, rowa, rowb, rowc);
  else {
    if (c == add(a,b)[0])
      return 1.0;
    else
      return 0.0;
  }
}

  
}
