
#include "nevpt2_info.h"
#include "assert.h"
#include <vector>
#include <boost/serialization/serialization.hpp>
#include "nevpt2_operators.h"
#include "nevpt2_util.h"

namespace SpinAdapted{
  //============================================================================
  // get the orbital window
  //    OrbWin[0] => first internal orbital
  //    OrbWin[1] => last internal orbital
  //    OrbWin[2] => first active orbital
  //    OrbWin[3] => last active orbital
  //    OrbWin[4] => first virtual orbital
  //    OrbWin[5] => last virtual orbital
  //============================================================================
  void NEVPT2Info::GetOrbWin(int *OW){
    if (!HaveData) ReadData();
    OW[0] = OrbWin[0];
    OW[1] = OrbWin[1];
    OW[2] = OrbWin[2];
    OW[3] = OrbWin[3];
    OW[4] = OrbWin[4];
    OW[5] = OrbWin[5];
    OW[6] = OrbWin[6];
  }

  //============================================================================
  //get the nuclear energy
  //============================================================================
  double NEVPT2Info::GetENuc(){
    if (!HaveData) ReadData();
    return ENuc;
  }
  
  //============================================================================
  //get the effective one-electron matrices
  //============================================================================
  void NEVPT2Info::GetH(vector<Matrix> &h)const{
    if (h.size()!=H.size()) h.resize(H.size());
    for (int i=0;i<H.size();i++){
      h[i] = H[i];
    }//i
  }
  
  //============================================================================
  //set the effective one-electron matrices
  //============================================================================
  void NEVPT2Info::AddH(Matrix &h){
    H.push_back(h);
  }
  
  //============================================================================
  //get the time required for a given perturbation class
  //============================================================================
  void NEVPT2Info::GetTimePerClass(vector<double> &TPC)const{
    int size = TimePerClass.size();
    //make sure the output vector has the right size  
    if (TPC.size()!=8) {
      TPC.clear();
      TPC.resize(size);
    }
    //copy the vector
    for (int i=0;i<size;i++){
      TPC[i] = TimePerClass[i];
    }
  }
  
  //============================================================================
  //get the energy for a given root
  //============================================================================
  void NEVPT2Info::getE(int iroot, vector<double>& e)const{
    assert((E.size()>iroot)&&(iroot>=0));
    //check the size of the output vector
    if (e.size()!=E[iroot].size()) e.resize(E[iroot].size());
    for (int iclass=0;iclass<e.size();iclass++){
      e[iclass] = E[iroot][iclass];
    }//iclass
  }
  
  //============================================================================
  //get all avaliable energies
  //============================================================================
  void NEVPT2Info::getE(vector<vector<double> >& e)const{
    //check the size of the vector
    if (e.size()!=E.size())e.resize(E.size()); 
    for (int iroot=0;iroot<E.size();iroot++){
      //check the size of the vector
      if (e[iroot].size()!=E[iroot].size())e[iroot].resize(E[iroot].size());
      for (int iclass=0;iclass<E[iroot].size();iclass++){
        e[iroot][iclass] = E[iroot][iclass];
      }//iclass
    }//iroot
  }

  //============================================================================
  //get all zero'th order valence energies
  //============================================================================
  void NEVPT2Info::getE0(vector<double> &e0)const{
    //check the size of the output vector
    if (e0.size()!=E0.size()) e0.resize(E0.size());
    for (int iroot=0;iroot<e0.size();iroot++){
      e0[iroot] = E0[iroot];
    }//iroot
  }

  
  //============================================================================
  //get the conventional overlap flag
  //============================================================================
  bool  NEVPT2Info::ConventionalOverlap(){
    if (!HaveData) ReadData();
    return ConvOverlap;
  }
  
  //============================================================================
  //set the number of roots. Note: the size fo the energy vectors is immediately
  // changed accordingly
  //============================================================================
  void NEVPT2Info::SetNRoots(int nr){
      nroots = nr;
      E.resize(nroots);
      for (int i=0;i<nroots;i++){
        E[i].resize(8);
      }
      E0.resize(nroots);
      EVal.resize(nroots);
    }
  
  //set the energy for a given root for all eight classes
  void NEVPT2Info::SetE(int iroot, const vector<double> &val){
    assert(E.size()>=iroot);
    if (val.size()!=E[iroot].size())E[iroot].resize(val.size());
    for (int iclass=0;iclass<E[iroot].size();iclass++){
      E[iroot][iclass] = val[iclass];
    }
  }
  
  //============================================================================
  //get the orbital ordering
  //============================================================================
  void NEVPT2Info::GetOrbOrder(vector<int> &orborder){
    if (!HaveOrbOrder)ReadOrbOrder();
    //check the size of the vector
    if (orborder.size()!=OrbOrder.size()) orborder.resize(OrbOrder.size());
    //copy the vector
    for (int iorb=0;iorb<OrbOrder.size();iorb++){
      orborder[iorb] = OrbOrder[iorb];
    }//iorb
  }
  
  //============================================================================
  //get the orbital energies
  //============================================================================
  void NEVPT2Info::GetOrbEnergies(vector<double> &orbenergies){
    if (!HaveOrbEnergies) ReadOrbEnergies();
    //check the size of the output vector
    if (orbenergies.size()!=EOrb.size()) orbenergies.resize(EOrb.size());
    //copy the vector
    for  (int iorb=0;iorb<EOrb.size();iorb++){
      orbenergies[iorb] = EOrb[iorb];
    }//iorb
  }
  
  
  //============================================================================
  //determine the boundaries for the two NEVPT2 sweeps
  //============================================================================
  void NEVPT2Info::CalcMaxBlockIter(int n_iters){
    if (n_iters%2==0){
      MaxBlockIter[0] = n_iters/2-1;
      MaxBlockIter[1] = n_iters/2-1;
    }
    else{
      MaxBlockIter[0] = n_iters/2;
      MaxBlockIter[1] = n_iters/2-1;
    }
  }
  
  //============================================================================
  //the function that reads the orbital order (if there is any given)
  //============================================================================
  void NEVPT2Info::ReadOrbOrder(){
    int NActive;
    char s[512];
    char msg[512];
    bool ReOrder = false;
    int res=1;
    
    //get the number of active orbitals
    if (HaveData){
      NActive = OrbWin[3]-OrbWin[2]+1;
    }
    else{
      ReadData();
      NActive = OrbWin[3]-OrbWin[2]+1;
    }
    
    OrbOrder.resize(NActive);
    for (int p=0;p<NActive;p++){
      OrbOrder[p] = dmrginp.reorder_vector()[p]+1;
    }
    /*
    //try to open the Fiedler vector file
    FILE *FiedlerFile;
    
    sprintf(msg, "%s%s", dmrginp.save_prefix().c_str(), "/fiedler_reorder.dat");
    FiedlerFile = fopen(msg,"r");
    if (FiedlerFile!=0){
      ReOrder = true;
      while(res==1){
        int i;
        res = fscanf(FiedlerFile,"%i",&i);
        if (res!=1) break;
        OrbOrder.push_back(i);  
      }
      fclose(FiedlerFile);
    }
    else{
      FILE *GeneticFile;
      sprintf(msg, "%s%s", dmrginp.save_prefix().c_str(), "/genetic_reorder.dat");
      GeneticFile = fopen(msg,"r");
      if (GeneticFile !=0){
        ReOrder = true;
        while (res==1){ 
          int i;
          res = fscanf(GeneticFile,"%i",&i);
          if (res!=1) break;
          OrbOrder.push_back(i);  
        }
        pout << s;
        fclose(GeneticFile);
      }
    }
    if (ReOrder){
      if (OrbOrder.size()!=NActive){
        sprintf(msg,"\nERROR BLOCK-NEVPT2: Number of reordered orbitals does not equal number of active orbitals");
        pout << msg;
        exit(0);
      }
    }//ReOrder
    else{
      OrbOrder.resize(NActive);
      for (int i=0;i<OrbOrder.size();i++){
        OrbOrder[i] = i+1;
      }
    }//use given order
     */
    PrintVector(OrbOrder,"OrbOrder.tmp");
    HaveOrbOrder = true;
  }
  
  //============================================================================
  //the function that reads in the nevpt2 input file
  //============================================================================
  void NEVPT2Info::ReadData(){
    FILE *f;
    int res;
    f = fopen("dmrg.nevpt2.inp","r");
    res = fscanf(f,"%s",BaseName);
    res = fscanf(f,"%i",&OrbWin[0]);
    res = fscanf(f,"%i",&OrbWin[1]);
    res = fscanf(f,"%i",&OrbWin[2]);
    res = fscanf(f,"%i",&OrbWin[3]);
    res = fscanf(f,"%i",&OrbWin[4]);
    res = fscanf(f,"%i",&OrbWin[5]);
    res = fscanf(f,"%lf",&ENuc);
    res = fscanf(f,"%i",&MaxCore);
    res = fscanf(f,"%lf",&RefDensWeight);
    int convOverlap=0;
    res = fscanf(f,"%i",&convOverlap);
    ConvOverlap = (bool) ConvOverlap;
    ConventionalNevpt2 = ((dmrginp.nevpt2())&&(dmrginp.read_higherpdm()));
    fclose(f);
    HaveData=true;
  }
  
  //============================================================================
  //the copy function
  //============================================================================
  NEVPT2Info NEVPT2Info::operator =(NEVPT2Info &Info){
      nroots = Info.getNRoots();
      for (int i =0;i<_NEV_TIME_DIM_;i++){
        Time[i] = Info.getTime(i);
      }
      Info.GetTimePerClass(TimePerClass);
      Info.getE(E);
      Info.getE0(E0);
      EVal.resize(nroots);
      for (int i=0;i<nroots;i++){
        EVal[i] = Info.getEVal(i);
      }
      ENuc   = Info.GetENuc();
      Info.GetOrbWin(OrbWin);
      Info.GetOrbEnergies(EOrb);
      ConvOverlap = Info.ConventionalOverlap();
      Info.GetOrbOrder(OrbOrder);
      Info.getBaseName(BaseName);
      HaveOrbOrder = Info.OrbOrderAvailable();
      HaveData = Info.DataAvailable();
      HaveOrbEnergies = Info.OrbEnergiesAvailable();
      return *this;
  }
  
  //============================================================================
  //Read the orbital energies form disk
  //============================================================================
  void NEVPT2Info::ReadOrbEnergies(){
    FILE *f;
    char msg[512];
    //make sure we have all the required information
    if (!HaveData) ReadData();
    if (!HaveOrbOrder) ReadOrbOrder();
    //the orbital spaces
    int i0 = OrbWin[0];//core
    int i1 = OrbWin[1];//core
    int t0 = OrbWin[2];//core
    int t1 = OrbWin[3];//core
    int a0 = OrbWin[4];//core
    int a1 = OrbWin[5];//core
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    int OrbDim = NInternal+NActive+NExternal;
    int i;
    //create a reorder vector for the complete set of orbitals
    vector <int> reorder;
    reorder.resize(a1-i0+1,0.0);
    //the internal part
    for (i=0;i<NInternal;i++){
      reorder[i] = i;
    }
    //the active part
    for (i=NInternal;i<NInternal+NActive;i++){
      reorder[i] = OrbOrder[i-NInternal]+NInternal-1;
    }
    //the external part
    for (i=NInternal+NActive;i<NInternal+NActive+NExternal;i++){
      reorder[i] = i;
    }
    //actually read the orbital energies
    EOrb.resize(OrbDim);
    sprintf(msg,"%s.EOrb.tmp",BaseName);
    f = fopen(msg,"r");
    double val = 0.0;
    for (i=0;i<EOrb.size();i++){
      int res = fscanf(f,"%lf",&val);
      EOrb[reorder[i]] = val;
    }
    fclose(f);
    //set the flag
    HaveOrbEnergies=true;
  }
  
  
}



