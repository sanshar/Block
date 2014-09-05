/* 
 * File:   nevpt2_info.h
 * Author: roemelt
 *
 * Created on January 17, 2014, 4:55 PM
 */

#ifndef NEVPT2_INFO_H
#define	NEVPT2_INFO_H

#include <vector>
#include <newmat.h>
#include <stdio.h>

#define _NEV_TIME_DIM_ 5

namespace SpinAdapted{
  class NEVPT2Info{
  private:
    //calculation infos
    int nroots;
    int OrbWin[6];
    bool ConvOverlap;
    bool ConventionalNevpt2;
    char BaseName[512];
    int MaxBlockIter[2];
    int NevSweep;
    double RefDensWeight;
    
    //important quantities
    double Time[_NEV_TIME_DIM_];
    vector<double> TimePerClass; 
    vector<vector<double> > E;
    vector<double> E0;
    vector<double> EVal;
    vector<double> EOrb;
    double ENuc;
    vector<int> OrbOrder;
    vector<Matrix> H;
    int MaxCore;
    
    //Data availablity flags
    bool HaveOrbOrder;
    bool HaveData;
    bool HaveOrbEnergies;
    
  public:
    NEVPT2Info(){
      SetNRoots(1);
      sprintf(BaseName,"dmrg");
      ENuc    = 0.0;
      for (int i =0;i<_NEV_TIME_DIM_;i++){
        Time[i] = 0.0;
      }
      TimePerClass.resize(8);
      OrbWin[0] =-1;
      OrbWin[0] =-1;
      OrbWin[0] =-1;
      OrbWin[0] =-1;
      OrbWin[0] =-1;
      OrbWin[0] =-1;
      MaxBlockIter[0] =-1;
      MaxBlockIter[1] =-1;
      HaveOrbOrder = false;
      HaveData     = false;
      ConvOverlap = false;
      HaveOrbEnergies = false;
      ConventionalNevpt2 = false;
      MaxCore = 1024;
      NevSweep=0;
      RefDensWeight=0.5;
    };
    //--------------------------
    //Constructor and Dsetructor
    //--------------------------
    NEVPT2Info(NEVPT2Info &Info){
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
    };
    ~NEVPT2Info(){
      E.clear();
      E0.clear();
      EVal.clear();
      EOrb.clear();
      OrbOrder.clear();
      H.clear();
    };
    
    //-----------
    //the getters
    //-----------
    int getNRoots()const{return nroots;}
    double getTime(int i)const{return Time[i];}
    void GetTimePerClass(vector<double> &TPC)const;
    double GetTimePerClass(int iclass)const{if (iclass<8)return TimePerClass[iclass];else return 0.0;}
    double getE(int iroot, int iclass)const{return E[iroot][iclass];}
    void getE(int iroot, vector<double> &e)const;
    void getE(vector<vector<double> > &e)const;
    double getE0(int iroot)const{return E0[iroot];}
    void getE0(vector<double> &e0)const;
    double getEVal(int iroot)const{return EVal[iroot];}
    void getBaseName(char *BN)const{sprintf(BN,"%s",BaseName);}
    void GetOrbOrder(vector<int> &orborder);
    void GetOrbEnergies(vector<double> &orbenergies);
    double GetENuc();
    void GetOrbWin(int *OW);
    bool GetConvOverlap();
    bool OrbOrderAvail()const{return HaveOrbOrder;}
    bool OrbEnergiesAvail()const{return HaveOrbEnergies;}
    bool DataAvail()const{return HaveData;}
    void GetH(int i, Matrix &h)const{h=H[i];}
    void GetH(vector<Matrix> &h)const;
    bool ConventionalOverlap();
    bool ConvNevPT2() const{return ConventionalNevpt2;}
    bool OrbOrderAvailable() const{return HaveOrbOrder;}
    bool OrbEnergiesAvailable() const{return HaveOrbEnergies;}
    bool DataAvailable() const{return HaveData;}
    int GetMaxBlockIter(int i) const{if (i<2)return MaxBlockIter[i];else return -1;}
    int GetMaxCore() const{return MaxCore;}
    int GetNevSweep() const{return NevSweep;}
    double GetRefDensWeight()const{return RefDensWeight;}
    
    //-----------
    //the setters
    //-----------
    void SetNRoots(int nr);
    void SetTime(int i, double t){Time[i] = t;}
    void AddTime(int i, double t){Time[i]+= t;}
    void SetTimePerClass(int i, double t){if (i<8)TimePerClass[i] = t;}
    void AddTimePerClass(int i, double t){if (i<8)TimePerClass[i]+= t;}
    void SetE(int iroot, int iclass, double val){E[iroot][iclass] = val;}
    void SetE(int iroot, const vector<double> &val);
    void SetE0(int iroot, double val){E0[iroot] = val;}
    void SetEVal(int iroot, double val){EVal[iroot] = val;}
    void SetENuc(double enuc){ENuc = enuc;}
    void SetBaseName(char *BN){sprintf(BaseName,"%s",BN);}
    void SetConvOverlap(bool co){ConvOverlap = co;}
    void AddH(Matrix &h);
    void CalcMaxBlockIter(int n_iters);
    void SetMaxCore(int MC){MaxCore = MC;}
    void ResetNevSweep(){NevSweep=0;}
    void IncrementNevSweep(){NevSweep++;}
    void SetRefDensWeight(double rdw){RefDensWeight=rdw;}
    
    //-------------------------------------------------
    //the functions that read the information from disk
    //-------------------------------------------------
    void ReadData();
    void ReadOrbOrder();
    void ReadOrbEnergies();
    
    //-----------------
    //the copy function
    //-----------------
    NEVPT2Info operator=(NEVPT2Info &Info);
    
  };
  

}

#endif	/* NEVPT2_INFO_H */

