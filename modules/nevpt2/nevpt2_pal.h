/* 
 * File:   nevpt2_pal.h
 * Author: mroemelt
 *
 * Created on July 17, 2014, 9:38 AM
 */

#ifndef NEVPT2_PAL_H
#define	NEVPT2_PAL_H

#include "wavefunction.h"
#include "nevpt2_operators.h"
#include <boost/serialization/serialization.hpp>

namespace SpinAdapted{
  //============================================================================
  // Add all contributions from different processes to a wavefunction
  //============================================================================
  void AddPalWavefunction(Wavefunction &WF);

  //============================================================================
  // Add all contributions from different processes to a set of wavefunctions
  //============================================================================
  void AddPalWavefunctions(vector<boost::shared_ptr<Wavefunction> > &VTrip, 
                          vector<boost::shared_ptr<Wavefunction> > &SigmaTrip, 
                           boost::shared_ptr<Wavefunction> &VSing,
                           boost::shared_ptr<Wavefunction> &SigmaSing);
  
  //============================================================================
  // Add all contributions from different processes to a set of wavefunctions
  // In this version we have an entire batch of sets of wavefunctions
  //============================================================================
  void AddPalWavefunctions(vector<vector<boost::shared_ptr<Wavefunction> > >&VTrip, 
                           vector<vector<boost::shared_ptr<Wavefunction> > >&SigmaTrip, 
                           vector<boost::shared_ptr<Wavefunction> >&VSing,
                           vector<boost::shared_ptr<Wavefunction> >&SigmaSing);

  //============================================================================
  // Send a local function V(tuv) and integral matrix around in a circle to use 
  // them in each process to evaluate local V(a) functions
  //            V(a) += Ktv(u,a) * Vtuv
  // so u must be u + NInternal
  //============================================================================
  void SendAroundVtuv_a(Wavefunction &Vtuv, Matrix &Ktv, WavefunctionArray &Ta, int u, bool Dummy);
  void SendAroundVtuv_a(Wavefunction &Vtuv, Wavefunction &Vtvu, Matrix &Kuv, Matrix &Kvu,
                        WavefunctionArray &Ta, int t, bool Dummy);
  void SendAroundVtuv_a(Wavefunction &Vtvu, Wavefunction &Vvtu, Wavefunction &Vtuv, 
                        Matrix &Ktv, Matrix &Kvt, Matrix &Ktu, 
                        WavefunctionArray &Ta, int u, int v, bool Dummy);

  //============================================================================
  // Send a local function V(tuv) and integral matrix around in a circle to use 
  // them in each process to evaluate local V(i) functions
  //            V(i) += Kut(v,i) * Vtuv
  // so u must be v + NInternal
  //============================================================================
  void SendAroundVtuv_i(Wavefunction &Vtuv, Wavefunction &Vutv, Matrix &Ktu, 
                        WavefunctionArray &Ti, int v, bool Dummy);
  void SendAroundVtuv_i(Wavefunction &Vtvu, Wavefunction &VtvuS, Wavefunction &Vvtu, 
                        Wavefunction &VvtuS, Matrix &Ktv, Matrix &Kut, 
                        WavefunctionArray &Ti, int u, int v, bool Dummy, double SingFac);
  void SendAroundVtuv_i(Wavefunction &Vtuv, Matrix &Kut, 
                        WavefunctionArray &Ti, int v, bool Dummy);

  
  //============================================================================
  //here's a class that facilitates the sending and receiving of multiple
  //Wavefunctions and Matrices via MPI
  //============================================================================
  class Packet{
  private:
      vector<boost::shared_ptr<Wavefunction> > WF;
      vector<boost::shared_ptr<Matrix> >M;
      bool Dummy;
      int t;
      int u;
      int v;
      
      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
        ar & WF;
        ar & M;
        ar & Dummy;
        ar & t;
        ar & u;
        ar & v;
      }
    
      
  public:
      //--------------------------
      //constructor and destructor
      //--------------------------
      Packet(){
        Dummy = false;
        t=-1;
        u=-1;
        v=-1;
      }
      Packet(const Packet &P){
          clear();
          for (int i=0;i<P.SizeWF();i++){
              WF.push_back(boost::make_shared<Wavefunction>(*(P.GetWF(i))));
          }
          for (int i=0;i<P.SizeM();i++){
              M.push_back(boost::make_shared<Matrix>(*(P.GetM(i))));
          }
          Dummy = P.dummy();
          t = P.get_t();
          u = P.get_u();
          v = P.get_v();
      }
      Packet(const Wavefunction &wf1, const Wavefunction &wf2, const Wavefunction &wf3, 
             const Matrix &m1, const Matrix &m2, const Matrix &m3,
             const int &T, const int &U, const int V, bool d){
          clear();
          WF.push_back(boost::make_shared<Wavefunction>(wf1));
          WF.push_back(boost::make_shared<Wavefunction>(wf2));
          WF.push_back(boost::make_shared<Wavefunction>(wf3));
          M.push_back(boost::make_shared<Matrix>(m1));
          M.push_back(boost::make_shared<Matrix>(m2));
          M.push_back(boost::make_shared<Matrix>(m3));
          t = T;
          u = U;
          v = V;
          Dummy = d;
      }
      Packet(const Wavefunction &wf1, const Wavefunction &wf2, const Wavefunction &wf3, 
             const Wavefunction &wf4, const Matrix &m1, const Matrix &m2,
             const int &U, const int V, bool d){
          clear();
          WF.push_back(boost::make_shared<Wavefunction>(wf1));
          WF.push_back(boost::make_shared<Wavefunction>(wf2));
          WF.push_back(boost::make_shared<Wavefunction>(wf3));
          WF.push_back(boost::make_shared<Wavefunction>(wf4));
          M.push_back(boost::make_shared<Matrix>(m1));
          M.push_back(boost::make_shared<Matrix>(m2));
          u = U;
          v = V;
          Dummy = d;
      }
      Packet(const Wavefunction &wf1, const Wavefunction &wf2,
             const Matrix &m1,
             const int &V, bool d){
          clear();
          WF.push_back(boost::make_shared<Wavefunction>(wf1));
          WF.push_back(boost::make_shared<Wavefunction>(wf2));
          M.push_back(boost::make_shared<Matrix>(m1));
          v = V;
          Dummy = d;
      }
      Packet(const Wavefunction &wf1, const Wavefunction &wf2,
             const Matrix &m1, const Matrix &m2,
             const int &T, const int &U, bool d){
          clear();
          WF.push_back(boost::make_shared<Wavefunction>(wf1));
          WF.push_back(boost::make_shared<Wavefunction>(wf2));
          M.push_back(boost::make_shared<Matrix>(m1));
          M.push_back(boost::make_shared<Matrix>(m2));
          t = T;
          u = U;
          Dummy = d;
      }
      Packet(const Wavefunction &wf1,
             const Matrix &m1,
             const int &T, bool d){
          clear();
          WF.push_back(boost::make_shared<Wavefunction>(wf1));
          M.push_back(boost::make_shared<Matrix>(m1));
          t = T;
          Dummy = d;
      }
      ~Packet(){
          for (int i=0;i<WF.size();i++){
              WF[i].reset();
          }//i
          WF.clear();
          for (int i=0;i<M.size();i++){
              M[i].reset();
          }//i
          M.clear();
      }
      //-----------
      //the getters
      //-----------
      boost::shared_ptr<Wavefunction> GetWF(int i) const {return WF[i];}
      boost::shared_ptr<Matrix> GetM(int i) const {return M[i];}
      int SizeM() const {return M.size();}
      int SizeWF() const {return WF.size();}
      bool dummy() const {return Dummy;}
      int get_t()const{return t;}
      int get_u()const{return u;}
      int get_v()const{return v;}
      //the copy function
      void operator=(const Packet &P){
          clear();
          for (int i=0;i<P.SizeWF();i++){
              WF.push_back(boost::make_shared<Wavefunction>(*(P.GetWF(i))));
          }
          for (int i=0;i<P.SizeM();i++){
              M.push_back(boost::make_shared<Matrix>(*(P.GetM(i))));
          }
          Dummy = P.dummy();
          t = P.get_t();
          u = P.get_u();
          v = P.get_v();
      }
      //-----------------
      //clear the vectors
      //-----------------
      void clear(){
          for (int i=0;i<WF.size();i++){
              WF[i].reset();
          }//i
          WF.clear();
          for (int i=0;i<M.size();i++){
              M[i].reset();
          }//i
          M.clear();
      }
      //-------------------
      //print the data size
      //-------------------
      void Print(int step=-1){
          char msg[512];
          sprintf(msg,"\nSize:%i  Rank:%i  step=%i",(int)WF.size(),mpi_rank(),step);pout << msg;
          for (int i=0;i<WF.size();i++){
              int dim1 = WF[i]->nrows();
              int dim2 = WF[i]->ncols();
              int dim3 = M[i]->Nrows();
              int dim4 = M[i]->Ncols();
              sprintf(msg,"  DimWF(%i,%i)  DimM(%i,%i)",dim1,dim2,dim3,dim4);pout << msg;
          }
      }
  };

  
  //============================================================================
  //This is a class that simulates the communication between different porcesses
  //via MPI by using a shared disk. This is here just for debugging reasons
  //============================================================================
  class PostMan{
  private:
      char BaseName[512];
      vector< vector<int> > IndexComps;
  public:
      //constructor
      PostMan(){
          sprintf(BaseName,"Message");
      }
      //destructor
      ~PostMan(){
          for (int i=0;i<IndexComps.size();i++){
              IndexComps[i].clear();
          }
          IndexComps.clear();
      };
      //set the Name
      void SetName(char *bn){sprintf(BaseName,"%s",bn);}
      //send a packet
      void send(int dest, int tag, const Packet &P){
          char tmp[512];
          int source = mpi_rank();
          sprintf(tmp,"%s.%i.%i.%i.tmp",BaseName,dest,source,tag);
          ofstream f;
          f.open(tmp,std::ios::binary);
          boost::archive::binary_oarchive write_message(f);
          write_message << P;
          f.close();
          vector<int> indexes;
          indexes.push_back(dest);indexes.push_back(source);indexes.push_back(tag);
          IndexComps.push_back(indexes);
      }
      //receive a packet
      void recv(int source, int tag, Packet &P){
          char tmp[512];
          int dest = mpi_rank();
          sprintf(tmp,"%s.%i.%i.%i.tmp",BaseName,dest,source,tag);
          ifstream f;
          f.open(tmp,std::ios::binary);
          boost::archive::binary_iarchive read_message(f);
          read_message >> P;
          f.close();
      }
      //send a vector of Wavefunctions
      void send(int dest, int tag, const vector<boost::shared_ptr<Wavefunction> > &P){
          char tmp[512];
          int source = mpi_rank();
          sprintf(tmp,"%s.%i.%i.%i.tmp",BaseName,dest,source,tag);
          ofstream f;
          f.open(tmp,std::ios::binary);
          boost::archive::binary_oarchive write_message(f);
          write_message << P;
          f.close();
          vector<int> indexes;
          indexes.push_back(dest);indexes.push_back(source);indexes.push_back(tag);
          IndexComps.push_back(indexes);
      }
      //receive a vector of Wavefunctions
      void recv(int source, int tag, vector <boost::shared_ptr<Wavefunction> > &P){
          char tmp[512];
          int dest = mpi_rank();
          sprintf(tmp,"%s.%i.%i.%i.tmp",BaseName,dest,source,tag);
          ifstream f;
          f.open(tmp,std::ios::binary);
          boost::archive::binary_iarchive read_message(f);
          read_message >> P;
          f.close();
      }
      //clear the buffer
      void ClearBuffer(){
          char msg[512];
          for (int i=0;i<IndexComps.size();i++){
              int dest = IndexComps[i][0];
              int source = IndexComps[i][1];
              int tag = IndexComps[i][2];
              sprintf(msg,"%s.%i.%i.%i.tmp",BaseName,dest,source,tag);
              remove(msg);
          }
      }
  };
  
  
}


#endif	/* NEVPT2_PAL_H */

