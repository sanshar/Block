#include "nevpt2_pal.h"
#include "nevpt2_mpi.h"
#include "MatrixBLAS.h"
#include <boost/serialization/serialization.hpp>
#include "nevpt2_operators.h"
#include "distribute.h"
#include "density.h"

#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

namespace SpinAdapted{
  
    //============================================================================
  // Add all contributions from different processes to a wavefunction
  //============================================================================
  void AddPalWavefunction(Wavefunction &WF){
#ifndef SERIAL
    
    //the communicator
    mpi::communicator world;
    int WorldSize=world.size();
    if (WorldSize==1) return;
    int rank = world.rank();
    
    //determine the send-recv structure
    std::vector<int> sendlist;
    makesendlist(sendlist);

    //receive and accumulate the data
    for (int isender=0;isender<sendlist.size();isender++){
      boost::shared_ptr<Wavefunction> tmp_recv (new Wavefunction());
      //receive
      world.recv(sendlist[isender],0,*tmp_recv);
      //accumulate
      ScaleAdd(1.0,*tmp_recv,WF);
    }//isender
      
    //send the data
    if (rank!=0) world.send(receivefrom(),0,WF);

    //broadcast the accumulated result
    boost::mpi::broadcast(world,WF,0);
#endif
  }

  //============================================================================
  // Add all contributions from different processes to a wavefunction
  //============================================================================
  void AddPalDensity(DensityMatrix &D){
#ifndef SERIAL
    
    //the communicator
    mpi::communicator world;
    int WorldSize=world.size();
    if (WorldSize==1) return;
    int rank = world.rank();
    
    //determine the send-recv structure
    std::vector<int> sendlist;
    makesendlist(sendlist);

    //receive and accumulate the data
    for (int isender=0;isender<sendlist.size();isender++){
      boost::shared_ptr<DensityMatrix> tmp_recv (new DensityMatrix());
      //DensityMatrix part;
      //receive
      world.recv(sendlist[isender],0,*tmp_recv);
      //accumulate
      ScaleAdd(1.0,*tmp_recv,D);
      //clear memory
      //part.Clear();
    }//isender

    //send the data
    if (rank!=0) world.send(receivefrom(),0,D);

    //broadcast the accumulated result
    boost::mpi::broadcast(world,D,0);
#endif
  }

  
  //============================================================================
  // Add all contributions from different processes to a set of wavefunctions
  //============================================================================
  void AddPalWavefunctions(vector<boost::shared_ptr<Wavefunction> > &VTrip, 
                           vector<boost::shared_ptr<Wavefunction> > &SigmaTrip, 
                           boost::shared_ptr<Wavefunction> &VSing,
                           boost::shared_ptr<Wavefunction> &SigmaSing){
#ifndef SERIAL
    int VTripSize = VTrip.size();
    int STripSize = SigmaTrip.size();
    
    //the send and recv packets
    vector<boost::shared_ptr<Wavefunction> >tmp_recv;
    vector<boost::shared_ptr<Wavefunction> >tmp_send;
    
    //the communicator
    mpi::communicator world;
    int WorldSize=world.size();
    if (WorldSize==1) return;
    int rank = world.rank();

    //determine the send-recv structure
    std::vector<int> sendlist;
    makesendlist(sendlist);

    //receive and accumulate the data
    for (int isender=0;isender<sendlist.size();isender++){
      //receive
      world.recv(sendlist[isender],0,tmp_recv);
      //accumulate
      ScaleAdd(1.0,*(tmp_recv[0]),*VSing);
      ScaleAdd(1.0,*(tmp_recv[1]),*SigmaSing);
      for (int i=0;i<VTripSize;i++){
        ScaleAdd(1.0,*(tmp_recv[i+2]),*(VTrip[i]));
      }
      for (int i=0;i<STripSize;i++){
        ScaleAdd(1.0,*(tmp_recv[i+VTripSize+2]),*(SigmaTrip[i]));
      }
    }//isender
      
    //generate a packet to be send
    tmp_send.push_back(boost::make_shared<Wavefunction>(*VSing));
    tmp_send.push_back(boost::make_shared<Wavefunction> (*SigmaSing));
    for (int i=0;i<VTrip.size();i++){
      tmp_send.push_back(boost::make_shared<Wavefunction> (*VTrip[i]));
    }
    for (int i=0;i<SigmaTrip.size();i++){
      tmp_send.push_back(boost::make_shared<Wavefunction> (*SigmaTrip[i]));
    }

    //send the data
    if (rank!=0) world.send(receivefrom(),0,tmp_send);

    //broadcast the accumulated result
    boost::mpi::broadcast(world,tmp_send,0);

    //clear the memory
    tmp_recv.clear();
    tmp_send.clear();
    
#endif
  }
  
  //============================================================================
  // Add all contributions from different processes to a set of wavefunctions
  // In this version we have an entire batch of sets of wavefunctions
  //============================================================================
  void AddPalWavefunctions(vector<vector<boost::shared_ptr<Wavefunction> > >&VTrip, 
                           vector<vector<boost::shared_ptr<Wavefunction> > >&SigmaTrip, 
                           vector<boost::shared_ptr<Wavefunction> >&VSing,
                           vector<boost::shared_ptr<Wavefunction> >&SigmaSing){
#ifndef SERIAL
    if (VSing.size()>0){
      int batchsize = VTrip.size(); 
      int VTripSize = VTrip[0].size();
      int STripSize = SigmaTrip[0].size();

      //the send and recv packets
      vector<vector<boost::shared_ptr<Wavefunction> > >tmp_recv;
      vector<vector<boost::shared_ptr<Wavefunction> > >tmp_send;

      //the communicator
      mpi::communicator world;
      int WorldSize=world.size();
      if (WorldSize==1) return;
      int rank = world.rank();
      
      //determine the send-recv structure
      std::vector<int> sendlist;
      makesendlist(sendlist);

      //receive and accumulate the data
      for (int isender=0;isender<sendlist.size();isender++){
        //receive
        world.recv(sendlist[isender],0,tmp_recv);
        //accumulate
        for (int b=0;b<batchsize;b++){
          ScaleAdd(1.0,*(tmp_recv[b][0]),*(VSing[b]));
          ScaleAdd(1.0,*(tmp_recv[b][1]),*(SigmaSing[b]));
          for (int i=0;i<VTripSize;i++){
            ScaleAdd(1.0,*(tmp_recv[b][i+2]),*(VTrip[b][i]));
          }
          for (int i=0;i<STripSize;i++){
            ScaleAdd(1.0,*(tmp_recv[b][i+VTripSize+2]),*(SigmaTrip[b][i]));
          }
        }//b

      }//isender

      //generate a packet to be send
      tmp_send.resize(batchsize);
      for (int b=0;b<batchsize;b++){
        tmp_send[b].push_back(boost::make_shared<Wavefunction>(*VSing[b]));
        tmp_send[b].push_back(boost::make_shared<Wavefunction> (*SigmaSing[b]));
        for (int i=0;i<VTrip[b].size();i++){
          tmp_send[b].push_back(boost::make_shared<Wavefunction> (*VTrip[b][i]));
        }
        for (int i=0;i<SigmaTrip[b].size();i++){
          tmp_send[b].push_back(boost::make_shared<Wavefunction> (*SigmaTrip[b][i]));
        }
      }//b
      
      //send the data
      if (rank!=0) world.send(receivefrom(),0,tmp_send);
      //mpi_barrier();  
      
      //broadcast the accumulated result
      boost::mpi::broadcast(world,tmp_send,0);
      
      //clear the memory
      tmp_recv.clear();
      tmp_send.clear();
    }//VSing.size()>0
#endif
  }

  

  
  //============================================================================
  // Send a local function V(tuv) and integral matrix around in a circle to use 
  // them in each process to evaluate local V(a) functions
  //============================================================================
  void SendAroundVtuv_a(Wavefunction &Vtuv, Matrix &Ktv, WavefunctionArray &Ta, int u, bool Dummy){
#ifndef SERIAL
    //copy the original function and matrix to the receive buffer
    Packet tmp_recv(Vtuv,Ktv,u,Dummy);
      
    //the communicator
    mpi::communicator world;
    int WorldSize=world.size();
    if (WorldSize==1) return;
    int rank = world.rank();
    int next = (rank+1)%WorldSize;
    int prev = (rank+WorldSize-1)%WorldSize;
    
    //loop over steps in the cycle
    for (int istep=0;istep<WorldSize-1;istep++){
      //copy the recv buffer into the send buffer
      Packet tmp_send(tmp_recv);
      
      //clear the recv buffer
      tmp_recv.clear();
      
      //send and receive the data
      mpi::request reqs[2];
      reqs[0] = world.isend(next,rank,tmp_send);
      reqs[1] = world.irecv(prev,prev,tmp_recv);
      mpi::wait_all(reqs,reqs+2);
      
      //if the received data are not dummy data, store them
      if (!tmp_recv.dummy()){
        Ta.ResetBuffer();
        bool EndOfInnerArray = false;
        int a,dummy;
        int u_=tmp_recv.get_t();
        while (!EndOfInnerArray){
          boost::shared_ptr<Wavefunction> Va = Ta.GetOpFromBuffer(dummy,a,EndOfInnerArray);
          if (!EndOfInnerArray){
            ScaleAdd(tmp_recv.GetM(0)->element(u_,a),*(tmp_recv.GetWF(0)),*Va);
            Ta.ReplaceOpInBuffer(Va,dummy,a);
          }
        }//a
      }//Dummy
      tmp_send.clear();
      world.barrier();
    }//istep
    
    //clear the memory
    tmp_recv.clear();
#endif
  }

  
  //============================================================================
  // Send a local function V(tuv) and integral matrix around in a circle to use 
  // them in each process to evaluate local V(a) functions
  //============================================================================
  void SendAroundVtuv_a(Wavefunction &Vtuv, Wavefunction &Vtvu, Matrix &Kuv, Matrix &Kvu,
                        WavefunctionArray &Ta, int t, bool Dummy){
#ifndef SERIAL
    //copy the original function and matrix to the receive buffer
    Packet tmp_recv(Vtuv,Vtvu,Kuv,Kvu,t,-1,Dummy);
    
    //the communicator
    mpi::communicator world;
    int WorldSize=world.size();
    if (WorldSize==1) return;
    int rank = world.rank();
    int next = (rank+1)%WorldSize;
    int prev = (rank+WorldSize-1)%WorldSize;
    
    //loop over steps in the cycle
    for (int istep=0;istep<WorldSize-1;istep++){
      //copy the recv buffer into the send buffer
      Packet tmp_send(tmp_recv);

      //clear the recv buffer
      tmp_recv.clear();
      
      //send and receive the data
      mpi::request reqs[2];
      reqs[0] = world.isend(next,rank,tmp_send);
      reqs[1] = world.irecv(prev,prev,tmp_recv);
      mpi::wait_all(reqs,reqs+2);

      //if the received data are not dummy data, store them
      if (!tmp_recv.dummy()){
        Ta.ResetBuffer();
        bool EndOfInnerArray = false;
        int a,dummy;
        int t_ = tmp_recv.get_t();
        while (!EndOfInnerArray){
          boost::shared_ptr<Wavefunction> Va = Ta.GetOpFromBuffer(dummy,a,EndOfInnerArray);
          if (!EndOfInnerArray){
            ScaleAdd(tmp_recv.GetM(0)->element(t_,a),*(tmp_recv.GetWF(0)),*Va);
            ScaleAdd(tmp_recv.GetM(1)->element(t_,a),*(tmp_recv.GetWF(1)),*Va);
            Ta.ReplaceOpInBuffer(Va,dummy,a);
          }
        }//a
      }//!Dummy
      tmp_send.clear();
      world.barrier();
    }//istep
    
    //clear the memory
    tmp_recv.clear();
#endif
  }

  //============================================================================
  // Send a local function V(tuv) and integral matrix around in a circle to use 
  // them in each process to evaluate local V(a) functions
  //============================================================================
  void SendAroundVtuv_a(Wavefunction &Vtvu, Wavefunction &Vvtu, Wavefunction &Vtuv, 
                        Matrix &Ktv, Matrix &Kvt, Matrix &Ktu, 
                        WavefunctionArray &Ta, int u, int v, bool Dummy){
#ifndef SERIAL
    //copy the original function and matrix to the receive buffer
    Packet tmp_recv(Vtvu,Vvtu,Vtuv,Ktv,Kvt,Ktu,-1,u,v,Dummy);

    //the communicator
    mpi::communicator world;
    int WorldSize=world.size();
    if (WorldSize==1) return;
    int rank = world.rank();
    int next = (rank+1)%WorldSize;
    int prev = (rank+WorldSize-1)%WorldSize;
    
    //loop over steps in the cycle
    for (int istep=0;istep<WorldSize-1;istep++){
      //copy the recv buffer into the send buffer
      Packet tmp_send(tmp_recv);
      
      //clear the recv buffer
      tmp_recv.clear();
      
      //send and receive the data
      mpi::request reqs[2];
      reqs[0] = world.isend(next,rank,tmp_send);
      reqs[1] = world.irecv(prev,prev,tmp_recv);
      mpi::wait_all(reqs,reqs+2);
      
      //if the received data are not dummy data, store them
      Ta.ResetBuffer();
      bool EndOfInnerArray = false;
      int a,dummy;
      int u_ = tmp_recv.get_u();
      int v_ = tmp_recv.get_v();
      if (!tmp_recv.dummy()){
        while (!EndOfInnerArray){
          boost::shared_ptr<Wavefunction> Va = Ta.GetOpFromBuffer(dummy,a,EndOfInnerArray);
          if (!EndOfInnerArray){
            ScaleAdd(tmp_recv.GetM(0)->element(u_,a),*tmp_recv.GetWF(0),*Va);
            ScaleAdd(tmp_recv.GetM(1)->element(u_,a),*tmp_recv.GetWF(1),*Va);
            ScaleAdd(tmp_recv.GetM(2)->element(v_,a),*tmp_recv.GetWF(2),*Va);
            Ta.ReplaceOpInBuffer(Va,dummy,a);
          }
        }//a
      }//!Dummy
      tmp_send.clear();
      world.barrier();
    }//istep
    
    //clear the memory
    tmp_recv.clear();
    
    
#endif
  }

  //============================================================================
  // Send a local function V(tuv) and integral matrix around in a circle to use 
  // them in each process to evaluate local V(i) functions
  //============================================================================
  void SendAroundVtuv_i(Wavefunction &Vtuv, Matrix &Kut, 
                        WavefunctionArray &Ti, int v, bool Dummy){
#ifndef SERIAL
    //copy the original function and matrix to the receive buffer
    Packet tmp_recv(Vtuv,Kut,v,Dummy);
      
    //the communicator
    mpi::communicator world;
    int WorldSize=world.size();
    if (WorldSize==1) return;
    int rank = world.rank();
    int next = (rank+1)%WorldSize;
    int prev = (rank+WorldSize-1)%WorldSize;
    
    //loop over steps in the cycle
    for (int istep=0;istep<WorldSize-1;istep++){
      //copy the recv buffer into the send buffer
      Packet tmp_send(tmp_recv);
      
      //clear the recv buffer
      tmp_recv.clear();
      
      //send and receive the data
      mpi::request reqs[2];
      reqs[0] = world.isend(next,rank,tmp_send);
      reqs[1] = world.irecv(prev,prev,tmp_recv);
      mpi::wait_all(reqs,reqs+2);
      
      //if the received data are not dummy data, store them
      if (!tmp_recv.dummy()){
        Ti.ResetBuffer();
        bool EndOfInnerArray = false;
        int i,dummy;
        int v_=tmp_recv.get_t();
        while (!EndOfInnerArray){
          boost::shared_ptr<Wavefunction> Vi = Ti.GetOpFromBuffer(dummy,i,EndOfInnerArray);
          if (!EndOfInnerArray){
            ScaleAdd(tmp_recv.GetM(0)->element(v_,i),*(tmp_recv.GetWF(0)),*Vi);
            Ti.ReplaceOpInBuffer(Vi,dummy,i);
          }
        }//a
      }//Dummy
      tmp_send.clear();
      world.barrier();
    }//istep
    
    //clear the memory
    tmp_recv.clear();
#endif
  }

  
  //============================================================================
  // Send a local function V(tuv) and integral matrix around in a circle to use 
  // them in each process to evaluate local V(i) functions
  //============================================================================
  void SendAroundVtuv_i(Wavefunction &Vtuv, Wavefunction &Vutv, Matrix &Ktu, 
                        WavefunctionArray &Ti, int v, bool Dummy){
#ifndef SERIAL
    //copy the original function and matrix to the receive buffer
    Packet tmp_recv(Vtuv,Vutv,Ktu,v,Dummy);
      
    //the communicator
    mpi::communicator world;
    int WorldSize=world.size();
    if (WorldSize==1) return;
    int rank = world.rank();
    int next = (rank+1)%WorldSize;
    int prev = (rank+WorldSize-1)%WorldSize;
    
    //loop over steps in the cycle
    for (int istep=0;istep<WorldSize-1;istep++){
      //copy the recv buffer into the send buffer
      Packet tmp_send(tmp_recv);
      
      //clear the recv buffer
      tmp_recv.clear();
      
      //send and receive the data
      mpi::request reqs[2];
      reqs[0] = world.isend(next,rank,tmp_send);
      reqs[1] = world.irecv(prev,prev,tmp_recv);
      mpi::wait_all(reqs,reqs+2);
      
      //if the received data are not dummy data, store them
      if (!tmp_recv.dummy()){
        Ti.ResetBuffer();
        bool EndOfInnerArray = false;
        int i,dummy;
        int v_=tmp_recv.get_v();
        while (!EndOfInnerArray){
          boost::shared_ptr<Wavefunction> Vi = Ti.GetOpFromBuffer(dummy,i,EndOfInnerArray);
          if (!EndOfInnerArray){
            ScaleAdd(tmp_recv.GetM(0)->element(i,v_),*(tmp_recv.GetWF(0)),*Vi);
            ScaleAdd(tmp_recv.GetM(0)->element(v_,i),*(tmp_recv.GetWF(1)),*Vi);
            Ti.ReplaceOpInBuffer(Vi,dummy,i);
          }
        }//a
      }//Dummy
      tmp_send.clear();
      world.barrier();
    }//istep
    
    //clear the memory
    tmp_recv.clear();
#endif
  }


  //============================================================================
  // Send a local function V(tuv) and integral matrix around in a circle to use 
  // them in each process to evaluate local V(i) functions
  //============================================================================
  void SendAroundVtuv_i(Wavefunction &Vtvu, Wavefunction &VtvuS, Wavefunction &Vvtu, 
                        Wavefunction &VvtuS, Matrix &Ktv, Matrix &Kut, 
                        WavefunctionArray &Ti, int u, int v, bool Dummy, double SingFac){
#ifndef SERIAL
    //copy the original function and matrix to the receive buffer
    Packet tmp_recv(Vtvu,VtvuS,Vvtu,VvtuS,Ktv,Kut,u,v,Dummy);
      
    //the communicator
    mpi::communicator world;
    int WorldSize=world.size();
    if (WorldSize==1) return;
    int rank = world.rank();
    int next = (rank+1)%WorldSize;
    int prev = (rank+WorldSize-1)%WorldSize;
    
    //loop over steps in the cycle
    for (int istep=0;istep<WorldSize-1;istep++){
      //copy the recv buffer into the send buffer
      Packet tmp_send(tmp_recv);
      
      //clear the recv buffer
      tmp_recv.clear();
      
      //send and receive the data
      mpi::request reqs[2];
      reqs[0] = world.isend(next,rank,tmp_send);
      reqs[1] = world.irecv(prev,prev,tmp_recv);
      mpi::wait_all(reqs,reqs+2);
      
      //if the received data are not dummy data, store them
      if (!tmp_recv.dummy()){
        Ti.ResetBuffer();
        bool EndOfInnerArray = false;
        int i,dummy;
        int u_=tmp_recv.get_u();
        int v_=tmp_recv.get_v();
        while (!EndOfInnerArray){
          boost::shared_ptr<Wavefunction> Vi = Ti.GetOpFromBuffer(dummy,i,EndOfInnerArray);
          if (!EndOfInnerArray){
            ScaleAdd(tmp_recv.GetM(0)->element(i,u_),*(tmp_recv.GetWF(0)),*Vi);
            ScaleAdd(tmp_recv.GetM(1)->element(i,v_)*SingFac,*(tmp_recv.GetWF(1)),*Vi);
            ScaleAdd(tmp_recv.GetM(0)->element(u_,i),*(tmp_recv.GetWF(2)),*Vi);
            ScaleAdd(tmp_recv.GetM(1)->element(i,v_)*SingFac,*(tmp_recv.GetWF(3)),*Vi);
            Ti.ReplaceOpInBuffer(Vi,dummy,i);
          }
        }//a
      }//Dummy
      tmp_send.clear();
      world.barrier();
    }//istep
    
    //clear the memory
    tmp_recv.clear();
#endif
  }
  
}


