
#include "nevpt2_operators.h"
#include "spinblock.h"
#include "wavefunction.h"
#include "guess_wavefunction.h"
#include "BaseOperator.h"
#include <vector>
#include <newmat.h>
#include "nevpt2_mpi.h"
#include "nevpt2_util.h"
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

namespace SpinAdapted {
  //=================================
  // Set the StateInfo
  //=================================
  template <class T>
  void OperatorArray<T>::SetStateInfo(const StateInfo &SI){
    char file[512];
    sprintf(file,"%s.StateInfo.tmp",FileName);
    if (mpi_rank()==0){
      remove(file);
      std::ofstream outfs(file, std::ios::binary);
      boost::archive::binary_oarchive save_info(outfs);
      save_info  << SI << *SI.leftStateInfo << *(SI.leftStateInfo->leftStateInfo);
      save_info << *(SI.leftStateInfo->rightStateInfo) << *SI.rightStateInfo;
      outfs.close();
    }
    HaveStateInfo = true;
  }
  
  //=================================
  // Get the current StateInfo
  //=================================
  template <class T>
  StateInfo OperatorArray<T>::GetStateInfo()const{
    StateInfo SI;
    SI.Allocate();
    char file[512];
    if (HaveStateInfo){
      if (mpi_rank()==0){
        sprintf(file,"%s.StateInfo.tmp",FileName);
        std::ifstream inpfs(file, std::ios::binary);
        boost::archive::binary_iarchive load_info(inpfs);
        load_info >> SI >> *SI.leftStateInfo >> *(SI.leftStateInfo->leftStateInfo)
                  >> *(SI.leftStateInfo->rightStateInfo) >> *SI.rightStateInfo;
        inpfs.close();
      }
    }
#ifndef SERIAL
      mpi::communicator world;
      mpi::broadcast(world,SI,0);
#endif
    return SI;
  }
  
  
  //=================================
  // open the file to read data
  //=================================
  template <class T>
  void OperatorArray<T>::OpenFileRead(int fn){
    if (!InCore){
      char msg[512];
      sprintf(msg,"%s%i",FileName,fn);
      ifs.open(msg,std::ios::binary);
      ActFile = fn;
    }
  }
  
  //=================================
  // open the file to write data
  //=================================
  template <class T>
  void OperatorArray<T>::OpenFileWrite(){
    if (!InCore){
      int fn;
      fn = mpi_rank();
      char msg[512];
      sprintf(msg,"%s%i",FileName,fn);
      ofs.open(msg,std::ios::binary|std::ios::ate|std::ios::app);
      ActFile = fn;
    }
  }
  
  //=================================
  // Close the file
  //=================================
  template <class T>
  void OperatorArray<T>::CloseFileRead(){
    ifs.close();
  }
  template <class T>
  void OperatorArray<T>::CloseFileWrite(){
    ofs.close();
  }
  
  //=================================
  // Set the FileName
  //=================================
  template <class T>
  void OperatorArray<T>::SetFileName(const char *name){
    sprintf(FileName,"%s",name);
  }
  
  //=================================
  // Rename the file
  //=================================
  template <class T>
  void OperatorArray<T>::ReNameFile(const char* NewName){
    char source[512];
    char target[512];
    //the data files (every process is responsible for renaming its own file)
    sprintf(source,"%s%i",FileName,mpi_rank());
    sprintf(target,"%s%i",NewName,mpi_rank());
    rename(source,target);
    sprintf(FileName,"%s",NewName);
    //the StateInfo files
    if (mpi_rank()==0){
      sprintf(source,"%s.StateInfo.tmp",FileName);
      sprintf(target,"%s.StateInfo.tmp",NewName);
      rename(source,target);
      sprintf(FileName,"%s",NewName);
    }//irank
  }
  
  //======================
  //store all data on disk
  //======================
  template <class T>
  void OperatorArray<T>::Store(bool sync){
    int nprocs = mpi_world_size();
    int iproc = mpi_rank();
    char FName[512];
    ofstream f;
    char msg[512];

    //------------------------
    //the operators themselves
    //------------------------
    //open the file 
    ofstream data;
    if (InCore){
      sprintf(FName,"%s%i",FileName,iproc);
      data.open(FName,std::ios::binary);
      for (int ielement=0;ielement<Operators.size();ielement++){
        //get the position
        long  pos = (long) data.tellp();
        //store the address
        int i = Operators[ielement].Geti();
        int j = Operators[ielement].Getj();
        if (IJAddr(i,j).first!= iproc){
          sprintf(msg,"\nERROR: Process %i stores an operator that was originally owned by process %i",iproc,IJAddr(i,j).first);
          pout << msg;
        }
        //store the element on disk
        Operators[ielement].write(data);
      }//elements
      //close the file
      data.close();
    }//InCore

    //-----------------------
    //store the overhead data
    //-----------------------
    //open the file
    sprintf(FName,"%s.%i.header.tmp",FileName,iproc);
    f.open(FName,std::ios::binary);
    //the address matrix
    //if requested, synchronize the address matrix. This is usually recommended
    if (sync) Synchronize();
    int nrows = IJAddr.Nrows();
    int ncols = IJAddr.Ncols();
    f.write(reinterpret_cast<char*>(&nrows), sizeof(int));
    f.write(reinterpret_cast<char*>(&ncols), sizeof(int));
    for (int i=0;i<nrows;i++){
      for (int j=0;j<ncols;j++){
        int  proc_addr = IJAddr(i,j).first;
        long file_addr = IJAddr(i,j).second;
        f.write(reinterpret_cast<char*>(&proc_addr), sizeof(int));
        f.write(reinterpret_cast<char*>(&file_addr), sizeof(long));
      }//j
    }//i

    f.write(reinterpret_cast<char*>(&InCore), sizeof(bool));//the Incore flag
    f.write(reinterpret_cast<char*>(&length), sizeof(int));//the number of Operators stored
    f.write(reinterpret_cast<char*>(&FileName), 512);//the data filename 
    f.write(reinterpret_cast<char*>(&NumFiles), sizeof(int));//the number of files
    f.write(reinterpret_cast<char*>(&BufferSize), sizeof(int));//the size of the buffer
    
    //the Spin Quantum of the stored operators and the Basis operator
    boost::archive::binary_oarchive save_data(f);
    save_data << Q;
    
    //close the header file
    f.close();
    
  }
  
  void ThreeIndOpArray::Store(bool sync){
    //----------------------
    //call the parent method
    //----------------------
    OpArray::Store(sync);
    //---------------------------------
    //Append the Three Index Dimensions
    //---------------------------------
    //get the FileName
    char fn[512];
    OpArray::GetFileName(fn);
    //open the file and go to the end of it
    ofstream f;
    char FName[512];
    int iproc = mpi_rank();
    sprintf(FName,"%s.%i.header.tmp",fn,iproc);
    f.open(FName,std::ios::binary | std::ios::app);
    
    //actually append the Dimensions
    f.write(reinterpret_cast<char*>(&Dimensions[0]), sizeof(int));
    f.write(reinterpret_cast<char*>(&Dimensions[1]), sizeof(int));
    f.write(reinterpret_cast<char*>(&Dimensions[2]), sizeof(int));
    //close the file
    f.close();
  }
  
  //==========================
  //Receive all data from disk
  //==========================
  template <class T>
  void OperatorArray<T>::Retrieve(const char *name, bool clear){
    char FName[512];
    int iproc = mpi_rank();
    int nproc = mpi_world_size();
    char msg[512];
    sprintf(FName,"%s.%i.header.tmp",name,iproc);
    //--------------------
    //read the header file
    //--------------------
    ifstream input;
    input.open(FName,std::ios::binary);
    //the dimension of the address matrix
    int nrows,ncols;
    input.read((char*)&nrows,sizeof(int));
    input.read((char*)&ncols,sizeof(int));
    //resize the address matrix
    IJAddr.ReSize(nrows,ncols);
    //read the addresses
    for (int i=0;i<nrows;i++){
      for (int j=0;j<ncols;j++){
        int  proc_addr=0;
        long file_addr=0;
        input.read((char*)&proc_addr,sizeof(int));
        input.read((char*)&file_addr,sizeof(long));
        IJAddr(i,j).first  = proc_addr;
        IJAddr(i,j).second = file_addr;
      }//j
    }//i
    
    input.read(reinterpret_cast<char*>(&InCore), sizeof(bool));//the Incore flag
    input.read(reinterpret_cast<char*>(&length), sizeof(int));//the number of Operators stored
    input.read(reinterpret_cast<char*>(&FileName), 512);//the data filename 
    input.read(reinterpret_cast<char*>(&NumFiles), sizeof(int));//the number of files
    input.read(reinterpret_cast<char*>(&BufferSize), sizeof(int));//the number of files

    //the Spin Quantum of the stored operators and the Basis operator
    boost::archive::binary_iarchive load_data(input);
    load_data >> Q;
    
    //close the file
    input.close();
    
    //if required, clear the data
    if (clear) remove(FName);
    
    //-------------
    //read the data
    //-------------
    bool FileAtEnd=false;
    if (InCore){
      //open the file
      sprintf(FName,"%s%i",FileName,iproc);
      input.open(FName,std::ios::binary);
      while(!FileAtEnd){
        //read an element from disk
        BufferElement<T> tmp;
        tmp.read(input);
        int i = tmp.Geti();
        int j = tmp.Getj();
        bool NeedTensorTrace = tmp.NeedTensorTrace();
        //add it to the array
        AppendOperator(tmp.GetElement(),i,j,NeedTensorTrace,false);
        //check whether the file is at the end
        input.peek();
        FileAtEnd = input.eof(); 
      }
      //close the file 
      input.close();
      //remove the file
      remove(FName);
    }//InCore
  }
  void ThreeIndOpArray::Retrieve(const char *name){
    //----------------------
    //call the parent method
    //----------------------
    OpArray::Retrieve(name, false);
    //-------------------------------
    //Read the Three Index Dimensions
    //-------------------------------
    char FName[512];
    char fn[512];
    int iproc = mpi_rank();
    //get the filename
    OpArray::GetFileName(fn);
    //open the file
    sprintf(FName,"%s.%i.header.tmp",name,iproc);
    ifstream input;
    input.open(FName,std::ios::binary);
    //go to the end of the file and then three steps back
    input.seekg (-3*sizeof(int), input.end);
    //now actually read the three index dimensions
    input.read(reinterpret_cast<char*>(&Dimensions[0]), sizeof(int));//the number of files
    input.read(reinterpret_cast<char*>(&Dimensions[1]), sizeof(int));//the number of files
    input.read(reinterpret_cast<char*>(&Dimensions[2]), sizeof(int));//the number of files
    //clear the disk space
    remove(FName);
  }
  
  //================
  // get an operator 
  //================
  template <class T>
  boost::shared_ptr<T> OperatorArray<T>::GetOperator(int i, int j){
    if (InCore){
      long PairIndex;
      //check if operator is available to the current process
      boost::shared_ptr<T> OP;
      if (mpi_rank()==IJAddr(i,j).first){
        PairIndex = IJAddr(i, j).second;
        OP=Operators[PairIndex].GetElement();
      }
#ifndef SERIAL
      mpi::communicator world;
      mpi::broadcast(world,*OP,IJAddr(i,j).first);
#endif
      return OP;
    }
    else{
      //get the address of operator(i,j)
      long address = IJAddr(i,j).second;
      if (address==-1){
        char msg[512];
        sprintf(msg,"ERROR: Operator could not be found!");
        boost::shared_ptr<T> null;
        null.reset();
        return null;
      }
      //check if the correct file is open
      char msg[512];
      int ifile = IJAddr(i,j).first;
      if (ifile!=ActFile){
        CloseFileRead();
        OpenFileRead(ifile);
      } 
      //go to the address
      ifs.seekg(address,ios::beg);
      //read the BufferElement from disk
      BufferElement<T> tmp;
      tmp.read(ifs);
      return tmp.GetElement();
    }
  };
  
  //Note: This routine opens the file and reads the matrix only if it the calling
  //process has also stored it.
  template <class T>
  boost::shared_ptr<T> OperatorArray<T>::GetOpPal(int i, int j){
    if (InCore){
      long PairIndex;
      //check if operator is available to the current process
      boost::shared_ptr<T> OP (new T);
      if (mpi_rank()==IJAddr(i,j).first){
        PairIndex = IJAddr(i, j).second;
        OP=Operators[PairIndex].GetElement();
      }
      if (IJAddr(i,j).second==-1){
        char msg[512];
        sprintf(msg,"ERROR: Operator could not be found!");
        boost::shared_ptr<T> null;
        null.reset();
        return null;
      }
#ifndef SERIAL
      mpi::communicator world;
      mpi::broadcast(world,*OP,IJAddr(i,j).first);
#endif
      return OP;
    }
    else{
      boost::shared_ptr<T> Wave (new T);
      //get the address of the operator(i,j)
      long address = IJAddr(i,j).second;
      if (address==-1){
        char msg[512];
        sprintf(msg,"ERROR: Operator could not be found!");
        boost::shared_ptr<T> null;
        null.reset();
        return null;
      }
      //check, who stored the requested operator
      int ifile = IJAddr(i,j).first;
#ifndef SERIAL
      mpi::communicator world;
      //only the process that stored the operator reads it
      if (ifile==world.rank()){
        //check if the correct file is open
        if (ActFile!=ifile){
          CloseFileRead();
          OpenFileRead(ifile);
        }
        //go to the address
        ifs.seekg(address,ios::beg);
        //read the BufferElement from disk
        BufferElement<T> tmp;
        tmp.read(ifs);
        return tmp.GetElement();
      }
      //send the operator to the other processes
      mpi::broadcast(world,*Wave,ifile);
      return Wave;
#else
      if (ifile!=ActFile){
        CloseFileRead();
        OpenFileRead(ifile);
      } 
      //go to the address
      ifs.seekg(address,ios::beg);
      //read the BufferElement from disk
      BufferElement<T> tmp;
      tmp.read(ifs);
      return tmp.GetElement();
#endif    
    }
  }
  
  //get one operator (only returns operator if it is available to this process, 
  //otherwise returns a zero pointer)
  template <class T>
  boost::shared_ptr<T> OperatorArray<T>::GetLocalOp(int i, int j){
    if (InCore){
      long PairIndex;
      //check if operator is available to the current process
      if (mpi_rank()==IJAddr(i,j).first){
        PairIndex = IJAddr(i, j).second;
        return Operators[PairIndex].GetElement();
      }
      else{
        boost::shared_ptr<T> null;
        null.reset();
        return null;
      }
    }
    else{
      boost::shared_ptr<T> Wave (new T);
      //get the address of the operator(i,j)
      long address = IJAddr(i,j).second;
      if (address==-1){
        char msg[512];
        sprintf(msg,"ERROR: Operator could not be found!");
        boost::shared_ptr<T> null;
        null.reset();
        return null;
      }
      //check, who stored the requested operator
      int ifile = IJAddr(i,j).first;
#ifndef SERIAL
      mpi::communicator world;
      //only the process that stored the operator reads it
      if (ifile==world.rank()){
        //check if the correct file is open
        if (ActFile!=ifile){
          CloseFileRead();
          OpenFileRead(ifile);
        }
        //go to the address
        ifs.seekg(address,ios::beg);
        //read the BufferElement from disk
        BufferElement<T> tmp;
        tmp.read(ifs);
        return tmp.GetElement();
      }
      else{
        boost::shared_ptr<T> null;
        null.reset();
        return null;
      }  
#else
      if (ifile!=ActFile){
        CloseFileRead();
        OpenFileRead(ifile);
      } 
      //go to the address
      ifs.seekg(address,ios::beg);
      //read the BufferElement from disk
      BufferElement<T> tmp;
      tmp.read(ifs);
      return tmp.GetElement();
#endif    
    }
  }
  
  
  //=======================================================
  //append one operator to the list and store its pairindex
  //=======================================================
  template <class T>
  void OperatorArray<T>::AppendOperator(boost::shared_ptr <T> A, int i, int j, bool TensorTrace, bool Print) {
    if (InCore){
      BufferElement<T> A_(A,i,j,TensorTrace);
      Operators.push_back(A_);
      int PairIndex = Operators.size()-1;
      //give message if an operator is stored that is already available
      if (IJAddr(i,j).first!=-1){
        if (Print){
          char msg[512];
          sprintf(msg,"\n!!!NOTE: Storing operator T(%i,%i) although it has been stored previously!!! Length = %i",i,j,length);
          pout << msg;
        }
        length--;
      }
      IJAddr(i,j).first = mpi_rank();
      IJAddr(i,j).second = PairIndex;
      length++;
    }
    else{
      //give message if an operator is stored that is already available
      if ((IJAddr(i,j).first!=-1)&&(Print)){
        char msg[512];
        sprintf(msg,"\n!!!NOTE: Storing operator T(%i,%i) although it has been stored previously!!!",i,j);
        pout << msg;
        length--;
      }
      //go to the end of the file
      ofs.seekp(0,ios::end);
      //get the position
      long  pos = (long) ofs.tellp();
      //store the address
      IJAddr(i,j).first = ActFile;
      IJAddr(i,j).second = pos;
      //store the element on disk
      BufferElement<T> A_(A,i,j,TensorTrace);
      A_.write(ofs);
      length++;
    }
  };
  
  //================================
  // Add an entire array to this one
  //================================
  //Note: we do not care about different processes here
  //Also: if an operator with the same address has been stored, we keep the old one
  template <class T>
  void OperatorArray<T>::AppendArray(OperatorArray<T> &OpArray) {
    //-------------------
    //open the two arrays
    //-------------------
    if (!ofs.is_open()){
      OpenFileWrite();
    }
    OpArray.OpenFileRead();
    //----------------------------------------------------------------------
    //get the stored operators from the other array and add them to this one
    //----------------------------------------------------------------------
    int i,j;
    bool EndOfArray=false;
    bool NeedTensorTrace=false;
    OpArray.ResetBuffer();
    //loop over operators from other array
    while (!EndOfArray){
      boost::shared_ptr<T> tmp = OpArray.GetOpFromBuffer(i,j,NeedTensorTrace,EndOfArray);
      //if not stored previously, store the operator
      if (IJAddr(i,j).first==-1){
        AppendOperator(tmp,i,j,NeedTensorTrace);
      }
    }
    //----------------
    //close the arrays
    //----------------
    OpArray.CloseFileRead();
    CloseFileWrite();
  }
  
  //=====================================
  // duplicate the address of an operator
  //=====================================
  template <class T>
  void OperatorArray<T>::DuplicateAddress(int originali, int originalj, 
                                               int duplicatei, int duplicatej){
    IJAddr(duplicatei,duplicatej) = IJAddr(originali,originalj);
  };
  
  
  //==============================
  // resize the orbital pair space
  //==============================
  template <class T>
  void OperatorArray<T>::ResizeOrbSpace(int m, bool init){
    IJAddr.resize(m,m);
    if (init){
      int i, j;
      for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
          IJAddr(i, j).first  = -1;
          IJAddr(i, j).second = -1;
        }//cols
      }//rows
    }//initialize?
  };
  template <class T>
  void OperatorArray<T>::ResizeOrbSpace(int m, int n, bool init){
    IJAddr.resize(m,n);
    int i, j;
    if (init){
      for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
          IJAddr(i, j).first  = -1;
          IJAddr(i, j).second = -1;
        }//cols
      }//rows
    }//initialize?
  };
  
  
  //==========================
  // delete all data from disk
  //==========================
  template <class T>
  void OperatorArray<T>::Clear(){
    int i;
    for(i=0;i<Operators.size();i++){
      //set the pointer to zero. This will (hopefully) free its allocated memory
      Operators[i].clear();
    }//i
    Operators.clear();
    Buffer.clear();
    length = 0;
    ActFile = 0;
    IJAddr.Clear();
    if (Basis)Basis.reset();
    char file[512];
    if (mpi_rank()==0){
      sprintf(file,"%s.StateInfo.tmp",FileName);
      remove(file);
    }
    //delete the file
    if (!InCore){
      char msg[512];
#ifdef SERIAL      
      for (i=0;i<NumFiles;i++){
        sprintf(msg,"%s%i",FileName,i);
        remove(msg);
        //if NumFiles is uninitialized(it should not be) try to delete all files anyways
        if (NumFiles==0){
          sprintf(msg,"%s%i",FileName,0);
          remove(msg);
        }//NumFiles==0

      }//i
#else
      sprintf(msg,"%s%i",FileName,mpi_rank());
      remove (msg);
      /*
      mpi::communicator world;
      if (world.rank()==0){
        for (i=0;i<NumFiles;i++){
          sprintf(msg,"%s%i",FileName,i);
          remove(msg);
        }//i
        //if NumFiles is uninitialized(it should not be) try to delete all files anyways
        if (NumFiles==0){
          for (i=0;i<world.size();i++){
            sprintf(msg,"%s%i",FileName,i);
            remove(msg);
          }//i
        }//NumFiles==0
      }//if master
       */
#endif
    }//!InCore?
  };
  
  //=============================
  //synchronize the parallel data
  //=============================
  template <class T>
  void OperatorArray<T>::Synchronize(){
#ifndef SERIAL
    int i;
    char msg[512];
    int foo=0;
    mpi::communicator world;
    if (world.size()>1){
      //--------------------------------------------
      //receive all address matrices and gather them
      //--------------------------------------------
      if (world.rank()==0){
        ObjectMatrix<std::pair<int,long> > addr;
        //loop over all slave processes
        for (i=1;i<world.size();i++){
          //receive the address matrices
          world.recv(i, i, addr);
          //add received addresses
          for (int k=0;k<IJAddr.Nrows();k++){
            for (int l=0;l<IJAddr.Ncols();l++){
              if ((IJAddr(k,l).first!=-1)&&(addr(k,l).first!=-1)){
                sprintf(msg,"\nERROR during OperatorArray::Synchronize():Two Processes own the same operator");
                pout << msg;
              }
              if ((IJAddr(k,l).first==-1)&&(addr(k,l).first!=-1)){
                IJAddr(k,l).first = addr(k,l).first;
                IJAddr(k,l).second = addr(k,l).second;
              }
            }//l
          }//k
        }//i
      }//master
      else{
        //loop over slave processes and send the address matrices
        world.send(0,world.rank(),IJAddr);
        
      }//slaves
      //---------------------------------------------------------
      //make the global address matrix available to all processes
      //---------------------------------------------------------
      mpi::broadcast(world,IJAddr,0);
      //also broadcast the updated length of the array
      /*
      int NewLength = 0;
      for (int i=0;i<IJAddr.Nrows();i++){
        for (int j=0;j<IJAddr.Ncols();j++){
          if ((IJAddr(i,j).first!=-1)&&(IJAddr(i,j).second!=-1)) NewLength++;
        }//j
      }//i
      length = NewLength;
      mpi::broadcast(world,length,0);
      */
    }//more than 1 process?
#endif 
  };
  
  
  //============================================================================
  //Evaluate the number of elements in the buffer according to max. available 
  //memory and M
  //============================================================================
  template <class T>
  void OperatorArray<T>::CalcBufferSize(int MaxMemSize, int M){
    long ElementSize = M * M * sizeof(double);
    double ElementSize_ = (double) M * M * sizeof(double);
    double tmp = (double) MaxMemSize;
    tmp *= 1024.0 * 1024.0;
    tmp /= ElementSize_;
    BufferSize = (long) tmp;
    if (BufferSize==0) BufferSize = 1;
  } 
  
  //============================================================================
  //get one operator from the buffer (if buffer is empty or we are at the end of
  // the buffer read the next chunk from the file).
  //Note: This function is closely related to and should be used together with the 
  //ReplaceOpInBuffer function.
  //============================================================================
  template <class T>
  boost::shared_ptr<T> OperatorArray<T>::GetOpFromBuffer(int &i, int &j,bool &EndOfFile){
    bool TensorTrace = false;
    return GetOpFromBuffer(i,j,TensorTrace,EndOfFile);
  }

  //the same function as above but here we also return the NeedTensorTrace flag
  template <class T>
  boost::shared_ptr<T> OperatorArray<T>::GetOpFromBuffer(int &i, int &j,bool &TensorTrace, bool &EndOfFile){
    if (InCore){
      //if necessary, set the end-of-file flag, initialize the CurrentElement and
      //return a null
      char msg[512];
      if (CurrentElement>=Operators.size()){
        EndOfFile = true;
        CurrentElement=0;
        //return a null
        boost::shared_ptr<T> Null;
        i = -1;
        j = -1;
        return Null;
      }
      else{
        //get the element and its indices and return it
        i = Operators[CurrentElement].Geti();
        j = Operators[CurrentElement].Getj();
        TensorTrace = Operators[CurrentElement].NeedTensorTrace();
        boost::shared_ptr<T> Wave = Operators[CurrentElement].GetElement();
        CurrentElement++;
        return Wave;
      }
    }//InCore
    else{
      //check if we are at the end of the buffer
      if (CurrentElement>=Buffer.size()){
        //if yes, free the buffer and read the next chunk from the file
        Buffer.clear();
        int iElements=0;
        //check if we are at the end of the file
        EndOfFile = false;
        ifs.peek();
        BufferPosition = ifs.tellg();
        bool FileAtEnd = ((ifs.eof())||(BufferPosition==-1));
        // if the buffer and the file are at end, set the flag and return a null
        if (FileAtEnd){
          //set the end of file flag
          EndOfFile = true;
          //rewind the input stream
          ifs.clear();
          ifs.seekg(0,ios::beg);
          BufferPosition = ifs.tellg();
          //return a null
          boost::shared_ptr<T> Null;
          i = -1;
          j = -1;
          return Null;
        }
        while((iElements<BufferSize)&&(!FileAtEnd)){
          //read an element from disk
          BufferElement<T> tmp;
          tmp.read(ifs);
          //add it to the buffer
          Buffer.push_back(tmp);
          //increase the counter fro the already read elements
          iElements++;
          //check whether the file is at the end
          ifs.peek();
          FileAtEnd = ifs.eof(); 
        }
        //set the current element to the beginning
        CurrentElement = 0;
      }//read new chunk
      BufferElement<T> TMP(Buffer[CurrentElement]);
      //get the orbital indices
      i = TMP.Geti();
      j = TMP.Getj();
      //get the TensorTrace Flag
      TensorTrace = TMP.NeedTensorTrace();
      //increase the current element counter
      CurrentElement++;
      //return the operator itself
      return TMP.GetElement();
    }//!InCore
  }

  //============================================================================
  //replace one operator with another one using the buffer vector
  //============================================================================
  //Note: This function is closely related to and should be used together with the 
  //GetOpFromoBuffer function.
  //!!!NOTE: They need to have the exact same size or the whole file is corrupted!!!
  template <class T>
  void OperatorArray<T>::ReplaceOpInBuffer(boost::shared_ptr<T> A,int i,int j, bool TensorTrace){
    if (InCore){
      long PairIndex;
      PairIndex = IJAddr(i, j).second;
      //if the operator has not yet been stored, append it
      if (IJAddr(i,j).first==-1){
        AppendOperator(A,i,j,TensorTrace);
      }
      //else replace the existing operator
      else{
        BufferElement<T> A_(A,i,j,TensorTrace);
        Operators[PairIndex].clear();
        Operators[PairIndex] = A_;
      }
    }
    else{
      //replace the operator
      BufferElement<T> Element(A,i,j,TensorTrace);
      Buffer[CurrentElement-1]=Element;//note: Currentelement gets incremented before the call of this function
      //if at the end of the buffer, write to disk
      if (CurrentElement==Buffer.size()){
        //keep the position in mind
        long pos = ifs.tellg();
        CloseFileRead();
        OpenFileWrite();
        ofs.seekp(BufferPosition,ios::beg);
        //Note: we don't have to store the addresses because all operators should 
        //be stored where they were stored originally
        for (int i=0;i<Buffer.size();i++){
          Buffer[i].write(ofs);
        }//elements in buffer
        //open the file for reading again and go to the original position
        CloseFileWrite();
        OpenFileRead();
        ifs.seekg(pos,ios::beg);
      }//end of buffer
    }
  } 

  //============================================================================
  //reset the buffer and rewind the input stream
  //============================================================================
  template <class T>
  void OperatorArray<T>::ResetBuffer(){
    Buffer.clear();
    BufferPosition=0;ifs.clear();
    ifs.seekg(0,ios::beg);
    BufferPosition = ifs.tellg();
    CurrentElement=0;
  }
  
  //============================================================================
  //Add to the Basis Operator
  //============================================================================
  template <class T>
  void OperatorArray<T>::AddToBasisOperator(boost::shared_ptr<T> &op,double fac){
    ScaleAdd(fac,*op,*Basis);
  }

  
  //============================================================================
  // Transform the stored operators to a new block configuration
  // i.e.       [S.][E] -> [S'.][E']
  // with S' = [S.]
  //============================================================================
  /*void WavefunctionArray::Transform(vector<Matrix> leftRotMatrix, const vector<Matrix> &rightRotMatrix, SpinBlock &big,
                                    const StateInfo &OldStateInfo, const StateInfo &NewStateInfo, 
                                    int root, bool Transpose){
  
    int i,j;
    int dim_i,dim_j;
    dim_i = GetAddressMatrix().Nrows();
    dim_j = GetAddressMatrix().Ncols();
    char msg[512];
    //generate the new array
    sprintf(msg,"TempArray.%i.%i.tmp",root,mpi_rank());
    WavefunctionArray NewArray;
    NewArray.Initialize(dim_i,msg,mpi_world_size());
    NewArray.SetBufferSize(GetBufferSize());
    NewArray.SetInCore(GetInCore());
    
    //Inverse the left rotation matrix
    for (int q = 0; q < leftRotMatrix.size (); ++q){
      if (leftRotMatrix [q].Nrows () > 0){
  *****************************************************************************************************
  FIXME:
  Left rotation matrix should be transposed rather than inverted. Fortunately, this bug doesn't affect
  to any results, since rotation matrix has singular values which are all equal to 1, meaning that
  the pseudo inverse of a rotation matrix is indeed, the transposition of it. From the same reason,
  it's not necessary to take the pseudo inverse of right rotation matrix.
  *****************************************************************************************************

  //    try
  //    {
  //      svd(inverseLeftRotationMatrix[q], D, U, V);
  //    }
  //    catch (Exception)
  //    {
  //      pout << Exception::what() << endl;
  //      pout << D << endl;
  //      pout << U << endl;
  //      pout << V << endl;
  //      abort();
  //    }
  //    Matrix vd = V;
  //    vd *= D.i();
  //    inverseLeftRotationMatrix[q].ReSize(V.Nrows(), U.Nrows());
  //    SpinAdapted::Clear(inverseLeftRotationMatrix[q]);
  //    MatrixMultiply(vd, 'n', U, 't', inverseLeftRotationMatrix[q], 1.);
        leftRotMatrix[q] = leftRotMatrix[q].t();
      }//NRows!=0
    }//q

    //open the files
    OpenFileRead();
    NewArray.OpenFileWrite();

    //prepare the buffer
    bool EndOfArray = false;
    ResetBuffer();
    //loop over all elements in array
    while (!EndOfArray){
      boost::shared_ptr<Wavefunction> oldWF = GetOpFromBuffer(i,j,EndOfArray);
      if (!EndOfArray){
        boost::shared_ptr<Wavefunction> newWF (new Wavefunction(oldWF->get_deltaQuantum(0),&big,true));
        //sprintf(msg,"\nNorm before=%lf  projection on WF=%lf",DotProduct(*oldWF,*oldWF),DotProduct(WF,*oldWF));pout << msg;
        //pout << OldStateInfo;
        //pout << "\nnewStateInfo";
        //pout << NewStateInfo;
        double oldNorm=DotProduct(*oldWF,*oldWF);
        //Scale(1.0/sqrt(oldNorm),*oldWF);
        GuessWave::onedot_transform_wavefunction_Q(OldStateInfo,NewStateInfo,*oldWF,leftRotMatrix,rightRotMatrix,*newWF,Transpose,oldWF->get_deltaQuantum(0));
        double newNorm=DotProduct(*newWF,*newWF);
        //Scale(sqrt(oldNorm/newNorm),*newWF);
        NewArray.AppendOperator(newWF,i,j);
      }//!EndOfarray
    }//ij
    //also transform the basis operator
    boost::shared_ptr<Wavefunction> Basis = GetBasisOperator();
    if (Basis.get()!=0){
      boost::shared_ptr<Wavefunction> NewBasis (new Wavefunction(Basis->get_deltaQuantum(0),&big,true));
      GuessWave::onedot_transform_wavefunction_Q(OldStateInfo,NewStateInfo,*Basis,leftRotMatrix,rightRotMatrix,*NewBasis,Transpose,Basis->get_deltaQuantum(0));
      NewArray.SetBasisOperator(NewBasis);
    }
    //close the files
    CloseFileRead();
    NewArray.CloseFileWrite();
    NewArray.SetQuanta(GetQuanta());
    
    //replace the old array by the new one
    GetFileName(msg);
    Clear();
    (*this) = NewArray;
    ReNameFile(msg);
    NewArray.Clear();
    
  }
  */
  //============================================================================
  // Transpose the stored operators 
  //============================================================================
  void WavefunctionArray::Transpose(const StateInfo &oldStateInfo, const StateInfo &newStateInfo,SpinBlock &big, int root){

    int i,j;
    int dim_i,dim_j;
    dim_i = GetAddressMatrix().Nrows();
    dim_j = GetAddressMatrix().Ncols();
    char msg[512];
    //generate the new array
    sprintf(msg,"TempArray.%i.%i.tmp",root,mpi_rank());
    WavefunctionArray NewArray;
    NewArray.Initialize(dim_i,msg,mpi_world_size());
    NewArray.SetBufferSize(GetBufferSize());
    NewArray.SetInCore(GetInCore());
    
    //open the files
    OpenFileRead();
    NewArray.OpenFileWrite();
  
    //prepare the buffer
    bool EndOfArray = false;
    ResetBuffer();
    //loop over all elements in array
    while (!EndOfArray){
      boost::shared_ptr<Wavefunction> oldWF = GetOpFromBuffer(i,j,EndOfArray);
      if (!EndOfArray){
        boost::shared_ptr<Wavefunction> newWF (new Wavefunction(oldWF->get_deltaQuantum(0),&big,true));
        GuessWave::onedot_transpose_wavefunction(oldStateInfo,newStateInfo,*oldWF,*newWF);
        NewArray.AppendOperator(newWF,i,j);
      }//!EndOfarray
    }//ij
    //also transform the basis operator
    boost::shared_ptr<Wavefunction> Basis = GetBasisOperator();
    if (Basis.get()!=0){
      boost::shared_ptr<Wavefunction> NewBasis (new Wavefunction(Basis->get_deltaQuantum(0),&big,true));
      GuessWave::onedot_transpose_wavefunction(oldStateInfo,newStateInfo,*Basis,*NewBasis);
      NewArray.SetBasisOperator(NewBasis);
    }
    //close the files
    CloseFileRead();
    NewArray.CloseFileWrite();
    NewArray.SetQuanta(GetQuanta());
    
    //replace the old array by the new one
    GetFileName(msg);
    Clear();
    (*this) = NewArray;
    ReNameFile(msg);
    NewArray.Clear();

  }

  //============================================================================
  //Initialize the basis operator in a given representation
  //============================================================================
  void WavefunctionArray::ResetBasisOperator(SpinBlock &b,SpinQuantum &Q){
    boost::shared_ptr<Wavefunction> EmptyBasis (new Wavefunction(Q,&b,true));
    SetBasisOperator(EmptyBasis);
  }
  
  //============================================================================
  // TensorTrace the stored operators
  //============================================================================
  void OpArray::TensorTrace(SpinBlock &big, const StateInfo &NewStateInfo){
    int i,j;
    bool NeedTensorTrace=true;
    int dim_i = GetAddressMatrix().nrows();
    int dim_j = GetAddressMatrix().ncols();
    char msg[512];
    
    //the two blocks
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    SpinBlock *oldSystem  = leftBlock->get_leftBlock();
    //generate the new array
    sprintf(msg,"TempArray.%i.tmp",mpi_rank());
    OpArray NewArray;
    NewArray.Initialize(dim_i,dim_j,msg,mpi_world_size());
    NewArray.SetBufferSize(GetBufferSize());
    NewArray.SetInCore(GetInCore());
    
    //open the files
    OpenFileRead();
    NewArray.OpenFileWrite();

    //prepare the buffer
    bool EndOfArray = false;
    ResetBuffer();
    //loop over all elements in array
    while (!EndOfArray){
      boost::shared_ptr<Cre> Op = GetOpFromBuffer(i,j,NeedTensorTrace,EndOfArray);
      if (!EndOfArray){
        //if the operator was not built in the current step of the sweep, apply tensortrace
        if (NeedTensorTrace){
          //generate the new operator
          boost::shared_ptr<Cre> NewOp(new Cre);
          NewOp->set_orbs() = Op->get_orbs();
          NewOp->set_initialised() = true;
          NewOp->resize_deltaQuantum(1);
          NewOp->set_deltaQuantum(0) = Op->get_deltaQuantum(0);
          NewOp->allocate(NewStateInfo);
          operatorfunctions::TensorTrace(oldSystem,*Op,leftBlock,&NewStateInfo,*NewOp);
          //store the operator in the new Array
          NewArray.AppendOperator(Op,i,j,true);
        }//TensorTrace required
        else{
          //store the operator in the new Array
          NewArray.AppendOperator(Op,i,j);
        }//no TensorTrace required
      }//!EndOfarray
    }//ij
    //also transform the basis operator
    boost::shared_ptr<Cre> Basis = GetBasisOperator();
    if (Basis.get()!=0){
      NewArray.SetBasisOperator(Basis);
    }
    //close the files
    CloseFileRead();
    NewArray.CloseFileWrite();
    NewArray.SetQuanta(GetQuanta());
    
    //replace the old array by the new one
    GetFileName(msg);
    Clear();
    (*this) = NewArray;
    ReNameFile(msg);
    NewArray.Clear();
  }

  
  //============================================================================
  // Renormalize the stored operators 
  //============================================================================
  void OpArray::RenormalizeTransform(SpinBlock &big, const StateInfo &NewStateInfo, const vector<Matrix> &RotationMatrix){
    int i,j;
    bool NeedTensorTrace=true;
    int dim_i = GetAddressMatrix().nrows();
    int dim_j = GetAddressMatrix().ncols();
    char msg[512];
    
    //the two blocks
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    SpinBlock *oldSystem  = leftBlock->get_leftBlock();
    //generate the new array
    sprintf(msg,"TempArray.%i.tmp",mpi_rank());
    OpArray NewArray;
    NewArray.Initialize(dim_i,dim_j,msg,mpi_world_size());
    NewArray.SetBufferSize(GetBufferSize());
    NewArray.SetInCore(GetInCore());
    
    //open the files
    OpenFileRead();
    NewArray.OpenFileWrite();

    //prepare the buffer
    bool EndOfArray = false;
    ResetBuffer();
    //loop over all elements in array
    while (!EndOfArray){
      boost::shared_ptr<Cre> Op = GetOpFromBuffer(i,j,NeedTensorTrace,EndOfArray);
      if (!EndOfArray){
        //actually do the transformation
        Op->renormalise_transform(RotationMatrix,&NewStateInfo);
        //store the operator in the new Array
        NewArray.AppendOperator(Op,i,j);
      }//!EndOfarray
    }//ij
    //also transform the basis operator
    boost::shared_ptr<Cre> Basis = GetBasisOperator();
    if (Basis.get()!=0){
      Basis->renormalise_transform(RotationMatrix,&NewStateInfo);
      NewArray.SetBasisOperator(Basis);
    }
    //close the files
    CloseFileRead();
    NewArray.CloseFileWrite();
    NewArray.SetQuanta(GetQuanta());
    
    //replace the old array by the new one
    GetFileName(msg);
    Clear();
    (*this) = NewArray;
    ReNameFile(msg);
    NewArray.Clear();

  }

  //============================================================================
  // Renormalize the stored operators (three index version)
  //============================================================================
  void ThreeIndOpArray::TensorTrace(SpinBlock &big, const StateInfo &NewStateInfo){
    int i,j,k;
    bool NeedTensorTrace=true;
    int Dim[3];
    GetThreeIndexDims(Dim);
    int dim_i = Dim[0];
    int dim_j = Dim[1];
    int dim_k = Dim[2];
    char msg[512];

    //the two blocks
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    SpinBlock *oldSystem  = leftBlock->get_leftBlock();
    //generate the new array
    sprintf(msg,"TempArray.%i.tmp",mpi_rank());
    ThreeIndOpArray NewArray;
    NewArray.Initialize(dim_i,dim_j,dim_k,msg,mpi_world_size());
    NewArray.SetBufferSize(GetBufferSize());
    NewArray.SetInCore(GetInCore());
    
    //open the files
    OpenFileRead();
    NewArray.OpenFileWrite();

    //prepare the buffer
    bool EndOfArray = false;
    ResetBuffer();
    //loop over all elements in array
    while (!EndOfArray){
      boost::shared_ptr<Cre> Op = GetOpFromBuffer(i,j,k,NeedTensorTrace,EndOfArray);
      if (!EndOfArray){
        //if the operator was not built in the current step of the sweep, apply tensortrace
        if (NeedTensorTrace){
          //generate the new operator
          boost::shared_ptr<Cre> NewOp(new Cre);
          NewOp->set_orbs() = Op->get_orbs();
          NewOp->set_initialised() = true;
          NewOp->set_deltaQuantum(1,Op->get_deltaQuantum(0));
          NewOp->allocate(leftBlock->get_stateInfo());
          operatorfunctions::TensorTrace(oldSystem,*Op,leftBlock,&(leftBlock->get_stateInfo()),*NewOp);
          //store the operator in the new Array
          NewArray.AppendOperator(NewOp,i,j,k,true);
        }//TensorTrace required
        else{
          //store the operator in the new Array
          NewArray.AppendOperator(Op,i,j,k,true);
        }//no TensorTrace required
      }//!EndOfarray
    }//ij
    //also transform the basis operator
    boost::shared_ptr<Cre> Basis = GetBasisOperator();
    if (Basis.get()!=0){
      NewArray.SetBasisOperator(Basis);
    }
    //close the files
    CloseFileRead();
    NewArray.CloseFileWrite();
    NewArray.SetQuanta(GetQuanta());

    //replace the old array by the new one
    GetFileName(msg);
    Clear();
    (*this) = NewArray;
    ReNameFile(msg);
    NewArray.Clear();
    mpi_barrier();
  }
  
  //============================================================================
  // Renormalize the stored operators (three index version)
  //============================================================================
  void ThreeIndOpArray::RenormalizeTransform(SpinBlock &big, const StateInfo &NewStateInfo, const vector<Matrix> &RotationMatrix){
    int i,j,k;
    bool NeedTensorTrace=true;
    int Dim[3];
    GetThreeIndexDims(Dim);
    int dim_i = Dim[0];
    int dim_j = Dim[1];
    int dim_k = Dim[2];
    char msg[512];

    //the two blocks
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    SpinBlock *oldSystem  = leftBlock->get_leftBlock();
    //generate the new array
    sprintf(msg,"TempArray.%i.tmp",mpi_rank());
    ThreeIndOpArray NewArray;
    NewArray.Initialize(dim_i,dim_j,dim_k,msg,mpi_world_size());
    NewArray.SetBufferSize(GetBufferSize());
    NewArray.SetInCore(GetInCore());

    //open the files
    OpenFileRead();
    NewArray.OpenFileWrite();

    //prepare the buffer
    bool EndOfArray = false;
    ResetBuffer();
    //loop over all elements in array
    while (!EndOfArray){
      boost::shared_ptr<Cre> Op = GetOpFromBuffer(i,j,k,NeedTensorTrace,EndOfArray);
      if (!EndOfArray){
        //actually do the transformation
        Op->renormalise_transform(RotationMatrix,&NewStateInfo);
        //store the operator in the new Array
        NewArray.AppendOperator(Op,i,j,k,true);
      }//!EndOfarray
    }//ij
    //also transform the basis operator
    boost::shared_ptr<Cre> Basis = GetBasisOperator();
    if (Basis.get()!=0){
      Basis->renormalise_transform(RotationMatrix,&NewStateInfo);
      NewArray.SetBasisOperator(Basis);
    }
    //close the files
    CloseFileRead();
    NewArray.CloseFileWrite();
    NewArray.SetQuanta(GetQuanta());
    
    //replace the old array by the new one
    GetFileName(msg);
    Clear();
    (*this) = NewArray;
    ReNameFile(msg);
    NewArray.Clear();
    mpi_barrier();
  }
  

  //instantiate the required code
  template class OperatorArray<Wavefunction>;
  template class OperatorArray<Cre>;
  
  
  //============================================================================
  // a class that manages the storage and access of two-electron integrals
  //============================================================================
  
  //=================================
  // open the file to read data
  //=================================
  void IntegralContainer::OpenFileRead(){
    if (!InCore){
      char fn[512];
      sprintf(fn,"%s%i",FileName,mpi_rank());
      ifs.open(fn,std::ios::binary);
    }
  }
  
  //=================================
  // open the file to write data
  //=================================
  void IntegralContainer::OpenFileWrite(){
    if (!InCore){
      char fn[512];
      sprintf(fn,"%s%i",FileName,mpi_rank());
      ofs.open(fn,std::ios::binary|std::ios::ate|std::ios::app);
    }
  }
  
  //=========================
  // resize the orbital space
  //=========================
  void IntegralContainer::Resize(int i, int j, int k, int l){
    dim1 = i;dim2 = j;dim3 = k;dim4 = l;
    int N=0;
    switch (symm){
      case _NO_SYMM_:
        //resize the vector that holds the integral matrices
        rep.resize(dim1*dim2);
        //resize and fill the IndexMap
        IndexMap.ReSize(dim1,dim2);
        if (InCore){
          for (int i=0;i<dim1;i++){
            for (int j=0;j<dim2;j++){
              IndexMap(i,j) = (long) i * dim1 + j;
            }
          }
        }
        else{
          for (int i=0;i<dim1;i++){
            for (int j=0;j<dim2;j++){
              IndexMap(i,j) = -1;
            }
          }
        }
        break;

      case _COULOMB_:
        //resize the vector that holds the integral matrices
        N = (dim1+1)*(dim1)/2;
        rep.resize(N);
        //resize and fill the IndexMap
        IndexMap.ReSize(dim1,dim2);
        if (InCore){
          for (int i=0;i<dim1;i++){
            for (int j=0;j<=i;j++){
              IndexMap(i,j) = (long) i * dim1 + j;
              IndexMap(j,i) = (long) i * dim1 + j;
            }
          }
        }
        else{
          for (int i=0;i<dim1;i++){
            for (int j=0;j<dim2;j++){
              IndexMap(i,j) = -1;
            }
          }
        }
        break;

      case _EXCHANGE_:
        //resize the vector that holds the integral matrices
        N = (dim1+1)*(dim1)/2;
        rep.resize(N);
        //resize and fill the IndexMap
        IndexMap.ReSize(dim1,dim2);
        if (InCore){
          for (int i=0;i<dim1;i++){
            for (int j=0;j<=i;j++){
              IndexMap(i,j) = (long) i * dim1 + j;
              IndexMap(j,i) = (long) i * dim1 + j;
            }
          }
        }
        else{
          for (int i=0;i<dim1;i++){
            for (int j=0;j<dim2;j++){
              IndexMap(i,j) = -1;
            }
          }
        }
        break;
        
      default:
        //resize the vector that holds the integral matrices
        rep.resize(dim1*dim2);
        //resize and fill the IndexMap
        IndexMap.ReSize(dim1,dim2);
        if (InCore){
          for (int i=0;i<dim1;i++){
            for (int j=0;j<dim2;j++){
              IndexMap(i,j) = (long) i * dim1 + j;
            }
          }
        }
        else{
          for (int i=0;i<dim1;i++){
            for (int j=0;j<dim2;j++){
              IndexMap(i,j) = -1;
            }
          }
        }
        break;
    }//switch(symm)
    TransposeMap.Clear();
    TransposeMap.ReSize(dim1,dim2);
    for (int p=0;p<dim1;p++){
      for (int q=0;q<dim2;q++){
        TransposeMap(p,q)=0;
      }//q
    }//p
  }
 
  //=========================
  // store an integral matrix
  //=========================
  void IntegralContainer::SetMatrix(boost::shared_ptr<Matrix> M, const int& i, const int& j){
    if (!InCore){
      //go to the end of the file
      ofs.seekp(0,ios::end);
      //get the position
      long  pos = (long) ofs.tellp();
      //store the address
      switch(symm){
        case _NO_SYMM_:
          IndexMap(i,j) = pos;
          break;
        case _COULOMB_:
          IndexMap(i,j) = pos;
          IndexMap(j,i) = pos;
          break;
        case _EXCHANGE_:
          IndexMap(i,j) = pos;
          IndexMap(j,i) = pos;
          if (i!=j){
            TransposeMap(i,j) = 0;
            TransposeMap(j,i) = 1;
          }
          break;
      }
      //store the wavefunction
      boost::archive::binary_oarchive save_mat(ofs);
      save_mat << (*M);
    }
    else{
      long Index = IndexMap(i,j);
      rep[Index] = M;
      if ((i!=j)&&(symm==_EXCHANGE_)){
        TransposeMap(i,j) = 0;
        TransposeMap(j,i) = 1;
      }
    }
  }
  
  
  //======================================
  // Get an integral matrix
  //
  // ON INPUT  orbital indices
  //           pointer to matrix
  //
  // ON OUTPUT pointer to integral matrix
  //======================================
  boost::shared_ptr<Matrix> IntegralContainer::GetMatrix(const int& i, const int& j){
    if (!InCore){
      //get the address of operator(i,j)
      long Address = IndexMap(i,j);
      if (Address==-1){
        char msg[512];
        sprintf(msg,"\nERROR: Integral (%i,%i) could not be found!",i,j);
        pout << msg;
        boost::shared_ptr<Matrix> null;
        null.reset();
        return null;
      }
      boost::shared_ptr<Matrix> Mat(new Matrix);
      //go to the address
      ifs.seekg(Address,ios::beg);
      //check if the file-stream is open
      if (!ifs.is_open()){
        char msg[512];
        sprintf(msg,"\nERROR: Integral File is not open");
        pout << msg;
        //exit(EXIT_FAILURE);
      }
      //load the wavefunction
      boost::archive::binary_iarchive load_mat(ifs);
      load_mat >> *Mat;
      if ((symm==_EXCHANGE_)&&(TransposeMap(i,j)==true)){
        char msg[512];
        Transpose(*Mat);
      }//exchange matrices
      return Mat;
    }//not InCore
    else{
      long Index = IndexMap(i,j);
      return rep[Index];
    }//InCore
  }
  
  //=====================================================
  // Invert index ordering of matrix container
  // 
  // ON INPUT  uninitialized pointer to matrix container
  //
  // ON OUTPUT pointer points to matrix container with 
  //           reversed index ordering
  //           e.g. IJAB -> ABIJ
  //=====================================================
  boost::shared_ptr<IntegralContainer> IntegralContainer::ReverseOrder(const char *NewName){
    boost::shared_ptr<IntegralContainer> ABIJ(new IntegralContainer(dim3,dim4,dim1,dim2,symm));
    char msg[512];
    sprintf(msg,"%s",NewName);
    ABIJ->SetFileName(msg);
    int a,b,i,j;
    if (mpi_rank()==0){
      boost::shared_ptr<Matrix> Kij;
      //open the container
      OpenFileRead();
      ABIJ->OpenFileWrite();
      for (a=0;a<dim3;a++){
        //the half reversed integrals
        std::vector<boost::shared_ptr<Matrix> > Jia(dim1);
        for (i=0;i<dim1;i++){
          //allocate the required memory
          boost::shared_ptr<Matrix> M(new Matrix(dim2,dim4));
          Jia[i] = M;
          for (j=0;j<dim2;j++){
            //get the matrix in original order
            Kij = GetMatrix(i,j);
            for (b=0;b<dim4;b++){
              //define the halfway reversed integrals
              Jia[i]->element(j,b) = Kij->element(a,b);
            }//b
          }//j
        }//i
        switch (symm){
          case (_NO_SYMM_):  
            for (b=0;b<dim4;b++){
              boost::shared_ptr<Matrix> Kab(new Matrix(dim1,dim2));
              for (i=0;i<dim1;i++){
                for (j=0;j<dim2;j++){
                  //define the fully reversed integrals
                  Kab->element(i,j) = Jia[i]->element(j,b);
                }//j
              }//i
              //store the fully reversed integrals
              ABIJ->SetMatrix(Kab,a,b);
            }//b
          break;
          case (_COULOMB_):
          case (_EXCHANGE_):
            for (b=0;b<=a;b++){
              boost::shared_ptr<Matrix> Kab(new Matrix(dim1,dim2));
              for (i=0;i<dim1;i++){
                for (j=0;j<dim2;j++){
                  //define the fully reversed integrals
                  Kab->element(i,j) = Jia[i]->element(j,b);
                }//j
              }//i
              //store the fully reversed integrals
              ABIJ->SetMatrix(Kab,a,b);
            }//b
          break;
        }//symm

        //clear the memory
        Jia.clear();
      }//a
      //close the container
      CloseFileRead();
      ABIJ->CloseFileWrite();
    }
    
    //broadcast the rearranged container
    ABIJ->Broadcast(0);
    
    return ABIJ;
  }
  
  //==============================
  //Broadcast the stored integrals
  //==============================
  void IntegralContainer::Broadcast(int sender){
#ifndef SERIAL
    mpi::communicator world;
    //open the sender file
    if (world.rank()==sender){
      OpenFileRead();
    }
    else{
      //remove all integrals before receiving
      Clear();
      //resize the index and transpose map
      Resize(dim1,dim2,dim3,dim4);
      //open the receiver file
      OpenFileWrite();
    }
    //now broadcast the integrals
    switch (symm){
      case _NO_SYMM_:
        for (int i=0;i<dim1;i++){
          for (int j=0;j<dim2;j++){
            boost::shared_ptr<Matrix> Mij (new Matrix(dim3,dim4));
            //if you are the sender, read the integral matrix
            if (world.rank()==sender){
              Mij = GetMatrix(i,j);
            }
            if (Mij) mpi::broadcast(world,*Mij,0);
            //if you are not the sender, store the integral matrix
            if (world.rank()!=sender){
              SetMatrix(Mij,i,j);
            }
          }//j
        }//i
        break;
      case _COULOMB_:
      case _EXCHANGE_:
        for (int i=0;i<dim1;i++){
          for (int j=0;j<=i;j++){
            boost::shared_ptr<Matrix> Mij (new Matrix(dim3,dim4));
            //if you are the sender, read the integral matrix
            if (world.rank()==sender){
              Mij = GetMatrix(i,j);
            }
            if (Mij) mpi::broadcast(world,*Mij,0);
            //if you are not the sender, store the integral matrix
            if (world.rank()!=sender){
              SetMatrix(Mij,i,j);
            }
          }//j
        }//i
        break;
    }//symm
    mpi_barrier();
    //close the files
    if (world.rank()==sender){
      CloseFileRead();
    }
    else{
      CloseFileWrite();
    }
#endif
  }
  
        
}
