/* 
 * File:   ripdm_operators.h
 * Author: roemelt
 *
 * Created on April 4, 2013, 11:55 AM
 */

#ifndef NEVPT2_OPERATORS_H
#define	NEVPT2_OPERATORS_H

#define _3PDM_NORMAL_ 0
#define _3PDM_A_ 1
#define _3PDM_A_PRIME_ 2


#include "spinblock.h"
#include "wavefunction.h"
#include "BaseOperator.h"
#include "nevpt2_mpi.h"
#include <vector>
#include <newmat.h>
#include "para_array.h"
#include <boost/serialization/serialization.hpp>
#include <cassert>
#include <boost/make_shared.hpp>


namespace SpinAdapted {
  
  //============================================================================
  // An element of the buffer for OperatorArrays: Operator and  two indices
  //============================================================================
  template <class T> class BufferElement : private boost::shared_ptr<T> {
  private:
    int i;
    int j;
    bool TensorTrace;//does this operator need to be traced
        
  public:
    //-----------
    //Constructor
    //-----------
    BufferElement(){
      i = -1;
      j = -1;
      TensorTrace = true;
    };
    BufferElement(boost::shared_ptr<T> Op,int first, int second):boost::shared_ptr<T>(Op){
      i = first;
      j = second;
      TensorTrace = true;
    };
    BufferElement(boost::shared_ptr<T> Op,int first, int second, bool tt):boost::shared_ptr<T>(Op){
      i = first;
      j = second;
      TensorTrace = tt;
    };
    BufferElement(const BufferElement &E):boost::shared_ptr<T>(E){
      i           = E.Geti();
      j           = E.Getj();
      TensorTrace = E.NeedTensorTrace();
    };
    //----------
    //Destructor
    //----------
    ~BufferElement(){clear();};
    
    //-----------
    //the getters
    //-----------
    boost::shared_ptr<T> GetElement() const{return (*this);}            //get the element
    int Geti() const{return i;}                                         //get the indices
    int Getj() const{return j;}                                         //get the indices
    bool NeedTensorTrace() const{return TensorTrace;}                   //get the tensor trace flag
    
    //-----------
    //the setters
    //------------
    void SetTensorTrace(bool tt){TensorTrace = tt;}                         //set the tensor trace flag

    //-------------------
    //assignment function
    //-------------------
    BufferElement<T> operator=(const BufferElement<T> &E){
      boost::shared_ptr<T>::operator =(E);
      i           = E.Geti();
      j           = E.Getj();
      TensorTrace = E.NeedTensorTrace();
      return *this;
    }
    
    //-----------------------
    //disk handling functions
    //-----------------------
    //write the element to disk (open output file stream as input required)
    void write(ofstream &ofs){
      //store the indices
      ofs.write(reinterpret_cast<char*>(&i), sizeof(int));
      ofs.write(reinterpret_cast<char*>(&j), sizeof(int));
      //store the tensortrace flag
      ofs.write(reinterpret_cast<char*>(&TensorTrace), sizeof(bool));
      //store the wavefunction
      boost::archive::binary_oarchive save_wave(ofs);
      save_wave << *(*this);
    }
    
    //write the element to disk (open input file stream as input required)
    void read(ifstream &ifs){
      //allocate some memory
      this->reset(new T);
      //read the indices
      ifs.read((char*)&i,sizeof(int));
      ifs.read((char*)&j,sizeof(int));
      //read the TensorTrace flag
      ifs.read((char*)&TensorTrace,sizeof(bool));
      //load the operator
      boost::archive::binary_iarchive load_wave(ifs);
      load_wave >> *(*this);
    }
    //-----------------
    //clear the element
    //-----------------
    void clear(){
      this->reset();
    }
  };
  
  
  //============================================================================
  // the class template that manages the access to operators 
  //============================================================================
  template <class T> class OperatorArray {
    private:
        std::vector<BufferElement<T> > Operators;       //the vector that holds the matrices in the InCore case
        ObjectMatrix<std::pair<int,long> >IJAddr;       //this matrix holds the address for any pair of orbital indices
        bool InCore;                                    //should the data be store in memory
        int length;                                    //the number of Operators stored. It is also used for creating pairindices
        char FileName[512];                             //the FileName
        int NumFiles;                                   //the number of files used
        int ActFile;                                    //the number of the open file
        ofstream ofs;                                   //the output file stream
        ifstream ifs;                                   //the input file stream
        std::vector<BufferElement<T> > Buffer;          //InCore buffer for operators read from disk
        long BufferSize;                                //number of element in Buffer
        int CurrentElement;                             //the current element read from the buffer
        long BufferPosition;                            //the position of the buffer in the file
        bool HaveStateInfo;                             //do we have a StateInfo file?
        SpinQuantum Q;                                  //the quanta of the stored operators
        boost::shared_ptr<T> Basis;                     //the basis operator without the integrals
        bool BasisOnly;                                 //a flag that indicated whether only the basis operator should be altered
        friend class boost::serialization::access;

        template<class Archive> void serialize(Archive & ar, const unsigned int version) {
            ar & IJAddr & length & FileName & InCore & NumFiles & ActFile & BufferSize;
            ar & BufferPosition & HaveStateInfo & Q;
        }

    public:
        //---------------------------
        // Constructor and Destructor
        //---------------------------
        OperatorArray() {
          length = 0; 
          InCore = true;
          sprintf(FileName,"TerminalOperators.tmp");
          NumFiles = mpi_world_size();
          ActFile  = 0;
          BufferSize = 0;
          CurrentElement=0;
          BufferPosition = 0;
          HaveStateInfo = false;
          BasisOnly = false;
        };
        OperatorArray(const OperatorArray<T> &O){
          Operators     = O.GetInCoreData();
          IJAddr        = O.GetAddressMatrix();
          InCore        = O.GetInCore();
          length        = O.GetLength();
          O.GetFileName(FileName);
          NumFiles      = O.GetNumFiles();
          BufferSize    = O.GetBufferSize();
          Buffer        = O.GetBuffer();
          BufferPosition= O.GetBufferPosition();
          Q             = O.GetQuanta();
          Basis         = O.GetBasisOperator();
          BasisOnly     = O.GetBasisOnly();
          if (O.GetHaveStateInfo()){SetStateInfo(O.GetStateInfo());}
        }
        ~OperatorArray() {
          //Note that the files on disk are left untouched by this routine
          //In order to delete all files from disk, call clear() function)
          int i;
          for(i=0;i<Operators.size();i++){
            Operators[i].clear();
          }//i
          Operators.clear();
          Buffer.clear();
          length = 0;
          NumFiles = 0;
          ActFile = 0;
          IJAddr.Clear();
          if (!Basis) Basis.reset();
          if (mpi_rank()==0){
            char file[512];
            sprintf(file,"%s.StateInfo.tmp",FileName);
            remove(file);
          }
          BasisOnly = false;
        };
        OperatorArray<T> operator=(const OperatorArray<T> &O){
          Operators     = O.GetInCoreData();
          IJAddr        = O.GetAddressMatrix();
          InCore        = O.GetInCore();
          length        = O.GetLength();
          O.GetFileName(FileName);
          NumFiles      = O.GetNumFiles();
          BufferSize    = O.GetBufferSize();
          Buffer        = O.GetBuffer();
          BufferPosition= O.GetBufferPosition();
          HaveStateInfo = O.GetHaveStateInfo();
          if (HaveStateInfo){SetStateInfo(O.GetStateInfo());}
          Q             = O.GetQuanta();
          Basis         = O.GetBasisOperator();
          BasisOnly     = O.GetBasisOnly();
          return *this;
        }
        //-------------
        //File Handling
        //-------------
        void OpenFileRead(int fn=mpi_rank());                                                            //Open the file to read data
        void OpenFileWrite();                                                                   //Open the file to write data
        void CloseFileRead();                                                                   //Close the file
        void CloseFileWrite();                                                                  //Close the file
        void SetFileName(const char *Name);                                                           //set the FileName
        void ReNameFile(const char *NewName);                                                         //renames the file
        void SetInCore(bool ic){InCore = ic;};                                                  //decide whether the data are stored in memory or on disk
        bool GetInCore() const{return InCore;};                                                 //tell if the data are stored in core 
        void GetFileName(char *Name) const{sprintf(Name,"%s",FileName);};                       //Get the FileName
        int GetNumFiles() const{return NumFiles;};                                              //get the number of files used
        void SetNumFiles(int NF){NumFiles = NF;};                                               //set the number of files used
        void Store(bool sync=true);                                                             //Store ALL the Data on disk
        void Retrieve(const char *name, bool clear=true);                                             //Retrieve the data stored on disk
        
        //-----------
        //the getters
        //-----------
        int GetLength() const{return length;};                                                  //get the number of operators stored
        std::vector<BufferElement<T> > GetInCoreData()const{return Operators;};                      //get all data if InCore
        boost::shared_ptr<T> GetOperator(int i, int j);                                         //get one operator
        boost::shared_ptr<T> GetOpPal(int i, int j);                                            //get one operator (one proc reads and broadcasts the operator) 
        boost::shared_ptr<T> GetLocalOp(int i, int j);                                          //get one operator (only returns operator if it is available to this process, otherwise returns a zero pointer)
        boost::shared_ptr<T> GetOpFromBuffer(int &i, int &j,bool &EndOfFile);                   //get one operator from the buffer (if buffer is empty or we are at the end of the buffer read the next chunk from the file)
        boost::shared_ptr<T> GetOpFromBuffer(int &i, int &j,bool &TensorTrace,bool &EndOfFile); //get one operator from the buffer (if buffer is empty or we are at the end of the buffer read the next chunk from the file)
        std::pair<int,long> GetAdress(int i, int j) const {return IJAddr(i, j);};               //get the address of a pair of orbital indices
        int GetNumberOfOperators() const {return length;};                                      //get the number of stored operators
        ObjectMatrix<std::pair<int,long> > GetAddressMatrix() const{return IJAddr;};            //get the entire address matrix
        StateInfo GetStateInfo()const;                                                          //get the Stateinfo corresponding to the representation of the stored operators
        bool GetHaveStateInfo()const{return HaveStateInfo;}                                     //get the StateInfo flag
        SpinQuantum GetQuanta()const{return Q;}                                                 //get the operator quanta
        boost::shared_ptr<T> GetBasisOperator()const{return Basis;}                             //get the basis operator
        bool GetBasisOnly()const{return BasisOnly;}                                             //get the BasisOnly flag
        
        //----------
        //Initialize
        //----------
        void Initialize(int dim, const char *name,int nf){                                            //initialize the array
          ResizeOrbSpace(dim,true);
          SetFileName(name);
          NumFiles = nf;
        }
        void Initialize(int dim, int dim_, const char *name,int nf){                                  //initialize the array
          ResizeOrbSpace(dim,dim_,true);
          SetFileName(name);
          NumFiles = nf;
        }
        
        //-----------
        //the setters
        //-----------
        void DuplicateAddress(int originali, int originalj,int duplicatei, int duplicatej);     //duplicate the address of an operator 
        void AppendOperator(boost::shared_ptr<T> A,int i,int j, bool TensorTrace=false, bool Print=true);//append one operator
        void ReplaceOpInBuffer(boost::shared_ptr<T> A,int i,int j,bool TensorTrace=true);       //replace one operator with another one using the buffer vector!!!NOTE: They need to have the exact same size or the whole file is corrupted!!!
        void SetStateInfo(const StateInfo &SI);                                                 //set the Stateinfo corresponding to the representation of the stored operators
        void ResizeOrbSpace(int m, bool init);                                                  //resize the address matrix
        void ResizeOrbSpace(int m, int n, bool init);                                           //resize the address matrix
        void SetQuanta(SpinQuantum q){Q=q;}                                                     //Set the operator quanta
        void SetBasisOperator(boost::shared_ptr<T> &op){Basis = op;}                            //Set the basis operator
        void AddToBasisOperator(boost::shared_ptr<T> &op,double fac = 1.0);                     //Add to the basis operator
        void SetBasisOnly(bool bo){BasisOnly=bo;}                                               //Set the BasisOnly flag
        void AppendArray(OperatorArray<T> &OpArray);                                            //append an entire array to this one
        
        //-----------------------------
        //the Buffer Handling functions
        //-----------------------------
        std::vector<BufferElement<T> > GetBuffer() const{return Buffer;};                       //get the buffer
        int GetBufferSize() const{return BufferSize;}                                           //get the size of the buffer
        void SetBufferSize(int BS){BufferSize=BS;}                                              //Set the size of the buffer
        void CalcBufferSize(int MaxMemSize, int M);                                             //Evaluate the number of elements in the buffer according to max. available memory and M
        int GetBufferPosition() const{return BufferPosition;};                                  //get the postion of the buffer in the file
        void ResetBuffer();                                                                     //reset the buffer and rewind the input stream
        
        //------------------------------------------------------------
        // Transform the stored operators to a new block configuration
        //------------------------------------------------------------
        virtual void Transform(){};
        //transpose the stored operators
        virtual void Transpose(){};
        
        //-------------
        //Clean Up data
        //-------------
        void Clear();                                                                           //deletes all data from disk
        
        //------------------------------------------------
        //Synchronize the available data between processes
        //------------------------------------------------
        void Synchronize();                                                                     //gather all data from parallel processes
        
        
    };

    //==========================================================================
    // An Array of Wavefunctions
    //==========================================================================
    class WavefunctionArray : public OperatorArray<Wavefunction>{
    public:
        // Transform the stored operators to a new block configuration
        void Transform(vector<Matrix> leftRotMatrix, const vector<Matrix> &rightRotMatrix, SpinBlock &big,
                                const StateInfo &OldStateInfo, const StateInfo &NewStateInfo, 
                                int root, bool Transpose);
        //transpose the stored operators
        void Transpose(const StateInfo &oldStateInfo, const StateInfo &newStateInfo,SpinBlock &big, int root);
        //Initialize the basis operator in a given representation
        void ResetBasisOperator(SpinBlock &b,SpinQuantum &Q);
    };
    
    //==========================================================================
    // An Array of Operators
    //==========================================================================
    //Note: The operators are not actually Cre operators, this is just a dummy 
    //choice
    class OpArray: public OperatorArray<Cre>{
    public:
        // Renormalize the stored operators 
        void RenormalizeTransform(SpinBlock &big, const StateInfo &NewStateInfo, const vector<Matrix> &RotationMatrix);
        // TensorTrace the stored operators
        void TensorTrace(SpinBlock &big, const StateInfo &NewStateInfo);
        //transpose the stored operators
        void Transpose(){};
    };
    
    //==========================================================================
    // An Array of Operators for three index operators
    //==========================================================================
    class ThreeIndOpArray: public OpArray{
    private:
        int Dimensions[3];//the dimensions when three indices are used
        
        //-------------------------------------------------------------------
        //make a unique mapping from a three index space to a two index space
        //-------------------------------------------------------------------
        pair<int, int> GenerateIndexPair(const int &i, const int &j, const int k){
          return make_pair(i,j*Dimensions[2]+k);
        }
        vector<int> GenerateIndexTriple(const int &i, const int &j){
          vector<int> tmp;
          tmp.push_back(i);
          int j_ = j/Dimensions[2];tmp.push_back(j_);
          int k_ = j%Dimensions[2];tmp.push_back(k_);
          return tmp;
        }
        
    public:
        //--------------------------
        //Constructor and Destructor
        //--------------------------
        ThreeIndOpArray(){
          for (int i=0;i<3;i++){
            Dimensions[i]=-1;
          }
        };
        ThreeIndOpArray(const ThreeIndOpArray &T):OpArray(T){
          T.GetThreeIndexDims(Dimensions);
        };
        ~ThreeIndOpArray(){};
        //------------------
        //Assigment function
        //------------------
        ThreeIndOpArray operator=(const ThreeIndOpArray &T){
          OpArray::operator =(T);
          T.GetThreeIndexDims(Dimensions);
          return (*this);
        }
        //--------------------------------------
        //get the dimensions in the 3 index case
        //--------------------------------------
        void GetThreeIndexDims(int *dim) const {        
          for (int i=0;i<3;i++){
            dim[i] = Dimensions[i];
          }
        }
        //--------------------------------
        //override the initialize function
        //--------------------------------
        void Initialize(int dim, int dim_, int dim__, const char *name,int nf){
          ResizeOrbSpace(dim,dim_*dim__,true);
          SetFileName(name);
          SetNumFiles(nf);
          Dimensions[0]=dim;Dimensions[1]=dim_;Dimensions[2]=dim__;
        }
        void Initialize(int dim, const char *name,int nf){
          ResizeOrbSpace(dim,dim*dim,true);
          SetFileName(name);
          SetNumFiles(nf);
          Dimensions[0]=dim;
          Dimensions[1]=dim;
          Dimensions[2]=dim;
        }
        
        //--------------------------------
        //the Get Operator functions
        //--------------------------------
        //get one operator
        boost::shared_ptr<Cre> GetOperator(int i, int j, int k){                                  
          pair<int, int> tmp = GenerateIndexPair(i,j,k);
          int i_ = tmp.first;
          int j_ = tmp.second;
          return OpArray::GetOperator(i_,j_);
        }
        //get one operator (one proc reads and broadcasts the operator) 
        boost::shared_ptr<Cre> GetOpPal(int i, int j, int k){                                     
          pair<int, int> tmp = GenerateIndexPair(i,j,k);
          int i_ = tmp.first;
          int j_ = tmp.second;
          return OpArray::GetOpPal(i_,j_);
        }
        //get one operator (only returns operator if it is available to this process, otherwise returns a zero pointer)
        boost::shared_ptr<Cre> GetLocalOp(int i, int j, int k){                                   
          pair<int, int> tmp = GenerateIndexPair(i,j,k);
          int i_ = tmp.first;
          int j_ = tmp.second;
          return OpArray::GetLocalOp(i_,j_);
        }
        //get one operator from the buffer (if buffer is empty or we are at the end of the buffer read the next chunk from the file)
        boost::shared_ptr<Cre> GetOpFromBuffer(int &i, int &j,int &k, bool &EndOfFile){
          int i_,j_;
          boost::shared_ptr<Cre> tmp = OpArray::GetOpFromBuffer(i_,j_,EndOfFile);
          vector<int> indices = GenerateIndexTriple(i_,j_);
          i = indices[0];
          j = indices[1];
          k = indices[2];
          return tmp;
        }
        //get one operator from the buffer (if buffer is empty or we are at the end of the buffer read the next chunk from the file)
        boost::shared_ptr<Cre> GetOpFromBuffer(int &i, int &j,int &k, bool &TensorTrace,bool &EndOfFile){
          int i_,j_;
          boost::shared_ptr<Cre> tmp = OpArray::GetOpFromBuffer(i_,j_,TensorTrace,EndOfFile);
          vector<int> indices = GenerateIndexTriple(i_,j_);
          i = indices[0];
          j = indices[1];
          k = indices[2];
          return tmp;
        }
        
        //--------------------------------
        //the Set Operator functions
        //--------------------------------
        void AppendOperator(boost::shared_ptr<Cre> A,int i,int j, int k, bool TensorTrace=false){    //append one operator
          pair<int, int> tmp = GenerateIndexPair(i,j,k);
          int i_ = tmp.first;
          int j_ = tmp.second;
          OpArray::AppendOperator(A,i_,j_,TensorTrace);
        }
        void ReplaceOpInBuffer(boost::shared_ptr<Cre> A,int i,int j,int k,bool TensorTrace=true){       //replace one operator with another one using the buffer vector
          pair<int, int> tmp = GenerateIndexPair(i,j,k);
          int i_ = tmp.first;
          int j_ = tmp.second;
          OpArray::ReplaceOpInBuffer(A,i_,j_,TensorTrace);
        }
        
        //-------------------------------------------------
        // Renormalize and TensorTrace the stored operators 
        //--------------------------------------------------
        void RenormalizeTransform(SpinBlock &big, const StateInfo &NewStateInfo, const vector<Matrix> &RotationMatrix);
        void TensorTrace(SpinBlock &big, const StateInfo &NewStateInfo);
        
        //--------------------------------
        //the Store and Retrieve functions
        //--------------------------------
        void Store(bool sync=true);
        void Retrieve(const char *name);                                                              


    };
    
      
    
    //==============================================================================
    // A class that stores three particle density matrices
    //==============================================================================
  class array_6d {
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
      ar & dim1 & dim2 & dim3 & dim4 & dim5 & dim6 & d2_x_d3_x_d4_x_d5_x_d6 & d3_x_d4_x_d5_x_d6
              & d4_x_d5_x_d6 & d5_x_d6;
      ar & data;
    }
    std::vector<double> data; // the vector that holds the actual data
    int dim1;
    int dim2;
    int dim3;
    int dim4;
    int dim5;
    int dim6;
    int d2_x_d3_x_d4_x_d5_x_d6;
    int d3_x_d4_x_d5_x_d6;
    int d4_x_d5_x_d6;
    int d5_x_d6;

  public:
    //------------
    // Constructor
    //------------
    array_6d() : dim1(0), dim2(0), dim3(0), dim4(0), dim5(0), dim6(0), d2_x_d3_x_d4_x_d5_x_d6(0),
    d3_x_d4_x_d5_x_d6(0), d4_x_d5_x_d6(0), d5_x_d6(0) {};
    
    array_6d(int i, int j, int k, int l, int m, int n) {
      dim1 = i;
      dim2 = j;
      dim3 = k;
      dim4 = l;
      dim5 = m;
      dim6 = n;
      d2_x_d3_x_d4_x_d5_x_d6 = dim2 * dim3 * dim4 * dim5*dim6;
      d3_x_d4_x_d5_x_d6 = dim3 * dim4 * dim5*dim6;
      d4_x_d5_x_d6 = dim4 * dim5*dim6;
      d5_x_d6 = dim5*dim6;
      data.resize(dim1 * dim2 * dim3 * dim4 * dim5*dim6,0.0);
      //boost::shared_ptr<std::vector<double> > new_data(new vector<double>(dim1 * dim2 * dim3 * dim4 * dim5*dim6, 0.0));
      //data = new_data;
    }

    array_6d(int i) {
      dim1 = i;
      dim2 = i;
      dim3 = i;
      dim4 = i;
      dim5 = i;
      dim6 = i;
      d2_x_d3_x_d4_x_d5_x_d6 = dim2 * dim3 * dim4 * dim5*dim6;
      d3_x_d4_x_d5_x_d6 = dim3 * dim4 * dim5*dim6;
      d4_x_d5_x_d6 = dim4 * dim5*dim6;
      d5_x_d6 = dim5*dim6;
      data.resize(dim1 * dim2 * dim3 * dim4 * dim5*dim6,0.0);
      //boost::shared_ptr<std::vector<double> > new_data(new vector<double>(dim1 * dim2 * dim3 * dim4 * dim5*dim6, 0.0));
      //data = new_data;
    }
    
    array_6d(const array_6d &a){
      dim1 = a.get_dim1();
      dim2 = a.get_dim2();
      dim3 = a.get_dim3();
      dim4 = a.get_dim4();
      dim5 = a.get_dim5();
      dim6 = a.get_dim6();
      d2_x_d3_x_d4_x_d5_x_d6 = dim2 * dim3 * dim4 * dim5*dim6;
      d3_x_d4_x_d5_x_d6 = dim3 * dim4 * dim5*dim6;
      d4_x_d5_x_d6 = dim4 * dim5*dim6;
      d5_x_d6 = dim5*dim6;
      data.resize(dim1 * dim2 * dim3 * dim4 * dim5*dim6,0.0);
      for (int i=0;i<data.size();i++){
        data[i] = a(i);
      }
    }
    
    ~array_6d(){
      clear();
    }
    
    array_6d operator=(const array_6d &a){
      dim1 = a.get_dim1();
      dim2 = a.get_dim2();
      dim3 = a.get_dim3();
      dim4 = a.get_dim4();
      dim5 = a.get_dim5();
      dim6 = a.get_dim6();
      d2_x_d3_x_d4_x_d5_x_d6 = dim2 * dim3 * dim4 * dim5*dim6;
      d3_x_d4_x_d5_x_d6 = dim3 * dim4 * dim5*dim6;
      d4_x_d5_x_d6 = dim4 * dim5*dim6;
      d5_x_d6 = dim5*dim6;
      data.resize(dim1 * dim2 * dim3 * dim4 * dim5*dim6,0.0);
      for (int i=0;i<data.size();i++){
        data[i] = a(i);
      }
      return *this;
    }
    
    //------------------
    //get the dimensions
    //------------------
    int get_dim1() const {return dim1;}
    int get_dim2() const {return dim2;}
    int get_dim3() const {return dim3;}
    int get_dim4() const {return dim4;}
    int get_dim5() const {return dim5;}
    int get_dim6() const {return dim6;}
    int get_size() const {return data.size();}

    //-------------------
    // access to the data
    //-------------------
    double& operator()(int i, int j, int k, int l, int m, int n) {
      assert((0 <= i)&&(i < dim1));
      assert((0 <= j)&&(j < dim2));
      assert((0 <= k)&&(k < dim3));
      assert((0 <= l)&&(l < dim4));
      assert((0 <= m)&&(m < dim5));
      assert((0 <= n)&&(n < dim6));
      return data[i * d2_x_d3_x_d4_x_d5_x_d6 + j * d3_x_d4_x_d5_x_d6 + k * d4_x_d5_x_d6 + l * d5_x_d6 + m * dim6 + n];
    }

    double operator()(int i, int j, int k, int l, int m, int n) const {
      assert((0 <= i)&&(i < dim1));
      assert((0 <= j)&&(j < dim2));
      assert((0 <= k)&&(k < dim3));
      assert((0 <= l)&&(l < dim4));
      assert((0 <= m)&&(m < dim5));
      assert((0 <= n)&&(n < dim6));
      return data[i * d2_x_d3_x_d4_x_d5_x_d6 + j * d3_x_d4_x_d5_x_d6 + k * d4_x_d5_x_d6 + l * d5_x_d6 + m * dim6 + n];
    }

    double& operator() (int i) {
      return data[i];
    }

    double operator() (int i) const {
      return data[i];
    }
    
    //-----------------
    // resize the array
    //-----------------
    void resize(int i, int j, int k, int l, int m, int n) {
      clear();
      dim1 = i;
      dim2 = j;
      dim3 = k;
      dim4 = l;
      dim5 = m;
      dim6 = n;
      d2_x_d3_x_d4_x_d5_x_d6 = dim2 * dim3 * dim4 * dim5*dim6;
      d3_x_d4_x_d5_x_d6 = dim3 * dim4 * dim5*dim6;
      d4_x_d5_x_d6 = dim4 * dim5*dim6;
      d5_x_d6 = dim5*dim6;
      data.resize(dim1 * dim2 * dim3 * dim4 * dim5*dim6,0.0);
      //boost::shared_ptr<std::vector<double> > new_data(new vector<double>(dim1 * dim2 * dim3 * dim4 * dim5*dim6, 0.0));
      //data = new_data;
      }

    void resize(int i) {
      clear();
      dim1 = i;
      dim2 = i;
      dim3 = i;
      dim4 = i;
      dim5 = i;
      dim6 = i;
      d2_x_d3_x_d4_x_d5_x_d6 = dim2 * dim3 * dim4 * dim5*dim6;
      d3_x_d4_x_d5_x_d6 = dim3 * dim4 * dim5*dim6;
      d4_x_d5_x_d6 = dim4 * dim5*dim6;
      d5_x_d6 = dim5*dim6;
      data.resize(dim1 * dim2 * dim3 * dim4 * dim5*dim6,0.0);
      //boost::shared_ptr<std::vector<double> > new_data(new vector<double>(dim1 * dim2 * dim3 * dim4 * dim5*dim6, 0.0));
      //data = new_data;
    }

    //-----------------
    // Inititialize
    //-----------------
    void initialize(){
      for (int i=0;i<data.size();i++) data[i] = 0.0;
    }
    
    //------------------
    // add another array
    //------------------
    array_6d& operator+=(const array_6d &A) {
      assert(A.get_dim1() == get_dim1());
      assert(A.get_dim2() == get_dim2());
      assert(A.get_dim3() == get_dim3());
      assert(A.get_dim4() == get_dim4());
      assert(A.get_dim5() == get_dim5());
      assert(A.get_dim6() == get_dim6());
      for (int i = 0; i < get_size(); i++) {
        data[i] += A(i);
      }
      return *this;
    }
    
    //--------------------
    // Clean up the memory
    //--------------------
    void clear(){
      data.clear();
      dim1 = 0;
      dim2 = 0;
      dim3 = 0;
      dim4 = 0;
      dim5 = 0;
      dim6 = 0;
      d2_x_d3_x_d4_x_d5_x_d6 = 0;
      d3_x_d4_x_d5_x_d6 = 0;
      d4_x_d5_x_d6 = 0;
      d5_x_d6 = 0;
      
    }
    
  };
  
  
#define _NO_SYMM_ 0
#define _COULOMB_ 1
#define _EXCHANGE_ 2  

  //============================================================================
  // a class that manages the storage and access of two-electron integrals
  //============================================================================
  class IntegralContainer{
  private:
    std::vector<boost::shared_ptr<Matrix> >  rep;
    int dim1;
    int dim2;
    int dim3;
    int dim4;
    int symm;
    ObjectMatrix<long> IndexMap;
    ObjectMatrix<int> TransposeMap;
    char FileName[512];
    ofstream ofs;                                   
    ifstream ifs;
    bool InCore;
        
    
  public:
    //--------------------------
    //Constructor and Destructor
    //--------------------------
    IntegralContainer(){
      dim1 = 0;dim2 = 0;dim3 = 0;dim4 = 0;
      symm = _NO_SYMM_;
      rep.resize(dim1*dim2);
      IndexMap.ReSize(0,0);
      TransposeMap.ReSize(0,0);
      InCore = false;
      sprintf(FileName,"IJKL.tmp");
    }
    IntegralContainer(int i, int j, int k, int l, int s){
      symm = s;
      InCore = false;
      Resize(i,j,k,l);
      IndexMap.resize(i,j);
      TransposeMap.resize(i,j);
      for (int p=0;p<i;p++){
        for (int q=0;q<j;q++){
          TransposeMap(p,q) = 0;
        }//q
      }//p
      sprintf(FileName,"IJKL.tmp");
      dim1 = i;dim2 = j;dim3 = k;dim4 = l;
    }
    IntegralContainer(const IntegralContainer &IC){
      Resize(IC.GetDim1(),IC.GetDim2(),IC.GetDim3(),IC.GetDim4());
      symm = IC.GetSymm();
      IC.GetFileName(FileName);
      IndexMap = IC.GetIndexMap();
      TransposeMap = IC.GetTransposeMap();
      InCore = IC.GetInCore();
      if (InCore){
        IC.GetRepository(rep);
      }
    }
    ~IntegralContainer(){
      Clear();
    }

    //-----------------------------
    //Resize the Integral Container
    //-----------------------------
    void Resize(int i, int j, int k, int l);
    
    //------------------
    //access to the data
    //------------------
    boost::shared_ptr<Matrix> GetMatrix(const int &i, const int &j);
    
    void SetMatrix(boost::shared_ptr<Matrix> M, const int &i, const  int &j);
    
    //-----------------------
    //open and close the file
    //-----------------------
    void OpenFileRead();        //Open the file to read data
    void OpenFileWrite();       //Open the file to write data
    void CloseFileRead(){ifs.close();}  //Close the file
    void CloseFileWrite(){ofs.close();} //Close the file

    //-----------------
    // General features
    //-----------------
    void SetFileName(const char *Name){sprintf(FileName,"%s",Name);}       //set the FileName
    void SetInCore(bool IC){InCore = IC;}                       //set the InCore flag
    bool GetInCore() const{return InCore;}                      //get the InCore flag
    int GetDim1() const{return dim1;}                           //get the dimension of the orbital space
    int GetDim2() const{return dim2;}                           //get the dimension of the orbital space
    int GetDim3() const{return dim3;}                           //get the dimension of the orbital space
    int GetDim4() const{return dim4;}                           //get the dimension of the orbital space
    int GetSymm() const{return symm;}                           //get the symmetry of the Container
    void GetFileName(char *Name) const{                         //get the FileName
      sprintf(Name,"%s",FileName);
    }
    ObjectMatrix<long> GetIndexMap() const{return IndexMap;}          //get the address matrix
    ObjectMatrix<int> GetTransposeMap() const{return TransposeMap;}   //get the transpose matrix
    int GetTranspose(int i,int j) const{return TransposeMap(i,j);}    //does matrix (i,j) have to be transposed
    void GetRepository(vector<boost::shared_ptr<Matrix> > &V) const{  //get the data (if InCore)
      V.clear();
      V.resize(rep.size());
      for (int i=0;i<rep.size();i++){
        V[i] = rep[i];
      }
    }
    IntegralContainer operator =(const IntegralContainer &IC){
      Clear();
      Resize(IC.GetDim1(),IC.GetDim2(),IC.GetDim3(),IC.GetDim4());
      symm = IC.GetSymm();
      IC.GetFileName(FileName);
      IndexMap = IC.GetIndexMap();
      TransposeMap = IC.GetTransposeMap();
      InCore = IC.GetInCore();
      if (InCore){
        IC.GetRepository(rep);
      }
      return *this;
    }
    
    //------------------------------------
    // reverse ordering of orbital indices
    //------------------------------------
    boost::shared_ptr<IntegralContainer> ReverseOrder(const char* NewName);
    
    //-------------------------------
    // broadcast the stored integrals
    //-------------------------------
    void Broadcast(int sender);
    
    //------------------------------
    //Free the memory and disk space
    //------------------------------
    void Clear(){
      IndexMap.Clear();
      TransposeMap.Clear();
      if (!InCore){
        char fn[512];
        sprintf(fn,"%s%i",FileName,mpi_rank());
        remove(fn);
      }
      else{
        rep.clear();
      }
    }
    
    
    
  };
  
}




#endif	/* RIPDM_OPERATORS_H */

