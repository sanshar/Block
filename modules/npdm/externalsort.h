#ifndef EXTERNAL_SORT_H
#define EXTERNAL_SORT_H

#include "global.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <vector> 
#include <fstream>
#include <iostream>
#define Buff_SIZE 1024*1024*8
#define merge_buff_size (Buff_SIZE*100) // during the merge of file, only one processor works. More momery available.

namespace SpinAdapted{
namespace Npdm{

  using namespace std;

template<class T>
class has_index{
  template<class U>
  static std::true_type __test(typename U::indextype);
  template<class >
  static std::false_type __test(...);
public:
  static constexpr const bool value = std::is_same<std::true_type, decltype(__test<T>(0))>::value;
};

  class index_element{
    public:
    typedef long indextype ;
    typedef double valuetype;
    long index;
    double element;
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & index;
      ar &element;
    }
    friend ostream& operator<<(ostream& os, const index_element& onepiece) {
      os << onepiece.index<< "\t\t\t" << onepiece.element;
      return os;
    }
    bool operator<(const index_element& onepiece2) const{
      if(index < onepiece2.index) return true;
      if(index > onepiece2.index) return false;
      if(index == onepiece2.index){
        cout << "there are elements with same index"<<endl;
        if(element !=onepiece2.element) 
        {
          cout << " and their values are different, abort"<<endl;
          abort();
        }

      }
    }
  };

//std::pair cannot be sent by boost mpi
//I need a bool to tell receivers that this is the last data, prepare to close the receiving. 
//Otherwise, the time that ending signal arrives cannot be made sure of. 
//The receivers will always wait. Or checking the ending signal is needed.
template< class T> 
struct info_pair
{
  std::vector<T> first;
  bool second;

  info_pair() = default;
  info_pair(const info_pair& x) = default;
  info_pair& operator=(const info_pair& x) = default;
  ~info_pair() = default;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & this->first;
    ar & this->second;
  }

};

template<class T>
info_pair<T> make_infopair(const std::vector<T>& x, bool finished){
  info_pair<T> pair;
  pair.first=x;
  pair.second=finished;
  return pair;
};

template< class T> 
class cache
{
  public:
  typedef T valuetype;

  cache(){}
  cache(const cache& cache0): 
    finishedread    ( cache0.finishedread), 
    current_pointer ( cache0.current_pointer),
    begin_pointer   ( cache0.begin_pointer),
    end_pointer     ( cache0.end_pointer),
    cache_size      ( cache0.cache_size),
    inputfile       ( cache0.inputfile)
  {
  }
  cache& operator=(const cache& cache0){
    finishedread    = cache0.finishedread;
    current_pointer = cache0.current_pointer;
    begin_pointer   = cache0.begin_pointer;
    end_pointer     = cache0.end_pointer;
    cache_size      = cache0.cache_size;
    inputfile       = cache0.inputfile;

  }
    
  //seek_set of file makes it possible to seperate one file into different parts.
  //Now, it is not implented, because there is no option to detect the end of pieces of one file.
  cache(char* filename, size_t buffer_size, size_t seek_set_of_file =0) {
    cache_size=buffer_size;
    begin_pointer=(valuetype*) malloc(sizeof(valuetype)*cache_size);
    if(begin_pointer==NULL){
      cout << "cannot allocate memory for cache, abort\n";
      abort();
    }

    inputfile=fopen(filename,"rb"); 
    if(inputfile==NULL){
      current_pointer = NULL;
      cout << "cannot open :"<<filename<<endl;
      return;
    }
    //fseek(inputfile,sizeof(valuetype)*seek_set_of_file,SEEK_SET);
    size_t real_size=fread(begin_pointer,sizeof(valuetype),cache_size,inputfile);
    if(real_size ==0) {
      fclose(inputfile); 
      current_pointer=NULL;
      return;
    }
    else if ( real_size< cache_size){
      fclose(inputfile);
      finishedread=true;
    }
    end_pointer = begin_pointer+real_size;
    current_pointer = begin_pointer;

  }
  ~cache(){
    //FIXME
    //std::vector:push_back creat and delete objects serveral times. Why?
    //cout << "deleted"<<endl;
    //FIXME, if it is clear here, there are problems in the copy of cache.
    //free(begin_pointer);
    //fclose(inputfile);
  }
  void clear(){
    free(begin_pointer);
    fclose(inputfile);
  }
  bool forward(){
    if(current_pointer!=end_pointer-1){
      current_pointer++;
      return true;
    }
    else{
      if(finishedread) return false;
      size_t real_size=fread(begin_pointer,sizeof(valuetype),cache_size,inputfile);
      if(real_size ==0) 
      {
        fclose(inputfile); 
        return false;
      }
      if ( real_size< cache_size){
        fclose(inputfile);
        finishedread=true;
      }
      end_pointer = begin_pointer+real_size;
      current_pointer = begin_pointer;
      return true;

    }
  }
  valuetype* current_position(){
    return current_pointer;
  }
  valuetype& value(){
    return *current_pointer;
  }

    //size_t cache_size;
    int cache_size;
    bool finishedread;
    valuetype* current_pointer;
    valuetype* begin_pointer;
    valuetype* end_pointer;
    FILE* inputfile;
  private:

  

};

class batch_index{
  public:
  batch_index(){}
  batch_index(long index, long pos, int size, int Nmpi):
    element_index(index),position(pos),batch_size(size),mpi_number(Nmpi){}
  batch_index(const batch_index& x):
    element_index(x.element_index),position(x.position),batch_size(x.batch_size),mpi_number(x.mpi_number){}
  batch_index& operator=(const batch_index& x){
    element_index=x.element_index;
    position=x.position;
    batch_size=x.batch_size;
    mpi_number=x.mpi_number;
  }
  long element_index;
  long position;
  int batch_size;
  int mpi_number;
  bool operator<(const batch_index& x)const{
    if(element_index< x.element_index) return true;
    else if(element_index> x.element_index) return false;
    else{
      cout << " Warning: one nonspinbatch are calculated twice."<<endl;
      if(batch_size!= x.batch_size){
        cout << " And their size are different, abort."<<endl;
        abort();
      }
    }

    return true;

  }

  friend ostream& operator<<(ostream& os, const batch_index& x){
    os << x.element_index<<"\t\t"<<x.position<<"\t\t"<<x.batch_size<<"\t\t"<<x.mpi_number<<endl;
    return os;
  }

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & element_index;
    ar & position;
    ar & batch_size;
    ar & mpi_number;
  }
};

template<class T,class = typename std::enable_if<has_index<T>::value>::type> 
void partition_data(long number_of_data, char* inputfilename, char* outputfilename);

void mergefile(char* filename);

//outfilename and inputfilename can be the same. There is no conflict.
template<class T, bool amdzero=true > 
void externalsort(char* inputfilename, char* outputfilename, long partition_index);

void parallel_external_sort( char* filename);

}
}

#endif
