#ifndef EXTERNAL_SORT_H
#define EXTERNAL_SORT_H

#include "global.h"
#include <boost/format.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/filesystem.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <vector> 
#include <fstream>
#include <iostream>
#define Buff_SIZE 1024*1024*8
#define merge_buff_size (Buff_SIZE*100) // during the merge of file, only one processor works. More momery available.

namespace SpinAdapted{
namespace Npdm{
namespace Sortpdm{

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
      if(index < onepiece2.index)
	{ return true; }
      else if(index > onepiece2.index)
	{ return false; }
      else // this means the indices are the same 
	//      if(index == onepiece2.index)
	{
	  pout << "there are elements with same index"<<endl;
	  if(element !=onepiece2.element) 
	    {
	      pout << " and their values are different, abort"<<endl;
	      abort();
	    }
	  return false; 
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
    return *this;
  }
    
  //seek_set of file makes it possible to seperate one file into different parts.
  //Now, it is not implented, because there is no option to detect the end of pieces of one file.
  cache(char* filename, size_t buffer_size, size_t seek_set_of_file =0) {
    cache_size=buffer_size;
    begin_pointer=(valuetype*) malloc(sizeof(valuetype)*cache_size);
    if(begin_pointer==NULL){
      pout << "cannot allocate memory for cache, abort\n";
      abort();
    }

    inputfile=fopen(filename,"rb"); 
    if(inputfile==NULL){
      current_pointer = NULL;
      pout << "cannot open :"<<filename<<endl;
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
    //pout << "deleted"<<endl;
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
    return *this;
  }
  long element_index;
  long position;
  int batch_size;
  int mpi_number;
  bool operator<(const batch_index& x)const{
    if(element_index< x.element_index) return true;
    else if(element_index> x.element_index) return false;
    else{
      pout << " Warning: one nonspinbatch are calculated twice."<<endl;
      if(batch_size!= x.batch_size){
        pout << " And their size are different, abort."<<endl;
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

template<class T, bool addzero = true > 
void externalsort(char* inputfilename, char* outputfilename, long number_of_data){
#ifndef SERIAL
  boost::mpi::communicator world;
  long partition_index=number_of_data/world.size();
#endif
  FILE* datafile;
  FILE* outputfile;
  long read_size;
  bool read_success =true;
  datafile=fopen(inputfilename,"rb");
  int piecesnumber=0;
  // sort within different pieces;
  while(read_success){
    char tmpfile[100];
#ifndef SERIAL
    sprintf (tmpfile, "%s%d%s%d","tmpsortingfile", piecesnumber,".",world.rank());
#else
    sprintf (tmpfile, "%s%d","tmpsortingfile", piecesnumber);
#endif
    outputfile=fopen(tmpfile,"wb");
    //index_element onepiece[Buff_SIZE];
    T onepiece[Buff_SIZE];
    read_size=fread(onepiece,sizeof(T),Buff_SIZE,datafile);
    if(read_size == 0) break;
    else if( read_size != Buff_SIZE) read_success =false;

    std::vector<index_element> sortedpiece(onepiece,onepiece+read_size);
    boost::sort(sortedpiece);
    std::copy(sortedpiece.begin(),sortedpiece.end(),onepiece);
  //  loaddata >> onepiece;
    //file >> onepiece.index;
    //file >> onepiece.element;
    //for(int i =0 ; i< Buff_SIZE;i++)
    //  pout << onepiece[i]<<std::endl;
    //data.push_back(onepiece[0]);
    fwrite(onepiece,sizeof(T),read_size,outputfile);
    piecesnumber++;
    fclose(outputfile);
    //pout << "piecenumber: "<<piecesnumber<<std::endl;

  };
  fclose(datafile);
  // Merge these pieces;
  
  std::vector<cache<T>> filecache;
  char tmpfile[100];
  for(int i =0; i< piecesnumber; i++){
#ifndef SERIAL
    sprintf (tmpfile, "%s%d%s%d","tmpsortingfile", i,".",world.rank());
#else
    sprintf (tmpfile, "%s%d","tmpsortingfile", i);
#endif

    //cache<T> tmpcache( tmpfile, Buff_SIZE/piecesnumber);
    cache<T> tmpcache( tmpfile, Buff_SIZE/20);
    //cache<T> tmpcache( tmpfile, Buff_SIZE);
    if(tmpcache.current_position()!= NULL)
      filecache.push_back(tmpcache);
      //filecache.emplace_back(tmpcache);

  }
  //pout <<"cachesize"<< filecache.size()<<std::endl;
  char sortedfile[100];

#ifndef SERIAL
  typename T::indextype previous_index= partition_index*world.rank()-1;
#else
  typename T::indextype previous_index= -1;
#endif
  outputfile = fopen(outputfilename,"wb");
  setvbuf(outputfile,NULL,_IOFBF,3*Buff_SIZE);

  for(;;){
    // select the smallest one in the current positions of different caches.
    int smallest = 0;
    for(int i=1 ; i< filecache.size(); i++){
      if(filecache[i].value() < filecache[smallest].value())
        smallest =i;
    }
    T current_value = filecache[smallest].value();
    //outputbuff[outputbuff_position++]=filecache[smallest].value();

    if(addzero){ 
      int n=(current_value.index - previous_index);
      for(int i=1; i<n;i++){
        typename T::valuetype tmp=0.;
        fwrite(&tmp,sizeof(typename T::valuetype),1,outputfile);
      }
      fwrite(&(current_value.element),sizeof(typename T::valuetype),1,outputfile);
      previous_index=current_value.index;
    }
    else
      fwrite(&(filecache[smallest].value()),sizeof(T),1,outputfile);

    if(!filecache[smallest].forward()){
      filecache.erase(filecache.begin()+smallest);
      if (filecache.size()==0) {
        if(addzero){
          int n;
#ifndef SERIAL
          if(world.rank()== world.size()-1)
            n =number_of_data-current_value.index;
          else
            n =(partition_index*(world.rank()+1)-current_value.index);
#else
          n =number_of_data-current_value.index;
#endif
          
          for(int i =1; i<n; i++)
          {
            typename T::valuetype tmp=0.;
            fwrite(&tmp,sizeof(typename T::valuetype),1,outputfile);
          }
        }
      break;
      }
    }
  }


  fclose(outputfile);
  for(int i=0; i< piecesnumber;i++)
  {
    char tmpfile[100];
#ifndef SERIAL
    sprintf (tmpfile, "%s%d%s%d","tmpsortingfile", i,".",world.rank());
#else
    sprintf (tmpfile, "%s%d","tmpsortingfile", i);
#endif

    boost::filesystem::remove(tmpfile);
  }
}

#ifndef SERIAL
template<class T,class= typename std::enable_if<has_index<T>::value>::type>
void partition_data(long number_of_data, char* inputfilename, char* outputfilename)
{
  boost::mpi::communicator world;
  long partition_index=number_of_data/world.size();
  //std::pout << "begin partition\n";

  std::vector<info_pair<T>> recv_buff(world.size());
  std::vector<std::vector<T>> send_buff(world.size());
  //std::vector<T> output_buff;
  FILE* outputfile = fopen(outputfilename,"wb");
  setvbuf(outputfile,NULL,_IOFBF,Buff_SIZE*3);
  //std::vector<bool> openrank(world.size,true)
  std::vector<int> activeworld; activeworld.resize(world.size());
  for(int i=0; i< world.size();i++)
  {
    activeworld[i]=i;
  }


  //std::vector<info_pair<T>> send_buff;
  FILE* inputfile = fopen(inputfilename,"rb");
  if(inputfile==NULL){
    pout <<"processor "<<world.rank()<< " cannot open "<<inputfilename<<std::endl;
    abort();
  }
  T inputbuff[Buff_SIZE];
  bool send_end= false;
  long realsize=fread(inputbuff,sizeof(T),Buff_SIZE,inputfile);
  if(realsize==0) {
    fclose(inputfile);
    std::vector<boost::mpi::request> req(world.size());
    for(int i=0; i< world.size();i++)
      req[i]=world.isend(i,0,make_infopair(send_buff[i],true));
    send_end= true;
  }
  bool finished= realsize==Buff_SIZE? false: true;// finished means this is the last piece of data.
 // pout <<"partition_index "<< partition_index<<std::endl;
  std::vector<boost::mpi::request> sendreq(world.size());
  std::vector<boost::mpi::request> recvreq(activeworld.size());
  for(;;){
   //send data
    if(!send_end)
    {
      if(finished) send_end=true;
  for(int i=0; i< world.size();i++)
    sendreq[i].wait();
      for(int i=0; i< world.size();i++)
        send_buff[i].clear();
      for(int i=0; i< realsize; i++)
      {
        if(inputbuff[i].index>= partition_index*world.size()) pout << " too big index"<<std::endl;
        send_buff[(int) floor(inputbuff[i].index/partition_index)].push_back(inputbuff[i]);
      }

      if(!finished){
        realsize=fread(inputbuff,sizeof(T),Buff_SIZE,inputfile);
        if(realsize==0){
          finished=true;
          send_end=true;
        }
      }

      for(int i=0; i< world.size();i++)
        sendreq[i]=world.isend(i,0,make_infopair(send_buff[i],finished));

      if(realsize < Buff_SIZE){

      finished= true;
      }
    }
   //receive data
    {
      for(int i=0; i< activeworld.size();i++)
      {
        recvreq[activeworld[i]]=world.irecv(activeworld[i],0,recv_buff[activeworld[i]]);
      }
        //FIXME
        //Different receiver should be more independent.
        //req[i].wait();
      for(int i=0; i< activeworld.size();){
        recvreq[activeworld[i]].wait();
        if(recv_buff[activeworld[i]].first.size()!=0)
        {
          fwrite(&(recv_buff[activeworld[i]].first[0]),sizeof(T),recv_buff[activeworld[i]].first.size(),outputfile);
          fflush(outputfile);
          recv_buff[activeworld[i]].first.clear();

        }
        if(recv_buff[activeworld[i]].second == true)
        {
          activeworld.erase(activeworld.begin()+i);
        }
        else 
        {
          i++;
        }
      }
    }
  if(activeworld.size() == 0)  break;
 }

  
  for(int i=0; i< world.size();i++)
    sendreq[i].wait();

  fclose(inputfile);
  fclose(outputfile);
  
}

void mergefile(char* filename);

void parallel_external_sort( char* filename);

#endif

}
}
}

#endif
