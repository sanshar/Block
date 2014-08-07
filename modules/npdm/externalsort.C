#include <boost/format.hpp>
#include <boost/range/algorithm.hpp>
//#include <mpi.h>
#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <utility>
#include <vector> 
#include <fstream>
#include <iostream>
#include <thread>
#include "externalsort.h"

namespace SpinAdapted{
namespace Npdm{

using namespace std;

template<class T,class > 
void partition_data(long number_of_data, char* inputfilename, char* outputfilename)
{
  boost::mpi::communicator world;
  long partition_index=number_of_data/world.size();
  std::cout << "begin partition\n";

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
    std::cout <<"processor "<<world.rank()<< " cannot open "<<inputfilename<<std::endl;
    abort();
  }
  T inputbuff[Buff_SIZE];
  bool send_end= false;
  long realsize=fread(inputbuff,sizeof(T),Buff_SIZE,inputfile);
  if(realsize==0) {
    fclose(inputfile);
    boost::mpi::request req[world.size()];
    for(int i=0; i< world.size();i++)
      req[i]=world.isend(i,0,make_infopair(send_buff[i],true));
    send_end= true;
  }
  bool finished= realsize==Buff_SIZE? false: true;// finished means this is the last piece of data.
 // std::cout <<"partition_index "<< partition_index<<std::endl;
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
        if(inputbuff[i].index>= partition_index*world.size()) std::cout << " too big index"<<std::endl;
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

void mergefile(char* filename){
  boost::mpi::communicator world;
  if(world.rank()!=0){
    FILE* inputfile=fopen(filename,"rb");
    char buff[merge_buff_size];
    std::vector<char> sendvector(merge_buff_size);
    for(;;){
      long realsize=fread(buff,sizeof(char),merge_buff_size,inputfile);
      if( realsize==0){
        sendvector.resize(0);
        world.send(0,0,sendvector);
        break;
      }
      sendvector.resize(realsize);
      std::move(buff,buff+realsize,sendvector.begin());

      //sendvector.assign(buff,buff+realsize);
      world.send(0,0,sendvector);
      if(realsize!=merge_buff_size) break;
    }
  fclose(inputfile);
  boost::filesystem::remove(filename);
  }
  else{
    FILE* outputfile=fopen(filename,"ab");
    std::vector<char> buff;
    buff.reserve(merge_buff_size);
    for(int i=1; i< world.size();i++){
      for(;;){
        world.recv(i,0,buff);
        if(buff.size()==0) break;
        fwrite(&(buff[0]),sizeof(char),buff.size(),outputfile);
        if(buff.size()< merge_buff_size) break;
      }

    }
  fclose(outputfile);
  }


}

//outfilename and inputfilename can be the same. There is no conflict.
template<class T, bool addzero > 
void externalsort(char* inputfilename, char* outputfilename, long number_of_data){
  boost::mpi::communicator world;
  long partition_index=number_of_data/world.size();
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
    //  std::cout << onepiece[i]<<std::endl;
    //data.push_back(onepiece[0]);
    fwrite(onepiece,sizeof(T),read_size,outputfile);
    piecesnumber++;
    fclose(outputfile);
    //std::cout << "piecenumber: "<<piecesnumber<<std::endl;

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
  //std::cout <<"cachesize"<< filecache.size()<<std::endl;
  char sortedfile[100];

  typename T::indextype previous_index= partition_index*world.rank()-1;

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
          if(world.rank()== world.size()-1)
            n =number_of_data-current_value.index;
          else
            n =(partition_index*(world.rank()+1)-current_value.index);
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

void parallel_external_sort( char* filename){

  //int provided;
  //MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  //if (provided < MPI_THREAD_MULTIPLE)
  //{
  //      printf("ERROR: The MPI library does not have full thread support\n");
  //          MPI_Abort(MPI_COMM_WORLD, 1);
  //}
  boost::mpi::communicator world;
  char inputfilename[200];
  sprintf (inputfilename, "%s%d.%d.%d%s", filename,0,0,world.rank(),".bin");
  char tmpfilename[200];
  sprintf (tmpfilename, "%s%d.%d.%d%s", filename,0,0,world.rank(),".tmp");
  char sortedfilename[200];
  sprintf (sortedfilename, "%s%d.%d.%d%s", filename,0,0,world.rank(),".sorted");
  partition_data<index_element>((long)pow(6,6),inputfilename,tmpfilename);
  std::cout << "finished partition_data" << world.rank()<<std::endl;
  externalsort<index_element>(tmpfilename,sortedfilename,(long)pow(6,6));
  mergefile(sortedfilename);

}

#if 0

//FIXME
//multithread and mpi are combined to do parallel external sort. 
//However, this version of openmpi does not support MPI_THREAD_MULTIPLE. They could run, but they had segmentation faults.
//And I did not test it on another version of openmpi.
//I just put them here, wishing that they should be useful later.

template< class T,class = typename std::enable_if<has_index<T>::value>::type > 
void partition_data_sendthread(char* inputfilename, long partition_index){
  mpi::communicator world;

  cout << world.size()<<endl;
  std::vector<std::vector<T>> send_buff;
  send_buff.resize(world.size());
  //std::vector<info_pair<T>> send_buff;
  FILE* inputfile = fopen(inputfilename,"rb");
  if(inputfile==NULL){
    cout << " cannot open "<<inputfilename<<endl;
    abort();
  }
  T inputbuff[Buff_SIZE];
  long realsize=fread(inputbuff,sizeof(T),Buff_SIZE,inputfile);
  if(realsize==0) {
    fclose(inputfile);
    mpi::request req[world.size()];
    for(int i=0; i< world.size();i++)
      req[i]=world.isend(i,0,make_infopair(send_buff[i],true));
    return;
  }
  bool finished= realsize==Buff_SIZE? false: true;// finished means this is the last piece of data.
 // cout <<"partition_index "<< partition_index<<endl;
  for(;;){
    for(int i=0; i< realsize; i++)
    {
//      cout <<"send elements: "<<i<<"  "<< inputbuff[i]<<endl;

      if(inputbuff[i].index>= partition_index*world.size()) cout << " too big index"<<endl;
      send_buff[(int) floor(inputbuff[i].index/partition_index)].push_back(inputbuff[i]);
    }

    if(!finished){
      realsize=fread(inputbuff,sizeof(T),Buff_SIZE,inputfile);
      if(realsize==0){

      finished=true;
      }
    }

    std::vector<mpi::request> req(world.size());
    for(int i=0; i< world.size();i++)
      req[i]=world.isend(i,0,make_infopair(send_buff[i],finished));
      //world.send(i,0,make_infopair(send_buff[i],finished));

    for(int i=0; i< world.size();i++)
      req[i].wait();
    //mpi::wait_all(req,req+world.size());
    if(finished) break;


    for(int i=0; i< world.size();i++)
      send_buff[i].clear();
    if(realsize < Buff_SIZE){

    finished= true;
    }
  }
  fclose(inputfile);

  //for(int i=0; i< world.size();i++)
  //  world.isend(i,1,false)


}

template< class T,class = typename std::enable_if<has_index<T>::value>::type > 
void partition_data_recvthread(char* outputfilename){
  mpi::communicator world;
  cout << world.size()<<endl;
  
  std::vector<info_pair<T>> recv_buff;
  recv_buff.resize(world.size());
  //std::vector<T> output_buff;
  FILE* outputfile = fopen(outputfilename,"wb");
  setvbuf(outputfile,NULL,_IOFBF,Buff_SIZE*3);
  //std::vector<bool> openrank(world.size,true)
  std::vector<int> activeworld; activeworld.resize(world.size());
  for(int i=0; i< world.size();i++)
  {
    activeworld[i]=i;
  }


  for(;;)
  {
    if(activeworld.size() == 0) break;
    std::vector<mpi::request> req(world.size());
    //mpi::wait_all(req,req+activeworld.size());
    for(int i=0; i< activeworld.size();)
    {
      //req[i]=world.irecv(activeworld[i],0,recv_buff[activeworld[i]]);
      mpi::request req=world.irecv(activeworld[i],0,recv_buff[activeworld[i]]);
      
      //world.recv(activeworld[i],0,recv_buff[activeworld[i]]);
      //std::move(recv_buff[i].first.begin(),recv_buff[i].first.end(),std::back_inserter(output_buff));
      //fwrite(&output_buff[0],sizeof(T),output_buff.size(),outputfile);
      //FIXME
      //Different receiver should be more independent.
      //req[i].wait();
      req.wait();
      if(recv_buff[activeworld[i]].first.size()!=0)
      {
        fwrite(&(recv_buff[activeworld[i]].first[0]),sizeof(T),recv_buff[activeworld[i]].first.size(),outputfile);
        fflush(outputfile);
        recv_buff[activeworld[i]].first.clear();

      }
    world.barrier();
      if(recv_buff[activeworld[i]].second == true)
      {
        activeworld.erase(activeworld.begin()+i);
        //recv_buff.erase(recv_buff.begin()+i);
      }
      else 
      {
        i++;
      }
    }
  }

  fclose(outputfile);
  
//
//  for(;;){
//    for(int i=0; i< activeworld.size(); ){
//      if(openrank[tmp_activeworld[i]]== false)
//        activeworld.erase(activeworld.begin()+i); //      else //        i++;
//    }
//    if(activeworl.size()== 0) break;
//
//    mpi::request recreq[activeworld.size()];
//    for(int i=0; i< activeworld.size();i++){
//        req[i]=world.irecv(i,0,recv_buff[i]);
//    }
//    mpi::wait_all(recreq,recreq+activeworld.size());
//    for(int i=0; i< activeworld.size();i++)
//      std::move(recv_buff[i].begin(),recv_buff.end(),std::back_inserter(output_buff));
//    fwrite(&outputbuff[0],sizeof(T),output_buff.size(),outputfile);
//  }
}

template< class T,class = typename std::enable_if<has_index<T>::value>::type > 
void partition_data_multithread(long number_of_data, char* inputfilename, char* outputfilename)
{
  mpi::communicator world;
  long partition_index = number_of_data/world.size();
  cout << "begin partition\n";
  std::thread datasend(partition_data_sendthread<T>,inputfilename,partition_index);
  std::thread datarecv(partition_data_recvthread<T>,outputfilename);

  datasend.join();
  datarecv.join();

  world.barrier();

}

#endif

}
}
