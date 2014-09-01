#include <boost/format.hpp>
#include <boost/range/algorithm.hpp>
#include "global.h"
//#include <boost/mpi.hpp>
//#include <boost/mpi/environment.hpp>
//#include <boost/mpi/communicator.hpp>
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
namespace Sortpdm{



#ifndef SERIAL 

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
  //pout << "finished partition_data" << world.rank()<<std::endl;
  externalsort<index_element>(tmpfilename,sortedfilename,(long)pow(6,6));
  mergefile(sortedfilename);

}

#endif

#if 0

//FIXME
//multithread and mpi are combined to do parallel external sort. 
//However, this version of openmpi does not support MPI_THREAD_MULTIPLE. They could run, but they had segmentation faults.
//And I did not test it on another version of openmpi.
//I just put them here, wishing that they should be useful later.

template< class T,class = typename std::enable_if<has_index<T>::value>::type > 
void partition_data_sendthread(char* inputfilename, long partition_index){
  mpi::communicator world;

  pout << world.size()<<endl;
  std::vector<std::vector<T>> send_buff;
  send_buff.resize(world.size());
  //std::vector<info_pair<T>> send_buff;
  FILE* inputfile = fopen(inputfilename,"rb");
  if(inputfile==NULL){
    pout << " cannot open "<<inputfilename<<endl;
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
 // pout <<"partition_index "<< partition_index<<endl;
  for(;;){
    for(int i=0; i< realsize; i++)
    {
//      pout <<"send elements: "<<i<<"  "<< inputbuff[i]<<endl;

      if(inputbuff[i].index>= partition_index*world.size()) pout << " too big index"<<endl;
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
  pout << world.size()<<endl;
  
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
//}

template< class T,class = typename std::enable_if<has_index<T>::value>::type > 
void partition_data_multithread(long number_of_data, char* inputfilename, char* outputfilename)
{
  mpi::communicator world;
  long partition_index = number_of_data/world.size();
  pout << "begin partition\n";
  std::thread datasend(partition_data_sendthread<T>,inputfilename,partition_index);
  std::thread datarecv(partition_data_recvthread<T>,outputfilename);

  datasend.join();
  datarecv.join();

  world.barrier();

}

#endif

}
}
}
