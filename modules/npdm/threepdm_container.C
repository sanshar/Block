/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include <boost/format.hpp>
#include "threepdm_container.h"
#include "npdm_permutations.h"
#include <boost/filesystem.hpp>
#include <boost/range/algorithm.hpp>

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Threepdm_container::Threepdm_container( int sites )
{
  if ( dmrginp.store_spinpdm() ) {
    if(dmrginp.spinAdapted())
      threepdm.resize(2*sites,2*sites,2*sites,2*sites,2*sites,2*sites);
    else
      threepdm.resize(sites,sites,sites,sites,sites,sites);
    threepdm.Clear();
  } 
  if ( !dmrginp.spatpdm_disk_dump() ) {
    if(dmrginp.spinAdapted())
      spatial_threepdm.resize(sites,sites,sites,sites,sites,sites);
    else
      spatial_threepdm.resize(sites/2,sites/2,sites/2,sites/2,sites/2,sites/2);
    spatial_threepdm.Clear();
  }

  else{
    elements_stride_.resize(6);
    for(int i=0; i<6; i++ )
      elements_stride_[i]=pow(sites,5-i);

    char file[5000];
    sprintf (file, "%s%s%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",mpigetrank(),".tmp");
    //std::ofstream spatpdm_disk(file, std::ios::binary);
    spatpdm_disk=fopen(file,"w");
    setvbuf(spatpdm_disk,NULL,_IOFBF,1024*1024*32);
    if(mpigetrank()==0)
      fprintf(spatpdm_disk,"%d\n",sites);
    //32M buffer;
    //spatpdm_disk.open(file, std::ios::binary);
    //spatpdm_disk.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_npdms(const int& i, const int& j)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  world.barrier();
#endif
  Timer timer;
  if ( dmrginp.store_spinpdm() ) {
    accumulate_npdm();
    save_npdm_binary(i, j);
    save_npdm_text(i, j);
  }
  if ( dmrginp.spatpdm_disk_dump() ) {
    merge_diskfile(i,j);
  }
  else
  {
    accumulate_spatial_npdm();
    save_spatial_npdm_text(i, j);
    save_spatial_npdm_binary(i, j);
  }

#ifndef SERIAL
  world.barrier();
#endif
  ecpu = timer.elapsedcputime();ewall=timer.elapsedwalltime();
  p3out << "3PDM save full array time " << ewall << " " << ecpu << endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/threepdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << threepdm.dim1() << endl;

    double trace = 0.0;
    for(int i=0; i<threepdm.dim1(); ++i)
      for(int j=0; j<threepdm.dim2(); ++j)
        for(int k=0; k<threepdm.dim3(); ++k)
          for(int l=0; l<threepdm.dim4(); ++l)
            for(int m=0; m<threepdm.dim5(); ++m)
              for(int n=0; n<threepdm.dim6(); ++n) {
                if ( abs(threepdm(i,j,k,l,m,n)) > NUMERICAL_ZERO ) {
                  ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % threepdm(i,j,k,l,m,n);
                  if ( (i==n) && (j==m) && (k==l) ) trace += threepdm(i,j,k,l,m,n);
                }
              }
    ofs.close();
    pout << "Spin-orbital 3PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_spatial_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << spatial_threepdm.dim1() << endl;

    double trace = 0.0;
    for(int i=0; i<spatial_threepdm.dim1(); ++i)
      for(int j=0; j<spatial_threepdm.dim2(); ++j)
        for(int k=0; k<spatial_threepdm.dim3(); ++k)
          for(int l=0; l<spatial_threepdm.dim4(); ++l)
            for(int m=0; m<spatial_threepdm.dim5(); ++m)
              for(int n=0; n<spatial_threepdm.dim6(); ++n) {
                if ( abs(spatial_threepdm(i,j,k,l,m,n)) > NUMERICAL_ZERO ) {
                  ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % spatial_threepdm(i,j,k,l,m,n);
                  if ( (i==n) && (j==m) && (k==l) ) trace += spatial_threepdm(i,j,k,l,m,n);
                }
              }
    ofs.close();
    pout << "Spatial      3PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_npdm_binary(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/threepdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << threepdm;
    ofs.close();
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

/*
void Threepdm_container::external_sort_index(const int &i, const int &j)
{

  boost::sort(nonspin_batch);
#ifndef SERIAL
  boost::mpi::communicator world;
  if(mpigetrank() != 0){
      world.send(0,0, nonspin_batch);
  }
  else{
    //Store index from different processors on disk of root node.
    for(int p=0; p< world.size();p++){
      if(p!=0){
        world.recv(p,0, nonspin_batch);
      }
      char file[5000];
      //batch_index tmpbuffer[1000000];
      //std::copy(nonspin_batch.begin(),nonspin_batch.end(),tmpbuffer);
      sprintf (file, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm_index.", i, j,p,".bin");
      FILE* inputfile=fopen(file,"wb");
      fwrite(&nonspin_batch[0],sizeof(Sortpdm::batch_index),nonspin_batch.size(),inputfile);
      fclose(inputfile);
      nonspin_batch.clear();
      if(world.size() == 0) return;
    }
    //external sort nonspin_batch
    //TODO
    //It is not parallel.
    char outfilename[5000];
    sprintf (outfilename, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm_index.",i,j,".bin");
    FILE* outputfile = fopen(outfilename,"wb");
    long sorting_buff= 1024*1024*(32/world.size());
    //For batch_index, the sorting buff is about 96M/world.size();
    std::vector<Sortpdm::cache<Sortpdm::batch_index>> filecache;
    for(int p=0; p< world.size();p++){
      char file[5000];
      sprintf (file, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm_index.", i, j,p,".bin");
      Sortpdm::cache<Sortpdm::batch_index> tmpcache( file, sorting_buff);
      filecache.push_back(tmpcache);
    }
    long outputbuffsize=sorting_buff*4;
    long outputbuff_position = 0;
    Sortpdm::batch_index outputbuff[outputbuffsize];
    for(;;){
      // select the smallest one in the current positions of different caches.
      int smallest = 0;
      for(int i=1 ; i< filecache.size(); i++){
        if(filecache[i].value() < filecache[smallest].value())
          smallest =i;
      }
      outputbuff[outputbuff_position++]=filecache[smallest].value();

      if(!filecache[smallest].forward()){
        filecache[smallest].clear();
        filecache.erase(filecache.begin()+smallest);
        if (filecache.size()==0) {
        fwrite(outputbuff,sizeof(Sortpdm::batch_index),outputbuff_position,outputfile);
        break;
        }
      }

      if(outputbuff_position == outputbuffsize){
        fwrite(outputbuff,sizeof(Sortpdm::batch_index),outputbuffsize,outputfile);
        outputbuff_position=0;
      }

      }
      fclose(outputfile);
      //Finish external sort of index.

      //Clean up.
      for(int p=0; p< world.size();p++){
        char file[5000];
        sprintf (file, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm_index.", i, j,p,".bin");
        boost::filesystem::remove(file);
      }
    }
#else
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/threepdm_index.", i, j,".bin");
    FILE* outfile=fopen(file,"wb");
    fwrite(&nonspin_batch[0],sizeof(Sortpdm::batch_index),nonspin_batch.size(),outfile);
    nonspin_batch.clear();
#endif
}
*/

void Threepdm_container::merge_diskfile(const int &i, const int &j)
{
  fclose(spatpdm_disk);
  if(mpigetrank()==0)
  {
    char oldfile[5000];
    sprintf (oldfile, "%s%s%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",mpigetrank(),".tmp");
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.", i, j,".txt");
    std::rename(oldfile,file);
  }

#ifndef SERIAL
  boost::mpi::communicator world;
  world.barrier();

  int chunksize = 100*1024*1024; // 100M buffer for reading.
  char* buffer = new char[chunksize];

  if(mpigetrank()==0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.", i, j,".txt");
    FILE* merged_file=fopen(file,"ab");
    for (int rank = 1; rank<mpigetsize();rank++)
    {
      long size;
      while(true)
      {
        MPI_Recv(&size, 1, MPI_LONG, rank, 0, world,  MPI_STATUS_IGNORE);
        if(size==-1)
          break;
        MPI_Recv(buffer, size, MPI_CHAR, rank, 0, world,  MPI_STATUS_IGNORE);
        fwrite(buffer,1,size,merged_file);
      }


    }
    fclose(merged_file);

  }
  else
    {
    char file[5000];
    sprintf (file, "%s%s%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",mpigetrank(),".tmp");
    //std::ofstream spatpdm_disk(file, std::ios::binary);
    FILE* local_file=fopen(file,"rb");
    fseek(local_file,0,SEEK_END);
    long file_size = ftell(local_file);
    fseek(local_file,0,SEEK_SET);
    do
    {
      long size=fread(buffer, 1, chunksize, local_file);
      MPI_Send(&size, 1, MPI_LONG, 0, 0, world);
      MPI_Send(buffer, size, MPI_CHAR, 0, 0, world);
    }
    while(ftell(local_file)<file_size);

    long size=-1;

    MPI_Send(&size, 1, MPI_LONG, 0, 0, world);
    fclose(local_file);
    std::remove(file);
    }
  delete[] buffer;
#endif

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_spatial_npdm_binary(const int &i, const int &j)
{ 
  if(!dmrginp.spatpdm_disk_dump())
  {
    if( mpigetrank() == 0)
    {
      char file[5000];
      sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",i,j,".bin");
#ifndef MOLCAS
      std::ofstream ofs(file, std::ios::binary);
      boost::archive::binary_oarchive save(ofs);
      save << spatial_threepdm;
      ofs.close();
#else
      FILE* f = fopen(file,"wb");
      fwrite(spatial_threepdm.data(),sizeof(double),spatial_threepdm.size(),f);
      fclose(f);
#endif
    }
  }
  /*
  else{
    fclose(spatpdm_disk);
    //When spatpdm_disk is opened, the state numbers, i and j, are not known. 
    char oldfile[5000];
    sprintf (oldfile, "%s%s%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",mpigetrank(),".tmp");
    char newfile[5000];
    sprintf (newfile, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",i,j,mpigetrank(),".bin");
    boost::filesystem::rename(oldfile,newfile);
    if(dmrginp.pdm_unsorted())
      external_sort_index(i,j);
    
    else{
      char file[5000];
      sprintf (file, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",i,j,mpigetrank(),".bin");
      char finalfile[5000];
      sprintf (finalfile, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",i,j,".bin");
#ifndef SERIAL
      boost::mpi::communicator world;
      world.barrier();
      Timer timer1;
      char tmpfile[5000];
      sprintf (tmpfile, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",i,j,mpigetrank(),".tmp");
      char sortedfile[5000];
      sprintf (sortedfile, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",i,j,mpigetrank(),".bin");
      Sortpdm::partition_data<Sortpdm::index_element>((long)pow(dmrginp.last_site(),6),file,tmpfile);
      //TODO
      //tmpfile and sortedfile can be the same file. Because tmpfile is divided into many small files. It can be overwritten by sortedfile.
      //However, when they are different, it is a little faster. Maybe compiler can do some optimizations. 
      Sortpdm::externalsort<Sortpdm::index_element>(tmpfile,sortedfile,(long)pow(dmrginp.last_site(),6));
      world.barrier();
      ecpu = timer1.elapsedcputime();ewall=timer1.elapsedwalltime();
      p3out << "3PDM parallel external sort time " << ewall << " " << ecpu << endl;
      Timer timer;
      Sortpdm::mergefile(sortedfile);
      world.barrier();
      if(mpigetrank()==0) boost::filesystem::rename(sortedfile,finalfile);
      boost::filesystem::remove(tmpfile);
      ecpu = timer.elapsedcputime();ewall=timer.elapsedwalltime();
      p3out << "3PDM merge sorted file time " << ewall << " " << ecpu << endl;
#else
      Timer timer2;
      Sortpdm::externalsort<Sortpdm::index_element>(file,finalfile,(long)pow(dmrginp.last_site(),6));
      boost::filesystem::remove(file);
      ecpu = timer2.elapsedcputime();ewall=timer2.elapsedwalltime();
      p3out << "3PDM external sort time " << ewall << " " << ecpu << endl;
#endif
    }
  }
  */
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::load_npdm_binary(const int &i, const int &j) { abort(); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::accumulate_npdm()
{
#ifndef SERIAL
  mpi::communicator world;
#ifdef NDEBUG
  if( mpigetrank() == 0)
    MPI_Reduce(MPI_IN_PLACE, spatial_threepdm.data(),  spatial_threepdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
  else
    MPI_Reduce(spatial_threepdm.data(), spatial_threepdm.data(),  spatial_threepdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
#else
  array_6d<double> tmp_recv;
  if( mpigetrank() == 0)
  {
    for(int p=1; p<world.size(); ++p) {
      world.recv(p, p, tmp_recv);
      for(int i=0; i<threepdm.dim1(); ++i)
        for(int j=0; j<threepdm.dim2(); ++j)
          for(int k=0; k<threepdm.dim3(); ++k)
            for(int l=0; l<threepdm.dim4(); ++l)
              for(int m=0; m<threepdm.dim5(); ++m)
                for(int n=0; n<threepdm.dim6(); ++n) {
                  if ( abs(tmp_recv(i,j,k,l,m,n)) > NUMERICAL_ZERO ) {
                    // Test if any duplicate elements built on different processors
                    assert ( abs(threepdm(i,j,k,l,m,n)) < NUMERICAL_ZERO );
                    threepdm(i,j,k,l,m,n) = tmp_recv(i,j,k,l,m,n);
                  }
                }
    }
  }
  else 
  {
    world.send(0, mpigetrank(), threepdm);
  }
#endif
#endif
} 
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::accumulate_spatial_npdm()
{
#ifndef SERIAL
  mpi::communicator world;
#ifdef NDEBUG
  if( mpigetrank() == 0)
    MPI_Reduce(MPI_IN_PLACE, spatial_threepdm.data(),  spatial_threepdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
  else
    MPI_Reduce(spatial_threepdm.data(), spatial_threepdm.data(),  spatial_threepdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
#else
  array_6d<double> tmp_recv;

  if( mpigetrank() == 0)
  {
    for(int p=1; p<world.size(); ++p) {
      world.recv(p, p, tmp_recv);
      for(int i=0; i<spatial_threepdm.dim1(); ++i)
        for(int j=0; j<spatial_threepdm.dim2(); ++j)
          for(int k=0; k<spatial_threepdm.dim3(); ++k)
            for(int l=0; l<spatial_threepdm.dim4(); ++l)
              for(int m=0; m<spatial_threepdm.dim5(); ++m)
                for(int n=0; n<spatial_threepdm.dim6(); ++n) {
                  if( abs(tmp_recv(i,j,k,l,m,n)) > NUMERICAL_ZERO ) {
                    // Test if any duplicate elements built on different processors
                    assert(abs(spatial_threepdm(i,j,k,l,m,n)) < NUMERICAL_ZERO );
                    spatial_threepdm(i,j,k,l,m,n) = tmp_recv(i,j,k,l,m,n);
                  }
                }
    }
  }
  else 
  {
    world.send(0, mpigetrank(), spatial_threepdm);
  }
#endif

#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::update_full_spin_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{

  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    double val = it->second;
    if ( abs(val) < NUMERICAL_ZERO ) continue;

    assert( (it->first).size() == 6 );
    int i0 = (it->first)[0];
    int j0 = (it->first)[1];
    int k0 = (it->first)[2];
    int l0 = (it->first)[3];
    int m0 = (it->first)[4];
    int n0 = (it->first)[5];
    int i = ro.at(i0/2)*2 + i0%2;
    int j = ro.at(j0/2)*2 + j0%2;
    int k = ro.at(k0/2)*2 + k0%2;
    int l = ro.at(l0/2)*2 + l0%2;
    int m = ro.at(m0/2)*2 + m0%2;
    int n = ro.at(n0/2)*2 + n0%2;

    //if ( abs(val) > 1e-8 ) {
    //  pout << "so-threepdm val: i,j,k,l,m,n = " 
    //       << i << "," << j << "," << k << "," << l << "," << m << "," << n
    //       << "\t\t" << val << endl;
    //}

    // Test for duplicates
    threepdm(i,j,k,l,m,n) = val;
  }

}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// This routine assumes that no spin-orbital indices are generated more than once

void Threepdm_container::update_full_spatial_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 6 );

    // Store significant elements only
    if ( abs(it->second) > NUMERICAL_ZERO ) {
      // Spin indices
      int i = (it->first)[0];
      int j = (it->first)[1];
      int k = (it->first)[2];
      int l = (it->first)[3];
      int m = (it->first)[4];
      int n = (it->first)[5];

//      if ( i%2 != n%2 ) continue;
//      if ( j%2 != m%2 ) continue;
//      if ( k%2 != l%2 ) continue;

//      spatial_threepdm( ro.at(i/2), ro.at(j/2), ro.at(k/2), ro.at(l/2), ro.at(m/2), ro.at(n/2) ) += it->second;
      spatial_threepdm( ro.at(i), ro.at(j), ro.at(k), ro.at(l), ro.at(m), ro.at(n) ) = it->second;
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

long Threepdm_container::oneindex_spin(const std::vector<int> & orbital_element_index)
{
  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  assert( orbital_element_index.size() == 6);
  long linearindex=0;
  for(int i=0; i< 6; i++)
    linearindex+=(ro.at(orbital_element_index[i]))*elements_stride_[i];
  return linearindex;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------


void Threepdm_container::dump_binary_to_disk(std::vector< std::pair< std::vector<int>, double > > & spin_batch)
{

  const std::vector<int>& ro = dmrginp.reorder_vector();
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 6 );

    // Store significant elements only
    if ( abs(it->second) > NUMERICAL_ZERO ) {
      std::vector<int> sites;
      sites.resize(6);
      for (int i=0;i<6;i++)
        sites[i] = ro.at((it->first)[i]);

      fwrite(sites.data(),sizeof(int),6,spatpdm_disk);
      fwrite(&(it->second),sizeof(double),1,spatpdm_disk);
    }
  }
}

void Threepdm_container::dump_text_to_disk(std::vector< std::pair< std::vector<int>, double > > & spin_batch)
{

  const std::vector<int>& ro = dmrginp.reorder_vector();
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 6 );

    // Store significant elements only
    if ( abs(it->second) > NUMERICAL_ZERO ) {

      int i = (it->first)[0];
      int j = (it->first)[1];
      int k = (it->first)[2];
      int l = (it->first)[3];
      int m = (it->first)[4];
      int n = (it->first)[5];

      fprintf(spatpdm_disk,"%d %d %d %d %d %d %20.14e\n", ro.at(i), ro.at(j), ro.at(k), ro.at(l), ro.at(m), ro.at(n) ,it->second);
    }
  }
}


  //TODO
  //The order of rdm elements does not seem useful. 
  //Now, we just need to dump element one by one.
  /*
void Threepdm_container::dump_to_disk(std::vector< std::pair< std::vector<int>, double > > & spin_batch)
{

  long spatpdm_disk_position= ftell(spatpdm_disk);
  std::map < long, double>  index_and_elements;
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 6 );

    // Store significant elements only
    if ( abs(it->second) > NUMERICAL_ZERO ) {
      long linearindex = oneindex_spin(it->first);
      index_and_elements.insert(std::pair<long,double>(linearindex,it->second));
    }
  }
  if(index_and_elements.size()==0) return;
  Sortpdm::index_element index_elements[72];
  assert(72>= index_and_elements.size());
  int i=0;
  for(auto it = index_and_elements.begin(); it!=index_and_elements.end();it++){
    index_elements[i].index=it->first;
    index_elements[i].element=it->second;
    i++;
  }
  fwrite(index_elements,sizeof(Sortpdm::index_element),index_and_elements.size(),spatpdm_disk);
  if(dmrginp.pdm_unsorted()){
  Sortpdm::batch_index onerecord(index_and_elements.begin()->first,spatpdm_disk_position/sizeof(Sortpdm::index_element),index_and_elements.size(),mpigetrank());
  nonspin_batch.push_back(onerecord);
  }


  //char file[5000];
  //sprintf (file, "%s%s%d%s", dmrginp.save_prefix().c_str(),"/testspatial_threepdm.",mpigetrank(),".bin");
  //std::ofstream spatpdm_disk(file, std::ios::ate| std::ios::binary);
  //spatpdm_disk.open(file, std::ios::binary);
  //boost::archive::binary_oarchive save(spatpdm_disk);
  //for(auto it= index_elements.begin(); it!=index_elements.end(); it++)
  //{
  //  //pout <<"element:  "<<it-> index<< "\t\t"<<it->element<<endl;
  //  save << *it;
  //  //spatpdm_disk << it->index;
  //  //spatpdm_disk << it->element;

  //}
  //spatpdm_disk.close();
  
#if 0
#ifndef SERIAL
  //mpi 
  boost::mpi::communicator world;
  std::vector<std::map < long, double>> indexelement_array;

  // dump std::map of index and elements into disk;
  if(mpigetrank()==0){
    boost::mpi::gather(world, index_and_elements, indexelement_array, 0);
    std::map < long, double> elements;
    for(auto it = indexelement_array.begin(); it!=indexelement_array.end();it++)
      elements.insert(it->begin(),it->end());

    for(auto it = elements.begin(); it!=elements.end();it++)
    {
      spatialpdm_disk.seekp(it->first*sizeof(double),ios::beg);
        spatialpdm_disk << it->second;
    }
  }
  else
    boost::mpi::gather(world, index_and_elements, 0);

#else
    for(auto it = index_and_elements.begin(); it!=index_and_elements.end();it++)
    {
     // spatialpdm_disk.seekp(it->first*sizeof(double),ios_base::beg);
      spatialpdm_disk.seekp((it->first)*sizeof(double));
      spatialpdm_disk << it->second;
    }
#endif
#endif
}
  */


void Threepdm_container::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
  Threepdm_permutations perm;
  if(dmrginp.spinAdapted())
  {

    std::vector< std::pair< std::vector<int>, double > > spatial_batch;
    perm.get_spatial_batch(new_spin_orbital_elements,spatial_batch);
    if(dmrginp.spatpdm_disk_dump() ){
      dump_text_to_disk(spatial_batch);
    }
    else update_full_spatial_array(spatial_batch);
    if( dmrginp.store_spinpdm())
    {
      std::vector< std::pair< std::vector<int>, double > > spin_batch;
      perm.process_new_elements( new_spin_orbital_elements, nonredundant_elements, spin_batch );
      update_full_spin_array( spin_batch );
    }
  }
  else{
    std::vector< std::pair< std::vector<int>, double > > spin_batch;
    perm.process_new_elements( new_spin_orbital_elements, nonredundant_elements, spin_batch );
    update_full_spin_array( spin_batch );
  }
}

//===========================================================================================================================================================

}
}


