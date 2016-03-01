/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include <boost/format.hpp>
#include "fourpdm_container.h"
#include "npdm_permutations.h"
#include <boost/range/algorithm.hpp>
#include <boost/filesystem.hpp>

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Fourpdm_container::Fourpdm_container( int sites )
{
  elements_stride_.resize(8);
  for(int i=0; i<8; i++ )
    elements_stride_[i]=pow(sites,7-i);

  if ( dmrginp.store_spinpdm() ) {
    if(dmrginp.spinAdapted())
      fourpdm.resize(2*sites,2*sites,2*sites,2*sites,2*sites,2*sites,2*sites,2*sites);
    else
      fourpdm.resize(sites,sites,sites,sites,sites,sites,sites,sites);
    fourpdm.Clear();
  } 
  if ( !dmrginp.spatpdm_disk_dump() ){
    if(dmrginp.spinAdapted())
      spatial_fourpdm.resize(sites,sites,sites,sites,sites,sites,sites,sites);
    else
      spatial_fourpdm.resize(sites/2,sites/2,sites/2,sites/2,sites/2,sites/2,sites/2,sites/2);
    spatial_fourpdm.Clear();
  } 

  else{
    char file[5000];
    sprintf (file, "%s%s%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.",mpigetrank(),".tmp");
    //std::ofstream spatpdm_disk(file, std::ios::binary);
    spatpdm_disk=fopen(file,"w");
    //32M buffer;
    setvbuf(spatpdm_disk,NULL,_IOFBF,1024*1024*32);
    if(mpigetrank()==0)
      fprintf(spatpdm_disk,"%d\n",sites);
    //spatpdm_disk.open(file, std::ios::binary);
    //spatpdm_disk.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_npdms(const int& i, const int& j)
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
  p3out << "4PDM save full array time " << ewall << " " << ecpu << endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/fourpdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << fourpdm.dim1() << endl;

    double trace = 0.0;
    for(int i=0; i<fourpdm.dim1(); ++i)
      for(int j=0; j<fourpdm.dim2(); ++j)
        for(int k=0; k<fourpdm.dim3(); ++k)
          for(int l=0; l<fourpdm.dim4(); ++l)
            for(int m=0; m<fourpdm.dim5(); ++m)
              for(int n=0; n<fourpdm.dim6(); ++n)
                for(int p=0; p<fourpdm.dim7(); ++p)
                  for(int q=0; q<fourpdm.dim8(); ++q) {
                    if ( abs(fourpdm(i,j,k,l,m,n,p,q)) > NUMERICAL_ZERO ) {
                      ofs << boost::format("%d %d %d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % p % q % fourpdm(i,j,k,l,m,n,p,q);
                      if ( (i==q) && (j==p) && (k==n) && (l==m) ) trace += fourpdm(i,j,k,l,m,n,p,q);
                    }
                  }
    ofs.close();
    pout << "Spin-orbital 4PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_spatial_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << spatial_fourpdm.dim1() << endl;

    double trace = 0.0;
    for(int i=0; i<spatial_fourpdm.dim1(); ++i)
      for(int j=0; j<spatial_fourpdm.dim2(); ++j)
        for(int k=0; k<spatial_fourpdm.dim3(); ++k)
          for(int l=0; l<spatial_fourpdm.dim4(); ++l)
            for(int m=0; m<spatial_fourpdm.dim5(); ++m)
              for(int n=0; n<spatial_fourpdm.dim6(); ++n)
                for(int p=0; p<spatial_fourpdm.dim7(); ++p)
                  for(int q=0; q<spatial_fourpdm.dim8(); ++q) {
                    if ( abs(spatial_fourpdm(i,j,k,l,m,n,p,q)) > NUMERICAL_ZERO ) {
                      ofs << boost::format("%d %d %d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % p % q % spatial_fourpdm(i,j,k,l,m,n,p,q);
                      if ( (i==q) && (j==p) && (k==n) && (l==m) ) trace += spatial_fourpdm(i,j,k,l,m,n,p,q);
                    }
                  }
    ofs.close();
    pout << "Spatial      4PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_npdm_binary(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/fourpdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << fourpdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::external_sort_index(const int &i, const int &j)
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
      sprintf (file, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm_index.", i, j,p,".bin");
      FILE* inputfile=fopen(file,"wb");
      fwrite(&nonspin_batch[0],sizeof(Sortpdm::batch_index),nonspin_batch.size(),inputfile);
      fclose(inputfile);
      nonspin_batch.clear();
    }
    //external sort nonspin_batch
    //TODO
    //It is not parallel.
    char outfilename[5000];
    sprintf (outfilename, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm_index.",i,j,".bin");
    FILE* outputfile = fopen(outfilename,"wb");
    long sorting_buff= 1024*1024*(32/world.size());
    //For batch_index, the sorting buff is about 96M/world.size();
    std::vector<Sortpdm::cache<Sortpdm::batch_index>> filecache;
    for(int p=0; p< world.size();p++){
      char file[5000];
      sprintf (file, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm_index.", i, j,p,".bin");
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
        sprintf (file, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm_index.", i, j,p,".bin");
        boost::filesystem::remove(file);
      }
    }
#else
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/fourpdm_index.", i, j,".bin");
    FILE* outfile=fopen(file,"wb");
    fwrite(&nonspin_batch[0],sizeof(Sortpdm::batch_index),nonspin_batch.size(),outfile);
    nonspin_batch.clear();
#endif
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_spatial_npdm_binary(const int &i, const int &j)
{
  if(!dmrginp.spatpdm_disk_dump())
  {
    if( mpigetrank() == 0)
    {
      char file[5000];
      sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.", i, j,".bin");
      std::ofstream ofs(file, std::ios::binary);
      boost::archive::binary_oarchive save(ofs);
      save << spatial_fourpdm;
      ofs.close();
    }
  }
  else{
    fclose(spatpdm_disk);
    //When spatpdm_disk is opened, the state numbers, i and j, are not known. 
    char oldfile[5000];
    sprintf (oldfile, "%s%s%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.",mpigetrank(),".tmp");
    char newfile[5000];
    sprintf (newfile, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.",i,j,mpigetrank(),".bin");
    boost::filesystem::rename(oldfile,newfile);

    if(dmrginp.pdm_unsorted())
      external_sort_index(i,j);
    else{
      char file[5000];
      sprintf (file, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.",i,j,mpigetrank(),".bin");
      char finalfile[5000];
      sprintf (finalfile, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.",i,j,".bin");
#ifndef SERIAL
      boost::mpi::communicator world;
      world.barrier();
      Timer timer1;
      char tmpfile[5000];
      sprintf (tmpfile, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.",i,j,mpigetrank(),".tmp");
      char sortedfile[5000];
      sprintf (sortedfile, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.",i,j,mpigetrank(),".bin");

      Sortpdm::partition_data<Sortpdm::index_element>((long)pow(dmrginp.last_site(),8),file,tmpfile);
      Sortpdm::externalsort<Sortpdm::index_element>(tmpfile,sortedfile,(long)pow(dmrginp.last_site(),8));
      world.barrier();
      ecpu = timer1.elapsedcputime(); ewall = timer1.elapsedwalltime();
      p3out << "4PDM parallel external sort time " << ewall << " " << ecpu << endl;
      Timer timer;
      Sortpdm::mergefile(sortedfile);
      world.barrier();
      if(mpigetrank()==0) boost::filesystem::rename(sortedfile,finalfile);
      boost::filesystem::remove(tmpfile);
      ecpu = timer.elapsedcputime();ewall=timer.elapsedwalltime();
      p3out << "4PDM merge sorted file time " << ewall << " " << ecpu << endl;
#else
      Timer timer2;
      Sortpdm::externalsort<Sortpdm::index_element> (file,finalfile,(long)pow(dmrginp.last_site(),8));
      boost::filesystem::remove(file);
      ecpu = timer2.elapsedcputime();ewall=timer2.elapsedwalltime();
      p3out << "4PDM external sort time " << ewall << " " << ecpu << endl;
#endif
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::load_npdm_binary(const int &i, const int &j) { abort(); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::accumulate_npdm()
{
#ifndef SERIAL
  mpi::communicator world;
#ifdef NDEBUG
  if( mpigetrank() == 0)
    MPI_Reduce(MPI_IN_PLACE, fourpdm.data(),  fourpdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
  else
    MPI_Reduce(fourpdm.data(), fourpdm.data(),  fourpdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
#else
  array_8d<double> tmp_recv;
  if( mpigetrank() == 0)
  {
    for(int u=1; u<world.size(); ++u) {
      world.recv(u, u, tmp_recv);
      for(int i=0; i<fourpdm.dim1(); ++i)
        for(int j=0; j<fourpdm.dim2(); ++j)
          for(int k=0; k<fourpdm.dim3(); ++k)
            for(int l=0; l<fourpdm.dim4(); ++l)
              for(int m=0; m<fourpdm.dim5(); ++m)
                for(int n=0; n<fourpdm.dim6(); ++n) 
                  for(int p=0; p<fourpdm.dim7(); ++p) 
                    for(int q=0; q<fourpdm.dim8(); ++q) {
                      if ( abs(tmp_recv(i,j,k,l,m,n,p,q)) > NUMERICAL_ZERO ) {
                        // Test for duplicates
                        assert(abs(fourpdm(i,j,k,l,m,n,p,q)) < NUMERICAL_ZERO );
                        fourpdm(i,j,k,l,m,n,p,q) = tmp_recv(i,j,k,l,m,n,p,q);
                      }
                    }
    }
  }
  else 
  {
    world.send(0, mpigetrank(), fourpdm);
  }
#endif
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::accumulate_spatial_npdm()
{
#ifndef SERIAL
  mpi::communicator world;
#ifdef NDEBUG
  if( mpigetrank() == 0)
    MPI_Reduce(MPI_IN_PLACE, spatial_fourpdm.data(),  spatial_fourpdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
  else
    MPI_Reduce(spatial_fourpdm.data(), spatial_fourpdm.data(),  spatial_fourpdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
#else
  array_8d<double> tmp_recv;
  if( mpigetrank() == 0)
  {
    for(int u=1; u<world.size(); ++u) {
      world.recv(u, u, tmp_recv);
      for(int i=0; i<spatial_fourpdm.dim1(); ++i)
        for(int j=0; j<spatial_fourpdm.dim2(); ++j)
          for(int k=0; k<spatial_fourpdm.dim3(); ++k)
            for(int l=0; l<spatial_fourpdm.dim4(); ++l)
              for(int m=0; m<spatial_fourpdm.dim5(); ++m)
                for(int n=0; n<spatial_fourpdm.dim6(); ++n) 
                  for(int p=0; p<spatial_fourpdm.dim7(); ++p) 
                    for(int q=0; q<spatial_fourpdm.dim8(); ++q) {
                      if ( abs(tmp_recv(i,j,k,l,m,n,p,q)) > NUMERICAL_ZERO ) {
                        // Test for duplicates
                        assert(abs(spatial_fourpdm(i,j,k,l,m,n,p,q)) < NUMERICAL_ZERO );
                        spatial_fourpdm(i,j,k,l,m,n,p,q) = tmp_recv(i,j,k,l,m,n,p,q);
                      }
                    }
    }
  }
  else 
  {
    world.send(0, mpigetrank(), spatial_fourpdm);
  }
#endif
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::update_full_spin_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    double val = it->second;
    if ( abs(val) < NUMERICAL_ZERO ) continue;

    assert( (it->first).size() == 8 );
    // Spin indices
    int i0 = (it->first)[0];
    int j0 = (it->first)[1];
    int k0 = (it->first)[2];
    int l0 = (it->first)[3];
    int m0 = (it->first)[4];
    int n0 = (it->first)[5];
    int p0 = (it->first)[6];
    int q0 = (it->first)[7];
    int i = ro.at(i0/2)*2 + i0%2;
    int j = ro.at(j0/2)*2 + j0%2;
    int k = ro.at(k0/2)*2 + k0%2;
    int l = ro.at(l0/2)*2 + l0%2;
    int m = ro.at(m0/2)*2 + m0%2;
    int n = ro.at(n0/2)*2 + n0%2;
    int p = ro.at(p0/2)*2 + p0%2;
    int q = ro.at(q0/2)*2 + q0%2;

    //if ( abs(val) > 1e-8 ) {
    //  pout << "so-fourpdm val: i,j,k,l,m,n,p,q = " 
    //       << i << "," << j << "," << k << "," << l << "," << m << "," << n << "," << p << "," << q
    //      << "\t\t" << val << endl;
    //}

    fourpdm(i,j,k,l,m,n,p,q) = val;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// This routine assumes that no spin-orbital indices are generated more than once

void Fourpdm_container::update_full_spatial_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 8 );

    // Store significant elements only
    if ( abs(it->second) > NUMERICAL_ZERO ) {
      // Spin indices
      int i = (it->first)[0];
      int j = (it->first)[1];
      int k = (it->first)[2];
      int l = (it->first)[3];
      int m = (it->first)[4];
      int n = (it->first)[5];
      int p = (it->first)[6];
      int q = (it->first)[7];
  
//      if ( i%2 != q%2 ) continue;
//      if ( j%2 != p%2 ) continue;
//      if ( k%2 != n%2 ) continue;
//      if ( l%2 != m%2 ) continue;
//
//      spatial_fourpdm( ro.at(i/2), ro.at(j/2), ro.at(k/2), ro.at(l/2), ro.at(m/2), ro.at(n/2), ro.at(p/2), ro.at(q/2) ) += it->second;
      spatial_fourpdm( ro.at(i), ro.at(j), ro.at(k), ro.at(l), ro.at(m), ro.at(n), ro.at(p), ro.at(q) ) = it->second;
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void Fourpdm_container::build_spatial_elements( std::map< std::vector<int>, double >& spin_batch, 
//                                                std::map< std::vector<int>, double >& spatial_batch )
//{
//  double factor = 1.0;
//
//  for (auto it = spatial_batch.begin(); it != spatial_batch.end(); ++it) {
//    assert( (it->first).size() == 8 );
//    int i = (it->first)[0];
//    int j = (it->first)[1];
//    int k = (it->first)[2];
//    int l = (it->first)[3];
//    int m = (it->first)[4];
//    int n = (it->first)[5];
//    int p = (it->first)[6];
//    int q = (it->first)[7];
//    // Sum over spin indices
//    double val = 0.0;
//    for (int s=0; s<2; s++) {
//      for (int t=0; t<2; t++) {
//        for (int u=0; u<2; u++) {
//          for (int v=0; v<2; v++) {
//            std::vector<int> idx = { 2*i+s, 2*j+t, 2*k+u, 2*l+v, 2*m+v, 2*n+u, 2*p+t, 2*q+s };
//            val += spin_batch[ idx ];
//          }
//        }
//      }
//    }
//    // Store significant elements only
//    if ( abs(val) > 1e-14 ) {
//      if ( store_sparse_spatial_array_ ) sparse_spatial_pdm[ it->first ] = factor * val;
//      if ( store_full_spatial_array_ ) {
//        if ( abs( spatial_fourpdm(i,j,k,l,m,n,p,q) ) > 1e-14 ) {
//          pout << "repeated spatial indices!\n";
//          abort();
//        }
//        spatial_fourpdm(i,j,k,l,m,n,p,q) = factor * val;
//      }
//    }
//  }
//
//}
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void Fourpdm_container::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
//{
//  assert( new_spin_orbital_elements.size() == 70 );
//  // Temporary batches of npdm elements
//  std::map< std::vector<int>, double > spin_batch;
//  std::map< std::vector<int>, double > spatial_batch;
//
//  for (int idx=0; idx < new_spin_orbital_elements.size(); ++idx) {
//    // Get all spin-index permutations
//    Fourpdm_permutations p;
//    std::map< std::vector<int>, int > spin_indices = p.get_spin_permutations( new_spin_orbital_elements[idx].first );
//    double val = new_spin_orbital_elements[idx].second;
//    for (auto it = spin_indices.begin(); it != spin_indices.end(); ++it) {
//      // Initialize spatial indices
//      std::vector<int> vec;
//      for (int i=0; i < (it->first).size(); ++i)
//        vec.push_back( (it->first)[i]/2 );
//      spatial_batch[ vec ] = 0.0;
//      // Assign temporary batch of spin-orbital elements
//      spin_batch[ it->first ] = it->second * val;
//      if ( store_sparse_spin_array_ && (abs(val) > 1e-14) ) sparse_spin_pdm[ it->first ] = it->second * val;
//    }
//  }
//
//  // Build and store new spatial elements
//  build_spatial_elements( spin_batch, spatial_batch );
//  if ( store_full_spin_array_ ) update_full_spin_array( spin_batch );
//}
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

long Fourpdm_container::oneindex_spin(const std::vector<int> & orbital_element_index)
{
  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  assert( orbital_element_index.size() == 8);
  long linearindex=0;
  for(int i=0; i< 8; i++)
    linearindex+=(ro.at(orbital_element_index[i]/2))*elements_stride_[i];
  return linearindex;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::dump_binary_to_disk(std::vector< std::pair< std::vector<int>, double > > & spin_batch)
{

  const std::vector<int>& ro = dmrginp.reorder_vector();
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 8 );

    // Store significant elements only
    if ( abs(it->second) > NUMERICAL_ZERO ) {
      std::vector<int> sites;
      sites.resize(8);
      for (int i=0;i<8;i++)
        sites[i] = ro.at((it->first)[i]);

      fwrite(sites.data(),sizeof(int),8,spatpdm_disk);
      fwrite(&(it->second),sizeof(double),1,spatpdm_disk);
    }
  }
}

void Fourpdm_container::dump_text_to_disk(std::vector< std::pair< std::vector<int>, double > > & spin_batch)
{

  const std::vector<int>& ro = dmrginp.reorder_vector();
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 8 );

    // Store significant elements only
    if ( abs(it->second) > NUMERICAL_ZERO ) {

      int i = (it->first)[0];
      int j = (it->first)[1];
      int k = (it->first)[2];
      int l = (it->first)[3];
      int m = (it->first)[4];
      int n = (it->first)[5];
      int p = (it->first)[6];
      int q = (it->first)[7];

      fprintf(spatpdm_disk,"%d %d %d %d %d %d %d %d %20.14e\n", ro.at(i), ro.at(j), ro.at(k), ro.at(l), ro.at(m), ro.at(n), ro.at(p), ro.at(q) ,it->second);
    }
  }
}

void Fourpdm_container::merge_diskfile(const int &i, const int &j)
{
  fclose(spatpdm_disk);
  if(mpigetrank()==0)
  {
    char oldfile[5000];
    sprintf (oldfile, "%s%s%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.",mpigetrank(),".tmp");
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.", i, j,".txt");
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
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.", i, j,".txt");
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
    sprintf (file, "%s%s%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.",mpigetrank(),".tmp");
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
/*
void Fourpdm_container::dump_to_disk(std::vector< std::pair< std::vector<int>, double > > & spin_batch)
{
  long spatpdm_disk_position= ftell(spatpdm_disk);
  std::map < long, double>  index_and_elements;
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 8 );

    // Store significant elements only
    if ( abs(it->second) > NUMERICAL_ZERO ) {
      // Spin indices
      if ( it->first[0]%2 != it->first[7]%2 ) continue;
      if ( it->first[1]%2 != it->first[6]%2 ) continue;
      if ( it->first[2]%2 != it->first[5]%2 ) continue;
      if ( it->first[3]%2 != it->first[4]%2 ) continue;
      long linearindex = oneindex_spin(it->first);
      std::map < long, double>::iterator findit= index_and_elements.find(linearindex);
      if(findit == index_and_elements.end()){
        index_and_elements.insert(std::pair<long,double>(linearindex,it->second));
      }
      else
        findit->second += it->second;
    }
  }
  if(index_and_elements.size()==0) return;
  Sortpdm::index_element index_elements[1152];
  // number of permutation is (4*3*2)*(4*3*2) , and then transpose(x2) .
  assert(1152>= index_and_elements.size());
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
}
*/


//===========================================================================================================================================================
void Fourpdm_container::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
  assert( new_spin_orbital_elements.size() == 70 );
  Fourpdm_permutations perm;
  if(dmrginp.spinAdapted())
  {
    //std::vector< std::pair< std::vector<int>, double > > spin_batch;
    std::vector< std::pair< std::vector<int>, double > > spatial_batch;
    // Work with the non-redundant elements only, and get all unique spin-permutations as a by-product
    //perm.process_new_elements( new_spin_orbital_elements, nonredundant_elements, spin_batch );
    perm.get_spatial_batch(new_spin_orbital_elements,spatial_batch);

    //if ( dmrginp.store_spinpdm() ) update_full_spin_array( spin_batch );
//    if ( !dmrginp.spatpdm_disk_dump() ) update_full_spatial_array( spin_batch );
//    else{
//      if(dmrginp.store_nonredundant_pdm()) 
//        dump_to_disk(nonredundant_elements);
//      else
    if(dmrginp.spatpdm_disk_dump() ){
        dump_text_to_disk(spatial_batch);
    }
    else update_full_spatial_array(spatial_batch);
    //if ( ! dmrginp.store_nonredundant_pdm() || dmrginp.spatpdm_disk_dump() ) nonredundant_elements.clear();
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
