/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_PERMUTATIONS_H
#define NPDM_PERMUTATIONS_H

#include <vector>
#include <utility>

namespace SpinAdapted{

//===========================================================================================================================================================

class Npdm_permutations {
  public:
    std::vector<std::pair<std::vector<int>,int> > reorders;
    void get_permute(const std::pair<std::vector<int>,int>& origin, int start, int n, std::vector<std::pair<std::vector<int>,int> >& reorders);
    Npdm_permutations() {}
    virtual ~Npdm_permutations() {}

    void process_new_elements( const std::vector<std::pair<std::vector<int>, double > >& in, 
                               std::vector< std::pair< std::vector<int>, double > >& nonredundant_elements,
                               std::vector< std::pair< std::vector<int>, double > >& spin_batch );

//    virtual void get_spatial_batch( const std::vector< std::pair< std::vector<int>, double > >& in, 
//                                              std::vector< std::pair< std::vector<int>, double > >& spatial_perms ) =0;
  private:
    virtual void get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
                                        const std::vector<int>& indices, const double& val ) = 0;
};

//===========================================================================================================================================================

class Onepdm_permutations : public Npdm_permutations {
  public:
    Onepdm_permutations(){
      if(reorders.size()>0) return;
      std::vector<int> origin = {0};
      reorders.push_back(std::make_pair(origin,0));
    }
    void get_spatial_batch( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& spatial_perms );
  private:
    void get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, const std::vector<int>& indices, const double& val );
};

//===========================================================================================================================================================

class Twopdm_permutations : public Npdm_permutations {
  public:
    Twopdm_permutations(){
      if(reorders.size()>0) return;
      std::vector<int> origin_order = {0,1};
      std::pair<std::vector<int>,int> origin = std::make_pair(origin_order,1);
      get_permute(origin,0,2,reorders);
    }
    void get_spatial_batch( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& spatial_perms );
  private:
    void get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, const std::vector<int>& indices, const double& val );
};

//===========================================================================================================================================================

class Threepdm_permutations : public Npdm_permutations {
  public:
    Threepdm_permutations(){
      if(reorders.size()>0) return;
      std::vector<int> origin_order = {0,1,2};
      std::pair<std::vector<int>,int> origin = std::make_pair(origin_order,1);
      get_permute(origin,0,3,reorders);
    }
    void get_spatial_batch( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& spatial_perms );
  private:
    void get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, const std::vector<int>& indices, const double& val );
};

//===========================================================================================================================================================

class Fourpdm_permutations : public Npdm_permutations {
  public:
    Fourpdm_permutations(){
      if(reorders.size()>0) return;
      std::vector<int> origin_order = {0,1,2,3};
      std::pair<std::vector<int>,int> origin = std::make_pair(origin_order,1);
      get_permute(origin,0,4,reorders);
    }
    void get_spatial_batch( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& spatial_perms );
  private:
    void get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, const std::vector<int>& indices, const double& val );
//    void get_even_and_odd_perms( const std::vector<int> mnpq, std::vector<std::vector<int> > & even_perms, std::vector<std::vector<int> > & odd_perms );
//    void get_spin_permutations_general(std::vector<std::pair<std::vector<int>,double> >& spin_batch, const std::vector<int>& indices, const double& val);
};

//===========================================================================================================================================================

class Pairpdm_permutations : public Npdm_permutations {
  private:
    void get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, const std::vector<int>& indices, const double& val );
};
}

#endif
