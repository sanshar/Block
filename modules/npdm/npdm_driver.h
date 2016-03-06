#ifndef NPDM_DRIVER_HEADER_H
#define NPDM_DRIVER_HEADER_H

#include "npdm_spin_adaptation.h"
#include "npdm_expectations.h"
#include "onepdm_container.h"
#include "twopdm_container.h"
#include "threepdm_container.h"
#include "fourpdm_container.h"
#include "pairpdm_container.h"
#include "npdm.h"
#include "wavefunction.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

class Npdm_driver_base {
  public:
    Npdm_driver_base() {}
    virtual ~Npdm_driver_base() {}
    virtual void clear() = 0;
    virtual void save_data( const int i, const int j ) = 0;
    virtual void compute_npdm_elements( std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos ) = 0;
};

//===========================================================================================================================================================

class Npdm_driver {

  public:
		double diskread_time =0;
    double write_intermediate_time =0;
    explicit Npdm_driver(NpdmOrder order, Npdm_container& container) : npdm_order_(order), container_(container) {}
    ~Npdm_driver() {}
   void clear() { container_.clear();}
   void save_data( const int i, const int j ) { container_.save_npdms(i,j); }
   void compute_npdm_elements( std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos );

  private:
    NpdmOrder npdm_order_;
    Npdm_container& container_;
    Npdm_spin_adaptation spin_adaptation_;
		std::vector<boost::shared_ptr<NpdmSpinOps> > inner_Operators;
		std::vector<boost::shared_ptr<std::map<std::vector<int>, Wavefunction> > > inner_intermediate;

    void loop_over_operator_patterns( Npdm_patterns& patterns, Npdm_expectations& expectations, const SpinBlock& big );
    void do_inner_loop( const char inner, Npdm_expectations & npdm_expectations, 
                        NpdmSpinOps_base & outerOps, NpdmSpinOps & innerOps, NpdmSpinOps & dotOps );
    void do_inner_loop( const char inner, Npdm_expectations & npdm_expectations, 
                        NpdmSpinOps_base & outerOps, NpdmSpinOps & dotOps, std::map<std::vector<int>, Wavefunction>& waves );
    void loop_over_block_operators( Npdm_expectations& npdm_expectations, NpdmSpinOps& outerOps, NpdmSpinOps& innerOps, NpdmSpinOps& dotOps );

    void par_loop_over_block_operators( const char inner, Npdm_expectations & npdm_expectations, 
                                        NpdmSpinOps & lhsOps, NpdmSpinOps & outerOps, NpdmSpinOps & dotOps, bool innerdot );

    void do_parallel_lhs_loop( const char inner, Npdm_expectations & npdm_expectations,
                               NpdmSpinOps & outerOps, NpdmSpinOps & innerOps, NpdmSpinOps & dotOps, bool skip );
    void do_parallel_intermediate_loop( const char inner, Npdm_expectations & npdm_expectations,
                               NpdmSpinOps & outerOps, NpdmSpinOps & innerOps, NpdmSpinOps & dotOps, bool skip );
    void get_inner_Operators( const char inner, Npdm_expectations & npdm_expectations, boost::shared_ptr<NpdmSpinOps> lhsOps, boost::shared_ptr<NpdmSpinOps> dotOps, boost::shared_ptr<NpdmSpinOps> rhsOps) ;
    int get_mpi_max_size( int my_size );
    bool broadcast_lhs( int lhs_size, int rhs_size );
    bool skip_this_mpi_rank( NpdmSpinOps & outerOps, NpdmSpinOps & innerOps );
    bool skip_parallel( NpdmSpinOps & outerOps, NpdmSpinOps & innerOps, bool lhsrhsdot );
    bool screen(const std::vector<CD> &lhs_cd_type,const std::vector<CD> &dot_cd_type);
    void clear_npdm_intermediate(Npdm::Npdm_expectations& expectations);

};
  
//===========================================================================================================================================================
// Here we set up some convenience classes which group appropriate drivers and npdm containers

class Onepdm_driver : public Npdm_driver_base {
  public:
    explicit Onepdm_driver( int sites ) : container( Onepdm_container(sites) ), driver( Npdm_driver(NPDM_ONEPDM, container) ) {}
    void clear() { driver.clear(); }
    void save_data( const int i, const int j ) { driver.save_data(i,j); }
    void compute_npdm_elements( std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos ) 
      { driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos ); }
  private:
    Onepdm_container container;
    Npdm_driver driver;
};

//===========================================================================================================================================================

class Twopdm_driver : public Npdm_driver_base {
  public:
    explicit Twopdm_driver( int sites ) : container( Twopdm_container(sites) ), driver( Npdm_driver(NPDM_TWOPDM, container) ) {}
    void clear() { driver.clear(); }
    void save_data( const int i, const int j ) { driver.save_data(i,j); }
    void compute_npdm_elements( std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos ) 
      { driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos ); }
  private:
    Twopdm_container container;
    Npdm_driver driver;
};

//===========================================================================================================================================================

class Threepdm_driver : public Npdm_driver_base {
  public:
    explicit Threepdm_driver( int sites ) : container( Threepdm_container(sites) ), driver( Npdm_driver(NPDM_THREEPDM, container) ) {}
    void clear() { driver.clear(); container= Threepdm_container(dmrginp.last_site());/*Reinitialize container to reopen the temporary file.*/}
    void save_data( const int i, const int j ) { driver.save_data(i,j); }
    void compute_npdm_elements( std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos ) 
      { driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos ); }
  private:
    Threepdm_container container;
    Npdm_driver driver;
};

//===========================================================================================================================================================

class Fourpdm_driver : public Npdm_driver_base {
  public:
    explicit Fourpdm_driver( int sites ) : container( Fourpdm_container(sites) ), driver( Npdm_driver(NPDM_FOURPDM, container) ) {}
    void clear() { driver.clear(); container= Fourpdm_container(dmrginp.last_site()); /*Reinitialize container to reopen the temporary file.*/}
    void save_data( const int i, const int j ) { driver.save_data(i,j); }
    void compute_npdm_elements( std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos ) 
      { driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos ); }
  private:
    Fourpdm_container container;
    Npdm_driver driver;
};

//===========================================================================================================================================================

class Pairpdm_driver : public Npdm_driver_base {
  public:
    explicit Pairpdm_driver( int sites ) : container( Pairpdm_container(sites) ), driver( Npdm_driver(NPDM_PAIRMATRIX, container) ) {}
    void clear() { driver.clear(); }
    void save_data( const int i, const int j ) { driver.save_data(i,j); }
    void compute_npdm_elements( std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos ) 
      { driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos ); }
  private:
    Pairpdm_container container;
    Npdm_driver driver;
};

}
}

#endif

