/*                                                                           Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_SWEEP_PARAMS_HEADER
#define SPIN_SWEEP_PARAMS_HEADER
#include <boost/serialization/serialization.hpp>

namespace SpinAdapted{
enum guessWaveTypes {BASIC, TRANSFORM, TRANSPOSE};

class SweepParams
{
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & block_iter & sweep_iter & n_iters & forward_starting_size & backward_starting_size & keep_states & keep_states_ls;
    ar & keep_qstates & sys_add & env_add & noise & additional_noise & davidson_tol & lowest_energy & lowest_energy_spins & guesstype;
    ar & error & largest_dw & onedot & currentRoot;
  }

  int restart_iter;
  bool onedot;
  int block_iter;
  int sweep_iter;
  int n_iters;
  int forward_starting_size;
  int backward_starting_size;
  int keep_states;
  int keep_states_ls;
  int keep_qstates;
  int sys_add;
  int env_add;
  int currentRoot;
  double noise;
  double additional_noise;
  double davidson_tol;
  double error;
  double largest_dw;
  vector<double> lowest_energy;
  vector<double> lowest_energy_spins;
  guessWaveTypes guesstype;

public:
  SweepParams();
  vector<double> ls_dw;
  vector<double> ls_energy;
  void set_sweep_parameters();
  void savestate(const bool &forward, const int &size);
  void restorestate(bool &forward, int &size);
  void calc_niter();

  const int &current_root() const {return currentRoot;}
  const bool &get_onedot() const { return onedot; }
  const int &get_block_iter() const { return block_iter; }
  const int &get_sweep_iter() const { return sweep_iter; }
  const int &get_n_iters() const { return n_iters; }
  const int &get_forward_starting_size() const { return forward_starting_size; }
  const int &get_backward_starting_size() const { return backward_starting_size; }
  const int &get_keep_states() const { return keep_states; }
  const int &get_keep_states_ls() const { return keep_states_ls; }
  const int &get_keep_qstates() const { return keep_qstates; }
  const int &get_sys_add() const { return sys_add; }
  const int &get_env_add() const { return env_add; }
  const double &get_noise() const { return noise; }
  const double &get_additional_noise() const {return additional_noise;}
  const double &get_davidson_tol() const { return davidson_tol; }
  const double &get_lowest_error() const { return error; }
  const double &get_largest_dw() const { return largest_dw; }
  const vector<double> &get_lowest_energy() const { return lowest_energy; }
  const vector<double> &get_lowest_energy_spins() const { return lowest_energy_spins; }
  const guessWaveTypes &get_guesstype() const { return guesstype; }
  const int &get_restart_iter() const {return restart_iter;}

  int &current_root() {return currentRoot;}
  int &set_restart_iter() {return restart_iter;}
  bool &set_onedot() { return onedot; }
  int &set_block_iter() { return block_iter; }
  int &set_sweep_iter() { return sweep_iter; }
  int &set_n_iters() { return n_iters; }
  int &set_forward_starting_size() { return forward_starting_size; }
  int &set_backward_starting_size() { return backward_starting_size; }
  int &set_keep_states() { return keep_states; }
  int &set_keep_states_ls() { return keep_states_ls; }
  int &set_keep_qstates() { return keep_qstates; }
  int &set_sys_add() { return sys_add; }
  int &set_env_add() { return env_add; }
  double &set_noise() { return noise; }
  double &set_additional_noise() {return additional_noise;}
  double &set_davidson_tol() { return davidson_tol; }
  double &set_lowest_error() { return error; }
  double &set_largest_dw() { return largest_dw; }
  vector<double> &set_lowest_energy() { return lowest_energy; }
  vector<double> &set_lowest_energy_spins() { return lowest_energy_spins; }
  guessWaveTypes &set_guesstype() { return guesstype; }
};
}
#endif
