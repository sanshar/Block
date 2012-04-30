#include <iostream>
#include <fstream>
#include <communicate.h>
#include "Symmetry.h"
#include "global.h"
#include "MatrixBLAS.h"
#include "spinblock.h"
#include "couplingCoeffs.h"
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
using namespace std;

namespace SpinAdapted {
string sym;
}

void SpinAdapted::Input::initialize_defaults()
{
#ifndef SERIAL
  mpi::communicator world;
  m_thrds_per_node = vector<int>(world.size(),1);
#else
  m_thrds_per_node = vector<int>(1,1);
#endif
  m_ham_type = QUANTUM_CHEMISTRY;
  m_algorithm_type = TWODOT;
  m_noise_type = RANDOM;
  m_calc_type = DMRG;
  
  m_outputlevel = 0;
  m_nquanta = 2;
  m_sys_add = 1;
  m_env_add = 1;
  m_total_symmetry_number = IrrepSpace( 0 );
  m_total_spin = 0;
  m_guess_permutations = 25;

  m_direct = true;
  m_nroots = 1;
  m_weights.resize(m_nroots);
  m_weights[0] = 1.;

  m_deflation_min_size = 2;
  m_deflation_max_size = 20;

  m_add_noninteracting_orbs = true;
  m_no_transform = false;
  m_do_fci = false;
  m_do_cd = false;
  m_maxiter = 1000;
  m_screen_tol = 1.00.e-20;

  m_load_prefix = ".";
  m_save_prefix = ".";

  m_maxj = 15;
  m_ninej.init(m_maxj);
  m_set_Sz = false;

  m_sweep_tol = 1e-8;
  m_restart = false;
  m_fullrestart = false;
  m_restart_warm = false;
  m_reset_iterations = false;

  m_do_diis = false;
  m_diis_error = 1e-2;
  m_start_diis_iter = 8;
  m_diis_keep_states = 6;
  m_diis_error_tol = 1e-8;
  m_num_spatial_orbs = 0;
}

//SpinAdapted::Input::Input(const string& config_name) : m_ninej(ninejCoeffs::getinstance())
SpinAdapted::Input::Input(const string& config_name)
{
  std::string symmetry;
  if(mpigetrank() == 0)
  {
    initialize_defaults();

    ifstream conf_file(config_name.c_str());
    conf_file >> m_norbs >> m_alpha >> m_beta;
  
    m_num_spatial_orbs = m_norbs/2;
    m_spatial_to_spin.resize(m_norbs/2+1);
    m_spin_to_spatial.resize(m_norbs);
    for (int i=0; i<m_norbs; i++) {
      m_spatial_to_spin[i/2] = i - i%2;
      m_spin_to_spatial[i] = i/2;
    }
    m_spatial_to_spin[m_norbs/2] = m_norbs;
    m_molecule_quantum = SpinQuantum(m_alpha + m_beta, abs(m_alpha - m_beta), m_total_symmetry_number);
    m_total_spin = abs(m_alpha-m_beta);

    string order_type;
    conf_file >> order_type;
    m_spin_vector.resize(m_norbs);
    for (int i=0; i<m_norbs; i++)
      m_spin_vector[i] = (i&1) ? -1 : 1;

    /* for the moment disabled
    if (order_type == "picklowest")
    {
      //make the hf determinant
      std::vector<double> alpha_one_particle_energies(m_norbs/2);
      std::vector<double> beta_one_particle_energies(m_norbs/2); 
      for (int i = 0; i < m_norbs; i+=2)                         
	alpha_one_particle_energies[i/2] = v_1(i, i);                                                                                                     
      for (int i = 1; i < m_norbs; i+=2)                                                                                                                 
	beta_one_particle_energies[i/2] = v_1(i, i);                                                       
      std::vector<int> alpha_sorted_indices, beta_sorted_indices;  
      get_sorted_indices(alpha_one_particle_energies, alpha_sorted_indices);                                                                                  
      get_sorted_indices(beta_one_particle_energies, beta_sorted_indices);        
      m_hf_occupancy.resize(m_norbs);                                       
      for (int i = 0; i < m_alpha; ++i)                                     
	m_hf_occupancy[2*alpha_sorted_indices[i]] = 1;                      
      for (int i = 0; i < m_beta; ++i)                                      
	m_hf_occupancy[2*beta_sorted_indices[i]+1] = 1;                           
    }
    */
    if (order_type == "hfocc")
    {
      m_hf_occupancy.resize(m_norbs);
      for (int i = 0; i < total_particle_number(); ++i)
	{
	  int occi;
	  conf_file >> occi;
	  m_hf_occupancy[occi] = 1;
	}
      int beta = 0, alpha = 0;
      for (int i = 0; i < m_norbs; i+=2)
	alpha += m_hf_occupancy[i];
      for (int i = 1; i < m_norbs; i+=2)
	beta += m_hf_occupancy[i];
      if(m_alpha != alpha || m_beta != beta)
	{
	  cout << "Error :: hfocc beta and alpha particles do not match the state\n";
	  cout << "Found " << alpha << " alpha particles from hfocc and expected " << m_alpha << endl;
	  cout << "Found " << beta << " beta particles from hfocc and expected " << m_beta << endl;
	  abort();
	}
    }
    else if (order_type == "defaulthf")
    {
      m_spin_vector.resize(m_norbs);
      for (int i = 0; i < m_norbs; ++i)
	m_spin_vector[i] = (i & 1) ? -1 : 1;
      m_hf_occupancy.resize(m_norbs);
      for (int i = 0; i < m_alpha; ++i)
	m_hf_occupancy[2*i] = 1;
      for (int i = 0; i < m_beta; ++i)
	m_hf_occupancy[2*i+1] = 1;
    }

    else
    {
      cout << "Error :: Did not find a guess (e.g. hfocc, or defaulthf) " << endl;
      abort();
    }

    string schedule;
    conf_file >> schedule;
    if (schedule != "schedule") 
    { 
      cout << "Error :: schedule must follow the guess (e.g. hfocc, picklowest, defaulthf)" << endl;
      cout << "Instead keyword found: "<<schedule<<endl;

      abort();
    }
    int nentry;
    conf_file >> nentry;


    m_sweep_iter_schedule.resize(nentry);
    m_sweep_state_schedule.resize(nentry);
    m_sweep_qstate_schedule.resize(nentry);
    m_sweep_tol_schedule.resize(nentry);
    m_sweep_noise_schedule.resize(nentry);
    m_sweep_additional_noise_schedule.resize(nentry);

    for (int i = 0; i < nentry; ++i)
    {
      int iter_val;
      int state_val;
      double tol_val;
      double noise_val;
      conf_file >> m_sweep_iter_schedule[i] >> m_sweep_state_schedule[i] >> m_sweep_qstate_schedule[i] >> m_sweep_tol_schedule[i] >> m_sweep_noise_schedule[i] >> m_sweep_additional_noise_schedule[i];
      if (i>0 && m_sweep_additional_noise_schedule[i] != 0.0 && m_sweep_additional_noise_schedule[i] == 0.0) {
	pout <<"Cannot have a non-zero additional noise after additional noise has been made zero in previous sweeps"<<endl; exit(0);}
      pout << m_sweep_iter_schedule[i] << " " << m_sweep_state_schedule[i] << " " << m_sweep_qstate_schedule[i] << " " << m_sweep_tol_schedule[i] << " " << m_sweep_noise_schedule[i] << " "<<m_sweep_additional_noise_schedule[i]<<endl;
    }

    m_quantaToKeep.resize(m_sweep_state_schedule[0]);

    string next_entry, discard;
    m_spin_orbs_symmetry.resize(m_norbs);
    symmetry = "c1";
    sym = "c1";
    while(conf_file >> next_entry)
    {
      size_t pos;
      pos = next_entry.find("!");
      if (int(pos) == 0) {
	conf_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	continue;
      }
      if (pos != string::npos) {
	next_entry = next_entry.substr(0,int(pos));
	conf_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
      if (next_entry == "sym")
      {
	m_spatial_to_spin.clear();
	m_spin_to_spatial.clear();

	m_num_spatial_orbs = 0;
        //string sym;
        conf_file >> sym;
        symmetry = sym;
        Symmetry::InitialiseTable(sym);
        pout << "Symmetry group " << sym << endl;
	pout << "Irreps of spatial orbitals "<<endl;
        for (int i = 0; i < m_norbs; ++i)
        {
          int ir;
          conf_file >> ir;
          m_spin_orbs_symmetry[i] = ir;
	  pout << ir<<" ";
	  if (sym == "dinfh") {
	    if (ir < -1 && i%2 == 0) {
	      m_num_spatial_orbs ++;
	      m_spatial_to_spin.push_back(i);
	    }
	    else if( (ir == 0 || ir ==1) && i%2 == 0) {
	      m_num_spatial_orbs++;
	      m_spatial_to_spin.push_back(i);
	    }
	    m_spin_to_spatial[i] = m_num_spatial_orbs-1;
	  }
	  else {
	    if (i%2 == 0) {
	      m_num_spatial_orbs++;
	      m_spatial_to_spin.push_back(i);
	    }
	    m_spin_to_spatial[i] = m_num_spatial_orbs-1;
	  }
        }
	m_spatial_to_spin.push_back(m_norbs);
	m_spin_to_spatial.push_back(m_norbs);
	cout << endl;
	cout <<"# spatial orbitals "<<m_num_spatial_orbs<<endl;
      }
      else if (next_entry == "Sz")
      {
        std::string wave_sym_str;
        conf_file >> wave_sym_str;
        m_Sz = atoi (wave_sym_str.c_str()) ;
	m_set_Sz = true;
      }
      else if (next_entry == "wave_sym")
      {
        std::string wave_sym_str;
        conf_file >> wave_sym_str;
        m_total_symmetry_number = IrrepSpace( atoi (wave_sym_str.c_str()) );
        m_molecule_quantum = SpinQuantum(m_alpha + m_beta, m_alpha - m_beta, m_total_symmetry_number);
        pout << "Solving for wave with symmetry :: " << m_molecule_quantum << endl;
      }
      else if (next_entry == "hubbard")
	m_ham_type = HUBBARD;
      else if (next_entry == "dmrg")
	m_calc_type = DMRG;
      else if (next_entry == "maxj")
	conf_file >> m_maxj;
      else if (next_entry == "fci")
	m_calc_type = FCI;
      else if (next_entry == "genblock")
	m_calc_type = GENBLOCK;
      else if (next_entry == "onepdm")
	m_calc_type = ONEPDM;
      else if (next_entry == "twopdm")
	m_calc_type = TWOPDM;
      else if (next_entry == "smalldavid")
      {
        m_solve_type = SMALL_DAVIDSON;
      }
      else if (next_entry == "bigdavid")
      {
        m_solve_type = BIG_DAVIDSON;
      }
      else if(next_entry == "prefix")
      {
	conf_file >> m_load_prefix ;
	m_save_prefix = m_load_prefix;
      }
      else if (next_entry == "noisetype")
      {
	//we add perturbative or random and perturbative or excited and perturbative
	std::string noise;
	conf_file >> noise;
	if (noise == "random")
	  m_noise_type = RANDOM;
	else if (noise == "excitedstate")
	  m_noise_type = EXCITEDSTATE;
	else {
	  pout << "noise type can be either random or excitedstate"<<endl;exit(0);
	}
      }
      else if (next_entry == "add_noninteracting_orbs")
	m_add_noninteracting_orbs = true;
      else if (next_entry == "diis_start_iter")
	conf_file >> m_start_diis_iter; 
      else if (next_entry == "diis_keep_states")
	conf_file >> m_diis_keep_states; 
      else if (next_entry == "diis_start_error")
	conf_file >> m_diis_error;
      else if (next_entry == "do_diis")
	m_do_diis = true;
      else if(next_entry == "nroots")
      {
        std::string nroots_str;
        conf_file >> nroots_str;
        m_nroots = atoi(nroots_str.c_str());
        if(m_deflation_min_size < m_nroots)
          m_deflation_min_size = m_nroots;
        //next entry must be nweights
        std::string weights_str;
        conf_file >> weights_str;
        if(weights_str != "weights")
        {
          cout << "Error :: config file entry following nroots must be weights" << endl;
          cout << "entry found was " << weights_str << endl;
          abort();
        }
        m_weights.resize(m_nroots);
        double norm = 0.;
        for(int i=0;i<m_nroots;++i)
        {
          conf_file >> weights_str;
          m_weights[i] = atof(weights_str.c_str());  
          norm += m_weights[i];
        }  
        for(int i=0;i<m_nroots;++i)
          m_weights[i] /= norm;
        if(m_nroots > 1)
          m_solve_type = BIG_DAVIDSON;
        cout << "Using the following weightings..." << endl;
        for(int i=0;i<m_nroots;++i)
          cout << "\t State[" << i << "] :: " << m_weights[i] << endl;
      }
      else if (next_entry == "guess_permutations")
      {
	conf_file >> m_guess_permutations;
      }
      else if (next_entry == "nquanta")
      {
        std::string nquanta_str;
        conf_file >> nquanta_str;
        m_nquanta = atoi (nquanta_str.c_str());
        pout << "\t\t\t Maximum states per quanta for initial Slater guess :: " << m_nquanta << endl;
      }
      else if (next_entry == "sysadd")
      {
        std::string sys_add_str;
        conf_file >> sys_add_str;
        m_sys_add = atoi (sys_add_str.c_str());
        pout << " \t\t\t System dot size :: " << m_sys_add << endl;
      }
      else if (next_entry == "envadd")
      {
        std::string env_add_str;
        conf_file >> env_add_str;
        m_env_add = atoi (env_add_str.c_str());
        pout << "\t\t\t Environment dot size :: " << m_env_add << endl;
      }
      else if (next_entry == "dofci")
      {
        m_do_fci = true;
      }
      else if (next_entry == "docd")
      {
        m_do_cd = true;
      }
      else if(next_entry == "deflation_min_size")
      {
        std::string deflation_min_size_str;
        conf_file >> deflation_min_size_str;
        m_deflation_min_size = atoi(deflation_min_size_str.c_str());
      }
      else if(next_entry == "deflation_max_size")
      {
        std::string deflation_max_size_str;
        conf_file >> deflation_max_size_str;
        m_deflation_max_size = atoi(deflation_max_size_str.c_str());
      }
      else if(next_entry == "maxiter")
      {
        std::string maxiter_str;
        conf_file >> maxiter_str;
        m_maxiter = atoi(maxiter_str.c_str());
        cout << "Maximum sweep iterations :: " << m_maxiter << endl;
      }
      else if(next_entry == "thrds_per_node")
      {
        std::string thrds_per_node_str;
#ifndef SERIAL
        mpi::communicator world;
        for(int i=0;i<world.size();++i)
        {
          conf_file >> thrds_per_node_str;
          m_thrds_per_node[i] = atoi(thrds_per_node_str.c_str());
          cout << "Maximum threads for node[" << i << "] :: " << m_thrds_per_node[i] << endl;
        }
#else
        for(int i=0;i<1;++i)
        {
          conf_file >> thrds_per_node_str;
          m_thrds_per_node[i] = atoi(thrds_per_node_str.c_str());
          cout << "Maximum threads for node[" << i << "] :: " << m_thrds_per_node[i] << endl;
        }
#endif
      }
      else if(next_entry == "screen_tol")
      {
        std::string screen_tol_str;
        conf_file >> screen_tol_str;
        m_screen_tol = atof(screen_tol_str.c_str());
        cout << "Screen tolerance :: " << m_screen_tol << endl;
      }
      else if(next_entry == "onedot")
      {
        m_algorithm_type = ONEDOT;
        m_env_add = 0;
        cout << "Using one dot DMRG algorithm" << endl;
      }
      else if(next_entry == "twodot")
      {
        m_algorithm_type = TWODOT;
        m_env_add = 1;
        cout << "Using two dot DMRG algorithm" << endl;
      }
      else if(next_entry == "twodot_to_onedot")
      {
        m_algorithm_type = TWODOT_TO_ONEDOT;
        m_env_add = 1;
        std::string twodot_to_onedot_iter_str;
        conf_file >> twodot_to_onedot_iter_str;
        m_twodot_to_onedot_iter = atoi(twodot_to_onedot_iter_str.c_str());
        cout << "Starting with two dot DMRG algorithm and switching to one dot DMRG algorithm on iteration " << m_twodot_to_onedot_iter << endl;
      }
      else if (next_entry == "sweep_tol")
	{
	  conf_file >> m_sweep_tol;
	}
      else if (next_entry == "outputlevel")
	{
	  conf_file >> m_outputlevel;
	}
      else if (next_entry == "restart")
	{
	  m_restart = true;
	}
      else if (next_entry == "fullrestart")
	{
	  m_fullrestart = true;	  
	}
      else if (next_entry == "restartwarm")
	{
	  m_restart_warm = true;
	}
      else if (next_entry == "reset_iterations")
	{
	  m_reset_iterations = true;
	}
      else if (next_entry== "oneintegral")
	{
	  conf_file >> m_oneintegral;
	}
      else if (next_entry == "twointegral")
	{
	  conf_file >> m_twointegral;
	}
      else
      {
        cout << "Unrecognized option :: " << next_entry << endl;
        abort();
      }
    }
    conf_file.close();
  }

#ifndef SERIAL
  boost::mpi::communicator world;
  mpi::broadcast(world,*this,0);
  mpi::broadcast(world, sym, 0);
  if (sym != "c1")
    Symmetry::InitialiseTable(sym);
  mpi::broadcast(world,symmetry,0);
#endif
}

int SpinAdapted::Input::nroots(int sweep_iter) const
{
  int iter;
  int current = 0;
  for (iter = 0; iter < dmrginp.sweep_iter_schedule().size(); ++iter)
  {
     if (sweep_iter >= dmrginp.sweep_iter_schedule()[iter]) 
       current = iter; 
     
  }

  int nroots = m_nroots;
  if (m_noise_type == EXCITEDSTATE && m_sweep_additional_noise_schedule[current] != 0.0 && nroots == 1)
    nroots++;
  return nroots;
}

std::vector<double> SpinAdapted::Input::weights(int sweep_iter) const
{
  int iter;
  int current = 0;
  for (iter = 0; iter < dmrginp.sweep_iter_schedule().size(); ++iter)
  {
     if (sweep_iter >= dmrginp.sweep_iter_schedule()[iter]) 
       current = iter; 
     
  }

  int nroots = this->nroots(sweep_iter);
  std::vector<double> weights(nroots);
  for (int i=0; i< m_nroots; i++)
    weights[i] = m_weights[i];
  if (nroots == m_nroots+1)
    weights[nroots-1] = m_sweep_additional_noise_schedule[current];
  if (nroots != m_nroots+1 && nroots != m_nroots)
  {
    pout << "Something wrong with the number of weights of wavefunctions!"<<endl;exit(0);
  }
  return weights;
}

int SpinAdapted::Input::getHFQuanta(const SpinBlock& b) const
{
  const vector<int>& sites = b.get_sites();
  int n=0, s=0;
  for (int i=0; i<sites.size(); i++)
  {
    if (m_hf_occupancy[ 2*sites[i] ])
    {
      n++; s++;
    }
    if (m_hf_occupancy[ 2*sites[i]+1])
    {
      n++; s--;
    }
  }

  const StateInfo& sI = b.get_stateInfo();
  for(int i=0; i<sI.quanta.size(); i++)
    if (sI.quanta[i].get_n() == n && sI.quanta[i].get_s() == abs(s))
      return i;
  return 0;
}
