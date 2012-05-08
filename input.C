#include <iostream>
#include <fstream>
#include <communicate.h>
#include "Symmetry.h"
#include "global.h"
#include "MatrixBLAS.h"
#include "spinblock.h"
#include "couplingCoeffs.h"
#include <boost/tokenizer.hpp>
#include <string.h>
#include <ctype.h>
#include <boost/algorithm/string.hpp>
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
using namespace std;

namespace SpinAdapted {
string sym;
}

void SpinAdapted::Input::ReadMeaningfulLine(ifstream& input, string& msg, int msgsize)
{
  bool readmore = true;
  while (readmore && !input.eof()) {
    char msgctr[msgsize];
    input.getline(msgctr, msgsize+1);

    msg=string(msgctr);
    if(msg.size() == msgsize) {
      cerr << "in the process of reading line begining with "<<endl<<msg<<endl;
      cerr<< "this line is too long"<<endl;
      abort();
    }
    int pos = msg.find("!");
    msg = msg.substr(0, pos);
    trim(msg);
    if (msg.size() != 0)
      readmore = false;
  }
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

  m_norbs = 0;
  m_alpha = 0;
  m_beta = 0;

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
  m_screen_tol = 1.00e-20;

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

SpinAdapted::Input::Input(const string& config_name)
{
  //first collect all the data
  sym = "c1";
  if(mpigetrank() == 0)
  {
    pout << "Reading input file"<<endl;
    bool PROVIDED_WEIGHTS = false;

    initialize_defaults();

    ifstream input(config_name.c_str());

    string msg; int msgsize = 5000;
    ReadMeaningfulLine(input, msg, msgsize);
    while(msg.size() != 0) {
      vector<string> tok;
      boost::split(tok, msg, is_any_of(" \t"), token_compress_on);
      string keyword = *tok.begin();

      if (boost::iequals(keyword,  "norbs") || boost::iequals(keyword,  "norb") || boost::iequals(keyword,  "norbitals"))
      {
	if (tok.size() != 2) {
	  cerr << "keyword norbs should be followed by a single integer and then an end line"<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}
	m_norbs = 2*atoi(tok[1].c_str());;
	//set some default values 
	m_num_spatial_orbs = m_norbs/2;
	m_spatial_to_spin.resize(m_norbs/2+1);
	m_spin_to_spatial.resize(m_norbs+1);
	m_spin_orbs_symmetry.resize(m_norbs);
	for (int i=0; i<m_norbs; i++) {
	  m_spatial_to_spin[i/2] = i - i%2;
	  m_spin_to_spatial[i] = i/2;
	  m_spin_orbs_symmetry[i] = 0;
	}
	m_spatial_to_spin[m_norbs/2] = m_norbs;
	m_spin_to_spatial[m_norbs] = m_norbs;
      }


      else if (boost::iequals(keyword,  "schedule"))
      {
	if (tok.size() != 2) {
	  cerr << "keyword schedule should be followed by a single number and then an end line"<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}	
	int nentry = atoi(tok[1].c_str());

	m_sweep_iter_schedule.resize(nentry);
	m_sweep_state_schedule.resize(nentry);
	m_sweep_qstate_schedule.resize(nentry);
	m_sweep_tol_schedule.resize(nentry);
	m_sweep_noise_schedule.resize(nentry);
	m_sweep_additional_noise_schedule.resize(nentry);


	for (int i = 0; i < nentry; ++i)
	  {
	    ReadMeaningfulLine(input, msg, msgsize);
	    vector<string> schd_tok;
	    boost::split(schd_tok, msg, is_any_of(" \t"), token_compress_on);

	    if (schd_tok.size() != 4) {
	      cerr << "Each line of the schedule contain four entries sweep_iteration   #retained states   davidson tolerance     noise"<<endl;
	      cerr << "error found at the following line "<<endl;
	      cerr<< msg<<endl;
	      abort();
	    }

	    m_sweep_iter_schedule[i] = atoi(schd_tok[0].c_str());
	    m_sweep_state_schedule[i] = atoi(schd_tok[1].c_str());
	    m_sweep_qstate_schedule[i] = 0;  //DEPRECATED OPTION
	    m_sweep_tol_schedule[i] = atof(schd_tok[2].c_str());
	    m_sweep_noise_schedule[i] = atof(schd_tok[3].c_str());
	    m_sweep_additional_noise_schedule[i] = 0.0;  //DEPRECATED OPTION

	    if (i>0 && m_sweep_iter_schedule[i] <= m_sweep_iter_schedule[i-1]) {
	      cerr << "Sweep iteration at a given line should be higher than the previous sweep iteration"<<endl;
	      cerr << "this sweep iteration "<<m_sweep_iter_schedule[i] <<endl;
	      cerr << "previous sweep iteration "<<m_sweep_iter_schedule[i-1]<<endl;
	      cerr << "error found in the following line "<<endl;
	      cerr << msg<<endl;
	      abort();
	    }
	  }

      }


      else if (boost::iequals(keyword,  "sym") || boost::iequals(keyword,  "symmetry"))
      {
	m_spatial_to_spin.clear();
	m_spin_to_spatial.clear();

	m_num_spatial_orbs = 0;
        //string sym;
	if (tok.size() !=  2) {
	  cerr << "keyword sym should be followed by a single string and then an end line"<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}	
        sym = tok[1];
	boost::algorithm::to_lower(sym); //store as lower case string
        Symmetry::InitialiseTable(sym);

	ReadMeaningfulLine(input, msg, msgsize);
	boost::split(tok, msg, is_any_of(" \t"), token_compress_on);

	if (boost::iequals(tok[0], "orbsym") || boost::iequals(tok[0], "orbitalsymmmetry") || boost::iequals(tok[0], "orbsymmetry") || 
	    boost::iequals(tok[0], "orb_sym") || boost::iequals(tok[0], "orbital_symmmetry") || boost::iequals(tok[0], "orb_symmetry"))
	{
	  if (m_norbs == 0) {
	    cerr << "Need to define the number of orbitals before specifying their irreducible representations"<<endl;
	    abort();
	  }
	  if (tok.size() != m_norbs/2+1) {
	    cerr << "keyword orbsym should be followed by "<<m_norbs<<" irreducible representations for each orbital and then an endline"<<endl;
	    cerr << "error found in the following line "<<endl;
	    cerr << msg<<endl;
	    abort();
	  }

	  for (int i=0; i < 2*(tok.size()-1); i+=2)
	  {
	    int ir = atoi(tok[i/2+1].c_str());
	    m_spin_orbs_symmetry[i] = ir;
	    m_spin_orbs_symmetry[i+1] = ir;

	    if (sym == "dinfh") {
	      if (ir < -1) {
		m_num_spatial_orbs ++;
		m_spatial_to_spin.push_back(i);
	      }
	      else if( (ir == 0 || ir ==1)) {
		m_num_spatial_orbs++;
		m_spatial_to_spin.push_back(i);
	      }
	      m_spin_to_spatial[i] = m_num_spatial_orbs-1;
	      m_spin_to_spatial[i+1] = m_num_spatial_orbs-1;
	    }
	    else {
	      m_num_spatial_orbs++;
	      m_spatial_to_spin.push_back(i);

	      m_spin_to_spatial[i] = m_num_spatial_orbs-1;
	      m_spin_to_spatial[i+1] = m_num_spatial_orbs-1;
	    }
	  }
	  m_spatial_to_spin.push_back(m_norbs);
	  m_spin_to_spatial.push_back(m_norbs);
	}
	else {
	  cerr << "Symmetry information of the molecule should follow by orbsym keyword and irreducible representations of all the orbitals"<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}		  
      }


      else if (boost::iequals(keyword,  "wave_sym") || boost::iequals(keyword,  "wave_symmetry"))
      {
	if (tok.size() < 2 || tok.size() >= 6) {
	  cerr << "keyword wave_sym should be followed by symmetry of the wavefunction in the form of nelec:spin:irrep and then an end line."<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}	
	string wavesymmetry = "";

	for (int i=1; i<tok.size(); i++)
	  wavesymmetry.append(tok[i]);

	vector<string > symtokens;
	boost::split(symtokens, wavesymmetry, is_any_of(": \t"), token_compress_on);

	if (symtokens.size() != 2 && symtokens.size() != 3) {
	  cerr<<"Symmetry string needs to have the form nelec:spin or nelec:spin:irrep"<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}
	string nelecstr = symtokens[0]; trim(nelecstr);
	int n_elec = atoi(nelecstr.c_str());

	string nspinstr = symtokens[1]; trim(nspinstr);
	int n_spin = atoi(nspinstr.c_str());

	if ( (n_elec-n_spin)%2 != 0) {
	  cerr<< "cannot have a spin of "<<n_spin<<"  with "<<n_elec<<" electrons "<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}
	if (symtokens.size() == 3) {
	  string nirrepstr = symtokens[2]; trim(nirrepstr);
	  m_total_symmetry_number = IrrepSpace(atoi(nirrepstr.c_str()));
	}
	else 
	  m_total_symmetry_number = IrrepSpace(0);
	m_alpha = (n_elec + n_spin)/2;
	m_beta = (n_elec - n_spin)/2;
        m_molecule_quantum = SpinQuantum(m_alpha + m_beta, m_alpha - m_beta, m_total_symmetry_number);
      }


      else if (boost::iequals(keyword,  "hubbard"))
	m_ham_type = HUBBARD;
      else if (boost::iequals(keyword,  "dmrg"))
	m_calc_type = DMRG;
      else if (boost::iequals(keyword,  "maxj")) {
	if (tok.size() !=  2) {
	  cerr << "keyword maxj should be followed by a single integer and then an end line."<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}	
        m_maxj = atoi(tok[1].c_str());
      }
      else if (boost::iequals(keyword,  "fci"))
	m_calc_type = FCI;
      else if (boost::iequals(keyword,  "genblock")|| boost::iequals(keyword,  "genblocks") || boost::iequals(keyword,  "generateblock"))
	m_calc_type = GENBLOCK;
      else if (boost::iequals(keyword,  "onepdm"))
	m_calc_type = ONEPDM;
      else if (boost::iequals(keyword,  "twopdm"))
	m_calc_type = TWOPDM;
      else if(boost::iequals(keyword,  "prefix") || boost::iequals(keyword,  "scratch"))
      {
	m_load_prefix = tok[1];
	m_save_prefix = m_load_prefix;
      }


      else if(boost::iequals(keyword,  "nroots"))
      {
        std::string nroots_str;
	if (tok.size() != 2) {
	  cerr << "keyword nroots should be followed by a single integer and then an end line."<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}	

        m_nroots = atoi(tok[1].c_str());
        if(m_deflation_min_size < m_nroots)
          m_deflation_min_size = m_nroots;

	ReadMeaningfulLine(input, msg, msgsize);
	vector<string> weighttoken;
	boost::split(weighttoken, msg, is_any_of(" \t"), token_compress_on);
	
	if (!boost::iequals(weighttoken[0],  "weights")) 
	  continue;
	PROVIDED_WEIGHTS = true;

        m_weights.resize(m_nroots);

	if (weighttoken.size() != m_nroots +1 ) {
	  cerr << "keyword weights should be followed by floating point numbers providing weights for "<<m_nroots<<" states."<<endl;
	  cerr << "You could chose to omit the keyworkd weights in which case the weights will be distributed uniformly between the different roots"<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}
        double norm = 0.;
	
        for(int i=1;i<weighttoken.size();i++)
        {
          m_weights[i-1] = atof(weighttoken[i].c_str());  
          norm += m_weights[i-1];
	  if (m_weights[i-1] <1e-10) {
	    cerr<< "Weight of a state cannot be less than 1e.0e-10"<<endl;
	    cerr << "error found in the following line "<<endl;
	    cerr << msg<<endl;
	    abort();
	  }
        }  
	if (norm <= 1.e-10) {
	  cerr<< "Weights should add up to approximately 1.0. Currently they add up to "<<norm<<endl;
	  abort();
	}

        for(int i=0;i<m_nroots;++i)
          m_weights[i] /= norm;
      }


      else if (boost::iequals(keyword,  "docd") || boost::iequals(keyword,  "do_cd"))
      {
        m_do_cd = true;
      }


      else if(boost::iequals(keyword,  "deflation_max_size") || boost::iequals(keyword,  "max_deflation_size"))
      {
	if (tok.size() !=  2) {
	  cerr << "keyword "<<keyword<<" should be followed by a single number and then an endline"<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}
        m_deflation_max_size = atoi(tok[1].c_str());
      }


      else if(boost::iequals(keyword,  "maxiter"))
      {
	if (tok.size() !=  2) {
	  cerr << "keyword maxiter should be followed by a single integer and then an endline"<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}
        m_maxiter = atoi(tok[1].c_str());
      }


      else if(boost::iequals(keyword,  "screen_tol") || boost::iequals(keyword,  "screen_tolerance") || boost::iequals(keyword,  "screening_tol") || boost::iequals(keyword,  "screening_tolerance"))
      {
	if (tok.size() != 2) {
	  cerr << "keyword screen_tol should be followed by a single number and then an endline"<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}
        m_screen_tol = atof(tok[1].c_str());
      }


      else if(boost::iequals(keyword,  "onedot"))
      {
        m_algorithm_type = ONEDOT;
        m_env_add = 0;
      }


      else if(boost::iequals(keyword,  "twodot"))
      {
        m_algorithm_type = TWODOT;
        m_env_add = 1;
      }


      else if(boost::iequals(keyword,  "twodot_to_onedot"))
      {
	if (tok.size() !=  2) {
	  cerr << "keyword twodot_to_onedot should be followed by a single number and then an endline"<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}
        m_algorithm_type = TWODOT_TO_ONEDOT;
        m_env_add = 1;
        m_twodot_to_onedot_iter = atoi(tok[1].c_str());
      }


      else if (boost::iequals(keyword,  "sweep_tol") || boost::iequals(keyword,  "sweep_tolerance"))
      {
	if (tok.size() !=  2) {
	  cerr << "keyword sweep_tol should be followed by a single number and then an endline"<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}
        m_sweep_tol = atof(tok[1].c_str());
      }


      else if (boost::iequals(keyword,  "outputlevel")) {
	if (tok.size() != 2) {
	  cerr << "keyword outputlevel should be followed by a single integer and then an endline"<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}
        m_outputlevel = atoi(tok[1].c_str());
      }


      else if (boost::iequals(keyword,  "restart")) {
	m_restart = true;
      }


      else if (boost::iequals(keyword,  "fullrestart")) {
	m_fullrestart = true;	  
      }


      else if (boost::iequals(keyword,  "reset_iterations") || boost::iequals(keyword,  "reset_iter") || boost::iequals(keyword,  "reset_iters")) {
	m_reset_iterations = true;
      }


      else if (boost::iequals(keyword,  "oneintegral")) {
	if (tok.size() !=  2) {
	  cerr << "keyword oneintegral should be followed by the name of the file containing one electron integrals."<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}
	m_oneintegral =tok[1];
      }


      else if (boost::iequals(keyword,  "twointegral")) {
	if (tok.size() !=  2) {
	  cerr << "keyword twointegral should be followed by the name of the file containing two electron integrals."<<endl;
	  cerr << "error found in the following line "<<endl;
	  cerr << msg<<endl;
	  abort();
	}
	m_twointegral = tok[1];
      }
      else
      {
        pout << "Unrecognized option :: " << keyword << endl;
	cerr << "error found in the following line "<<endl;
	cerr << msg<<endl;
        abort();
      }
      msg.resize(0);
      ReadMeaningfulLine(input, msg, msgsize);
      
    }

  }

#ifndef SERIAL
  boost::mpi::communicator world;
  mpi::broadcast(world, sym, 0);
  if (sym != "c1")
    Symmetry::InitialiseTable(sym);
  if (mpigetrank() == 0) {
#endif
    
    pout << "Checking Input for errors"<<endl;
    performSanityTest();
    pout << "Summary of Input"<<endl;
    pout << "----------------"<<endl;
    writeSummary();
    pout << endl;
#ifndef SERIAL
  }
  mpi::broadcast(world,*this,0);
#endif

}

void SpinAdapted::Input::writeSummary()
{
#ifndef SERIAL
  if (mpigetrank() == 0) {
#endif
  printf("%-50s :   %-i\n", "Total number of orbitals", m_norbs/2);
  printf("%-50s :   %-i:%-i:%-i\n", "Symmetry of the targetted wavefunctions",m_alpha + m_beta, m_alpha - m_beta, m_total_symmetry_number.getirrep());
  printf("%-50s :   %-i\n", "Number of wavefunctions targetted", m_nroots);
  if (m_nroots >1) {
    printf("%-50s :   ", "The weights of the wavefunctions");
    for (int i=0; i<m_nroots; i++) 
      printf("%-10.2e", m_weights[i]);
    printf("\n");
  }
  printf("%-50s :   ", "Symmetry of the molecule"); cout << sym<<endl;;  
  if (sym != "c1") {
    printf("%-50s :   ", "Irreducible representations of the orbitals");
    for (int i=0; i<m_spin_orbs_symmetry.size(); i+=2) 
      cout << Symmetry::stringOfIrrep(m_spin_orbs_symmetry[i])<<"  ";
    printf("\n");
  }

  printf("\nSchedule\n");
  printf("--------\n");
  printf("%-10s : %-20s  %-20s  %-20s\n", "Iter", "# States", "Davidson_tol",  "Random_noise");
  for (int i=0; i<m_sweep_iter_schedule.size(); i++)
    printf("%-10i : %-20i  %-20.4e  %-20.4e\n", m_sweep_iter_schedule[i], m_sweep_state_schedule[i], m_sweep_tol_schedule[i], m_sweep_noise_schedule[i]);
  if (m_algorithm_type == TWODOT_TO_ONEDOT) 
    printf("%-50s :   %-i\n", "Switching from twodot to onedot algorithm", m_twodot_to_onedot_iter);
    
  printf("%-50s :   %-i\n", "Maximum sweep iterations", m_maxiter);
#ifndef SERIAL
  }
#endif
  
}

void SpinAdapted::Input::performSanityTest()
{
#ifndef SERIAL
  if (mpigetrank() == 0) {
#endif
  if (m_norbs <= 0) {
    cerr << "total number of orbitals has to be a positive number"<<endl;
    abort();
  }
  if (m_norbs > 100) {
    cerr << "Number of orbitals cannot be greater than 100"<<endl;
    abort();
  }
  if (m_alpha+m_beta <= 0) {
    cerr << "Total number of electrons cannot be negative"<<endl;
    abort();
  }
  if (m_alpha < m_beta) {
    cerr << "DMRG requires the spin to a positive number and less than the total number of electrons"<<endl;
    abort();
  }
  if (m_norbs < m_alpha+m_beta) {
    cerr<< "No of spin orbitals has to be greater than total number of electrons"<<endl;
    abort();
  }
  for (int i=0; i<m_spin_orbs_symmetry.size(); i+=2) {
    Symmetry::irrepAllowed(m_spin_orbs_symmetry[i]);
  }
  if (m_algorithm_type == TWODOT_TO_ONEDOT && m_twodot_to_onedot_iter >= m_maxiter) {
    cerr << "Switch from twodot to onedot algorithm cannot happen after maxiter"<<endl;
    cerr << m_twodot_to_onedot_iter <<" < "<<m_maxiter<<endl;
    abort();
  }
  if (m_maxiter <= m_sweep_iter_schedule.back()) {
    cerr << "maximum iterations allowed is less than or equal to the last sweep iteration in your schedule."<<endl;
    cerr << m_maxiter <<" <= "<< (m_sweep_iter_schedule.back())<<endl;
    cerr << "either increase the max_iter or reduce the number of sweeps"<<endl;
    abort();
  }
  Symmetry::irrepAllowed(m_total_symmetry_number.getirrep());
  //this is important so the user cannot break the code
#ifndef SERIAL
  }
#endif
  if (m_norbs <= 6)
    m_calc_type = TINYCALC;

  //make some initial hf guess
  m_spin_vector.resize(m_norbs);
  for (int i = 0; i < m_norbs; ++i)
    m_spin_vector[i] = (i & 1) ? -1 : 1;
  m_hf_occupancy.resize(m_norbs);
  for (int i = 0; i < m_alpha; ++i)
    m_hf_occupancy[2*i] = 1;
  for (int i = 0; i < m_beta; ++i)
    m_hf_occupancy[2*i+1] = 1;

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
