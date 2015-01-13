/*                                                                           
Developed by Sandeep Sharma, Roberto Olivares-Amaya and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma, Roberto Olivares-Amaya and Garnet K.-L. Chan
*/


#include <iostream>
#include <fstream>
#include <communicate.h>
#include "Symmetry.h"
#include "global.h"
#include "MatrixBLAS.h"
#include "spinblock.h"
#include "couplingCoeffs.h"
#include "genetic/GAOptimize.h"
#include "genetic/ReadIntegral.h"
#include <boost/tokenizer.hpp>
#include <string.h>
#include <ctype.h>
#include <boost/algorithm/string.hpp>
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "fiedler.h"
#include "pario.h"

using namespace std;

namespace SpinAdapted {
string sym;
bool NonabelianSym;
}
void CheckFileExistence(string filename, string filetype);
void CheckFileInexistence(string filename, string filetype);

void SpinAdapted::Input::ReadMeaningfulLine(ifstream& input, string& msg, int msgsize)
{
  bool readmore = true;
  while (readmore && !input.eof()) {
    char msgctr[msgsize];
    input.getline(msgctr, msgsize+1);

    msg=string(msgctr);
    if(msg.size() == msgsize) {
      pout << "in the process of reading line beginning with "<<endl<<msg<<endl;
      pout<< "this line is too long"<<endl;
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
  m_algorithm_type = TWODOT_TO_ONEDOT;
  m_noise_type = RANDOM;
  m_calc_type = DMRG;
  m_solve_type = DAVIDSON;
  m_stateSpecific = false;
  m_implicitTranspose = true; //dont make DD just use CC^T to evaluate it
  m_occupied_orbitals = -1;
  m_num_Integrals = 1;
  v_2.resize(1, TwoElectronArray(TwoElectronArray::restrictedNonPermSymm));
  v_1.resize(1);
  coreEnergy.resize(1);
  m_baseState.resize(1,0);
  m_projectorState.resize(0);
  m_targetState = -1;
  
  m_spinAdapted = true;
  m_Bogoliubov = false;
  m_sys_add = 1;
  m_env_add = 1;
  
  m_twodot_to_onedot_iter = 0;
  m_integral_disk_storage_thresh = 1000;
  m_max_lanczos_dimension = 5000;

  m_norbs = 0;
  m_alpha = 0;
  m_beta = 0;
  m_hf_occ_user = "";

  m_outputlevel = 0;
  m_nquanta = 2;
  m_total_symmetry_number = IrrepSpace( 0 );
  m_total_spin = 0;
  m_guess_permutations = 10;

  m_direct = true;
  m_nroots = 1;
  m_weights.resize(m_nroots);
  m_weights[0] = 1.;

  m_deflation_min_size = 2;
  m_deflation_max_size = 20;

  m_add_noninteracting_orbs = true;
  m_no_transform = false;
  m_do_fci = false;
  m_do_npdm_ops = false;
  m_do_npdm_in_core = false;
  m_new_npdm_code = false;
  m_do_pdm = false;
  m_store_spinpdm = false;
  m_spatpdm_disk_dump = false;
  m_store_nonredundant_pdm =false;
  m_pdm_unsorted = false;
  m_npdm_intermediate= true;
  m_npdm_multinode= true;
 
  m_maxiter = 10;
  m_oneindex_screen_tol = NUMERICAL_ZERO;
  m_twoindex_screen_tol = NUMERICAL_ZERO;

  m_load_prefix = ".";
  m_save_prefix = ".";

  m_maxj = 15;
  m_ninej.init(m_maxj);
  m_set_Sz = false;

  n_twodot_noise = 0;
  m_twodot_noise = 1.0e-4;
  m_twodot_gamma = 3.0e-1;

  m_sweep_tol = 1.0e-5;
  m_restart = false;
  m_fullrestart = false;
  m_restart_warm = false;
  m_backward = false;
  m_reset_iterations = false;

  m_do_diis = false;
  m_diis_error = 1e-2;
  m_start_diis_iter = 8;
  m_diis_keep_states = 6;
  m_diis_error_tol = 1e-8;
  m_num_spatial_orbs = 0;
  m_schedule_type_default = false;
  m_schedule_type_backward = false;
  m_maxM = 0;
  m_lastM = 500;
  m_startM = 250;
  
  m_calc_ri_4pdm=false;
  m_store_ripdm_readable=false;
  m_nevpt2 = false;
  m_conventional_nevpt2 = false;
  m_kept_nevpt2_states = -1;
  NevPrint.first=false;
  NevPrint.second=0;

  //reorder options, by default it does fiedler
  m_reorderType = FIEDLER;
  m_reorderfile = "";
  m_gaconffile = "default";

  m_orbformat=MOLPROFORM;

  m_warmup = LOCAL0;
}

void SpinAdapted::Input::usedkey_error(string& key, string& line) {
  pout << "keyword "<<key<<" already used once"<<endl;
  pout << "error found at line "<<endl;
  pout << line<<endl;
  abort();
}

SpinAdapted::Input::Input(const string& config_name) {
  //first collect all the data
  std::vector<int> usedkey(NUMKEYWORDS, -1);
  int n_elec = -1;
  int n_spin = -1;
  sym = "c1";
  NonabelianSym = false;
  std::vector<string> orbitalfile;


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
      boost::split(tok, msg, is_any_of(", \t"), token_compress_on);
      string keyword = *tok.begin();

      if (boost::iequals(keyword,  "orbs") || boost::iequals(keyword,  "orbitals") || boost::iequals(keyword, " orbitals"))
      {
	if(usedkey[ORBS] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[ORBS] = 0;
	if (tok.size() < 2) {
	  pout << "keyword orbs should be followed by atleast a single filename and then an end line"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}
	m_num_Integrals = tok.size()-1;
	orbitalfile.resize(m_num_Integrals);
	for (int l=0; l<m_num_Integrals; l++) 	  
	  orbitalfile[l] = tok[l+1];
      }
      else if (boost::iequals(keyword, "maxM")) {
	if(usedkey[MAXM] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[MAXM] = 0;
	if (tok.size() != 2) {
	  pout << "keyword maxM should be followed by a single number and then an end line"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
	m_maxM = atoi(tok[1].c_str());

	if (m_maxM <= 0) {
	  pout << "maxM cannot be less than equal to 0"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
	  
      }
      else if (boost::iequals(keyword, "nonspinadapted")) {
	if(usedkey[NONSPINADAPTED] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[NONSPINADAPTED] = 0;
	if (tok.size() !=  1) {
	  pout << "keyword nonspinadated is a stand alone keyword"<<endl;
	  pout << msg<<endl;
	  abort();
	}
	m_add_noninteracting_orbs = false;
        m_spinAdapted = false;
      }
      else if (boost::iequals(keyword, "bogoliubov")) {
    if(usedkey[BOGOLIUBOV] == 0)
      usedkey_error(keyword, msg);
    usedkey[BOGOLIUBOV] = 0;
	if (tok.size() !=  1) {
	  pout << "keyword bogoliubov is a stand alone keyword"<<endl;
	  pout << msg<<endl;
	  abort();
	}
    m_Bogoliubov = true;
    m_ham_type = BCS;
      }
      else if (boost::iequals(keyword, "warmup")) {
        if (usedkey[WARMUP] == 0)
          usedkey_error(keyword, msg);
        usedkey[WARMUP] = 0;
        if (tok.size() != 2) {
          pout << "must specify warmup type with keyword warmup" << endl;
          pout << msg << endl;
          abort();
        }
        if (boost::iequals(tok[1], "wilson")) {
          m_warmup = WILSON; // default option select the lowest energy slater determinants
        } else if (boost::iequals(tok[1], "local_0site")) {
          m_warmup = LOCAL0;
        } else if (boost::iequals(tok[1], "local_2site")) {
          m_warmup = LOCAL2;
        } else if (boost::iequals(tok[1], "local_3site")) {
          m_warmup = LOCAL3;
        } else if (boost::iequals(tok[1], "local_4site")) {
          m_warmup = LOCAL4;
        } else {
          pout << "warm-up algorithm not defined" << endl;
          pout << tok[1] << endl;
          abort();
        }
      }
      else if (boost::iequals(keyword, "startM")) {
	if(usedkey[STARTM] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[STARTM] = 0;
	if (tok.size() != 2) {
	  pout << "keyword startM should be followed by a single number and then an end line"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
	m_startM = atoi(tok[1].c_str());

	if (m_startM < 50) {
	  pout << "startM cannot be less than 50"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
	  
      }
      else if (boost::iequals(keyword, "statespecific")) {
	pout<<"--------------------------------------------------------------------"<<endl;
	pout << "WARNING: THIS OPTION IMPLIES THAT A PREVIOUS DMRG CALCULATION HAS ALREADY BEEN PERFORMED"<<endl;
	pout << "THIS CALCULATION WILL TAKE THE PREVIOUS WAVEFUNCTIONS AND REFINE THEM"<<endl;
	pout<<"--------------------------------------------------------------------"<<endl;
	m_stateSpecific = true;
      }
      else if (boost::iequals(keyword, "donttranspose")) {
	pout<<"--------------------------------------------------------------------"<<endl;
	pout << "DONT USE THIS OPTION BECAUSE IT SLOWS DOWN THE CODE."<<endl;
	pout<<"--------------------------------------------------------------------"<<endl;
	m_implicitTranspose = false;
      }
      else if (boost::iequals(keyword, "occ")) {
	if (tok.size() != 2) {
	  pout << "keyword occ should be followed by a single number and then an end line"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
	m_occupied_orbitals = atoi(tok[1].c_str());
      }

      else if (boost::iequals(keyword, "lastM")) {
	if(usedkey[LASTM] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[LASTM] = 0;
	if (tok.size() != 2) {
	  pout << "keyword lastM should be followed by a single number and then an end line"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
	m_lastM = atoi(tok[1].c_str());
      }
      else if (boost::iequals(keyword,  "reorder")) {
	if(usedkey[REORDER] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[REORDER] = 0;
	if (tok.size() != 2) {
	  pout << "keyword reorder should be followed by the filename and then an end line"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
	m_reorderType = MANUAL;
	m_reorderfile = tok[1];
      }
      else if (boost::iequals(keyword,  "gaopt"))
      {
        if (tok.size() != 2) {
          perr << "keyword gaopt should be followed by the filename and then an end line"<<endl;
          perr << "error found in the following line "<<endl;
          perr << msg<<endl;
          abort();
        }       
	m_reorderType = GAOPT;
	m_gaconffile = tok[1];
      }
      else if (boost::iequals(keyword, "fiedler")){
        if (tok.size() != 1) {
          perr << "keyword fiedler should not be followed by anything"<<endl;
          perr << "error found in the following line "<<endl;
          perr << msg<<endl;
          abort();
        }       
	m_reorderType = FIEDLER;
      }
      else if (boost::iequals(keyword, "noreorder") || boost::iequals(keyword, "nofiedler")) {
	m_reorderType = NOREORDER; 
      }

      else if (boost::iequals(keyword,  "schedule"))
      {
	if(usedkey[SCHEDULE] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[SCHEDULE] = 0;

	if (tok.size() != 2 && tok.size() != 1) {
	  pout << "keyword schedule should be followed by the keyword \"default\" and then an end line or just an endline"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}
	
	msg.resize(0);
	ReadMeaningfulLine(input, msg, msgsize);
	vector<string> schd_tok;
	boost::split(schd_tok, msg, is_any_of(" \t"), token_compress_on);
	if (tok.size() != 1) {
	  if (boost::iequals(tok[1], "default")) { 
	    m_schedule_type_default = true;
	    continue;
	  }
	}
	
	m_sweep_iter_schedule.resize(0);
	m_sweep_state_schedule.resize(0);
	m_sweep_qstate_schedule.resize(0);
	m_sweep_tol_schedule.resize(0);
	m_sweep_noise_schedule.resize(0);
	m_sweep_additional_noise_schedule.resize(0);

	int i = 0;
	while(!boost::iequals(schd_tok[0], "END"))
	{
	  
	  if (schd_tok.size() != 4) {
	    pout << "Each line of the schedule contain four entries sweep_iteration   #retained states   Davidson tolerance     noise"<<endl;
	    pout << "error found at the following line "<<endl;
	    pout<< msg<<endl;
	    abort();
	  }
	  
	  m_sweep_iter_schedule.push_back( atoi(schd_tok[0].c_str()));
	  m_sweep_state_schedule.push_back( atoi(schd_tok[1].c_str()));
	  m_sweep_tol_schedule.push_back( atof(schd_tok[2].c_str()));
	  m_sweep_noise_schedule.push_back( atof(schd_tok[3].c_str()));
	  
	  if (m_sweep_state_schedule[i] <= 0) {
	    pout << "Number of retained states cannot be less than 0"<<endl;
	    pout << "error found in the following line "<<endl;
	    pout << msg<<endl;
	  }
	  if (i>0 && m_sweep_iter_schedule[i] <= m_sweep_iter_schedule[i-1]) {
	    pout << "Sweep iteration at a given line should be higher than the previous sweep iteration"<<endl;
	    pout << "this sweep iteration "<<m_sweep_iter_schedule[i] <<endl;
	    pout << "previous sweep iteration "<<m_sweep_iter_schedule[i-1]<<endl;
	    pout << "error found in the following line "<<endl;
	    pout << msg<<endl;
	    abort();
	  }
	  i++;
	  ReadMeaningfulLine(input, msg, msgsize);
	  boost::split(schd_tok, msg, is_any_of(" \t"), token_compress_on);
	}
	
      }


      else if (boost::iequals(keyword,  "sym") || boost::iequals(keyword,  "symmetry")|| boost::iequals(keyword,  "point_group")|| boost::iequals(keyword,  "pg"))
      {
	if(usedkey[SYM] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[SYM] = 0;
	m_spatial_to_spin.clear();
	m_spin_to_spatial.clear();

	m_num_spatial_orbs = 0;
        //string sym;
	if (tok.size() !=  2) {
	  pout << "keyword sym should be followed by a single string and then an end line"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
        sym = tok[1];
	boost::algorithm::to_lower(sym); //store as lower case string
        Symmetry::InitialiseTable(sym);
      }
      else if (boost::iequals(keyword, "davidson")) 
	m_solve_type = DAVIDSON;
      else if (boost::iequals(keyword, "lanczos")) 
	m_solve_type = LANCZOS;
      else if (boost::iequals(keyword, "thrds_per_node") || boost::iequals(keyword, "threads_per_node")) {

	int nprocs = 1;
#ifndef SERIAL       
	mpi::communicator world;
	nprocs = world.size();
#endif
	if (tok.size() !=  nprocs+1) {
	  pout << "keyword number of threads for each of the "<<nprocs<<" processes should be specified!"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
	for (int i=1; i<tok.size(); i++)
	  m_thrds_per_node[i-1] = atoi(tok[i].c_str());
      }
      else if (boost::iequals(keyword,  "nelecs") || boost::iequals(keyword,  "nelec")) {
	if(usedkey[NELECS] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[NELECS] = 0;
	if (tok.size() !=  2) {
	  pout << "keyword nelec should be followed by a single number and then an end line"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
	n_elec = atoi(tok[1].c_str());
	if(n_elec < 0) {
	  pout << "Number of electrons cannot be less than 0."<<endl;
	  pout << "See the manual for further details."<<endl;
	  abort();
	}
      }
      else if (boost::iequals(keyword,  "nspin") || boost::iequals(keyword,  "spin")) {
	if(usedkey[SPIN] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[SPIN] = 0;
	if (tok.size() !=  2) {
	  pout << "keyword spin should be followed by a single number and then an end line"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
	n_spin = atoi(tok[1].c_str());
	if(n_spin < 0) {
	  pout << "Spin of the wavefunction cannot be less than 0."<<endl;
	  pout << "See the manual for further details."<<endl;
	  abort();
	}
      }
      else if (boost::iequals(keyword,  "nirrep") || boost::iequals(keyword,  "irrep")) {
	if(usedkey[IRREP] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[IRREP] = 0;
  // When there are 2 irrep number, it means that calcultions of transition density matrix between wavefunctions with different irrep.
  if (tok.size()==2 )
	  m_total_symmetry_number = IrrepSpace(atoi(tok[1].c_str())-1);
  else if (tok.size()==3 ){
	  m_bra_symmetry_number = IrrepSpace(atoi(tok[1].c_str())-1);
	  m_total_symmetry_number = IrrepSpace(atoi(tok[2].c_str())-1);
    m_transition_diff_spatial_irrep=true;
  }
  else{
	  pout << "keyword irrep should be followed by one or two numbers and then an end line"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
      }
      else if (boost::iequals(keyword,  "hubbard"))
	m_ham_type = HUBBARD;
      else if (boost::iequals(keyword,  "heisenberg"))
	m_ham_type = HEISENBERG;
      else if (boost::iequals(keyword,  "dmrg"))
	m_calc_type = DMRG;
      else if (boost::iequals(keyword,  "calcoverlap"))
	m_calc_type = CALCOVERLAP;
      else if (boost::iequals(keyword,  "calchamiltonian"))
	m_calc_type = CALCHAMILTONIAN;
      else if (boost::iequals(keyword,  "response")) {
	if (tok.size() != 1) {
	  pout << "The keyword response should not be followed by anything!"<<endl;
	  abort();
	}
	m_calc_type = RESPONSE;
	m_solve_type = CONJUGATE_GRADIENT;
	if(m_targetState == -1)
	  m_targetState = 2;
	
      }
      else if (boost::iequals(keyword,  "compress")) {
	m_calc_type = COMPRESS;
	if (tok.size() !=  2) {
	  pout << "keyword compression should be followed by a single number and then an end line"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
	m_baseState.resize(1, atoi(tok[1].c_str()));

      }
      else if (boost::iequals(keyword,  "maxj")) {
	if (tok.size() !=  2) {
	  pout << "keyword maxj should be followed by a single integer and then an end line."<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
        m_maxj = atoi(tok[1].c_str());
      }
      else if (boost::iequals(keyword,  "fci"))
	m_calc_type = FCI;
      else if (boost::iequals(keyword,  "onepdm") || boost::iequals(keyword,  "onerdm") || boost::iequals(keyword,  "ordm"))
	m_calc_type = ONEPDM;
      else if (boost::iequals(keyword,  "twopdm") || boost::iequals(keyword,  "twordm") || boost::iequals(keyword,  "trdm"))
        m_calc_type = TWOPDM;
      else if (boost::iequals(keyword,  "threepdm"))
	m_calc_type = THREEPDM;
      else if (boost::iequals(keyword,  "fourpdm"))
	m_calc_type = FOURPDM;
      else if (boost::iequals(keyword,  "nevpt2_npdm"))
	m_calc_type = NEVPT2PDM;
      else if (boost::iequals(keyword,  "restart_onepdm") || boost::iequals(keyword,  "restart_onerdm") || boost::iequals(keyword,  "restart_ordm"))
	m_calc_type = RESTART_ONEPDM;
      else if (boost::iequals(keyword,  "restart_twopdm") || boost::iequals(keyword,  "restart_twordm") || boost::iequals(keyword,  "restart_trdm"))
	m_calc_type = RESTART_TWOPDM;
      else if (boost::iequals(keyword,  "restart_threepdm") || boost::iequals(keyword,  "restart_threerdm") )
	m_calc_type = RESTART_THREEPDM;
      else if (boost::iequals(keyword,  "restart_fourpdm") || boost::iequals(keyword,  "restart_fourrdm") )
	m_calc_type = RESTART_FOURPDM;
      else if (boost::iequals(keyword,  "transition_onepdm") || boost::iequals(keyword,  "transition_onerdm") || boost::iequals(keyword,  "tran_onepdm"))
	m_calc_type = TRANSITION_ONEPDM;
      else if (boost::iequals(keyword,  "transition_twopdm") || boost::iequals(keyword,  "transition_twordm") || boost::iequals(keyword,  "tran_twopdm"))
	m_calc_type = TRANSITION_TWOPDM;
      else if (boost::iequals(keyword,  "restart_tran_onepdm") || boost::iequals(keyword,  "restart_tran_onerdm") || boost::iequals(keyword,  "restart_tran_ordm"))
	m_calc_type = RESTART_T_ONEPDM;
      else if (boost::iequals(keyword,  "restart_tran_twopdm") || boost::iequals(keyword,  "restart_tran_twordm") || boost::iequals(keyword,  "restart_tran_trdm"))
	m_calc_type = RESTART_T_TWOPDM;
      else if (boost::iequals(keyword,  "nevpt2") || boost::iequals(keyword,  "pt2")){
	m_calc_type = NEVPT2;
        m_nevpt2 = true;
        m_transition_diff_spatial_irrep=false;
      }
      else if (boost::iequals(keyword,  "ripdm")){
	m_calc_type = NEVPT2;
        m_nevpt2 = false;
        m_transition_diff_spatial_irrep=false;
      }
      else if (boost::iequals(keyword,  "restart_nevpt2") || boost::iequals(keyword,  "restart_pt2")){
        m_calc_type = RESTART_NEVPT2;
        m_nevpt2 = true;
        m_transition_diff_spatial_irrep=false;
      }
      else if (boost::iequals(keyword,  "restart_ripdm")){
        m_calc_type = RESTART_NEVPT2;
        m_nevpt2 = false;
        m_transition_diff_spatial_irrep=false;
      }
      else if (boost::iequals(keyword,  "calc_ri4pdm") || boost::iequals(keyword,  "ri4pdm"))
      {
        m_calc_ri_4pdm = true;
      }
      else if (boost::iequals(keyword,  "ripdm_readable"))
      {
        m_store_ripdm_readable = true;
      }
      else if (boost::iequals(keyword, "conventional_nevpt2"))
      {
        m_conventional_nevpt2 = true;
      }
      else if (boost::iequals(keyword,  "M_nevpt2")){
	if (tok.size() != 2) {
	  pout << "keyword M_nevpt2 should be followed by a single float and then an end line."<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
          abort();
        }
        m_kept_nevpt2_states = atof(tok[1].c_str());
      }
      

      else if(boost::iequals(keyword,  "prefix") || boost::iequals(keyword,  "scratch"))
      {
	if(usedkey[PREFIX] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[PREFIX] = 0;
	m_load_prefix = tok[1];
	m_save_prefix = m_load_prefix;
      }

     // user-defined initial HF wave function occupancies, in spin orbital
      else if(boost::iequals(keyword, "hf_occ"))
      {
        if( usedkey[HF_OCC] == 0 ) 
           usedkey_error( keyword, msg );
        usedkey[HF_OCC] = 0;

        // the occupancies start from the second element of the tok string
        // ++it;
        vector<string> :: iterator it = ++tok.begin(); 
        // the occupancies start from the second element of the tok string
        if (tok.size() == 2 ) {
           m_hf_occ_user = (*it).c_str();
        }
        else{
           m_hf_occ_user = "manual";
           for( ; it != tok.end(); ++it ){
            int occ_i = atoi( (*it).c_str() );
            m_hf_occupancy.push_back( occ_i );
           }
        }
        pout << "m_hf_occ_user " << m_hf_occ_user << endl;

      }

      else if(boost::iequals(keyword,  "nroots"))
      {
	if(usedkey[NROOTS] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[NROOTS] = 0;
        std::string nroots_str;
	if (tok.size() != 2) {
	  pout << "keyword nroots should be followed by a single integer and then an end line."<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	

        m_nroots = atoi(tok[1].c_str());
        if(m_deflation_min_size < m_nroots) 
          m_deflation_min_size = m_nroots;
	  
	m_deflation_max_size = max(20, m_nroots+20);

	ReadMeaningfulLine(input, msg, msgsize);
	vector<string> weighttoken;
	boost::split(weighttoken, msg, is_any_of(" \t"), token_compress_on);
	
	if (!boost::iequals(weighttoken[0],  "weights")) {
	  //manually make all the weights equal
	  m_weights.resize(m_nroots);
	  for (int i=0; i<m_nroots; i++)
	    m_weights[i] = 1.0/m_nroots;
	  continue;
	}
	PROVIDED_WEIGHTS = true;

        m_weights.resize(m_nroots);

	if (weighttoken.size() != m_nroots +1 ) {
	  pout << "keyword weights should be followed by floating point numbers providing weights for "<<m_nroots<<" states."<<endl;
	  pout << "You could chose to omit the keyword weights in which case the weights will be distributed uniformly between the different roots"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}
        double norm = 0.;
	
        for(int i=1;i<weighttoken.size();i++)
        {
          m_weights[i-1] = atof(weighttoken[i].c_str());  
          norm += m_weights[i-1];
	  if (m_weights[i-1] <1e-10) {
	    pout<< "Weight of a state cannot be less than 1e.0e-10"<<endl;
	    pout << "error found in the following line "<<endl;
	    pout << msg<<endl;
	    abort();
	  }
        }  
	if (norm <= 1.e-10) {
	  pout<< "Weights should add up to approximately 1.0. Currently they add up to "<<norm<<endl;
	  abort();
	}

        for(int i=0;i<m_nroots;++i)
          m_weights[i] /= norm;
      }


      else if (boost::iequals(keyword,  "docd") || boost::iequals(keyword,  "do_npdm_ops"))
      {
        m_do_npdm_ops = true;
      }
      else if (boost::iequals(keyword,  "new_npdm_code"))
      {
        m_new_npdm_code = true;
      }
      else if (boost::iequals(keyword,  "do_npdm_in_core"))
      {
        m_do_npdm_in_core = true;
      }
      else if (boost::iequals(keyword, "store_spinpdm"))
      {
        m_store_spinpdm = true;
      }
      else if (boost::iequals(keyword, "disk_dump_pdm"))
      {
        m_spatpdm_disk_dump = true;
      }
      else if (boost::iequals(keyword, "nonredundant_pdm"))
      {
        m_store_nonredundant_pdm = true;
      }
      else if (boost::iequals(keyword, "pdm_unsorted"))
      {
        m_pdm_unsorted = true;
      }
      else if (boost::iequals(keyword, "npdm_intermediate"))
      {
        m_npdm_intermediate = true;
      }
      else if (boost::iequals(keyword, "npdm_no_intermediate"))
      {
        m_npdm_intermediate = false;
      }
      else if (boost::iequals(keyword, "npdm_multinode"))
      {
        m_npdm_multinode = true;
      }
      else if (boost::iequals(keyword, "npdm_no_multinode"))
      {
        m_npdm_multinode = false;
      }



      else if(boost::iequals(keyword,  "deflation_max_size") || boost::iequals(keyword,  "max_deflation_size"))
      {
	if (tok.size() !=  2) {
	  pout << "keyword "<<keyword<<" should be followed by a single number and then an endline"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}
        m_deflation_max_size = atoi(tok[1].c_str());
      }
      else if(boost::iequals(keyword,  "ProjectorStates") )
      {
	if (tok.size() < 2) {
	  pout << "keyword projectorstates should be followed by atleast a single number and then an endline"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}
	m_projectorState.resize(tok.size()-1);
	for (int l=0; l<tok.size()-1; l++)
	  m_projectorState[l] = atoi(tok[l+1].c_str());
      }
      else if(boost::iequals(keyword,  "TargetState") )
      {
	if (tok.size() !=  2) {
	  pout << "keyword "<<keyword<<" should be followed by a single number and then an endline"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}
        m_targetState = atoi(tok[1].c_str());
      }
      else if(boost::iequals(keyword,  "BaseStates") )
      {
	if (tok.size() < 2) {
	  pout << "keyword basestate should be followed by a single number and then an endline"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}
	m_baseState.resize(tok.size()-1);
	for (int l=0; l<tok.size()-1; l++)
	  m_baseState[l] = atoi(tok[l+1].c_str());
      }


      else if(boost::iequals(keyword,  "maxiter"))
      {
	if(usedkey[MAXITER] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[MAXITER] = 0;
	if (tok.size() !=  2) {
	  pout << "keyword maxiter should be followed by a single integer and then an endline"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}
        m_maxiter = atoi(tok[1].c_str());
      }

      else if(boost::iequals(keyword,  "screen_tol") || boost::iequals(keyword,  "screen_tolerance") || boost::iequals(keyword,  "screening_tol") || boost::iequals(keyword,  "screening_tolerance"))
      {
	if(usedkey[SCREEN_TOL] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[SCREEN_TOL] = 0;
	if (tok.size() != 2) {
	  pout << "keyword screen_tol should be followed by a single number and then an endline"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}
        m_oneindex_screen_tol = atof(tok[1].c_str());
        m_twoindex_screen_tol = atof(tok[1].c_str());
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
	  pout << "keyword twodot_to_onedot should be followed by a single number and then an endline"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}
        m_algorithm_type = TWODOT_TO_ONEDOT;
        m_env_add = 1;
        m_twodot_to_onedot_iter = atoi(tok[1].c_str());
      }

      else if(boost::iequals(keyword, "twodot_noise")) 
      {
	if(usedkey[TWODOT_NOISE] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[TWODOT_NOISE] = 0;
         if (tok.size() !=  3) {
           pout << "keyword twodot_noise should be followed by two single numbers and then an endline"<<endl;
           pout << "error found in the following line "<<endl;
           pout << msg<<endl;
           abort();
         }
         n_twodot_noise = 1;
         m_twodot_noise = atof(tok[1].c_str());
         m_twodot_gamma = atof(tok[2].c_str());
       }

      else if (boost::iequals(keyword,  "sweep_tol") || boost::iequals(keyword,  "sweep_tolerance"))
      {
	if(usedkey[SWEEP_TOL] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[SWEEP_TOL] = 0;
	if (tok.size() !=  2) {
	  pout << "keyword sweep_tol should be followed by a single number and then an endline"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}
	if (boost::iequals(tok[1], "loose"))
	  m_sweep_tol = 1.0e-5;
	else if (boost::iequals(tok[1], "tight"))
	  m_sweep_tol = 1.0e-7;
	else if (boost::iequals(tok[1], "verytight"))
	  m_sweep_tol = 1.0e-9;
	else
	  m_sweep_tol = atof(tok[1].c_str());
	if (m_sweep_tol < 1.0e-12) {
	  pout << "Sweep tolerance of less than 1.0e-12 may cause convergence problems. Changing it to from "<<m_sweep_tol<<" to 1.0e-12 instead."<<endl;
	  m_sweep_tol = 1.0e-12;
	}
      }


      else if (boost::iequals(keyword,  "outputlevel")) {
	if (tok.size() != 2) {
	  pout << "keyword outputlevel should be followed by a single integer and then an endline"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
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

      else if (boost::iequals(keyword,  "backward")) {
         m_backward = true;	  
         m_schedule_type_backward = true; 
         m_algorithm_type = TWODOT;
      }


      else if (boost::iequals(keyword,  "reset_iterations") || boost::iequals(keyword,  "reset_iter") || boost::iequals(keyword,  "reset_iters")) {
	m_reset_iterations = true;
      }

      else
      {
        pout << "Unrecognized option :: " << keyword << endl;
	pout << "error found in the following line "<<endl;
	pout << msg<<endl;
        abort();
      }
      msg.resize(0);
      ReadMeaningfulLine(input, msg, msgsize);
      
    }

    if (m_hf_occ_user==""){
      pout << "need to specify hf_occ. The recommended default setting is: hf_occ integral" << endl;
      abort();
    }

    if (n_elec == -1) {
      if (!m_Bogoliubov) {
        pout << "number of electrons has to be specified using the keyword nelec"<<endl;
        abort();
      } else {
        n_elec = m_norbs;
      }
    }

    if (n_spin == -1) {
      pout << "spin of the wavefunction has to be specified using the keyword spin"<<endl;
      abort();
    }

    m_alpha = (n_elec + n_spin)/2;
    m_beta = (n_elec - n_spin)/2;
    if (sym == "trans" || sym == "lzsym") 
      m_total_symmetry_number = IrrepSpace(m_total_symmetry_number.getirrep()+1); //in translational symmetry lowest irrep is 0 and not 1

    if (n_twodot_noise == 1) {
     if (m_algorithm_type == ONEDOT || fabs(m_twodot_noise*m_twodot_gamma) <= NUMERICAL_ZERO){
       pout << "twodot_noise is disabled using RANDOM noise" << endl;
       n_twodot_noise = 0;
       m_twodot_noise = 0.0;
       m_twodot_gamma = 0.0;
     }
    }
  }

#ifndef SERIAL
  boost::mpi::communicator world;
  mpi::broadcast(world, m_num_Integrals, 0);
  mpi::broadcast(world, sym, 0);
  mpi::broadcast(world, m_Bogoliubov, 0);
  mpi::broadcast(world, orbitalfile, 0);
#endif
  v_2.resize(m_num_Integrals, TwoElectronArray(TwoElectronArray::restrictedNonPermSymm));
  v_1.resize(m_num_Integrals);
  coreEnergy.resize(m_num_Integrals);
  if ( (m_calc_type==RESPONSE) && m_num_Integrals != m_baseState.size() + 1) {
    pout << "number of integrals should be 1 more than the number of base states"<<endl;
    pout << "about to exit"<<endl;
    abort();
  }

  if (m_Bogoliubov && m_num_Integrals >1 ) {
    pout << "Currently the response code does not work with non-particle number conserving hamiltonians!!";
    abort();
  }
    
//if (sym != "c1") // must be initialized even if c1 sym.
    Symmetry::InitialiseTable(sym);

  if (mpigetrank() == 0) {
    for (int l=0; l<orbitalfile.size(); l++)
      CheckFileExistence(orbitalfile[l], "Orbital file");
  }

  //read the orbitals
  for (int integral=0; integral < m_num_Integrals; integral++) {
    v_1[integral].rhf=true;
    v_2[integral].rhf=true;
    if (sym != "lzsym" && sym != "dinfh_abelian" && !NonabelianSym) {
      v_2[integral].permSymm = true;
    }
    else
      v_2[integral].permSymm = false;

    
    if (m_Bogoliubov) {
      v_2[integral].permSymm = false;
      v_cc.rhf = true;
      v_cccc.rhf = true;
      v_cccd.rhf = true;
    }
    
    // Kij-based ordering by GA opt.
#ifndef SERIAL
    mpi::broadcast(world,m_reorderType,0);
    mpi::broadcast(world,m_add_noninteracting_orbs,0);
#endif
    
    if (m_Bogoliubov) {
      v_cc.rhf=true;
      v_cccc.rhf=true;
      v_cccd.rhf=true;
      readorbitalsfile(orbitalfile[integral], v_1[integral], v_2[integral], coreEnergy[integral], v_cc, v_cccc, v_cccd);
      assert(!m_add_noninteracting_orbs);
    } 
    else {
      readorbitalsfile(orbitalfile[integral], v_1[integral], v_2[integral], coreEnergy[integral], integral);
    }
  }
  
  m_molecule_quantum = SpinQuantum(m_alpha + m_beta, SpinSpace(m_alpha - m_beta), m_total_symmetry_number);

  if (mpigetrank() == 0) {
    makeInitialHFGuess();

    pout << "Checking input for errors"<<endl;
    performSanityTest();
    generateDefaultSchedule();
    if (n_twodot_noise == 1) {
      pout << "\nScheduled random noise is disabled using twodot_noise \n" << endl;
      fill(m_sweep_noise_schedule.begin(),m_sweep_noise_schedule.end(),0.0);
    }
    //add twodot_toonedot(bla bla bla)
    pout << "Summary of input"<<endl;
    pout << "----------------"<<endl;
#ifndef MOLPRO
    writeSummary();
#else
    writeSummaryForMolpro();
#endif
    pout << endl;
  }

#ifndef SERIAL
  mpi::broadcast(world, *this,0);
  mpi::broadcast(world, v_1,0);
  mpi::broadcast(world, v_2,0);
  mpi::broadcast(world, coreEnergy,0);
  mpi::broadcast(world, v_cc,0);
  mpi::broadcast(world, v_cccc,0);
  mpi::broadcast(world, v_cccd,0);
  mpi::broadcast(world, NPROP, 0);
  mpi::broadcast(world, PROPBITLEN, 0);
#endif
}

void SpinAdapted::Input::readreorderfile(ifstream& dumpFile, std::vector<int>& oldtonew) {
  string msg; int msgsize = 5000;
  ReadMeaningfulLine(dumpFile, msg, msgsize);
  vector<string> tok;
  std::vector<int> diff;
  boost::split(tok, msg, is_any_of(", \t"), token_compress_on);
  while(msg.size() != 0) {
    for (int i=0; i<tok.size(); i++)
      if(tok[i].size() != 0) {
	if (find(oldtonew.begin(), oldtonew.end(), atoi(tok[i].c_str())-1) != oldtonew.end()) {
	  pout << "Orbital "<<atoi(tok[i].c_str())<<" appears twice in reorder file"<<endl;
	  pout << "error found at the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}
	//reorder is starting from 1 to n, but internally we store it from 0 to n-1
	oldtonew.push_back(atoi(tok[i].c_str())-1); 
	if (oldtonew.back() >m_norbs || oldtonew.back() < 0) {
	  pout << "Illegal orbital index "<<atoi(tok[i].c_str())<<" in reorder file"<<endl;
	  abort();
	}
      }
    msg.resize(0);
    ReadMeaningfulLine(dumpFile, msg, msgsize);
    boost::split(tok, msg, is_any_of(", \t"), token_compress_on);
  }

  //p_reorder contains old to new mapping
  if (oldtonew.size() != m_norbs/2) {
    pout << "Numbers or orbitals in reorder file should be equal to "<<m_norbs/2<<" instead "<<oldtonew.size()<<" found "<<endl;
    abort();
  }

}

void SpinAdapted::Input::readorbitalsfile(string& orbitalfile, OneElectronArray& v1, TwoElectronArray& v2, double& coreEnergy, int integralIndex)
{
  ifstream dumpFile; 
  dumpFile.open(orbitalfile.c_str(), ios::in);

  string msg; int msgsize = 5000;
  vector<string> tok;
  int offset = m_orbformat == DMRGFORM ? 0 : 1;
  std::ofstream ReorderFileOutput;
  std::ifstream ReorderFileInput;
  char ReorderFileName[5000];
  std::vector<int> reorder;

  if (mpigetrank() == 0) {
    ReadMeaningfulLine(dumpFile, msg, msgsize);
    boost::split(tok, msg, is_any_of("=, \t"), token_compress_on);
    
    if(offset != 0)
      m_norbs = 2*atoi(tok[2].c_str());
    else
      m_norbs = 2*atoi(tok[1].c_str());
    
    m_num_spatial_orbs = 0;
    m_spin_orbs_symmetry.resize(m_norbs);
    m_spin_to_spatial.resize(m_norbs);    

    //this is the file to which the reordering is written
    sprintf(ReorderFileName, "%s%s", save_prefix().c_str(), "/RestartReorder.dat");
  }
  boost::filesystem::path p(ReorderFileName);

#ifndef SERIAL
  boost::mpi::communicator world;
  mpi::broadcast(world,m_restart,0);
  mpi::broadcast(world,m_fullrestart,0);
#endif

  //do the reordering only if it is not a restart calculation
  //if it is then just read the reorder.dat from the scratch space
  if(get_restart() || get_fullrestart()) {
    if (mpigetrank() == 0) {
    ReorderFileInput.open(ReorderFileName);
    boost::filesystem::path ReorderFilePath(ReorderFileName);
    
    if(!boost::filesystem::exists(ReorderFilePath)) {
      pout << "==============="<<endl;
      pout << "This is a restart job and the reorder file "<<ReorderFileName<<" should be present"<<endl;
      abort();
    }
    else {
      pout << "================"<<endl;
      pout << "The Fiedler routine for finding the orbital ordering has already been run." << endl;
      pout << "Using the reorder file " << ReorderFileName << endl;
      pout << "================"<<endl;
      m_reorder.resize(m_norbs/2);
      for (int i=0; i<m_norbs/2; i++)
	ReorderFileInput >> m_reorder[i];
      ReorderFileInput.close();
    }
  }
  }
  else {
    if (mpigetrank() == 0) {
    ReorderFileOutput.open(ReorderFileName);
    }


    // read the reorder file or calculate the reordering using one of the many options  
    if (m_reorderType == FIEDLER) {
      
      if (mpigetrank() == 0){
        m_reorder=get_fiedler(orbitalfile);
      pout << "Fiedler-vector orbital ordering: ";     
      }
    }
    else if (m_reorderType == GAOPT) {

      ifstream gaconfFile;
      
      if (mpigetrank() == 0) {
      if(m_gaconffile != "default") 
         gaconfFile.open(m_gaconffile.c_str(), ios::in);
      m_reorder = get_fiedler(orbitalfile);      
      }
      //to provide as initial guess to gaopt
      m_reorder = getgaorder(gaconfFile, orbitalfile, m_reorder);      
      pout << "Genetic algorithm orbital ordering: ";
    }

    else if (m_reorderType == MANUAL) {
      if (mpigetrank() == 0) {
      ifstream reorderFile(m_reorderfile.c_str());
      CheckFileExistence(m_reorderfile, "Reorder file ");
      readreorderfile(reorderFile, m_reorder);
      pout << "Manually provided orbital ordering: ";
      }
    }
    else { //this is the no-reorder case 
      if (mpigetrank() == 0) {
      m_reorder.resize(m_norbs/2);
      for (int i=0; i<m_reorder.size(); i++)
	m_reorder.at(i) = i;
      pout << "No orbital reorder: ";
      }
    }

    if (mpigetrank() == 0) {
    for (int i=0; i<m_reorder.size(); i++){
      ReorderFileOutput << m_reorder[i]<<"  ";
  }
    ReorderFileOutput.close();    
    }
  }


  //the name reorder is confusing because it clases the with the m_reorder member of the input class and these two are inverse of each other.
  //the reorder here helps one to go from the unreordered matrices to the reordered matrices
  //and the m_reorder member of input helps one to go in the opposite direction
  //e.g. the user defined orbital order (m_reorder) could be
  // 2 3 1 4   which implies that O_reorder(1,2) -> O_unreordered(2,3)  (the former is stored internally and the latter is what is given during input)
  // but for this m_reorder the reorder vector below would be 3 1 2 4 and O_unreordered(1,2) -> O_reorder(3, 1)
  if (mpigetrank() == 0) {
  reorder.resize(m_norbs/2);
  for (int i=0; i<m_norbs/2; i++) {
    reorder.at(m_reorder[i]) = i;
    pout << m_reorder[i]+1 << " ";
  }
  pout << endl;
  pout << endl;



  int orbindex = 0;
  msg.resize(0);
  ReadMeaningfulLine(dumpFile, msg, msgsize);
  boost::split(tok, msg, is_any_of("=, \t"), token_compress_on);
  
  int readLine = 1, numRead = 1;
  bool RHF = true;
  while (!(boost::iequals(tok[0], "ISYM") || boost::iequals(tok[0], "&END"))) {
    for (int i=0; i<tok.size(); i++) {
      if (boost::iequals(tok[i], "ORBSYM") || tok[i].size() == 0) continue;

      int reorderOrbInd =  reorder.at(orbindex/2);

      int ir;
      if (atoi(tok[i].c_str()) >= 0 ) 
	ir = atoi(tok[i].c_str()) - offset;
      else if (atoi(tok[i].c_str()) < -1)
	ir = atoi(tok[i].c_str()) + offset;

      if (sym == "trans") ir += 1; //for translational symmetry the lowest irrep is 0
      if (sym == "lzsym") ir = atoi(tok[i].c_str());

      
      m_spin_orbs_symmetry[2*reorderOrbInd] = ir;
      m_spin_orbs_symmetry[2*reorderOrbInd+1] = ir;
      
      if (readLine == numRead) {
	m_num_spatial_orbs++;
	m_spatial_to_spin.push_back(orbindex);
	numRead = Symmetry::sizeofIrrep(ir);
	readLine = 0;
      }
      m_spin_to_spatial[orbindex] = m_num_spatial_orbs-1;
      m_spin_to_spatial[orbindex+1] = m_num_spatial_orbs-1;
      orbindex +=2;
      readLine++;
      

    }
    msg.resize(0);
    ReadMeaningfulLine(dumpFile, msg, msgsize);
    boost::split(tok, msg, is_any_of("=, \t"), token_compress_on);
    if(boost::iequals(tok[0], "IUHF")) RHF=false;
  }

  
  if(sym == "dinfh" ) {
    m_spatial_to_spin.clear();
    for (int i=0; i<m_spin_orbs_symmetry.size();) {
      int ir = m_spin_orbs_symmetry[i];
      m_spatial_to_spin.push_back(i);
      i += 2*Symmetry::sizeofIrrep(ir); 
    }
  }


  m_spatial_to_spin.push_back(m_norbs);
  m_spin_to_spatial.push_back(m_norbs);

  while (!((boost::iequals(tok[0], "&END")) || (boost::iequals(tok[0], "/")))) {
    int temp;
    if (boost::iequals(tok[0], "NPROP") ) {
      NPROP.push_back( atoi(tok[1].c_str()));
      NPROP.push_back( atoi(tok[2].c_str()));
      NPROP.push_back( atoi(tok[3].c_str()));
    }
    else if (boost::iequals(tok[0], "PROPBITLEN") ) {
      temp = atoi(tok[1].c_str());
    PROPBITLEN=1;
    for (int i=0; i<temp; i++)
      PROPBITLEN *= 2;
    }

    msg.resize(0);
    ReadMeaningfulLine(dumpFile, msg, msgsize);
    boost::split(tok, msg, is_any_of("=, \t"), token_compress_on);
    if(boost::iequals(tok[0], "IUHF")) RHF=false;
  }

  int AOrbOffset = 0, BOrbOffset = 0;
  if(!RHF) {
    v1.rhf = false;
    v2.rhf = false;    
  }
  v1.ReSize(m_norbs);  
  v2.ReSize(m_norbs);


  msg.resize(0);
  ReadMeaningfulLine(dumpFile, msg, msgsize); //this if the first line with integrals

  int i, j, k, l;
  double value;
  while(msg.size() != 0) {
    boost::split(tok, msg, is_any_of(" \t"), token_compress_on);
    if (tok.size() != 5) {
      pout << "error in reading orbital file"<<endl;
      pout << "error encountered at line "<<endl;
      pout << msg<<endl;
      abort();
    }
    value = atof(tok[0].c_str());
    i = atoi(tok[1].c_str())-offset;j = atoi(tok[2].c_str())-offset;k = atoi(tok[3].c_str())-offset;l = atoi(tok[4].c_str())-offset;

    if (i==-1 && j==-1 && k==-1 && l==-1) {
      coreEnergy += value;
      if (AOrbOffset == 0 && BOrbOffset == 0) //AA
	{AOrbOffset = 1; BOrbOffset = 1;} //got to BB}
      else if (AOrbOffset == 1 && BOrbOffset == 1) //BB
	{AOrbOffset = 1; BOrbOffset = 0;} //got to AB}
      else if (AOrbOffset == 1 && BOrbOffset == 0) //AB
	{AOrbOffset = 0; BOrbOffset = 0;} //got to AA}
    }
	
    else if (k==-1 && l==-1) { 
      v1(2*reorder.at(i)+AOrbOffset,2*reorder.at(j)+AOrbOffset) = value;  v1(2*reorder.at(j)+AOrbOffset,2*reorder.at(i)+AOrbOffset) = value;
    } 
    else {
      int I = reorder.at(i);int J = reorder.at(j);int K = reorder.at(k);int L = reorder.at(l);
      v2(2*I+BOrbOffset,2*K+AOrbOffset,2*J+BOrbOffset,2*L+AOrbOffset) = value;
    }
    msg.resize(0);
    ReadMeaningfulLine(dumpFile, msg, msgsize); //this if the first line with integrals
  }

  if (m_norbs/2 >= m_integral_disk_storage_thresh) //
  {
    for (int i=0; i<m_spatial_to_spin.size()-1; i++) {
      std::vector<int> orb;
      for (int j=m_spatial_to_spin[i]; j<m_spatial_to_spin[i+1]; j+=2) {
	orb.push_back(j/2);
      }
      PartialTwoElectronArray vpart(orb);
      vpart.populate(v2);
      vpart.Save(m_save_prefix, integralIndex);
    }
    v2.ReSize(0);
  }

  dumpFile.close();  
  }
}

void SpinAdapted::Input::readorbitalsfile(string& orbitalfile, OneElectronArray& v1, TwoElectronArray& v2, double& coreEnergy, PairArray& vcc, CCCCArray& vcccc, CCCDArray& vcccd) {

  ifstream dumpFile;
  dumpFile.open(orbitalfile.c_str(), ios::in);
  
  string msg; int msgsize = 5000;
  vector<string> tok;
  std::ofstream ReorderFileOutput;
  std::ifstream ReorderFileInput;
  char ReorderFileName[5000];
  std::vector<int> reorder;

  int offset = m_orbformat == DMRGFORM ? 0 : 1;

  if (mpigetrank() == 0) {
    ReadMeaningfulLine(dumpFile, msg, msgsize);
    boost::split(tok, msg, is_any_of("=, \t"), token_compress_on);

    if(offset != 0)
      m_norbs = 2*atoi(tok[2].c_str());
    else
      m_norbs = 2*atoi(tok[1].c_str());
  
    m_num_spatial_orbs = 0;
    m_spin_orbs_symmetry.resize(m_norbs);
    m_spin_to_spatial.resize(m_norbs);

    sprintf(ReorderFileName, "%s%s", save_prefix().c_str(), "/RestartReorder.dat");
  }
  boost::filesystem::path p(ReorderFileName);

#ifndef SERIAL
  boost::mpi::communicator world;
  mpi::broadcast(world,m_restart,0);
  mpi::broadcast(world,m_fullrestart,0);
#endif

  //do the reordering only if it is not a restart calculation
  //if it is then just read the reorder.dat from the scratch space
  if(get_restart() || get_fullrestart()) {
    if (mpigetrank() == 0) {
    ReorderFileInput.open(ReorderFileName);
    if(!boost::filesystem::exists(p)) {
      pout << "==============="<<endl;
      pout << "This is a restart job and the reorder file "<<ReorderFileName<<" should be present"<<endl;
      abort();
    } else {
      pout << "================"<<endl;
      pout << "The Fiedler routine for finding the orbital ordering has already been run." << endl;
      pout << "Using the reorder file " << ReorderFileName << endl;
      pout << "----------------"<<endl;
      m_reorder.resize(m_norbs/2);
      for (int i=0; i<m_norbs/2; i++)
	ReorderFileInput >> m_reorder[i];
      ReorderFileInput.close();
    }
    }
  } 
  else {
    if (mpigetrank() == 0) {
      ReorderFileOutput.open(ReorderFileName);
    }
    
    // read the reorder file or calculate the reordering using one of the many options  
    if (m_reorderType == FIEDLER) {
      if (mpigetrank() == 0) {
        m_reorder=get_fiedler_bcs(orbitalfile);
        pout << "Fiedler-vector orbital ordering: ";
      }
    } else if (m_reorderType == GAOPT) {
      if (mpigetrank() == 0) {
        pout << "GAOPT orbital ordering for BCS calculation not implemented";
      }
      abort();      
    } else if (m_reorderType == MANUAL) {
      if (mpigetrank() == 0) {
        ifstream reorderFile(m_reorderfile.c_str());
        CheckFileExistence(m_reorderfile, "Reorder file ");
        readreorderfile(reorderFile, m_reorder);
        pout << "Manually provided orbital ordering: ";
      }
    } else { //this is the no-reorder case
      if (mpigetrank() == 0) {
        m_reorder.resize(m_norbs/2);
        for (int i=0; i<m_reorder.size(); i++)
	      m_reorder.at(i) = i;
        pout << "No orbital reorder: ";
      }
    }
    if (mpigetrank() == 0) {
      // write the reorder file
      for (int i=0; i<m_reorder.size(); i++)
        ReorderFileOutput << m_reorder[i]<<"  ";
      ReorderFileOutput.close();
    }
  }

  //the name reorder is confusing because it clases the with the m_reorder member of the input class and these two are inverse of each other.
  //the reorder here helps one to go from the unreordered matrices to the reordered matrices
  //and the m_reorder member of input helps one to go in the opposite direction
  //e.g. the user defined orbital order (m_reorder) could be
  // 2 3 1 4   which implies that O_reorder(1,2) -> O_unreordered(2,3)  (the former is stored internally and the latter is what is given during input)
  // but for this m_reorder the reorder vector below would be 3 1 2 4 and O_unreordered(1,2) -> O_reorder(3, 1)
  if (mpigetrank() == 0) {
  reorder.resize(m_norbs/2);
  for (int i=0; i<m_norbs/2; i++) {
    reorder.at(m_reorder[i]) = i;
    pout << m_reorder[i]+1 << " ";
  }
  pout << endl;
  pout << endl;
  
  // orbital symmetry
  int orbindex = 0;
  msg.resize(0);
  ReadMeaningfulLine(dumpFile, msg, msgsize);
  boost::split(tok, msg, is_any_of("=, \t"), token_compress_on);
  
  int readLine = 1, numRead = 1;
  bool RHF = true;
  while (!(boost::iequals(tok[0], "ISYM") || boost::iequals(tok[0], "&END"))) {
    for (int i=0; i<tok.size(); i++) {
      if (boost::iequals(tok[i], "ORBSYM") || tok[i].size() == 0) continue;

      int reorderOrbInd =  reorder.at(orbindex/2);
      int ir;
      if (atoi(tok[i].c_str()) >= 0 ) 
	    ir = atoi(tok[i].c_str()) - offset;
      else if (atoi(tok[i].c_str()) < -1)
	    ir = atoi(tok[i].c_str()) + offset;

      if (sym == "trans") ir += 1; //for translational symmetry the lowest irrep is 0
      if (sym == "lzsym") ir = atoi(tok[i].c_str());
      
      m_spin_orbs_symmetry[2*reorderOrbInd] = ir;
      m_spin_orbs_symmetry[2*reorderOrbInd+1] = ir;
      
      if (readLine == numRead) {
    	m_num_spatial_orbs++;
    	m_spatial_to_spin.push_back(orbindex);
    	numRead = Symmetry::sizeofIrrep(ir);
    	readLine = 0;
      }
      m_spin_to_spatial[orbindex] = m_num_spatial_orbs-1;
      m_spin_to_spatial[orbindex+1] = m_num_spatial_orbs-1;
      orbindex +=2;
      readLine++;
    }
    msg.resize(0);
    ReadMeaningfulLine(dumpFile, msg, msgsize);
    boost::split(tok, msg, is_any_of("=, \t"), token_compress_on);
    if(boost::iequals(tok[0], "IUHF"))
      RHF=false;
  }

  if(sym == "dinfh" ) {
    m_spatial_to_spin.clear();
    for (int i=0; i<m_spin_orbs_symmetry.size(); i+=2) {
      int ir = m_spin_orbs_symmetry[i];
      if (ir < -1 || ir == 0 || ir == 1) 
	    m_spatial_to_spin.push_back(i);
    }
  }

  m_spatial_to_spin.push_back(m_norbs);
  m_spin_to_spatial.push_back(m_norbs);

  while (!((boost::iequals(tok[0], "&END")) || (boost::iequals(tok[0], "/")))) {
    int temp;
    if (boost::iequals(tok[0], "NPROP") ) {
      NPROP.push_back( atoi(tok[1].c_str()));
      NPROP.push_back( atoi(tok[2].c_str()));
      NPROP.push_back( atoi(tok[3].c_str()));
    } else if (boost::iequals(tok[0], "PROPBITLEN") ) {
      temp = atoi(tok[1].c_str());
      PROPBITLEN=1;
      for (int i=0; i<temp; i++)
        PROPBITLEN *= 2;
    }
    msg.resize(0);
    ReadMeaningfulLine(dumpFile, msg, msgsize);
    boost::split(tok, msg, is_any_of("=, \t"), token_compress_on);
    if(boost::iequals(tok[0], "IUHF")) RHF=false;
  }

  int section = 0;
  if(!RHF) {
    v1.rhf = false;
    v2.rhf = false; 
    vcc.rhf=false;
    vcccc.rhf=false;
    vcccd.rhf=false;
  }
  v2.ReSize(m_norbs);
  v1.ReSize(m_norbs);
  vcc.ReSize(m_norbs);
  vcccc.ReSize(m_norbs);
  vcccd.ReSize(m_norbs);
  
  msg.resize(0);
  ReadMeaningfulLine(dumpFile, msg, msgsize); //this if the first line with integrals
  
  int i, j, k, l;
  double value;
  while(msg.size() != 0) {
    boost::split(tok, msg, is_any_of(" \t"), token_compress_on);
    if (tok.size() != 5) {
      pout << "error in reading orbital file"<<endl;
      pout << "error encountered at line "<<endl;
      pout << msg<<endl;
      abort();
    }
    value = atof(tok[0].c_str());
    i = atoi(tok[1].c_str())-offset;j = atoi(tok[2].c_str())-offset;k = atoi(tok[3].c_str())-offset;l = atoi(tok[4].c_str())-offset;
    if (i==-1 && j==-1 && k==-1 && l==-1) {
      coreEnergy += value;
      section += 1;
    } else if (RHF) {
      if (section == 0) { // ccdd
        v2(2*reorder.at(i), 2*reorder.at(k), 2*reorder.at(j),2*reorder.at(l)) = value;
      } else if (section == 1) { // cccd
        vcccd.set(2*reorder.at(i), 2*reorder.at(j), 2*reorder.at(k)+1, 2*reorder.at(l), value);
      } else if (section == 2) { // cccc
        vcccc.set(2*reorder.at(i), 2*reorder.at(j), 2*reorder.at(l)+1, 2*reorder.at(k)+1, value);
      } else if (section == 3) { // cd
        if (!(k==-1 && l==-1)) {
          pout << "Orbital file error" << endl;
          abort();
        }
        v1(2*reorder.at(i), 2*reorder.at(j)) = value;
      } else if (section == 4) { // cc
        if (!(k==-1 && l==-1)) {
          pout << "Orbital file error" << endl;
          abort();
        }
        vcc(2*reorder.at(i), 2*reorder.at(j)+1) = value;
      } else {
        pout << "read orbital file error" << endl;
        abort();
      }
    } else {
      if (section == 0) { // ccdd_aa
        v2(2*reorder.at(i), 2*reorder.at(k), 2*reorder.at(j),2*reorder.at(l)) = value;
      } else if (section == 1) { // ccdd_bb
        v2(2*reorder.at(i)+1, 2*reorder.at(k)+1, 2*reorder.at(j)+1,2*reorder.at(l)+1) = value;
      } else if (section == 2) { // ccdd_ab
        v2(2*reorder.at(i), 2*reorder.at(k)+1, 2*reorder.at(j),2*reorder.at(l)+1) = value;
      } else if (section == 3) { // cccd_a
        vcccd.set(2*reorder.at(i), 2*reorder.at(j), 2*reorder.at(k)+1, 2*reorder.at(l), value);
      } else if (section == 4) { // cccd_b
        vcccd.set(2*reorder.at(i)+1, 2*reorder.at(j)+1, 2*reorder.at(k), 2*reorder.at(l)+1, value);
      } else if (section == 5) { // cccc  w_{ijkl}C_ia C_ja C_kb C_lb
        vcccc.set(2*reorder.at(i), 2*reorder.at(j), 2*reorder.at(l)+1, 2*reorder.at(k)+1, value);
      } else if (section == 6) { // cd_a
        if (!(k==-1 && l==-1)) {
          pout << "Orbital file error" << endl;
          abort();
        }
        v1(2*reorder.at(i), 2*reorder.at(j)) = value;
      } else if (section == 7) { // cd_b
        if (!(k==-1 && l==-1)) {
          pout << "Orbital file error" << endl;
          abort();
        }
        v1(2*reorder.at(i)+1, 2*reorder.at(j)+1) = value;
      } else if (section == 8) { // cc
        if (!(k==-1 && l==-1)) {
          pout << "Orbital file error" << endl;
          abort();
        }
        vcc(2*reorder.at(i), 2*reorder.at(j)+1) = value;
      } else {
        pout << "read orbital file error" << endl;
        abort();
      }
    }
    msg.resize(0);
    ReadMeaningfulLine(dumpFile, msg, msgsize); //this if the first line with integrals
  }
  if (m_norbs/2 >= m_integral_disk_storage_thresh) //
  {
    for (int i=0; i<m_spatial_to_spin.size()-1; i++) {
      std::vector<int> orb;
      for (int j=m_spatial_to_spin[i]; j<m_spatial_to_spin[i+1]; j+=2) {
	orb.push_back(j/2);
      }
      PartialTwoElectronArray vpart(orb);
      vpart.populate(v2);
      vpart.Save(m_save_prefix);
    }
    v2.ReSize(0);
  }
  dumpFile.close();
  }
}

std::vector<int> SpinAdapted::Input::get_fiedler(string& dumpname){
     Matrix fiedler; 
     ifstream dumpFile;
     dumpFile.open(dumpname.c_str(), ios::in);
     genetic::ReadIntegral(dumpFile, fiedler);
     dumpFile.close();
     SymmetricMatrix fiedler_sym;
     fiedler_sym << fiedler;
     std::vector<int> findices = fiedler_reorder(fiedler_sym);
     return findices;
}

std::vector<int> SpinAdapted::Input::get_fiedler_bcs(string& dumpname) {
  Matrix fiedler; 
  ifstream dumpFile;
  dumpFile.open(dumpname.c_str(), ios::in);
  genetic::ReadIntegral_BCS(dumpFile, fiedler);
  dumpFile.close();
  SymmetricMatrix fiedler_sym;
  fiedler_sym << fiedler;
  std::vector<int> findices = fiedler_reorder(fiedler_sym);
  return findices;
}

std::vector<int> SpinAdapted::Input::getgaorder(ifstream& gaconfFile, string& orbitalfile, std::vector<int>& fiedlerorder)
{
  ifstream dumpFile; dumpFile.open(orbitalfile.c_str());
   return genetic::gaordering(gaconfFile, dumpFile, fiedlerorder).Gen().Sequence();
}

#ifdef MOLPRO
void SpinAdapted::Input::writeSummaryForMolpro()
{
   pout << setw(50) << "Total number of orbitals : "  ;
   pout << m_norbs/2 << endl;
   pout << setw(50) << "Symmetry of targeted wavefunctions : " ;
   pout << m_alpha + m_beta << ":" << m_alpha-m_beta << ":" << m_total_symmetry_number.getirrep()+1 << endl;
   pout << setw(50) << "Number of wavefunctions targeted : " ;
   pout << m_nroots << endl;
   if (m_nroots >1) {
      pout << setw(50) << "The weights of the wavefunctions : ";
  for (int i=0; i<m_nroots; i++) 
     pout << setprecision(2) << m_weights[i];
  pout << endl;
   }
   pout << setw(50) << "Symmetry of the molecule : " ;
   pout << sym << endl;
   if (sym != "c1") {
     pout << setw(50) << "Irreducible representation of the orbitals : " ;
     for (int i=0; i<m_spin_orbs_symmetry.size(); i+=2) 
        pout << Symmetry::stringOfIrrep(m_spin_orbs_symmetry[i])<<"  "; 
     pout << endl;
   }


  pout << endl << "Schedule" << endl;
  pout << "========" << endl;
 // Need to add proper spacing here, with setw( n);
  pout << setw(10) << "Iter" ;
  pout << setw(20) <<  "# States" ;
  pout << setw(20) <<  "Davidson_tol" ;
  pout << setw(20) << "Random_noise" << endl;
  for (int i=0; i<m_sweep_iter_schedule.size(); i++) {
     pout << setw(10) << m_sweep_iter_schedule[i]; 
     pout << setw(20) << m_sweep_state_schedule[i];
     pout << setw(20) << setprecision(4) << m_sweep_tol_schedule[i] ;
     pout << setw(20) << scientific << setprecision(4) << m_sweep_noise_schedule[i] << endl;
  }
  if (m_algorithm_type == TWODOT_TO_ONEDOT) 
     pout << setw(50) << "Switching from twodot to onedot algorithm : " << m_twodot_to_onedot_iter << endl << endl;
  pout << setw(50) << "Maximum sweep iterations : " << m_maxiter << endl << endl;
}
#endif

void SpinAdapted::Input::writeSummary()
{
#ifndef SERIAL
  if (mpigetrank() == 0) {
#endif
//printf("%-50s :   %-i\n", "Total number of orbitals", m_norbs/2);
//if (m_Bogoliubov)
//  printf("%-50s :   even:%-i:%-i\n", "Symmetry of the targeted wavefunctions", m_alpha - m_beta, m_total_symmetry_number.getirrep()+1);    
//else
//  printf("%-50s :   %-i:%-i:%-i\n", "Symmetry of the targeted wavefunctions",m_alpha + m_beta, m_alpha - m_beta, m_total_symmetry_number.getirrep()+1);
//printf("%-50s :   %-i\n", "Number of wavefunctions targeted", m_nroots);
//if (m_nroots >1) {
//  printf("%-50s :   ", "The weights of the wavefunctions");
//  for (int i=0; i<m_nroots; i++) 
//    printf("%-10.2e", m_weights[i]);
//  printf("\n");
//}
//printf("%-50s :   %s\n", "Symmetry of the molecule", sym.c_str());
//if (sym != "c1") {
//  printf("%-50s :   ", "Irreducible representations of the orbitals");
//  for (int i=0; i<m_spin_orbs_symmetry.size(); i+=2) 
//    pout << Symmetry::stringOfIrrep(m_spin_orbs_symmetry[i])<<"  ";
//  printf("\n");
//}


//  printf("\nSchedule\n");
//  printf("--------\n");
//  printf("%-10s : %-20s  %-20s  %-20s\n", "Iter", "# States", "Davidson_tol",  "Random_noise");
//  for (int i=0; i<m_sweep_iter_schedule.size(); i++)
//    printf("%-10i : %-20i  %-20.4e  %-20.4e\n", m_sweep_iter_schedule[i], m_sweep_state_schedule[i], m_sweep_tol_schedule[i], m_sweep_noise_schedule[i]);
//  if (m_algorithm_type == TWODOT_TO_ONEDOT) 
//    printf("%-50s :   %-i\n", "Switching from twodot to onedot algorithm", m_twodot_to_onedot_iter);
//  
//  printf("%-50s :   %-i\n", "Maximum sweep iterations", m_maxiter);

  // removed printf dependencies
  pout << setw(50) << left << "Total number of orbitals"
       << " :   " << left << m_norbs/2 << endl;
  if (m_Bogoliubov)
  {
     pout << setw(50) << left << "Symmetry of the targeted wavefunctions"
          << " :   even:" << m_alpha-m_beta << ":" << m_total_symmetry_number.getirrep()+1 << endl;
  }
  else
  {
     pout << setw(50) << left << "Symmetry of the targeted wavefunctions"
          << " :   " << m_alpha+m_beta << ":" << m_alpha-m_beta << ":" << m_total_symmetry_number.getirrep()+1 << endl;
  }
  pout << setw(50) << left << "Number of wavefunctions targeted"
       << " :   " << m_nroots << endl;
  if (m_nroots > 1)
  {
     pout << setw(50) << left << "The weights of the wavefunctions" << " :   ";
     for (int i = 0; i < m_nroots; ++i) pout << left << setprecision(2) << setw(10) << scientific << m_weights[i];
     pout << endl;
  }
  pout << setw(50) << left << "Symmetry of the molecule" << " :   "
       << sym.c_str() << endl;
  if (sym != "c1")
  {
     pout << setw(50) << left << "Irreducible representations of the orbitals" << " :   ";
     for (int i = 0; i < m_spin_orbs_symmetry.size(); i+=2)
     {
        pout << Symmetry::stringOfIrrep(m_spin_orbs_symmetry[i]) << "  ";
     }
     pout << endl;
  }

  pout << endl << "Schedule" << endl;
  pout << "--------" << endl;
  pout << setw(10) << left << "Iter" << " : "
       << setw(20) << left << "# States" << "  "
       << setw(20) << left << "Davidson_tol" << "  "
       << setw(20) << left << "Random_noise" << endl;

  for (int i = 0; i < m_sweep_iter_schedule.size(); ++i)
  {
     pout << setw(10) << left << m_sweep_iter_schedule[i] << " : "
          << setw(20) << left << m_sweep_state_schedule[i] << "  "
          << setw(20) << left << scientific << m_sweep_tol_schedule[i] << "  "
          << setw(20) << left << scientific << m_sweep_noise_schedule[i] << endl;
  }
  if (m_algorithm_type == TWODOT_TO_ONEDOT) 
  {
     pout << setw(50) << left << "Switching from twodot to onedot algorithm" << " :   "
          << m_twodot_to_onedot_iter << endl;
  }
  pout << setw(50) << left << "Maximum sweep iterations" << " :   " << m_maxiter << endl;

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
    pout << "total number of orbitals has to be a positive number"<<endl;
    abort();
  }
  if (m_norbs/2 < 4 && m_calc_type == DMRG) {
    pout << "DMRG cannot be run with fewer than 4 orbitals"<<endl;
    abort();
  }
  if (m_norbs/2 > 200) {
    pout << "Number of orbitals cannot be greater than 130"<<endl;
    abort();
  }
  if (m_alpha+m_beta < 0 || (m_alpha+m_beta == 0 && !m_Bogoliubov)) {
    pout << "Total number of electrons cannot be negative"<<endl;
    abort();
  }
  if (m_alpha < m_beta) {
    pout << "DMRG requires the spin to a positive number and less than the total number of electrons"<<endl;
    abort();
  }
  if (m_norbs < m_alpha+m_beta) {
    pout<< "No of spin orbitals has to be greater than total number of electrons"<<endl;
    abort();
  }
  for (int i=0; i<m_spin_orbs_symmetry.size(); i+=2) {
    Symmetry::irrepAllowed(m_spin_orbs_symmetry[i]);
  }

  if (sym == "dinfh") {
    if (m_total_symmetry_number.getirrep() < 0) {
      pout << "Wavefunction irrep cannot be less than 0"<<endl;
      abort();
    }
  }
  else
    Symmetry::irrepAllowed(m_total_symmetry_number.getirrep());

  if (m_calc_type == RESPONSE && m_occupied_orbitals == -1) {
    pout << "For response type of calculation, number of occupied orbitals must be specified"<<endl;
    abort();
  }

  //this is important so the user cannot break the code
  if (m_schedule_type_default && !m_schedule_type_backward) {
    if (m_maxM == 0) {
      pout << "With default schedule a non-zero maxM has to be specified"<<endl;
      pout << "Current m_maxM = "<<m_maxM<<endl;
      abort();
    }
    if (m_maxM <= 0) {
      pout << "maxM cannot be less than 0"<<endl;
      abort();
    }
    if (m_maxM > 10000) {
      pout << "default schedule only works for maxM less than 10000"<<endl;
      abort();
    }
    if (m_maxM <= m_nroots) {
      pout << "maxM cannot be less than 0"<<endl;
      abort();
    }
    if (m_maxM < m_startM) {
       pout << "maxM is smaller than startM" << endl;
       pout << "Make sure you specify a maxM larger than " << m_startM << endl;
       pout << "or specify a startM smaller than " << m_maxM << endl;
       abort();
    }
  }
  else if (m_schedule_type_backward) {
    int lastM=0; //To be changed to a variable
    lastM = m_lastM;
    if (lastM==0) lastM=50;
    if (lastM > m_startM) {
       pout << "lastM is larger than startM" << endl;
       pout << "Make sure you specify a lastM larger than " << lastM << endl;
       pout << "or specify a larger startM " << endl;
       abort();
    }
  }
  else {
  // This needs to go at the end of sanity check, or move up schedule 
    if (m_maxM != 0) {
      pout << "With detailed schedule a non-zero maxM should not be specified"<<endl;
      abort();
    }
    if (m_sweep_iter_schedule.size() == 0) {
      pout << "Zero lines of schedule specified"<<endl;
      abort();
    }
  }

  //Still part of the schedule. Might need to move to Schedule.C
  //if(m_algorithm_type == TWODOT_TO_ONEDOT && m_twodot_to_onedot_iter == 0)
  //  m_twodot_to_onedot_iter = min(m_sweep_iter_schedule.back()+2, m_maxiter-1);


  if (m_algorithm_type == TWODOT_TO_ONEDOT && m_twodot_to_onedot_iter >= m_maxiter) {
    pout << "Switch from twodot to onedot algorithm cannot happen after maxiter"<<endl;
    pout << m_twodot_to_onedot_iter <<" < "<<m_maxiter<<endl;
    abort();
  }
  if (m_algorithm_type != ONEDOT && m_stateSpecific == true) {
    pout << "Only onedot algorithm is allowed with state specific calculation."<<endl;
    abort();
  }

#ifndef SERIAL
  }
#endif
  if (m_norbs <= 6)
    m_calc_type = TINYCALC;

}

void SpinAdapted::Input::makeInitialHFGuess() {

  std::vector<int> hf_occupancy_tmp(m_norbs,0);

  if (m_Bogoliubov) {
    // overwrite hf_occ_user option, since initial guess is always vacuum in BCS case
    m_hf_occupancy.assign(m_norbs, 0);
  }
  else if (m_hf_occ_user == "manual") {
    //check if n_orbs is correct and if n_elec is correct
    if (m_hf_occupancy.size() != m_norbs/2 ) {
      pout << "ERROR: The length of user-defined HF occupancies does not match the number of orbitals " << endl;
      pout << "Length of occupancies is: " << m_hf_occupancy.size() << ", and the number of orbitals is: " << m_norbs/2 << endl;
     abort();
    }
    int UserElectrons = 0;
    for (int i=0; i<m_hf_occupancy.size(); i++) {
      UserElectrons += m_hf_occupancy[i];
      if (m_hf_occupancy[i] == 2) { hf_occupancy_tmp[2*i] = 1; hf_occupancy_tmp[2*i+1] = 1;}
      else if (m_hf_occupancy[i] == 1) { hf_occupancy_tmp[2*i] = 1; hf_occupancy_tmp[2*i+1] = 0;}
      else if (m_hf_occupancy[i] == 0) { hf_occupancy_tmp[2*i] = 0; hf_occupancy_tmp[2*i+1] = 0;}
      else {pout << "the HF occupancy of orbital "<<i<<" should either be 2, 1, or 0. Instead "<<m_hf_occupancy[i]<<" is given."<<endl; abort();}
    }	

    if (UserElectrons != m_alpha+m_beta ) {
      pout << "ERROR: The total number of electrons in HF occupancy does not match the number given initially " << endl;
      pout << "#electron in HF occupancy: " << UserElectrons << ", and the number of electrons given initially: " << m_alpha+m_beta << endl;
     abort();
    }
    
  }
  else if (m_hf_occ_user == "canonical") {
    for (int i=0; i<m_alpha; i++) {
      hf_occupancy_tmp[2*i] = 1;
      if (i<m_beta) hf_occupancy_tmp[2*i+1] = 1;
    }
  }  
  else if (m_hf_occ_user == "integral") { 

    //arrange t(i,i) in a multimap and it will rearrange the ts in ascending order with i
//  multimap<double, int> ele_map;
//  for( int i = 0; i < m_norbs/2; ++i ){
//    ele_map.insert( pair<double, int>( v_1[0](2*i, 2*i), i ) );
//  }

//  multimap<double, int> :: iterator it_alpha = ele_map.begin();
//  for( int i = 0; i < m_alpha; ++i ){
//    int ia = it_alpha->second;
//    hf_occupancy_tmp.at( 2*ia ) = 1;
//    if (i < m_beta)
//    hf_occupancy_tmp.at( 2*ia+1 ) = 1;
//    ++it_alpha;
//  }
    hf_occupancy_tmp = this->hfOccGenerator_();
  }
  else {
    pout << "currently other options besides manual, integral  and canonical are not implemented."<<endl;
    abort();
  }


  //now reorder the hf_occupancy, 
  std::vector<int> reorder(m_norbs/2);
  for (int i=0; i<m_norbs/2; i++) {
    reorder.at(m_reorder[i]) = i;
  }

  m_hf_occupancy.resize(m_norbs,0);
  for (int i=0; i<m_norbs/2; i++) {
    m_hf_occupancy[2*reorder[i]] = hf_occupancy_tmp[2*i];
    m_hf_occupancy[2*reorder[i]+1] = hf_occupancy_tmp[2*i+1];
  }

  pout << "Initial HF occupancy guess: ";//  << endl;
  for( int i = 0; i < m_hf_occupancy.size(); ++i ){
   pout << m_hf_occupancy.at(i) << " " ;
  }
  pout << endl;


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

  if (setStateSpecific())
    return 1;

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

  if (setStateSpecific())
    return std::vector<double>(1,1.0);
  
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
    if (sI.quanta[i].get_n() == n && sI.quanta[i].get_s() == SpinSpace(s))
      return i;
  return 0;
}
