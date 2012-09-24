/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


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
#ifdef MOLPRO
#include "global/CxOutputStream.h"
#define pout if (mpigetrank() == 0) xout
#endif

using namespace std;

namespace SpinAdapted {
string sym;
}
void CheckFileExistance(string filename, string filetype);

void SpinAdapted::Input::ReadMeaningfulLine(ifstream& input, string& msg, int msgsize)
{
  bool readmore = true;
  while (readmore && !input.eof()) {
    char msgctr[msgsize];
    input.getline(msgctr, msgsize+1);

    msg=string(msgctr);
    if(msg.size() == msgsize) {
      pout << "in the process of reading line begining with "<<endl<<msg<<endl;
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
  m_twodot_to_onedot_iter = 0;
  m_integral_disk_storage_thresh = 100; //this is usually 100

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
  m_maxiter = 10;
  m_screen_tol = 1.00e-20;

  m_load_prefix = ".";
  m_save_prefix = ".";

  m_maxj = 15;
  m_ninej.init(m_maxj);
  m_set_Sz = false;

  m_sweep_tol = 1.0e-5;
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
  m_schedule_type_default = false;
  m_maxM = 0;
  m_core_energy = 0.0;
  m_reorder = false;
  m_reorderfile = "";

  m_orbformat=MOLPROFORM;
}

void SpinAdapted::Input::usedkey_error(string& key, string& line) {
  pout << "keyword "<<key<<" already used once"<<endl;
  pout << "error found at line "<<endl;
  pout << line<<endl;
  abort();
}

SpinAdapted::Input::Input(const string& config_name)
{
  //first collect all the data
  std::vector<int> usedkey(NUMKEYWORDS, -1);
  int n_elec = -1;
  int n_spin = -1;
  sym = "c1";
  string orbitalfile;
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

      if (boost::iequals(keyword,  "orbs") || boost::iequals(keyword,  "orbitals"))
      {
	if(usedkey[ORBS] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[ORBS] = 0;
	if (tok.size() != 2) {
	  pout << "keyword orbs should be followed by a single filename and then an end line"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}
	orbitalfile = tok[1];
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
      else if (boost::iequals(keyword,  "reorder"))
      {
	if(usedkey[REORDER] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[REORDER] = 0;
	if (tok.size() != 2) {
	  pout << "keyword reorder should be followed by the filename and then an end line"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
	m_reorder = true;
	m_reorderfile = tok[1];
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
	    pout << "Each line of the schedule contain four entries sweep_iteration   #retained states   davidson tolerance     noise"<<endl;
	    pout << "error found at the following line "<<endl;
	    pout<< msg<<endl;
	    abort();
	  }
	  
	  m_sweep_iter_schedule.push_back( atoi(schd_tok[0].c_str()));
	  m_sweep_state_schedule.push_back( atoi(schd_tok[1].c_str()));
	  m_sweep_qstate_schedule.push_back( 0);  //DEPRECATED OPTION
	  m_sweep_tol_schedule.push_back( atof(schd_tok[2].c_str()));
	  m_sweep_noise_schedule.push_back( atof(schd_tok[3].c_str()));
	  m_sweep_additional_noise_schedule.push_back( 0.0);  //DEPRECATED OPTION
	  
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
	if (tok.size() !=  2) {
	  pout << "keyword irrep should be followed by a single number and then an end line"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}	
	m_total_symmetry_number = IrrepSpace(atoi(tok[1].c_str())-1);
      }
      else if (boost::iequals(keyword,  "hubbard"))
	m_ham_type = HUBBARD;
      else if (boost::iequals(keyword,  "dmrg"))
	m_calc_type = DMRG;
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
      else if (boost::iequals(keyword,  "restart_onepdm") || boost::iequals(keyword,  "restart_onerdm") || boost::iequals(keyword,  "restart_ordm"))
	m_calc_type = RESTART_ONEPDM;
      else if (boost::iequals(keyword,  "restart_twopdm") || boost::iequals(keyword,  "restart_twordm") || boost::iequals(keyword,  "restart_trdm"))
	m_calc_type = RESTART_TWOPDM;
      else if(boost::iequals(keyword,  "prefix") || boost::iequals(keyword,  "scratch"))
      {
	if(usedkey[PREFIX] == 0) 
	  usedkey_error(keyword, msg);
	usedkey[PREFIX] = 0;
	m_load_prefix = tok[1];
	m_save_prefix = m_load_prefix;
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

	ReadMeaningfulLine(input, msg, msgsize);
	vector<string> weighttoken;
	boost::split(weighttoken, msg, is_any_of(" \t"), token_compress_on);
	
	if (!boost::iequals(weighttoken[0],  "weights")) 
	  continue;
	PROVIDED_WEIGHTS = true;

        m_weights.resize(m_nroots);

	if (weighttoken.size() != m_nroots +1 ) {
	  pout << "keyword weights should be followed by floating point numbers providing weights for "<<m_nroots<<" states."<<endl;
	  pout << "You could chose to omit the keyworkd weights in which case the weights will be distributed uniformly between the different roots"<<endl;
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


      else if (boost::iequals(keyword,  "docd") || boost::iequals(keyword,  "do_cd"))
      {
        m_do_cd = true;
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
	  pout << "keyword twodot_to_onedot should be followed by a single number and then an endline"<<endl;
	  pout << "error found in the following line "<<endl;
	  pout << msg<<endl;
	  abort();
	}
        m_algorithm_type = TWODOT_TO_ONEDOT;
        m_env_add = 1;
        m_twodot_to_onedot_iter = atoi(tok[1].c_str());
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
	if (m_sweep_tol <= 1.0e-10) {
	  pout << "Sweep tolerance of less than equal to 1.0e-10 may cause convergence problems. Changing it to from "<<m_sweep_tol<<" to 1.0e-9 instead."<<endl;
	  m_sweep_tol = 1.0e-9;
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

    if (n_elec == -1) {
      pout << "number of electrons has to be specified using the keyword nelec"<<endl;
      abort();
    }
    if (n_spin == -1) {
      pout << "spin of the wavefunction has to be specified using the keyword spin"<<endl;
      abort();
    }
    m_alpha = (n_elec + n_spin)/2;
    m_beta = (n_elec - n_spin)/2;
    if (sym == "trans") 
      m_total_symmetry_number = IrrepSpace(m_total_symmetry_number.getirrep()+1); //in translational symmetry lowest irrep is 0 and not 1
    m_molecule_quantum = SpinQuantum(m_alpha + m_beta, m_alpha - m_beta, m_total_symmetry_number);


  }


#ifndef SERIAL
  boost::mpi::communicator world;
  mpi::broadcast(world, sym, 0);
#endif
    
  if (sym != "c1")
    Symmetry::InitialiseTable(sym);

  ifstream orbitalFile;
  orbitalFile.open(orbitalfile.c_str(), ios::in);

  //read the orbitals
  v_1.rhf= true; 
  v_2.rhf=true;
  if (sym != "dinfh")
    v_2.permSymm = true;

  if (mpigetrank() == 0) {
    CheckFileExistance(orbitalfile, "Orbital file ");
    readorbitalsfile(orbitalFile, v_1, v_2);
    
//#ifndef MOLPRO
    pout << "Checking input for errors"<<endl;
    performSanityTest();
    pout << "Summary of input"<<endl;
    pout << "----------------"<<endl;
    writeSummary();
    pout << endl;
/*#else
    xout << "Checking input for errors"<<endl;
    performSanityTest();
    xout << "Summary of input"<<endl;
    xout << "----------------"<<endl;
    writeSummaryForMolpro();
    pout << endl;
#endif*/
  }
#ifndef SERIAL
  mpi::broadcast(world,*this,0);
  mpi::broadcast(world,v_1,0);
  mpi::broadcast(world,v_2,0);
  mpi::broadcast(world, NPROP, 0);
  mpi::broadcast(world, PROPBITLEN, 0);
#endif

}

void SpinAdapted::Input::readreorderfile(ifstream& dumpFile, std::vector<int>& preorder, std::vector<int>& oldtonew)
{
  string msg; int msgsize = 5000;
  ReadMeaningfulLine(dumpFile, msg, msgsize);
  vector<string> tok;
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
	oldtonew.push_back(atoi(tok[i].c_str())-1); //reorder is starting from 1 to n, but internally we store it from 0 to n
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
    pout << "Numbers or orbitals in reorder file should be equal to "<<m_norbs/2<<" instead "<<preorder.size()<<" found "<<endl;
    abort();
  }

  preorder.resize(m_norbs/2);
  for (int i=0; i<m_norbs/2; i++) {
    if (oldtonew[i] >= m_norbs/2 || oldtonew[i] < 0) {
      pout << "Orbital indices in reorder file should be in range 0 to "<<m_norbs/2<<endl;
      abort();
    }
    preorder.at(oldtonew[i]) = i;
  }
}

void SpinAdapted::Input::readorbitalsfile(ifstream& dumpFile, OneElectronArray& v1, TwoElectronArray& v2)
{
  string msg; int msgsize = 5000;
  ReadMeaningfulLine(dumpFile, msg, msgsize);
  vector<string> tok;
  boost::split(tok, msg, is_any_of("=, \t"), token_compress_on);

  int offset = m_orbformat == DMRGFORM ? 0 : 1;
  if(offset != 0)
    m_norbs = 2*atoi(tok[2].c_str());
  else
    m_norbs = 2*atoi(tok[1].c_str());

  m_num_spatial_orbs = 0;
  m_spin_orbs_symmetry.resize(m_norbs);
  m_spin_to_spatial.resize(m_norbs);

  v2.ReSize(m_norbs);
  v1.ReSize(m_norbs);

  std::vector<int> reorder, oldtonew;
  // read the reorder file
  if (m_reorder) {
    ifstream reorderFile(m_reorderfile.c_str());
    CheckFileExistance(m_reorderfile, "Reorder file ");
    readreorderfile(reorderFile, reorder, oldtonew);
  }

  int orbindex = 0;
  msg.resize(0);
  ReadMeaningfulLine(dumpFile, msg, msgsize);
  boost::split(tok, msg, is_any_of("=, \t"), token_compress_on);
  
  while (!(boost::iequals(tok[0], "ISYM") || boost::iequals(tok[0], "&END"))) {
    for (int i=0; i<tok.size(); i++) {
      if (boost::iequals(tok[i], "ORBSYM") || tok[i].size() == 0) continue;

      int reorderOrbInd = m_reorder ? reorder.at(orbindex/2) : orbindex/2;

      int ir;
      if (atoi(tok[i].c_str()) >= 0 ) 
	ir = atoi(tok[i].c_str()) - offset;
      else if (atoi(tok[i].c_str()) < -1)
	ir = atoi(tok[i].c_str()) + offset;

      if (sym == "trans") ir += 1; //for translational symmetry the lowest irrep is 0

      m_spin_orbs_symmetry[2*reorderOrbInd] = ir;
      m_spin_orbs_symmetry[2*reorderOrbInd+1] = ir;

      if (sym == "dinfh") {
	if (ir < -1 || ir == 0 || ir == 1) {
	  m_num_spatial_orbs ++;
	  m_spatial_to_spin.push_back(orbindex);
	}
	m_spin_to_spatial[orbindex] = m_num_spatial_orbs-1;
	m_spin_to_spatial[orbindex+1] = m_num_spatial_orbs-1;
      }
      else {
	m_num_spatial_orbs++;
	m_spatial_to_spin.push_back(orbindex);
	
	m_spin_to_spatial[orbindex] = m_num_spatial_orbs-1;
	m_spin_to_spatial[orbindex+1] = m_num_spatial_orbs-1;
      }
      orbindex += 2;
    }
    msg.resize(0);
    ReadMeaningfulLine(dumpFile, msg, msgsize);
    boost::split(tok, msg, is_any_of("=, \t"), token_compress_on);
  }

  
  if(sym == "dinfh" && m_reorder) {
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
    }
    else if (boost::iequals(tok[0], "PROPBITLEN") ) {
      temp = atoi(tok[1].c_str());
    }
    PROPBITLEN=1;
    for (int i=0; i<temp; i++)
      PROPBITLEN *= 2;

    msg.resize(0);
    ReadMeaningfulLine(dumpFile, msg, msgsize);
    boost::split(tok, msg, is_any_of("=, \t"), token_compress_on);
  }

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

    if (i==-1 && j==-1 && k==-1 && l==-1) m_core_energy = value;
    else if (k==-1 && l==-1) { 
      if(m_reorder){ v1(2*reorder.at(i),2*reorder.at(j)) = value;  v1(2*reorder.at(j),2*reorder.at(i)) = value;}
      else {v1(2*i,2*j) = value;v1(2*j,2*i) = value;}
    } 
    else {
      int I=i, J=j, K=k, L=l;
      if(m_reorder) {I = reorder.at(i);J = reorder.at(j);K = reorder.at(k);L = reorder.at(l);}
      v2(2*I,2*K,2*J,2*L) = value;
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
  
}

#ifdef MOLPRO
void SpinAdapted::Input::writeSummaryForMolpro()
{
#ifndef SERIAL
  if (mpigetrank() == 0) {
#endif
     xout << setw(50) << "Total number of orbitals : "  ;
     xout << m_norbs/2 << endl;
     xout << setw(50) << "Symmetry of tragetted wavefunctions : " ;
     xout << m_alpha + m_beta << ":" << m_alpha-m_beta << ":" << m_total_symmetry_number.getirrep()+1 << endl;
     xout << setw(50) << "Number of wavefunctions targetted : " ;
     xout << m_nroots << endl;
     if (m_nroots >1) {
        xout << setw(50) << "The weights of the wavefunctions : ";
    for (int i=0; i<m_nroots; i++) 
       xout << scientific << setprecision(2) << m_weights[i];
    xout << endl;
     }
     xout << setw(50) << "Symmetry of the molecule : " ;
     xout << sym << endl;
     if (sym != "c1") {
       xout << setw(50) << "Irreducible representation of the orbitals : " ;
       for (int i=0; i<m_spin_orbs_symmetry.size(); i+=2) 
          xout << Symmetry::stringOfIrrep(m_spin_orbs_symmetry[i])<<"  "; 
       xout << endl;
     }

  if (m_calc_type == DMRG) {
    xout << endl << "Schedule" << endl;
    xout << "--------" << endl;
   // Need to add proper spacing here, with setw( n);
    xout << setw(10) << "Iter" ;
    xout << setw(20) <<  "# States" ;
    xout << setw(20) <<  "Davidson_tol" ;
    xout << setw(20) << "Random_noise" << endl;
    for (int i=0; i<m_sweep_iter_schedule.size(); i++) {
       xout << setw(10) << m_sweep_iter_schedule[i]; 
       xout << setw(20) << m_sweep_state_schedule[i];
       xout << setw(20) << scientific << setprecision(4) << m_sweep_tol_schedule[i] ;
       xout << setw(20) << scientific << setprecision(4) << m_sweep_noise_schedule[i] << endl;
    }
    if (m_algorithm_type == TWODOT_TO_ONEDOT) 
       xout << setw(50) << "Switching from twodot to onedot algorithm : " << m_twodot_to_onedot_iter << endl << endl;
    xout << setw(50) << "Maximum sweep iterations : " << m_maxiter << endl << endl;
  }
#ifndef SERIAL
  }
#endif
}
#endif

void SpinAdapted::Input::writeSummary()
{
#ifndef SERIAL
  if (mpigetrank() == 0) {
#endif
  printf("%-50s :   %-i\n", "Total number of orbitals", m_norbs/2);
  printf("%-50s :   %-i:%-i:%-i\n", "Symmetry of the targetted wavefunctions",m_alpha + m_beta, m_alpha - m_beta, m_total_symmetry_number.getirrep()+1);
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

  if (m_calc_type == DMRG) {
    printf("\nSchedule\n");
    printf("--------\n");
    printf("%-10s : %-20s  %-20s  %-20s\n", "Iter", "# States", "Davidson_tol",  "Random_noise");
    for (int i=0; i<m_sweep_iter_schedule.size(); i++)
      printf("%-10i : %-20i  %-20.4e  %-20.4e\n", m_sweep_iter_schedule[i], m_sweep_state_schedule[i], m_sweep_tol_schedule[i], m_sweep_noise_schedule[i]);
    if (m_algorithm_type == TWODOT_TO_ONEDOT) 
      printf("%-50s :   %-i\n", "Switching from twodot to onedot algorithm", m_twodot_to_onedot_iter);
    
    printf("%-50s :   %-i\n", "Maximum sweep iterations", m_maxiter);
  }
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
  if (m_norbs/2 < 4) {
    pout << "DMRG cannot be run with fewer than 4 orbitals"<<endl;
    abort();
  }
  if (m_norbs/2 > 200) {
    pout << "Number of orbitals cannot be greater than 130"<<endl;
    abort();
  }
  if (m_alpha+m_beta <= 0) {
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

  //this is important so the user cannot break the code
  if (m_schedule_type_default) {
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

    if (m_sweep_tol <= 0.0) {
      pout << "Using the default tolerance sweep tolerance of 1.0e-5."<<endl;
      m_sweep_tol = 1.0e-5;
    }
    int nentry = 0;
    m_sweep_iter_schedule.resize(nentry);
    m_sweep_state_schedule.resize(nentry);
    m_sweep_qstate_schedule.resize(nentry,0);
    m_sweep_tol_schedule.resize(nentry);
    m_sweep_noise_schedule.resize(nentry);
    m_sweep_additional_noise_schedule.resize(nentry,0);
    
    double sweeptol = m_sweep_tol;
    if (m_maxM >= 50) {
      m_sweep_iter_schedule.push_back(0); m_sweep_state_schedule.push_back(50); m_sweep_tol_schedule.push_back(1.0e-5);  m_sweep_noise_schedule.push_back(1.0e-4);
    }
    else {
      m_sweep_iter_schedule.push_back(0); m_sweep_state_schedule.push_back(m_maxM); m_sweep_tol_schedule.push_back(1.0e-6);  m_sweep_noise_schedule.push_back(1.0e-4);
    }

    if (m_maxM >50) {
    if (m_maxM >= 100) {
      m_sweep_iter_schedule.push_back(4); m_sweep_state_schedule.push_back(100); m_sweep_tol_schedule.push_back(1.0e-5);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    else {
      m_sweep_iter_schedule.push_back(4); m_sweep_state_schedule.push_back(m_maxM); m_sweep_tol_schedule.push_back(1.0e-5);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    }

    if (m_maxM >100) {
    if (m_maxM >= 250) {
      m_sweep_iter_schedule.push_back(8); m_sweep_state_schedule.push_back(250); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    else {
      m_sweep_iter_schedule.push_back(8); m_sweep_state_schedule.push_back(m_maxM); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    }

    if (m_maxM >250) {
    if (m_maxM >= 500) {
      m_sweep_iter_schedule.push_back(12); m_sweep_state_schedule.push_back(500); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    else {
      m_sweep_iter_schedule.push_back(12); m_sweep_state_schedule.push_back(m_maxM); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    }

    if (m_maxM >500) {
    if (m_maxM >= 1000) {
      m_sweep_iter_schedule.push_back(16); m_sweep_state_schedule.push_back(1000); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    else {
      m_sweep_iter_schedule.push_back(16); m_sweep_state_schedule.push_back(m_maxM); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    }

    if (m_maxM >1000) {
    if (m_maxM >= 2000) {
      m_sweep_iter_schedule.push_back(19); m_sweep_state_schedule.push_back(2000); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    else {
      m_sweep_iter_schedule.push_back(19); m_sweep_state_schedule.push_back(m_maxM); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    }

    if (m_maxM >2000) {
    if (m_maxM >= 4000) {
      m_sweep_iter_schedule.push_back(22); m_sweep_state_schedule.push_back(4000); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    else {
      m_sweep_iter_schedule.push_back(22); m_sweep_state_schedule.push_back(m_maxM); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    }

    if (m_maxM >4000) {
    if (m_maxM >= 6000) {
      m_sweep_iter_schedule.push_back(25); m_sweep_state_schedule.push_back(6000); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    else {
      m_sweep_iter_schedule.push_back(25); m_sweep_state_schedule.push_back(m_maxM); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    }

    if (m_maxM >6000) {
    if (m_maxM >= 8000) {
      m_sweep_iter_schedule.push_back(28); m_sweep_state_schedule.push_back(8000); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    else {
      m_sweep_iter_schedule.push_back(28); m_sweep_state_schedule.push_back(m_maxM); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    }

    if (m_maxM >8000) {
    if (m_maxM >= 10000) {
      m_sweep_iter_schedule.push_back(31); m_sweep_state_schedule.push_back(10000); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    else {
      m_sweep_iter_schedule.push_back(31); m_sweep_state_schedule.push_back(m_maxM); m_sweep_tol_schedule.push_back(5.0e-6);  m_sweep_noise_schedule.push_back(5.0e-5);
    }
    }    

    int lastiter = m_sweep_iter_schedule.back();
    m_sweep_iter_schedule.push_back(lastiter+2); m_sweep_state_schedule.push_back(m_maxM); m_sweep_tol_schedule.push_back(sweeptol/10.0);  m_sweep_noise_schedule.push_back(0.0e-5);


    if (m_twodot_to_onedot_iter < 18 && m_algorithm_type == TWODOT_TO_ONEDOT) {
      if (m_twodot_to_onedot_iter <= 0)
	pout << "Sweep at which the switch from twodot to onedot will happen -> "<<lastiter+4<<endl;
      m_twodot_to_onedot_iter = lastiter+4;
    }
    if (m_maxiter <= m_sweep_iter_schedule.back()) {
      pout << "With the default schedule and maxM specified, maxiter has to be atleast "<<lastiter+6<<endl;
      pout << "changing maxiter to "<<lastiter+6<<endl;
      m_maxiter = lastiter+6;
    }
  }
  else {
    if (m_maxM != 0) {
      pout << "With detailed schedule a non-zero maxM should not be specified"<<endl;
      abort();
    }
    if (m_sweep_iter_schedule.size() == 0) {
      pout << "Zero lines of schedule specified"<<endl;
      abort();
    }
  }

  if(m_algorithm_type == TWODOT_TO_ONEDOT && m_twodot_to_onedot_iter == 0)
    m_twodot_to_onedot_iter = min(m_sweep_iter_schedule.back()+2, m_maxiter-1);

  if (m_algorithm_type == TWODOT_TO_ONEDOT && m_twodot_to_onedot_iter >= m_maxiter) {
    pout << "Switch from twodot to onedot algorithm cannot happen after maxiter"<<endl;
    pout << m_twodot_to_onedot_iter <<" < "<<m_maxiter<<endl;
    abort();
  }


  if (m_maxiter < m_sweep_iter_schedule.back()) {
    pout << "maximum iterations allowed is less than the last sweep iteration in your schedule."<<endl;
    pout << m_maxiter <<" < "<< (m_sweep_iter_schedule.back())<<endl;
    pout << "either increase the max_iter or reduce the number of sweeps"<<endl;
    abort();
  }

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
