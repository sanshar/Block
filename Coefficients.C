#include "global.h"
#include "wrapper.h"
#include "fciqmchelper.h"
#ifndef SERIAL
#include "mpi.h"
#endif
using namespace SpinAdapted;


void CheckFileExistence(string filename, string filetype);
vector<vector<bool>> ReadCoefFile() {
  CheckFileExistence("ComputeCoefficients", " Coefficient File ");
  ifstream file("ComputeCoefficients");
  string msg; int msgsize = 5000;
  Input::ReadMeaningfulLine(file, msg, msgsize);
  vector<string> schd_tok;
  boost::split(schd_tok, msg, is_any_of(" \t"), token_compress_on);
  int norbs = atoi(schd_tok[0].c_str());
  int ncoefs = atoi(schd_tok[1].c_str());
  vector<vector<bool>> coefs(ncoefs, vector<bool>(norbs, false));

  for (int i = 0; i < ncoefs; ++i) {
    msg.clear();    
    Input::ReadMeaningfulLine(file, msg, msgsize);
    schd_tok.clear();
    boost::split(schd_tok, msg, is_any_of(" \t"), token_compress_on);
    for (int j = 0; j < schd_tok.size()-1; ++j) {
      coefs[i][atoi(schd_tok[j].c_str())-1] = true;
    }
  }
  file.close();
  return coefs;
}

int main(int argc, char* argv[]) {
  int rank = 0, size = 1;
  #ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  #endif
  initBoostMPI(argc, argv);
  ReadInputFromC(argv[1], -1);
  readMPSFromDiskAndInitializeStaticVariables();  
  initializeGlobalMPS(0);  
  if (mpigetrank() == 0) {
    MPS state(0);
    auto coefs = ReadCoefFile();
    ofstream fout("Coefficients.0.0.txt");
    for (int i = 0; i < coefs.size(); ++i) {
      fout << state.get_coefficient(coefs[i]) << endl;
    }
    fout.close();
  }
#ifndef SERIAL
  MPI_Finalize();
#endif
  return 0;
}
