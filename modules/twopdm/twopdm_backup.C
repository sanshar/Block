#include "twopdm.h"
#include "MatrixBLAS.h"
#include <boost/format.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>


void compute_twopdm_sweep(const std::vector<Wavefunction>& wavefunctions, const SpinBlock& system, const SpinBlock& systemDot, const SpinBlock& newSystem, const SpinBlock& newEnvironment, const SpinBlock& big, const int numprocs)
{
  const int nroots = wavefunctions.size();
  array_4d<double> twopdm(big.size()*2, big.size()*2, big.size()*2, big.size()*2);
  cout <<"twopdm first dim "<< twopdm.dim1()<<endl;
  for(int i=0;i<nroots;++i)
    for(int j=0;j<=i;++j)
      {
	const Wavefunction &wavefunction1 = wavefunctions[i];
	const Wavefunction &wavefunction2 = wavefunctions[j];
	load_twopdm_binary(twopdm, i ,j);

	pout << "superblock_compute_onedot_two_pdm_sweep " << i << " " << j << endl;
	//const std::vector<int> distribute_work = distribute_procs(numprocs,4);

	//if(mpigetrank() == distribute_work[0])
	compute_two_pdm_2_2(wavefunction1, wavefunction2, big, twopdm);

	accumulate_twopdm(twopdm);
	save_twopdm_binary(twopdm, i ,j);
      }
}

void compute_two_pdm_2_2(const Wavefunction& wave1, const Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  // (1) CCDD 
  // cA cA
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();


  for (int ij = 0; ij < leftBlock->get_op_array(CRE_CRE_S2).get_size(); ++ij)
  {
    boost::shared_ptr<SparseMatrix> cicj2 = leftBlock->get_op_array(CRE_CRE_S2).get_local_element(ij).getworkingrepresentation(leftBlock);
    int ix = cicj2->get_orbs(0);
    int jx = cicj2->get_orbs(1);
    boost::shared_ptr<SparseMatrix> cicj0 = leftBlock->get_op_rep(CRE_CRE_S0, ix, jx);
    
    for (int kl =0; kl <rightBlock->get_op_array(CRE_CRE_S2).get_size(); kl++)
    {
      boost::shared_ptr<SparseMatrix> ckcl2 = rightBlock->get_op_array(CRE_CRE_S2).get_local_element(kl).getworkingrepresentation(rightBlock);
      int kx = ckcl2->get_orbs(0);
      int lx = ckcl2->get_orbs(1);
      boost::shared_ptr<SparseMatrix> ckcl0 = rightBlock->get_op_rep(CRE_CRE_S0, kx, lx);

      get_expectations_cccc(wave1, wave2, big, *cicj2, *cicj0, *ckcl2, *ckcl0, twopdm);

    }
  }      


  for (int ij = 0; ij < leftBlock->get_op_array(CRE_DES_S2).get_size(); ++ij)
  {
    boost::shared_ptr<SparseMatrix> cidj2 = leftBlock->get_op_array(CRE_DES_S2).get_local_element(ij).getworkingrepresentation(leftBlock);
    int ix = cidj2->get_orbs(0);
    int jx = cidj2->get_orbs(1);
    boost::shared_ptr<SparseMatrix> cidj0 = leftBlock->get_op_rep(CRE_DES_S0, ix, jx);
    
    for (int kl =0; kl <rightBlock->get_op_array(CRE_DES_S2).get_size(); kl++)
    {
      boost::shared_ptr<SparseMatrix> ckdl2 = rightBlock->get_op_array(CRE_DES_S2).get_local_element(kl).getworkingrepresentation(rightBlock);
      int kx = ckdl2->get_orbs(0);
      int lx = ckdl2->get_orbs(1);
      boost::shared_ptr<SparseMatrix> ckdl0 = rightBlock->get_op_rep(CRE_DES_S0, kx, lx);
      
      /*cout <<"wave1  " <<wave1<<endl;
      cout << "cidj2 "<<*cidj2<<endl;
      cout << "cidj0 "<<*cidj0<<endl;
      cout << "ckdl2 "<<*ckdl2<<endl;
      cout << "ckdl0 "<<*ckdl0<<endl;*/
      get_expectations_cdcd(wave1, wave2, big, *cidj2, *cidj0, *ckdl2, *ckdl0, twopdm);

    }
  }      

}

void get_expectations_cdcd(const Wavefunction& wave1, const Wavefunction& wave2, const SpinBlock& big, const SparseMatrix& cidj2, const SparseMatrix& cidj0, const SparseMatrix& ckdl2, const SparseMatrix& ckdl0, array_4d<double>& twopdm)
{
  int ix = 2*cidj2.get_orbs(0);
  int jx = 2*cidj2.get_orbs(1);
  int kx = 2*ckdl2.get_orbs(0);
  int lx = 2*ckdl2.get_orbs(1);


  double baab=0, abba=0, c2c2=0;
  double c2c0=0.0, c0c2=0.0, c0c0 =0.0;
  vector<Wavefunction> braij0(1), brakl0(1);
  braij0[0].initialise( (wave1.get_deltaQuantum()-cidj0.get_deltaQuantum())[0], &big, true);
  brakl0[0].initialise( (wave2.get_deltaQuantum()-ckdl0.get_deltaQuantum())[0], &big, true);
  vector<SpinQuantum> quanta = wave1.get_deltaQuantum()-cidj2.get_deltaQuantum();

  vector<Wavefunction> braij2(quanta.size()), brakl2(quanta.size());
  for (int i=0; i<quanta.size(); i++)
  {
    braij2[i].initialise(quanta[i], &big, true);
    brakl2[i].initialise(quanta[i], &big, true);
  }

  Opxwave(cidj2, wave1, braij2, big, true);
  Opxwave(cidj0, wave1, braij0, big, true);
  Opxwave(ckdl2, wave2, brakl2, big, false);
  Opxwave(ckdl0, wave2, brakl0, big, false);

  vector<int> Sz3(3); vector<int> Sz0(1);
  Sz0[0] = dmrginp.Sz();
  Sz3[0] = dmrginp.Sz()+2; Sz3[1] = dmrginp.Sz(); Sz3[2] = dmrginp.Sz()-2;
  vector<double> cidjcldk(3, 0.0);
  for (int j =0; j<quanta.size(); j++) {
    vector<double> temp = DotProduct(braij2[j], brakl2[j], Sz3, big);
    for (int i=0; i<3; i++) cidjcldk[i] += temp[i];
  }

  baab = cidjcldk[0];
  abba = cidjcldk[2];
  c2c2 = cidjcldk[1];

  assign_antisymmetric(twopdm, ix+1, lx, jx, kx+1, -baab); 
  assign_antisymmetric(twopdm, kx+1, jx, lx, ix+1, -baab); 
  assign_antisymmetric(twopdm, ix, lx+1, jx+1, kx, -abba); 
  assign_antisymmetric(twopdm, kx, jx+1, lx+1, ix, -abba); 


  vector<double> temp = DotProduct(braij0[0], brakl0[0], Sz0, big);
  c0c0 = temp[0];
  if (quanta.size() == 3)
  {
    cidjcldk = DotProduct(braij2[1], brakl0[0], Sz0, big); 
    c2c0 = cidjcldk[0];
    cidjcldk = DotProduct(braij0[0], brakl2[1], Sz0, big); 
    c0c2 = cidjcldk[0];
  }

  double aaaa = 1.0/2.0*(c0c0+c0c2+c2c0+c2c2);
  double bbaa = 1.0/2.0*(c0c0+c0c2-c2c0-c2c2);
  double aabb = 1.0/2.0*(c0c0-c0c2+c2c0-c2c2);
  double bbbb = 1.0/2.0*(c0c0-c0c2-c2c0+c2c2);

  assign_antisymmetric(twopdm, ix, lx, jx, kx, -aaaa); 
  assign_antisymmetric(twopdm, ix+1, lx, jx+1, kx, -bbaa); 
  assign_antisymmetric(twopdm, ix, lx+1, jx, kx+1, -aabb); 
  assign_antisymmetric(twopdm, ix+1, lx+1, jx+1, kx+1, -bbbb); 

  assign_antisymmetric(twopdm, kx, jx, lx, ix, -aaaa); 
  assign_antisymmetric(twopdm, kx, jx+1, lx, ix+1, -bbaa); 
  assign_antisymmetric(twopdm, kx+1, jx, lx+1, ix, -aabb); 
  assign_antisymmetric(twopdm, kx+1, jx+1, lx+1, ix+1, -bbbb); 

}




vector<double> DotProduct(const Wavefunction& w1, const Wavefunction& w2, vector<int>& Sz, const SpinBlock& big)
{
  int leftOpSz = big.get_leftBlock()->get_stateInfo().quanta.size ();
  int rightOpSz = big.get_rightBlock()->get_stateInfo().quanta.size ();
  const StateInfo* rS = big.get_stateInfo().rightStateInfo, *lS = big.get_stateInfo().leftStateInfo;

  vector<double> output(Sz.size(), 0.0);
  for (int lQ =0; lQ < leftOpSz; lQ++)
    for (int rQ = 0; rQ < rightOpSz; rQ++) {
      int lSp = lS->quanta[lQ].get_s(), rSp = rS->quanta[rQ].get_s();
      if (w1.allowed(lQ, rQ) && w2.allowed(lQ, rQ))
      {
	vector<double> coeff1(Sz.size(), 0.0);
	vector<double> coeff2(Sz.size(), 0.0);
	for (int lZ = -lSp; lZ <= lSp ; lZ=lZ+2)
	  for (int rZ = -rSp; rZ <= rSp ; rZ=rZ+2)
	    for (int sindex = 0; sindex < Sz.size(); sindex++)
	    {
	      if (lZ+rZ != Sz[sindex]) continue;
	      coeff1[sindex] += pow(cg(lSp, rSp, w1.get_spin(), lZ, rZ, Sz[sindex]),2);
	    }
	double b1b2 = MatrixDotProduct(w1(lQ, rQ), w2(lQ, rQ));
	for (int sindex = 0; sindex < Sz.size(); sindex++) {
	  output[sindex] += b1b2*coeff1[sindex];
	}
      }	
    }

  return output;
}


void writewave(const Wavefunction& w1, int Sz, const SpinBlock& big)
{
  int leftOpSz = big.get_leftBlock()->get_stateInfo().quanta.size ();
  int rightOpSz = big.get_rightBlock()->get_stateInfo().quanta.size ();
  const StateInfo* rS = big.get_stateInfo().rightStateInfo, *lS = big.get_stateInfo().leftStateInfo;

  for (int lQ =0; lQ < leftOpSz; lQ++)
    for (int rQ = 0; rQ < rightOpSz; rQ++) {
      int lSp = lS->quanta[lQ].get_s(), rSp = rS->quanta[rQ].get_s();
      if (w1.allowed(lQ, rQ))
      {
	for (int lZ = -lSp; lZ <= lSp ; lZ=lZ+2)
	  for (int rZ = -rSp; rZ <= rSp ; rZ=rZ+2)
	  {
	    if (lZ+rZ != Sz) continue;
	    double coeff = cg(lSp, rSp, w1.get_spin(), lZ, rZ, Sz);
	    Matrix m = w1(lQ, rQ);
	    MatrixScale(coeff, m);
	    cout << lS->quanta[lQ]<<" "<<lZ<<"      "<<rS->quanta[rQ]<<" "<<rZ<<"  coeff "<< coeff<<endl;
	    cout << m<<endl;
	  }
      }	
    }

}



void Opxwave(const SparseMatrix& op, const Wavefunction& wave, vector<Wavefunction>& waveout, const SpinBlock& big, bool left)
{
  int leftOpSz = big.get_leftBlock()->get_stateInfo().quanta.size ();
  int rightOpSz = big.get_rightBlock()->get_stateInfo().quanta.size ();
  const StateInfo* rS = big.get_stateInfo().rightStateInfo, *lS = big.get_stateInfo().leftStateInfo;

  ObjectMatrix<Matrix> bra;
  if (left)
    bra.ReSize(leftOpSz*leftOpSz, rightOpSz); 
  else
    bra.ReSize(leftOpSz, rightOpSz*rightOpSz);

  if (left)
  {
    for (int lQ =0; lQ < leftOpSz; lQ++)
      for (int rQ = 0; rQ < rightOpSz; rQ++)
	for (int lQPrime =0; lQPrime < leftOpSz; lQPrime++) 
	{
	  if (wave.allowed(lQPrime, rQ) && op.allowed(lQPrime, lQ))
	  {
	    double scale = op.get_scaling(lS->quanta[lQPrime].get_s(), lS->quanta[lQ].get_s());
	    int lindex = lQPrime*leftOpSz + lQ;
	    bra(lindex, rQ).ReSize(lS->getquantastates(lQ), rS->getquantastates(rQ));
	    bra(lindex, rQ) = 0.0;
	    MatrixMultiply (op.operator_element(lQPrime, lQ), TransposeOf(op.conjugacy()), wave(lQPrime, rQ), 'n',
			    bra(lindex, rQ), scale);
	  }
	}

    for (int waveindex = 0; waveindex<waveout.size(); waveindex++)
    {
      for (int lQ =0; lQ < leftOpSz; lQ++)
	for (int rQ = 0; rQ < rightOpSz; rQ++) {
	  if (waveout[waveindex].allowed(lQ,rQ)) {
	    int lSp = lS->quanta[lQ].get_s(), rSp = rS->quanta[rQ].get_s();
	    waveout[waveindex](lQ, rQ) = 0.0;
	    for (int lQPrime=0; lQPrime<leftOpSz; lQPrime++) {
	      if ( wave.allowed(lQPrime, rQ) && op.allowed(lQPrime, lQ))
	      {
		int lindex = lQPrime*leftOpSz + lQ;
		/*double scale = dmrginp.get_ninej()(lSp, rSp, waveout[waveindex].get_spin(),
						   op.get_spin(), 0 , op.get_spin(),
						   lS->quanta[lQPrime].get_s(), rSp, wave.get_spin());*/
		double scale = dmrginp.get_ninej()(lS->quanta[lQPrime].get_s(), rSp, wave.get_spin(),
						   op.get_spin(), 0 , op.get_spin(),
						   lS->quanta[lQ].get_s(), rSp, waveout[waveindex].get_spin());
		if (scale != 0.0)
		  MatrixScaleAdd(scale, bra(lindex, rQ), waveout[waveindex](lQ, rQ));
	      }
	    }
	    
	  }
	}
    }
  }
  else 
  {
    for (int lQ =0; lQ < leftOpSz; lQ++)
      for (int rQ = 0; rQ < rightOpSz; rQ++)
	for (int rQPrime =0; rQPrime < rightOpSz; rQPrime++)
	{
	  if (wave.allowed(lQ, rQPrime) && op.allowed(rQPrime, rQ))
	  {
	    double scale = op.get_scaling(rS->quanta[rQPrime].get_s(), rS->quanta[rQ].get_s());
	    int rindex = rQPrime*rightOpSz + rQ;
	    bra(lQ, rindex).ReSize(lS->getquantastates(lQ), rS->getquantastates(rQ));
	    bra(lQ, rindex) = 0.0;
	    MatrixMultiply (wave(lQ, rQPrime), 'n', op.operator_element(rQPrime, rQ), op.conjugacy(),
			    bra(lQ, rindex), scale);
	  }
	}

    for (int waveindex=0; waveindex<waveout.size(); waveindex++)
      for (int lQ =0; lQ < leftOpSz; lQ++)
	for (int rQ = 0; rQ < rightOpSz; rQ++) {
	  if (waveout[waveindex].allowed(lQ, rQ)) {
	    int lSp = lS->quanta[lQ].get_s(), rSp = rS->quanta[rQ].get_s();

	    for (int rQPrime = 0; rQPrime < rightOpSz; rQPrime++) { 
	      if (wave.allowed(lQ, rQPrime) && op.allowed(rQPrime, rQ))
	      {
		int rindex = rQPrime*rightOpSz+rQ;
		/*double scale = dmrginp.get_ninej()(lSp, rSp, waveout[waveindex].get_spin(),
						   0, op.get_spin(), op.get_spin(),
						   lSp, rS->quanta[rQPrime].get_s(), wave.get_spin());*/
		double scale = dmrginp.get_ninej()(lSp, rS->quanta[rQPrime].get_s(), wave.get_spin(),
						   0, op.get_spin(), op.get_spin(),
						   lSp, rS->quanta[rQ].get_s(), waveout[waveindex].get_spin());
		if (scale != 0.0)
		  MatrixScaleAdd(scale, bra(lQ, rindex), waveout[waveindex](lQ, rQ));
	      }
	    }
	  }
	}
  }


}


void get_expectations_cccc(const Wavefunction& wave1, const Wavefunction& wave2, const SpinBlock& big, const SparseMatrix& cicj2, const SparseMatrix& cicj0, const SparseMatrix& ckcl2, const SparseMatrix& ckcl0, array_4d<double>& twopdm)
{
  int ix = 2*cicj2.get_orbs(0);
  int jx = 2*cicj2.get_orbs(1);
  int kx = 2*ckcl2.get_orbs(0);
  int lx = 2*ckcl2.get_orbs(1);


  double bbbb=0, aaaa=0, bbaa =0, aabb =0, c2c2=0;
  double c2c0=0.0, c0c2=0.0, c0c0 =0.0;
  vector<Wavefunction> braij0(1), brakl0(1);
  braij0[0].initialise( (wave1.get_deltaQuantum()-cicj0.get_deltaQuantum())[0], &big, true);
  brakl0[0].initialise( (wave2.get_deltaQuantum()-ckcl0.get_deltaQuantum())[0], &big, true);
  vector<SpinQuantum> quanta = wave1.get_deltaQuantum()-cicj2.get_deltaQuantum();

  vector<Wavefunction> braij2(quanta.size()), brakl2(quanta.size());
  for (int i=0; i<quanta.size(); i++)
  {
    braij2[i].initialise(quanta[i], &big, true);
    brakl2[i].initialise(quanta[i], &big, true);
  }

  Opxwave(cicj2, wave1, braij2, big, true);
  Opxwave(cicj0, wave1, braij0, big, true);
  Opxwave(ckcl2, wave2, brakl2, big, false);
  Opxwave(ckcl0, wave2, brakl0, big, false);


  vector<int> Sz3(3); vector<int> Sz0(1);
  Sz0[0] = dmrginp.Sz();
  Sz3[0] = dmrginp.Sz()+2; Sz3[1] = dmrginp.Sz(); Sz3[2] = dmrginp.Sz()-2;
  vector<double> cicjckcl(3, 0.0);
  for (int j =0; j<quanta.size(); j++) {
    vector<double> temp = DotProduct(braij2[j], brakl2[j], Sz3, big);
    for (int i=0; i<3; i++) cicjckcl[i] += temp[i];
  }

  bbbb = cicjckcl[0];
  aaaa = cicjckcl[2];
  c2c2 = cicjckcl[1];

  assign_antisymmetric(twopdm, ix+1, jx+1, lx+1, kx+1, bbbb); 
  assign_antisymmetric(twopdm, ix, jx, lx, kx, aaaa); 
  assign_antisymmetric(twopdm, kx+1, lx+1, jx+1, ix+1, bbbb); 
  assign_antisymmetric(twopdm, kx, lx, jx, ix, aaaa); 


  cicjckcl = DotProduct(braij0[0], brakl0[0], Sz0, big);
  c0c0 = cicjckcl[0];
  if (quanta.size() == 3)
  {
    cicjckcl = DotProduct(braij2[1], brakl0[0], Sz0, big); 
    c2c0 = cicjckcl[0];
    cicjckcl = DotProduct(braij0[0], brakl2[1], Sz0, big); 
    c0c2 = cicjckcl[0];
  }
  double abab = 1.0/2.0*(c2c2+c2c0+c0c2+c0c0);
  double baab = 1.0/2.0*(c2c2+c2c0-c0c2-c0c0);
  double abba = 1.0/2.0*(c2c2-c2c0+c0c2-c0c0);
  double baba = 1.0/2.0*(c2c2-c2c0-c0c2+c0c0);

  assign_antisymmetric(twopdm, ix, jx+1, lx, kx+1, abab); 
  assign_antisymmetric(twopdm, ix+1, jx, lx, kx+1, baab); 
  assign_antisymmetric(twopdm, ix, jx+1, lx+1, kx, abba); 
  assign_antisymmetric(twopdm, ix+1, jx, lx+1, kx, baba); 

  assign_antisymmetric(twopdm, kx+1, lx, jx+1, ix, abab); 
  assign_antisymmetric(twopdm, kx+1, lx, jx, ix+1, baab); 
  assign_antisymmetric(twopdm, kx, lx+1, jx+1, ix, abba); 
  assign_antisymmetric(twopdm, kx, lx+1, jx, ix+1, baba); 

}


void save_twopdm_text(const array_4d<double>& twopdm, const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[30];
    sprintf (file, "%s%d.%d", "twopdm.", i, j);
    ofstream ofs(file);
    ofs << twopdm.dim1() << endl;
    for(int k=0;k<twopdm.dim1();++k)
      for(int l=0;l<twopdm.dim2();++l)
        for(int m=0;m<twopdm.dim3();++m)
          for(int n=0;n<twopdm.dim4();++n)
            ofs << boost::format("%d %d %d %d %20.14e\n") % k % l % m % n % twopdm(k,l,m,n);
    ofs.close();
  }
}

void save_twopdm_binary(const array_4d<double>& twopdm, const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[30];
    sprintf (file, "%s%d.%d", "twopdm.", i, j);
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << twopdm;
    ofs.close();
  }
}

void load_twopdm_binary(array_4d<double>& twopdm, const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[30];
    sprintf (file, "%s%d.%d", "twopdm.", i, j);
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load(ifs);
    load >> twopdm;
    ifs.close();
  }
  mpi::communicator world;
  mpi::broadcast(world,twopdm,0);
  if(mpigetrank())
    twopdm.Clear();
}

void save_averaged_twopdm(const int &nroots)
{
  if(!mpigetrank())
  {
    array_4d<double> twopdm;
    char file[30];
    for(int i=0;i<nroots;++i)
    {
      sprintf (file, "%s%d.%d", "twopdm.", i, i);
      ifstream ifs(file);
      int size = 0;
      ifs >> size;
      if(i==0)
        twopdm.resize(size,size,size,size);
      int k,l,m,n;
      double val;
      while(ifs >> k)
      {
        ifs >> l >> m >> n >> val;
        twopdm(k,l,m,n) += dmrginp.weights()[i]*val;
      }
    }
    sprintf (file, "%s", "twopdm");
    ofstream ofs(file);
    ofs << twopdm.dim1() << endl;
    for(int k=0;k<twopdm.dim1();++k)
      for(int l=0;l<twopdm.dim2();++l)
        for(int m=0;m<twopdm.dim3();++m)
          for(int n=0;n<twopdm.dim4();++n)
            ofs << boost::format("%d %d %d %d %20.14e\n") % k % l % m % n % twopdm(k,l,m,n);

    ofs.close();
  }
}

void accumulate_twopdm(array_4d<double>& twopdm)
{
  array_4d<double> tmp_recv;
  mpi::communicator world;
  if (!mpigetrank()) 
  {
    for(int i=1;i<world.size();++i)
    {
      world.recv(i, i, tmp_recv);
      for(int k=0;k<twopdm.dim1();++k)
        for(int l=0;l<twopdm.dim2();++l)
          for(int m=0;m<twopdm.dim3();++m)
            for(int n=0;n<twopdm.dim4();++n)
              if(tmp_recv(k,l,m,n) != 0.)
                twopdm(k,l,m,n) = tmp_recv(k,l,m,n);
    }
  }
  else
  {
    world.send(0, mpigetrank(), twopdm);
  }
}

void assign_antisymmetric(array_4d<double>& twopdm, const int i, const int j, const int k, const int l, const double val)
{
  twopdm(i, j, k, l) = val;
  twopdm(i, j, l, k) = -val;
  twopdm(j, i, k, l) = -val;
  twopdm(j, i, l, k) = val;
}

