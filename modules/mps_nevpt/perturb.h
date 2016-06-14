#ifndef PERTURB_HEAD
#define PERTURB_HEAD

#include "SpinQuantum.h"
#include "global.h"
#include "Symmetry.h"
#include <vector>
#include "wavefunction.h"

using namespace std;


namespace SpinAdapted{
//  namespace mps_nevpt{
//    class perturb
//    {
//      private: 
//      public:
//        ActivePerturbType type;
//        perturb() =delete;
//
//        //perturb(ActivePerturbType type_, const Wavefunction& w_, const PerturbArray& v_, const ActivePerturbArray& ActiveV_, const std::vector<SpinQuantum>& delta_): type(type_),w0(w),v(v_),ActiveV(ActiveV_),delta(delta_){
//        perturb(ActivePerturbType type_, const Wavefunction& w_, const PerturbArray& v_, const ActivePerturbArray& ActiveV_, const SpinQuantum& delta_): type(type_),w0(w),v(v_),ActiveV(ActiveV_),delta(delta_){
//        ketquanta = dmrginp.effective_molecule_quantum_vec(); 
//        braquanta = delta+ketquanta;
//        }
//
//        //perturb(const Wavefunction& w_,const PerturbArray& v_, const ActivePerturbArray& ActiveV_, std::vector<SpinQuantum> delta_): w0(w),v(v_),ActiveV(ActiveV_),delta(delta_){
//        perturb(const Wavefunction& w_,const PerturbArray& v_, const ActivePerturbArray& ActiveV_, SpinQuantum delta_): w0(w),v(v_),ActiveV(ActiveV_),delta(delta_){
//        ketquanta = dmrginp.effective_molecule_quantum_vec(); 
//        braquanta = delta+ketquanta;
//        }
//
//        const Wavefunction& w0;
//        const wavenumber(){ return type+1000;}
//        const PerturbArray& v;
//        const ActivePerturbArray& ActiveV;
//        std::vector<SpinQuantum> braquanta;
//        const std::vector<SpinQuantum> ketquanta;
//        //std::vector<SpinQuantum> delta;
//        SpinQuantum delta;
//        virtual run(const Wavefunction& w_, const SpinBlock& big)= 0;
//    }
//
//    class DDsinglePerturb : public perturb
//    {
//      private:
//      public:
//        ActivePerturbType type= DDsingle;
//        perturb(const Wavefunction& w_,const PerturbArray& v_, const ActivePerturbArray& ActiveV_): w0(w),v(v_),ActiveV(ActiveV_){
//        ketquanta = dmrginp.effective_molecule_quantum_vec(); 
//        //delta =std::vector(SpinQuantum(-2,0,0));
//        delta =SpinQuantum(-2,0,0);
//        braquanta = delta+ketquanta;
//        double CompFactor;
//        void run (const Wavefunction& w, const SpinBlock& big);
//        }
//    }
//
//
//
    //Perturber means a perturb function, not a type of perturb as definied before.
    
    class perturber
    {
      public:
        TwoPerturbType type_;
        bool initialized_ = false;
        const bool& initialized() const {return initialized_; }
        const TwoPerturbType& type() const { return type_; }
        vector<int> orbs;
        perturber(){;}
        static vector<double> ZeroEnergy;
        static vector<double> CoreEnergy;

        //perturber(ActivePerturbType type_, const Wavefunction& w_, const PerturbArray& v_, const ActivePerturbArray& ActiveV_, const std::vector<SpinQuantum>& delta_): type(type_),w0(w),v(v_),ActiveV(ActiveV_),delta(delta_){
        //perturber(ActivePerturbType type_, const Wavefunction& w_, const PerturbArray& v_, const ActivePerturbArray& ActiveV_, const SpinQuantum& delta_): type(type_),w0(w),v(v_),ActiveV(ActiveV_),delta(delta_){
        //ketquanta = dmrginp.effective_molecule_quantum_vec(); 
        //braquanta = delta+ketquanta;
        //}

        //perturber(const Wavefunction& w_,const PerturbArray& v_, const ActivePerturbArray& ActiveV_, std::vector<SpinQuantum> delta_): w0(w),v(v_),ActiveV(ActiveV_),delta(delta_){
        //perturber(const Wavefunction& w_,const PerturbArray& v_, const ActivePerturbArray& ActiveV_, SpinQuantum delta_): w0(w),v(v_),ActiveV(ActiveV_),delta(delta_){
        //ketquanta = dmrginp.effective_molecule_quantum_vec(); 
        //braquanta = delta+ketquanta;
        //}
        //virtual void init(const Wavefunction& w_, const vector<int>& orbs_) =0 ;
        //virtual void init(const Wavefunction& w_, int orb1, int orb2=-1, int orb3 =-1) =0 ;
//        virtual void init(int w_, const vector<int>& orbs_) =0 ;
//        virtual void init(int w_, int orb1, int orb2=-1, int orb3 =-1) =0 ;
//        virtual void init(const vector<int>& orbs_) =0 ;
//        virtual void init(int orb1, int orb2=-1, int orb3 =-1) =0 ;
        const vector<int>& orb() const { return orbs;}
        const int& orb(int i) const { return orbs[i];}
        vector<int>& orb() { return orbs;}
        int& orb(int i) { return orbs[i];}

        //const Wavefunction& w0;
        int w0;
        //virtual const int wavenumber(){ return static_cast<int>(type_)*10000;}
        virtual int wavenumber() const =0;
        //const PerturbArray& v;
        //const ActivePerturbArray& ActiveV;
        std::vector<SpinQuantum> braquanta;
        std::vector<SpinQuantum> ketquanta;
        //std::vector<SpinQuantum> delta;
        SpinQuantum delta;
//        virtual void run(const Wavefunction& w_, const SpinBlock& big)= 0;
    };

    class VaPerturber : public perturber
    {
      public: 
        VaPerturber(){
        type_= TwoPerturbType::Va;
        initialized_ = true;
        }
        //void init( const Wavefunction& w_, const vector<int>& orbs_)
        void init( int w_, const vector<int>& orbs_)
        {
          w0=w_;
          orbs= orbs_;
          assert(orbs.size()==1);
          //delta = vector<SpinQuantum>(1,SpinQuantum(-1,1,SymmetryOfSpatialOrb(orbs[0])));
          delta = -getSpinQuantum(orbs[0]);
          ketquanta = dmrginp.effective_molecule_quantum_vec();
          braquanta = delta+ketquanta;
        }
        
        void init( const vector<int>& orbs_)
        {
          orbs= orbs_;
          assert(orbs.size()==1);
          //delta = vector<SpinQuantum>(1,-SpinQuantum(1,1,SymmetryOfSpatialOrb(orbs[0])));
          delta = -getSpinQuantum( orbs[0] );
          //delta = -getSpinQuantum( orbs[0]);
          ketquanta = dmrginp.effective_molecule_quantum_vec();
          braquanta = delta+ketquanta;
        }
        
        void init(int orb)
        {
          //init(std::vector<int>(1,dmrginp.spinAdapted()? orb: 2*orb));
          init(std::vector<int>(1,orb));
        }

        int wavenumber() const { return static_cast<int>(type_)*10000+orbs[0]+100000;}
        
    };

    class ViPerturber : public perturber
    {
      public: 
        ViPerturber(){
        type_= TwoPerturbType::Vi;
        initialized_ = true;
        }
        //void init( const Wavefunction& w_, const vector<int>& orbs_)
        void init( int w_, const vector<int>& orbs_)
        {
          w0=w_;
          orbs= orbs_;
          assert(orbs.size()==1);
          //delta = vector<SpinQuantum>(1,SpinQuantum(-1,1,SymmetryOfSpatialOrb(orbs[0])));
          delta = getSpinQuantum(orbs[0]);
          ketquanta = dmrginp.effective_molecule_quantum_vec();
          braquanta = delta+ketquanta;
        }
        
        void init( const vector<int>& orbs_)
        {
          orbs= orbs_;
          assert(orbs.size()==1);
          //delta = vector<SpinQuantum>(1,-SpinQuantum(1,1,SymmetryOfSpatialOrb(orbs[0])));
          delta = getSpinQuantum( orbs[0] );
          //delta = -getSpinQuantum( orbs[0]);
          ketquanta = dmrginp.effective_molecule_quantum_vec();
          braquanta = delta+ketquanta;
        }
        
        void init(int orb)
        {
          //init(std::vector<int>(1,dmrginp.spinAdapted()? orb: 2*orb));
          init(std::vector<int>(1,orb));
        }

        int wavenumber() const { return static_cast<int>(type_)*10000+orbs[0]+100000;}
        
    };

//  }
}




#endif
