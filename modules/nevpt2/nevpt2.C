#include "nevpt2_operators.h"
#include "nevpt2_mpi.h"
#include "nevpt2.h"
#include "nevpt2_info.h"
#include "para_array.h"
#include "operatorfunctions.h"
#include "ripdm.h"
#include <boost/serialization/serialization.hpp>
#include "nevpt2_util.h"
#include "tensor_operator.h"
#include "davidson.h"
#include <string.h>
#include "nevpt2_pal.h"
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

#define _NEV_ZERO 1e-14

namespace SpinAdapted{
  //============================================================================
  // Evaluate the energy contribution from class V(1,i)
  //============================================================================
  double V_i(int *OrbWin, const vector<double> &EOrb, const SpinBlock &big,
             Wavefunction &WF, IntegralContainer &IKJL, 
             const array_6d &DC3, const array_4d<double> &DC2, const Matrix &D1, 
             const Matrix &D1_, array_6d &AuxA_, const Matrix &heff_, 
             const Matrix &heff, bool Conventional, bool ConventionalOverlap){
    
    char msg[512];
    int t,u,v;
    int i0,i1,t0,t1,a0,a1;
    int a,b;
    int tu;
    int i,j;
    double x,x_,y,y_,z,z_;
    double N=0.0;
    double D=0.0;
    double E=0.0;
    int t_,u_,v_;
    int at,au,av,at_,au_,av_;
    int d,e,f;
    int ad,ae,af;
    
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    
    //define the orbital spaces
    i0 = OrbWin[0];
    i1 = OrbWin[1];
    t0 = OrbWin[2];
    t1 = OrbWin[3];
    a0 = OrbWin[4];
    a1 = OrbWin[5];
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    
    if (!Conventional){
      
    }//!conventional
    else{
      double h__,h_,h;
      //-------------------------------
      //Generate the auxiliary A matrix
      //-------------------------------
      //open the integral container
      IKJL.OpenFileRead();
      for (d=0;d<NActive;d++){
        ad = d + NInternal;
        //get the integrals
        boost::shared_ptr<Matrix> Kdd = IKJL.GetMatrix(ad,ad);
        for (t=0;t<NActive;t++){
          h = heff_.element(d+t0,t+t0);
          at = t + NInternal;
          //get the integrals
          boost::shared_ptr<Matrix> Kdt = IKJL.GetMatrix(ad,at);
          for (t_=0;t_<NActive;t_++){
            at_ = t_ + NInternal;
            for (u_=0;u_<NActive;u_++){
              au_ = u_ + NInternal;
              for (v_=0;v_<NActive;v_++){
                av_ = v_ + NInternal;
                for (u=0;u<NActive;u++){
                  h__ = heff_.element(d+t0,u+t0);
                  au = u+ NInternal;
                  for (v=0;v<NActive;v++){
                    av = v + NInternal;
                    h_ = heff_.element(v+t0,d+t0);
                    //the one-electron part
                    AuxA_(t_,u_,v_,t,u,v) -= h  * DC3(v_,u,d,t_,u_,v);
                    AuxA_(t_,u_,v_,t,u,v) += h_ * DC3(v_,u,t,t_,u_,d);
                    AuxA_(t_,u_,v_,t,u,v) -= h__* DC3(v_,d,t,t_,u_,v);
                    for (e=0;e<NActive;e++){
                      ae = e + NInternal;
                      //the two-electron part
                      AuxA_(t_,u_,v_,t,u,v) -= Kdt->element(av,ae) * DC3(v_,u,e,t_,u_,d);
                      AuxA_(t_,u_,v_,t,u,v) += 0.5 * Kdt->element(ae,ae) * DC3(v_,u,d,t_,u_,v);
                      AuxA_(t_,u_,v_,t,u,v) += 0.5 * Kdd->element(av,ae) * DC3(v_,u,t,t_,u_,e);
                      AuxA_(t_,u_,v_,t,u,v) -= 0.5 * Kdd->element(au,ae) * DC3(v_,e,t,t_,u_,v);
                      
                      AuxA_(t_,u_,v_,t,u_,v)+= 2.0 * Kdt->element(au,ae) * DC3(v_,d,e,t_,u,v);
                      AuxA_(t_,u_,v_,v,u_,t)-= 2.0 * Kdt->element(au,ae) * DC3(v_,d,v,t_,u,e);
                      //NOTE: This actually disagrees with the paper but is in agreement with ORCA
                      AuxA_(t_,u_,v_,u,t,v) += 2.0 * Kdt->element(ae,au_)* DC3(v_,d,u,t_,e,v);
                    }//e
                    //the two-electron part (cont'd)
                    AuxA_(t_,u_,v_,t,u_,v) += 2.0 * Kdt->element(av,au) * DC2(v_,u,t_,d);
                    AuxA_(t_,u_,v_,t,u_,v) -= 1.0 * Kdt->element(au,au) * DC2(v_,d,t_,v);
                    AuxA_(t_,u_,v_,t,u_,v) -= 1.0 * Kdd->element(av,au) * DC2(v_,t,t_,u);
                    AuxA_(t_,u_,v_,t,u,v)  += 1.0 * Kdd->element(au,au_)* DC2(v_,t,t_,v);
                    
                  }//v
                  //the one-electron part (cont'd)
                  h_ = heff_.element(u+t0,d+t0);
                  AuxA_(t_,u_,v_,t,u_,u) += 2.0 * h  * DC2(v_,d,t_,u);
                  AuxA_(t_,u_,v_,t,u_,u) -= 2.0 * h_ * DC2(v_,t,t_,d);
                  AuxA_(t_,d,v_,t,u,u_)  += 2.0 * h__* DC2(v_,t,t_,u_);
                  
                }//u
              }//v_
            }//u_
          }//t_
        }//t
      }//d
      //-------------------------------
      //Generate the auxiliary B matrix
      //-------------------------------
      array_4d<double> AuxB_;
      AuxB_.resize(NActive,NActive,NActive,NActive);
      Initialize(AuxB_);
      for (d=0;d<NActive;d++){
        ad = d + NInternal;
        for (t=0;t<NActive;t++){
          at = t + NInternal;
          //get the integral
          boost::shared_ptr<Matrix> Kdt = IKJL.GetMatrix(ad,at);
          for (t_=0;t_<NActive;t_++){
            at_ = t_ + NInternal;
            for (u_=0;u_<NActive;u_++){
              au_ = u_ + NInternal;
              for (v_=0;v_<NActive;v_++){
                av_ = v_ + NInternal;
                //the one-electron part
                AuxB_(t_,u_,v_,t) -= heff.element(d+t0,t+t0) * DC2(v_,d,t_,u_);
                for (e=0;e<NActive;e++){
                  ae = e + NInternal;
                  for (f=0;f<NActive;f++){
                    af = f + NInternal;
                    AuxB_(t_,u_,v_,t) -= Kdt->element(af,ae) *  DC3(v_,e,d,t_,u_,f);
                  }//f
                  AuxB_(t_,u_,v_,t) += 2.0 * Kdt->element(ae,au_) *  DC2(v_,d,t_,e);
                }//e
              }//v_
              //the one-electron part (cont'd)
              AuxB_(t_,d,u_,t) += 2.0 * heff.element(d+t0,t+t0) * D1.element(u_,t_);
                
            }//u_
          }//t_
        }//t
      }//d
      //-------------------------------
      //Generate the auxiliary D matrix
      //-------------------------------
      Matrix AuxD_;
      AuxD_.ReSize(NActive,NActive);
      Initialize(AuxD_);
      for (d=0;d<NActive;d++){
        ad = d + NInternal;
        for (t=0;t<NActive;t++){
          at = t + NInternal;
          //get the integral
          boost::shared_ptr<Matrix> Kdt = IKJL.GetMatrix(ad,at);
          for (t_=0;t_<NActive;t_++){
            at_ = t_ + NInternal;
            //the one-electron part
            AuxD_.element(t_,t) += heff.element(at,ad) * D1_.element(d,t_);
            for (e=0;e<NActive;e++){
              ae = e + NInternal;
              for (f=0;f<NActive;f++){
                af = f + NInternal;
                //the two-electron part
                AuxD_.element(t_,t) -= Kdt->element(af,ae) * DC2(e,d,t_,f);
              }//f
              AuxD_.element(t_,t) += 2.0 * Kdt->element(ae,at_) * D1.element(d,e);
            }//e
          }//t_
        }//t
      }//d
      //---------------------------------
      // Evaluate the energy contribution
      //---------------------------------
      for (i=0;i<NInternal;i++){
        //reset the overlap and the denominator
        N = 0.0;
        D = 0.0;
        for (t_=0;t_<NActive;t_++){
          at_ = t_ + NInternal;
          h_ = heff_.element(t_+t0,i+i0);
          //get the integral
          boost::shared_ptr<Matrix> Kit_ = IKJL.GetMatrix(i,at_);
          for (t=0;t<NActive;t++){
            h = heff_.element(t+t0,i+i0);
            at = t + NInternal;
            //get the integral
            boost::shared_ptr<Matrix> Kit = IKJL.GetMatrix(i,at);
            for (u_=0;u_<NActive;u_++){
              au_ = u_ + NInternal;
              for (v_=0;v_<NActive;v_++){
                av_ = v_ + NInternal;
                x_ = Kit_->element(au_,av_);
                for (u=0;u<NActive;u++){
                  au = u+NInternal;
                  for (v=0;v<NActive;v++){
                    av = v + NInternal;
                    x = Kit->element(au,av);
                    
                    N -= x_ * x * DC3(v_,u,t,t_,u_,v);
                    D += x_ * x * AuxA_(t_,u_,v_,t,u,v);
                  }//v
                  x = Kit->element(au_,au);
                  N += 2.0 * x_ * x * DC2(v_,t,t_,u);
                }//u
                N -= 2.0 * x_ * h * DC2(v_,t,t_,u_);
                D += 2.0 * x_ * h * AuxB_(t_,u_,v_,t);
              }//v_
              x_ = Kit_->element(at,au_);
              N += 4.0 * x_ * h * D1.element(u_,t_);
            }//u_
            N += h_ * h * D1_.element(t_,t);
            D += h_ * h * AuxD_.element(t_,t);
          }//t
        }//t_
        if (N<_NEV_ZERO) continue;
        D *= 1/N;
        D -= EOrb[i+i0];
        E -= N/D;
      }//i
      //close the integral container
      IKJL.CloseFileRead();
      
    }//conventional
    
    return E;
    
  }
  
  //============================================================================
  // Evaluate the energy contribution from class V(-1,a)
  //============================================================================
  double V_a(int *OrbWin, const vector<double> &EOrb, const SpinBlock &big,
             Wavefunction &WF, IntegralContainer &IKJL, IntegralContainer & IKJA,
             boost::shared_ptr<IntegralContainer> KIAJ, const array_6d &DC3, 
             const array_4d<double> &DC2, const Matrix &D1, array_6d &AuxA, 
             const Matrix &heff_, const Matrix &heff, bool Conventional, bool ConventionalOverlap){
    
    char msg[512];
    int t,u,v;
    int i0,i1,t0,t1,a0,a1;
    int a,b;
    int tu;
    int i,j;
    double x,x_,y,y_,z,z_;
    double N=0.0;
    double D=0.0;
    double E=0.0;
    int t_,u_,v_;
    int at,au,av,at_,au_,av_;
    int d,e,f;
    int ad,ae,af;
    
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    
    //define the orbital spaces
    i0 = OrbWin[0];
    i1 = OrbWin[1];
    t0 = OrbWin[2];
    t1 = OrbWin[3];
    a0 = OrbWin[4];
    a1 = OrbWin[5];
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    
    if (!Conventional){
      
    }//!conventional
    else{
      //-------------------------------
      //Generate the auxiliary A matrix
      //-------------------------------
      //open the integral container
      IKJL.OpenFileRead();
      for (d=0;d<NActive;d++){
        ad = d + NInternal;
        //get the integrals
        boost::shared_ptr<Matrix> Kdd = IKJL.GetMatrix(ad,ad);
        for (t=0;t<NActive;t++){
          at = t + NInternal;
          //get the integrals
          boost::shared_ptr<Matrix> Kdt = IKJL.GetMatrix(ad,at);
          for (t_=0;t_<NActive;t_++){
            at_ = t_ + NInternal;
            for (u_=0;u_<NActive;u_++){
              au_ = u_ + NInternal;
              for (v_=0;v_<NActive;v_++){
                av_ = v_ + NInternal;
                for (u=0;u<NActive;u++){
                  au = u+ NInternal;
                  for (v=0;v<NActive;v++){
                    av = v + NInternal;
                    //the one-electron part
                    AuxA(t_,u_,v_,t,u,v) += heff_.element(d+t0,t+t0) * DC3(v_,u_,d,t_,u,v);
                    AuxA(t_,u_,v_,t,u,v) -= heff_.element(v+t0,d+t0) * DC3(v_,u_,t,t_,u,d);
                    AuxA(t_,u_,v_,t,u,v) -= heff_.element(u+t0,d+t0) * DC3(v_,u_,t,t_,d,v);
                    for (e=0;e<NActive;e++){
                      ae = e + NInternal;
                      x = Kdt->element(av,ae);
                      y = Kdt->element(ae,ae);
                      z = Kdd->element(av,ae);
                      //the two-electron part
                      //Note: the contracted four body terms are preevaluated elsewhere
                      AuxA(t_,u_,v_,t,u,v) += x * DC3(v_,u_,e,t_,u,d);
                      AuxA(t_,u_,v_,t,u,v) -= 0.5 * y * DC3(v_,u_,d,t_,u,v);
                      AuxA(t_,u_,v_,t,u,v) -= 0.5 * z * DC3(v_,u_,t,t_,u,e);
                      AuxA(t_,u_,v_,t,v,u) += 0.5 * z * DC3(v_,u_,t,t_,e,u);
                    }//e
                  }//v
                }//u
              }//v_
            }//u_
          }//t_
        }//t
      }//d
      //-------------------------------------
      //Generate the auxiliary B and C matrix
      //-------------------------------------
      array_4d<double> AuxB;
      AuxB.resize(NActive,NActive,NActive,NActive);
      //AuxC.resize(NActive,NActive,NActive,NActive);
      Initialize(AuxB);
      //Initialize(AuxC);
      for (d=0;d<NActive;d++){
        ad = d + NInternal;
        boost::shared_ptr<Matrix> Kdd = IKJL.GetMatrix(ad,ad);
        for (t=0;t<NActive;t++){
          at = t + NInternal;
          //get the integral
          boost::shared_ptr<Matrix> Kdt = IKJL.GetMatrix(ad,at);
          for (t_=0;t_<NActive;t_++){
            at_ = t_ + NInternal;
            for (u_=0;u_<NActive;u_++){
              au_ = u_ + NInternal;
              for (v_=0;v_<NActive;v_++){
                av_ = v_ + NInternal;
                //the one-electron part
                AuxB(t_,u_,v_,t) -= heff_.element(t+t0,d+t0) * DC2(v_,u_,t_,d);
                for (e=0;e<NActive;e++){
                  ae = e + NInternal;
                  //the two-electron part
                  AuxB(t_,u_,v_,t) += 0.5 * Kdd->element(at,ae) *  DC2(v_,u_,t_,e);
                  for (f=0;f<NActive;f++){
                    af = f + NInternal;
                    AuxB(t_,u_,v_,t) -= Kdt->element(af,ae) *  DC3(v_,u_,d,t_,e,f);
                  }//f
                }//e
              }//v_
            }//u_
          }//t_
        }//t
      }//d
      //copy the B matrix into the C matrix
      /*
       for (t_=0;t_<NActive;t_++){
        for (u_=0;u_<NActive;u_++){
          for (v_=0;v_<NActive;v_++){
            for (t=0;t<NActive;t++){
              //here we assume that we have a CAS-CI reference function
              AuxC(t,t_,u_,v_) = AuxB(t_,u_,v_,t);  
            }//t
          }//v_
        }//u_
      }//t_
       */
      //-------------------------------
      //Generate the auxiliary D matrix
      //-------------------------------
      Matrix AuxD;
      AuxD.ReSize(NActive,NActive);
      Initialize(AuxD);
      for (d=0;d<NActive;d++){
        ad = d + NInternal;
        //get the integral
        boost::shared_ptr<Matrix> Kdd = IKJL.GetMatrix(ad,ad);
        for (t=0;t<NActive;t++){
          at = t + NInternal;
          //get the integral
          boost::shared_ptr<Matrix> Kdt = IKJL.GetMatrix(ad,at);
          for (t_=0;t_<NActive;t_++){
            //the one-electron part
            AuxD.element(t_,t) -= heff.element(at,ad) * D1.element(d,t_);
            for (e=0;e<NActive;e++){
              ae = e + NInternal;
              AuxD.element(t_,t) += Kdd->element(at,ae) * D1.element(e,t_);
              for (f=0;f<NActive;f++){
                af = f + NInternal;
                AuxD.element(t_,t) -= Kdt->element(af,ae) * DC2(t_,d,e,f);
              }//f
            }//e
          }//t_
        }//t
      }//d
      //close the integral container
      IKJL.CloseFileRead();
      
      //---------------------------------
      // Evaluate the energy contribution
      //---------------------------------
      //open the integral container
      KIAJ->OpenFileRead();
      for (a=0;a<NExternal;a++){
        //reset the overlap and the denominator
        N = 0.0;
        D = 0.0;
        for (t_=0;t_<NActive;t_++){
          at_ = t_ + NInternal;
          //get the integrals
          boost::shared_ptr<Matrix> Kt_a = KIAJ->GetMatrix(at_,a);
          for (t=0;t<NActive;t++){
            at = t + NInternal;
            //get the integrals
            boost::shared_ptr<Matrix> Kta = KIAJ->GetMatrix(at,a);
            for (u_=0;u_<NActive;u_++){
              au_ = u_ + NInternal;
              for (v_=0;v_<NActive;v_++){
                av_ = v_ + NInternal;
                x_ = Kt_a->element(av_,au_);
                for (u=0;u<NActive;u++){
                  au = u+NInternal;
                  for (v=0;v<NActive;v++){
                    av = v + NInternal;
                    x = Kta->element(av,au);
                    
                    D += x_ * x * AuxA(t_,u_,v_,t,u,v);
                    N += x_ * x * DC3(v_,u_,t,t_,u,v);
                  }//v
                }//u
                //here we assume that we have a CAS-CI reference function
                D += 2.0 * x_ * heff_.element(a+a0,t+t0) * AuxB(t_,u_,v_,t);
                N += 2.0 * x_ * heff_.element(a+a0,t+t0) * DC2(v_,u_,t_,t);
              }//t
            }//v_
            D += heff_.element(a+a0,t_+t0) * heff_.element(a+a0,t+t0) * AuxD.element(t_,t);
            N += heff_.element(a+a0,t_+t0) * heff_.element(a+a0,t+t0) * D1.element(t_,t);
            
          }//u_
        }//t_
        if (N<_NEV_ZERO) continue;
        D *= 1/N;
        D += EOrb[a+a0];
        E -= N/D;
      }//a
      //close the integral conatiner
      KIAJ->CloseFileRead();
      
    }//conventional
    
    return E;
    
  }
  
  //============================================================================
  // Evaluate the energy contribution from class V(0,ijab)
  //============================================================================
  double V_ijab(int *OrbWin, vector<double> &EOrb, SpinBlock &big, Wavefunction &WF, IntegralContainer &IAJB){
    
    char msg[512];
    int i,j,a,b;
    int a_,b_;
    double E=0.0;
    double Nk=0.0;
    double D=0.0;
    double x,y;
    bool ij_eq,ab_eq;
    
    //the orbital spaces
    int i0 = OrbWin[0];
    int i1 = OrbWin[1];
    int t0 = OrbWin[2];
    int t1 = OrbWin[3];
    int a0 = OrbWin[4];
    int a1 = OrbWin[5];
    
    int NInternal = i1-i0+1;
    int NActive   = t1-t0+1;
    int NExternal = a1-a0+1;
    
    //open the integral container
    IAJB.OpenFileRead();
    
    for (i=i0;i<=i1;i++){
      for (j=i;j<=i1;j++){
        
        if (i==j){ ij_eq = true;}
        else {ij_eq = false;}
        
        //get the integrals
        boost::shared_ptr<Matrix> Kij = IAJB.GetMatrix(i-i0,j-i0);
        
        for (a=a0;a<=a1;a++){
          for (b=a;b<=a1;b++){
            a_ = a-a0;
            b_ = b-a0;
            if (a==b){ ab_eq = true;}
            else {ab_eq = false;}
            x = Kij->element(a_,b_);
            y = Kij->element(b_,a_);
            
            //Evaluate the norm
            if ((ab_eq)&&(ij_eq)){
              Nk = x * x;
            }//a==b && i==j
            else if (ij_eq){
              Nk = 2.0 * x * x;
            }
            else if (ab_eq){
              Nk = 2.0 * x * x;
            }
            else{
              Nk = 4 * (x*x+y*y-x*y);
            }
            if (Nk<_NEV_ZERO) continue;
            
            //Orb energies
            D = EOrb[a]+EOrb[b]-EOrb[i]-EOrb[j];
            //add mp2 energy contribution
            E -= Nk/D;
          }//b
        }//a
      }//j
    }//i
    //close the integral container
    IAJB.CloseFileRead();
    
    return E;
  }
  
  //============================================================================
  // Evaluate the energy contribution from class V(-1,iab)
  //============================================================================
  double V_iab(int *OrbWin, vector<double> &EOrb, SpinBlock &big, Wavefunction &WF, 
               IntegralContainer &IAJB, IntegralContainer &IKJL, Matrix &D1, 
               array_4d<double> &D2, Matrix &Heff){
    
    char msg[512];
    int i,j,a,b;
    int t,t_,u,v,w;
    double N=0.0;
    double D=0.0;
    double E=0.0;
    double x,x_,y,y_;
    double tmp;
    Matrix AuxK;
    
    //the orbital spaces
    int i0 = OrbWin[0];
    int i1 = OrbWin[1];
    int t0 = OrbWin[2];
    int t1 = OrbWin[3];
    int a0 = OrbWin[4];
    int a1 = OrbWin[5];
    
    int NInternal = i1-i0+1;
    int NActive   = t1-t0+1;
    int NExternal = a1-a0+1;

    //the buffer for the K integrals
    vector<boost::shared_ptr<Matrix> > Ki(0);
    
    //-------------------------------
    //Generate the auxiliary K matrix
    //-------------------------------
    AuxK.ReSize(NActive,NActive);
    Initialize(AuxK);
    //open the integral container
    IKJL.OpenFileRead();
    for (t_=t0;t_<=t1;t_++){
      for (t=t0;t<=t1;t++){
        for (u=t0;u<=t1;u++){
          //the one-electron term
          AuxK.element(t_-t0,t-t0) -= Heff.element(t,u) * D1.element(t_-t0,u-t0);
          //get the integral
          boost::shared_ptr<Matrix> Ktu = IKJL.GetMatrix(t-i0,u-i0);
          for (v=t0;v<=t1;v++){
            for (w=t0;w<=t1;w++){
              AuxK.element(t_-t0,t-t0) -= Ktu->element(v-i0,w-i0) * D2(t_-t0,u-t0,v-t0,w-t0);
            }//w
          }//q
        }//p
      }//t_
    }//t
    //close the integral container
    IKJL.CloseFileRead();
    
    //---------------------------------
    // Evaluate the energy contribution
    //---------------------------------
    //open the integral matrices
    IAJB.OpenFileRead();
    for (i=i0;i<=i1;i++){
      //get the integrals
      for (t=t0;t<=t1;t++){
        boost::shared_ptr<Matrix> Kit = IAJB.GetMatrix(i-i0,t-i0);
        Ki.push_back(Kit);
      }//t
      for (a=0;a<NExternal;a++){
        for (b=a;b<NExternal;b++){
          N=0.0;
          D=0.0;
          for (t_=0;t_<NActive;t_++){
            x_ = Ki[t_]->element(a,b);
            y_ = Ki[t_]->element(b,a);
            for (t=0;t<NActive;t++){
              x = Ki[t]->element(a,b);
              y = Ki[t]->element(b,a);
              if (a==b){
                tmp = x*x_;
              }//a==b
              else{
                tmp = 2.0*(x*x_+y*y_)-x_*y-x*y_;
              }//a!=b
              N += tmp * D1.element(t_,t);
              D += tmp * AuxK.element(t_,t);
            }//t_
          }//t
          if (N<_NEV_ZERO) continue;
          D *= 1.0/N;
          //the orbital energies part in the denominator
          D += EOrb[a+a0] + EOrb[b+a0] -EOrb[i];
          //Evaluate the energy contribution
          E -= N/D;
        }//b
      }//a
      //clear all integrals
      Ki.clear();
    }//i
    //close the integral container
    IAJB.CloseFileRead();
    
    return E;   
  }
  
  
  //============================================================================
  // Evaluate the energy contribution from class V(1,ija)
  //============================================================================
  double V_ija(int *OrbWin, vector<double> &EOrb, SpinBlock &big, Wavefunction &WF, 
               IntegralContainer &IKJA, IntegralContainer &IKJL, IntegralContainer &IJKL,
               Matrix &D1, Matrix &D1_, array_4d<double> &D2, Matrix &Heff){
    
    char msg[512];
    int i,j,a,b;
    int t,t_,u,v,w;
    double N=0.0;
    double D=0.0;
    double E=0.0;
    double x,x_,y,y_;
    double tmp;
    Matrix AuxK;
    int t00,t_0,u0,v0,w0;
    
    //the orbital spaces
    int i0 = OrbWin[0];
    int i1 = OrbWin[1];
    int t0 = OrbWin[2];
    int t1 = OrbWin[3];
    int a0 = OrbWin[4];
    int a1 = OrbWin[5];
    
    int NInternal = i1-i0+1;
    int NActive   = t1-t0+1;
    int NExternal = a1-a0+1;

    //the buffer for the K integrals
    vector<boost::shared_ptr<Matrix> > Ki(0);
    
    //-------------------------------
    //Generate the auxiliary K matrix
    //-------------------------------
    AuxK.ReSize(NActive,NActive);
    Initialize(AuxK);
    //open the integral container
    IKJL.OpenFileRead();
    IJKL.OpenFileRead();
    for (t_=t0;t_<=t1;t_++){
      t_0 = t_-t0;
      for (t=t0;t<=t1;t++){
        t00 = t-t0;
        //get the integrals
        boost::shared_ptr<Matrix> Jtt_ = IJKL.GetMatrix(t-i0,t_-i0);
        boost::shared_ptr<Matrix> Ktt_ = IKJL.GetMatrix(t-i0,t_-i0);
        for (u=t0;u<=t1;u++){
          u0 = u-t0;
          //the one-electron term
          AuxK.element(t_0,t00) += Heff.element(u,t) * D1_.element(u0,t_0);
          //get the integral
          boost::shared_ptr<Matrix> Ktu = IKJL.GetMatrix(t-i0,u-i0);
          for (v=t0;v<=t1;v++){
            v0=v-t0;
            AuxK.element(t_0,t00) += (2.0*Jtt_->element(u-i0,v-i0)-Ktt_->element(u-i0,v-i0))*D1.element(u0,v0);
            for (w=t0;w<=t1;w++){
              w0 = w-t0;
              AuxK.element(t_0,t00) -= Ktu->element(v-i0,w-i0) * D2(t_0,u0,v0,w0);      
            }//w
          }//q
        }//p
      }//t_
    }//t
    //close the integral container
    IKJL.CloseFileRead();
    IJKL.CloseFileRead();
    
    //---------------------------------
    // Evaluate the energy contribution
    //---------------------------------
    //open the integral container
    IKJA.OpenFileRead();
    for (i=i0;i<=i1;i++){
      for (j=i;j<=i1;j++){
        //get the integrals
        boost::shared_ptr<Matrix> Kij = IKJA.GetMatrix(i-i0,j-i0);
        boost::shared_ptr<Matrix> Kji = IKJA.GetMatrix(j-i0,i-i0);
        for (a=a0;a<=a1;a++){
          N = 0.0;
          D = 0.0;
          for (t_=t0;t_<=t1;t_++){
            t_0 = t_-t0;
            x_ = Kji->element(t_-i0,a-a0);
            y_ = Kij->element(t_-i0,a-a0);
            for (t=t0;t<=t1;t++){
              t00 = t-t0;
              x = Kji->element(t-i0,a-a0);
              y = Kij->element(t-i0,a-a0);
              
              if (i==j){
                tmp = x_*x;
              }//i==j
              else{
                tmp = 2.0 * (y_*y+x_*x) - y_*x - x_*y;
              }//i!=j
              N += tmp * D1_.element(t_0,t00);
              D += tmp * AuxK.element(t_0,t00);
            }//t
          }//t_
          if (N<_NEV_ZERO) continue;
          D *= 1.0/N;
          //the orbital energies part in the denominator
          D += EOrb[a]-EOrb[i]-EOrb[j];
          //Evaluate the energy contribution
          E -= N/D;
        }//a
      }//j
    }//i
    //close the integral container
    IKJA.CloseFileRead();
    
    return E;
  }
  
  //============================================================================
  // Evaluate the energy contribution from class V(-2,ab)
  //============================================================================
  double V_ab(int *OrbWin, const vector<double> &EOrb, const SpinBlock &big,
              Wavefunction &WF, IntegralContainer &IAJB, IntegralContainer &IKJL,
              const array_4d<double> &D2, const array_6d &D3, const Matrix &heff, 
              bool Conventional, char *BaseName, double E0, bool ConventionalOverlap,
              int MaxCore){
    
    char msg[512];
    int t,u;
    int i0,i1,t0,t1,a0,a1;
    int a,b;
    int tu;
    int i,j;
    double x,x_,y,y_;
    double N=0.0;
    double D=0.0;
    double E=0.0;
    int t_,u_;

    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    
    //define the orbital spaces
    i0 = OrbWin[0];
    i1 = OrbWin[1];
    t0 = OrbWin[2];
    t1 = OrbWin[3];
    a0 = OrbWin[4];
    a1 = OrbWin[5];
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;


    if (!Conventional){
    
      //get the quanta of the ground state Wavefunction
      SpinQuantum WFQ = WF.get_deltaQuantum(0);
      IrrepSpace Dummy(0);
      SpinQuantum TripOpQ(-2,SpinSpace(2),Dummy);
      SpinQuantum SingOpQ(-2,SpinSpace(0),Dummy);
      vector<SpinQuantum> VQTrip = WFQ + TripOpQ;
      vector<SpinQuantum> VQSing = WFQ + SingOpQ;
      
      //the container that hold the V(t,u) operators
      vector<boost::shared_ptr<WavefunctionArray> > VtuArray;
      boost::shared_ptr<WavefunctionArray> VtuSing(new WavefunctionArray);
      sprintf(msg, "%s%s", dmrginp.save_prefix().c_str(), "/VtuSing.tmp");
      VtuSing->Initialize(NActive,msg,mpi_world_size());
      VtuArray.push_back(VtuSing);
      for (int iquanta=0;iquanta<VQTrip.size();iquanta++){
        boost::shared_ptr<WavefunctionArray> VtuTrip(new WavefunctionArray);
        sprintf(msg, "%s%s.%i.tmp", dmrginp.save_prefix().c_str(), "/VtuTrip",iquanta);
        VtuTrip->Initialize(NActive,msg,mpi_world_size());
        VtuArray.push_back(VtuTrip);
      }
      //the matrix that holds the signs of the operators
      Matrix Sign(NActive,NActive);
      for (t=0;t<NActive;t++){
        for (u=0;u<NActive;u++){
          Sign.element(t,u) = 1.0;
        }//u
      }//t
      
      //evaluate the size of the buffers
      int MaxM = dmrginp.sweep_state_schedule()[dmrginp.sweep_qstate_schedule().size()];
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->CalcBufferSize(MaxCore,MaxM);
      }
      //--------------------------------------------------------------------------
      // generate V(t,u) = a(t)a(u)|psi> 
      //--------------------------------------------------------------------------
      //open the container
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->OpenFileWrite();
      }
      //--------------------------------------------------------
      //case 1: both indices are on the left side of the lattice
      //--------------------------------------------------------
      for (tu=0;tu<leftBlock->get_op_array(CRE_CRE).get_size();tu++){
        //------------------------
        //get the triplet operator
        //------------------------
        boost::shared_ptr<SparseMatrix> TripOp = leftBlock->get_op_array(CRE_CRE).get_local_element(tu)[1]->getworkingrepresentation(leftBlock);
        //get the orbital indices: Note: since we will always deal with the conjugated
        //operator, the orbital indices are reversed
        u = TripOp->get_orbs(0);
        t = TripOp->get_orbs(1);
        //----------------------------
        //generate the possible Quanta
        //----------------------------
        //get the Quanta of the Wavefunction
        SpinQuantum WFQ=WF.get_deltaQuantum(0);
        //the operator quanta
        SpinQuantum opQ=Transposeview(*TripOp).get_deltaQuantum(0);
        //the vector that holds all possible combinations of wavefunction plus operator
        std::vector<SpinQuantum> VQ  = opQ  + WFQ;
        //-------------------------------------------
        //Multiply the operator with the wavefunction
        //-------------------------------------------
        for (int iquanta=0;iquanta<VQ.size();iquanta++){
          //generate the wavefunction that holds the result of the multiplication
          boost::shared_ptr<Wavefunction> Vtu(new Wavefunction(VQ[iquanta],&big,true));
          //multiply the wavefunction with the triplet operator
          operatorfunctions::TensorMultiply(leftBlock,Transposeview(*TripOp),&big,WF,*Vtu,opQ,1.0);
          
          //------------------
          //store the operator
          //------------------
          VtuArray[iquanta+1]->AppendOperator(Vtu,t,u);
          if (t!=u){
            VtuArray[iquanta+1]->DuplicateAddress(t,u,u,t);
            Sign.element(u,t)=-1.0;
          }//t!=u
        }//iquanta
        VQ.clear();
        //------------------------
        //get the singlet operator
        //------------------------
        boost::shared_ptr<SparseMatrix> SingOp = leftBlock->get_op_array(CRE_CRE).get_local_element(tu)[0]->getworkingrepresentation(leftBlock);
        //get the orbital indices: Note: since we will always deal with the conjugated
        //operator, the orbital indices are reversed
        t = SingOp->get_orbs(1);
        u = SingOp->get_orbs(0);
        //the operator quanta
        opQ = Transposeview(*SingOp).get_deltaQuantum(0); 
        VQ = opQ + WFQ;
        //generate the wavefunction that holds the result of the multiplication
        boost::shared_ptr<Wavefunction> VtuSing(new Wavefunction(VQ[0],&big,true));
        //multiply the wavefunction with the singlet operator
        operatorfunctions::TensorMultiply(leftBlock,Transposeview(*SingOp),&big,WF,*VtuSing,opQ,-1.0);//the minus sign arises from the conjugation
        //------------------
        //store the operator
        //------------------
        VtuArray[0]->AppendOperator(VtuSing,t,u);
        if (t!=u) VtuArray[0]->DuplicateAddress(t,u,u,t);
      }//tu
      //--------------------------------------------------------
      //case 2: both indices are on the left side of the lattice
      //--------------------------------------------------------
      for (tu=0;tu<rightBlock->get_op_array(CRE_CRE).get_size();tu++){
        //------------------------
        //get the triplet operator
        //------------------------
        boost::shared_ptr<SparseMatrix> TripOp = rightBlock->get_op_array(CRE_CRE).get_local_element(tu)[1]->getworkingrepresentation(rightBlock);
        //get the orbital indices: Note: since we will always deal with the conjugated
        //operator, the orbital indices are reversed
        u = TripOp->get_orbs(0);
        t = TripOp->get_orbs(1);
        //----------------------------
        //generate the possible Quanta
        //----------------------------
        //get the Quanta of the Wavefunction
        SpinQuantum WFQ=WF.get_deltaQuantum(0);
        //the operator quanta
        SpinQuantum opQ=Transposeview(*TripOp).get_deltaQuantum(0);
        //the vector that holds all possible combinations of wavefunction plus operator
        std::vector<SpinQuantum> VQ = opQ + WFQ;
        //-------------------------------------------
        //Multiply the operator with the wavefunction
        //-------------------------------------------
        for (int iquanta=0;iquanta<VQ.size();iquanta++){
          //generate the wavefunction that holds the result of the multiplication
          boost::shared_ptr<Wavefunction> Vtu(new Wavefunction(VQ[iquanta],&big,true));
          //multiply the wavefunction with the triplet operator
          operatorfunctions::TensorMultiply(rightBlock,Transposeview(*TripOp),&big,WF,*Vtu,opQ,1.0);
          //------------------
          //store the operator
          //------------------
          VtuArray[iquanta+1]->AppendOperator(Vtu,t,u);
          if (t!=u){
            VtuArray[iquanta+1]->DuplicateAddress(t,u,u,t);
            Sign.element(u,t)=-1.0;
          }//t!=u
        }//iquanta
        VQ.clear();
        //------------------------
        //get the singlet operator
        //------------------------
        boost::shared_ptr<SparseMatrix> SingOp = rightBlock->get_op_array(CRE_CRE).get_local_element(tu)[0]->getworkingrepresentation(rightBlock);
        //get the orbital indices: Note: since we will always deal with the conjugated
        //operator, the orbital indices are reversed
        u = SingOp->get_orbs(0);
        t = SingOp->get_orbs(1);
        //the operator quanta
        opQ = Transposeview(*SingOp).get_deltaQuantum(0); 
        VQ = opQ + WFQ;
        //generate the wavefunction that holds the result of the multiplication
        boost::shared_ptr<Wavefunction> VtuSing(new Wavefunction(VQ[0],&big,true));
        //multiply the wavefunction with the singlet operator
        operatorfunctions::TensorMultiply(rightBlock,Transposeview(*SingOp),&big,WF,*VtuSing,opQ,-1.0);//the minus sign arises from the conjugation
        //------------------
        //store the operator
        //------------------
        VtuArray[0]->AppendOperator(VtuSing,t,u);
        if (t!=u) VtuArray[0]->DuplicateAddress(t,u,u,t);
      }//tu
      //--------------------------------------------------------
      //case 3: the indices are on the both sides of the lattice
      //--------------------------------------------------------
      int tstart,tstop;
      //divide the loop over the parallel processes
      PALDivideLoop(tstart,tstop,0,leftBlock->get_op_array(CRE).get_size());
      //loop over local elements on the left side
      for (int tcount=tstart;tcount<tstop;tcount++){
        //get the operator
        boost::shared_ptr<SparseMatrix> lop = leftBlock->get_op_array(CRE).get_local_element(tcount)[0]->getworkingrepresentation(leftBlock);
        //get the orbital index
        u = lop->get_orbs(0);
        for (int ucount=0;ucount<rightBlock->get_op_array(CRE).size();ucount++){
          //get the operator
          boost::shared_ptr<SparseMatrix> rop = rightBlock->get_op_array(CRE).get_local_element(ucount)[0]->getworkingrepresentation(rightBlock);
          //get the orbital index
          t = rop->get_orbs(0);
          //get the Quanta of the Wavefunction
          SpinQuantum WFQ=WF.get_deltaQuantum(0);
          //the possible operator quanta -singlet and triplet-
          std::vector<SpinQuantum> opQu = Transposeview(*lop).get_deltaQuantum(0) + Transposeview(*rop).get_deltaQuantum(0); 
          //----------------
          //triplet operator
          //----------------
          SpinQuantum opQ = opQu[1];
          //the vector that holds all possible combinations of wavefunction plus operator
          std::vector<SpinQuantum> VQ = opQ + WFQ;
          for (int iquanta=0;iquanta<VQ.size();iquanta++){
            //generate the wavefunction that holds the result of the multiplication
            boost::shared_ptr<Wavefunction> Vtu(new Wavefunction(VQ[iquanta],&big,true));
            //multiply the wavefunction with the triplet operator
            operatorfunctions::TensorMultiply(leftBlock,Transposeview(*lop),Transposeview(*rop),&big,WF,*Vtu,opQ,1.0);
            //------------------
            //store the operator
            //------------------
            VtuArray[iquanta+1]->AppendOperator(Vtu,t,u);
            if (t!=u){
              VtuArray[iquanta+1]->DuplicateAddress(t,u,u,t);
              Sign.element(u,t)=-1.0;
            }//t!=u
          }//iquanta
          VQ.clear();
          //----------------
          //singlet operator
          //----------------
          //the operator quanta
          opQ = opQu[0];
          VQ = opQ + WFQ;
          //generate the wavefunction that holds the result of the multiplication
          boost::shared_ptr<Wavefunction> VtuSing(new Wavefunction(VQ[0],&big,true));
          //multiply the wavefunction with the triplet operator
          operatorfunctions::TensorMultiply(leftBlock,Transposeview(*lop),Transposeview(*rop),&big,WF,*VtuSing,opQ,1.0);
          //------------------
          //store the operator
          //------------------
          VtuArray[0]->AppendOperator(VtuSing,t,u);
          VtuArray[0]->DuplicateAddress(t,u,u,t);
        }//u
      }//t
      
      //close the container
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->CloseFileWrite();
      }
      
      //invert the index order of the two-electron integrals
      sprintf(msg,"%s.block.AIBJ.tmp",BaseName);
      boost::shared_ptr<IntegralContainer> AIBJ=IAJB.ReverseOrder(msg);
      
      //------------------------------------------------------------------------
      // Construct the action of the Hamiltonian on the V(t,u) functions
      // sigma(t,u) = H * V(t,u)
      //------------------------------------------------------------------------
      //the container that hold the sigma(t,u) functions
      vector<boost::shared_ptr<WavefunctionArray> > SigmaArray;
      boost::shared_ptr<WavefunctionArray> SigmaSing(new WavefunctionArray);
      sprintf(msg, "%s%s", dmrginp.save_prefix().c_str(), "/VtuSing.tmp");
      SigmaSing->Initialize(NActive,msg,mpi_world_size());
      SigmaSing->CalcBufferSize(MaxCore,MaxM);
      SigmaArray.push_back(SigmaSing);
      for (int iquanta=0;iquanta<VQTrip.size();iquanta++){
        boost::shared_ptr<WavefunctionArray> SigmaTrip(new WavefunctionArray);
        sprintf(msg, "%s%s.%i.tmp", dmrginp.save_prefix().c_str(), "/SigmatuTrip",iquanta);
        SigmaTrip->Initialize(NActive,msg,mpi_world_size());
        SigmaTrip->CalcBufferSize(MaxCore,MaxM);
        SigmaArray.push_back(SigmaTrip);
      }
      //generate the sigma(t,u) functions
      GenerateActionOfH(WF,VtuArray,SigmaArray,TripOpQ,SingOpQ,big,OrbWin);
      
      //--------------------------------------------------------------------------
      // construct Vab = gamma(tu){fac * (au|bt) * Ttu}
      // and calculate the energy contribution
      //--------------------------------------------------------------------------
      //---------------
      // the prefactors
      //---------------
      vector<double> TripFac;
      for (i=0;i<VQTrip.size();i++){
        EvalCompFactors(TripFac, 2, VQTrip[i].get_s().getirrep(), WFQ.get_s().getirrep());
        /*
        double twoS_ = (double) VQTrip[i].get_s().getirrep();
        double twoS  = (double) WFQ.get_s().getirrep();
        double S_    = twoS_ / 2.0;
        double S     = twoS / 2.0;
        TripFac.push_back(sqrt((2.0*S_+1)/3.0*(2.0*S+1)));
         */
      }
      //open the integral and operator container
      AIBJ->OpenFileRead();
      //open the container
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->OpenFileRead();VtuArray[i]->ResetBuffer();
        SigmaArray[i]->OpenFileRead();SigmaArray[i]->ResetBuffer();
      }
      for (a=0;a<NExternal;a++){
        //for (b=a;b<NExternal;b++){
        //establish batching over second index b that runs from a to NExternal-1
        vector<pair<int,int> > bbatch;
        EstablishBatching(bbatch,MaxCore,MaxM,a,NExternal,VQTrip.size()+1);
        for (int ibatch=0;ibatch<bbatch.size();ibatch++){
          mpi_barrier();
          //generate the vector that holds Vsing,VTrip,SigmaSing and sigmaTrip
          //for the entire batch
          vector<boost::shared_ptr<Wavefunction> > Vab_S;
          vector<boost::shared_ptr<Wavefunction> > SigmaS;
          vector<vector<boost::shared_ptr<Wavefunction> > > Vab_T;
          vector<vector<boost::shared_ptr<Wavefunction> > > SigmaT;
          for (b=bbatch[ibatch].first;b<bbatch[ibatch].second;b++){
            //------------------------------
            //the singlet perturber function
            //------------------------------
            boost::shared_ptr<Wavefunction> VabSing(new Wavefunction(VQSing[0],&big,true));
            boost::shared_ptr<Wavefunction> SigmaSing(new Wavefunction(VQSing[0],&big,true));
            //-------------------------------
            //the triplet perturber functions
            //-------------------------------
            vector<boost::shared_ptr<Wavefunction> > VabTrip;
            vector<boost::shared_ptr<Wavefunction> > SigmaTrip;
            for (i=0;i<VQTrip.size();i++){
              boost::shared_ptr<Wavefunction> VabT(new Wavefunction(VQTrip[i],&big,true));
              VabTrip.push_back(VabT);
              boost::shared_ptr<Wavefunction> sigmaTrip(new Wavefunction(VQTrip[i],&big,true));
              SigmaTrip.push_back(sigmaTrip);
            }
            //----------------
            //get the integral
            //----------------
            boost::shared_ptr<Matrix> Kab = AIBJ->GetMatrix(a,b);

            //-------------
            //construct Vab
            //-------------
            /*
            for (t=0;t<NActive;t++){
              for (u=0;u<NActive;u++){
                //get the operators in active space
                boost::shared_ptr<Wavefunction> Ttusing  = T_int1.GetOperator(t,u);
                boost::shared_ptr<Wavefunction> Ttutrip1 = T_int2.GetOperator(t,u);
                boost::shared_ptr<Wavefunction> Ttutrip2 = T_int3.GetOperator(t,u);
                boost::shared_ptr<Wavefunction> Ttutrip3 = T_int4.GetOperator(t,u);

                //the integral (au|bt) times the sign of the operator
                x = Kab->element(u+NInternal,t+NInternal) * Sign.element(u,t);
                y = Kab->element(u+NInternal,t+NInternal);

                //add them to the according perturber functions
                ScaleAdd(y,*Ttusing,*VabSing);
                ScaleAdd(x*TripFac[0],*Ttutrip1,*(VabTrip[0]));
                if (VabTrip.size()>1) ScaleAdd(x*TripFac[1],*Ttutrip2,*(VabTrip[1]));
                if (VabTrip.size()>2) ScaleAdd(x*TripFac[2],*Ttutrip3,*(VabTrip[2]));
              }//u
            }//t
             */
            for (int iquanta=0;iquanta<VtuArray.size();iquanta++){
              bool EndOfArray=false;
              bool EndOfArray_=false;
              while (!EndOfArray){
                boost::shared_ptr<Wavefunction> Vtu  = VtuArray[iquanta]->GetOpFromBuffer(t,u,EndOfArray);
                boost::shared_ptr<Wavefunction> Sigmatu = SigmaArray[iquanta]->GetOpFromBuffer(t,u,EndOfArray_);
                if (EndOfArray!=EndOfArray_){sprintf(msg,"\n!!!Error NEVPT2 in V(a,b): VtuArray and SigmaArray are not synchronized!!!");pout << msg;}
                if (!EndOfArray){
                  //the integral (au|bt) times the sign of the operator
                  if (t!=u){
                    x = Kab->element(u+NInternal,t+NInternal) * Sign.element(u,t) + Kab->element(t+NInternal,u+NInternal) * Sign.element(t,u); 
                    y = Kab->element(u+NInternal,t+NInternal) + Kab->element(t+NInternal,u+NInternal);
                  }
                  else{
                    x = Kab->element(u+NInternal,t+NInternal) * Sign.element(u,t);
                    y = Kab->element(u+NInternal,t+NInternal);
                  }
                  if (iquanta==0){
                    ScaleAdd(y,*Vtu,*VabSing);
                    ScaleAdd(y,*Sigmatu,*SigmaSing);
                  }//singlet operator
                  else{
                    int itmp=iquanta-1;
                    ScaleAdd(x*TripFac[itmp],*Vtu,*(VabTrip[itmp]));
                    ScaleAdd(x*TripFac[itmp],*Sigmatu,*(SigmaTrip[itmp]));
                  }//triplet operator
                }
              }//tu
            }//iquanta
            //gather all functions
            Vab_S.push_back(VabSing);
            SigmaS.push_back(SigmaSing);
            Vab_T.push_back(VabTrip);
            SigmaT.push_back(SigmaTrip);
          }//b
          //-----------------------------------------------------
          //Add up all contributions from the different processes
          //-----------------------------------------------------
          AddPalWavefunctions(Vab_T,SigmaT,Vab_S,SigmaS);
          for (b=bbatch[ibatch].first;b<bbatch[ibatch].second;b++){
            //get the synchronized vectors
            int b0 = b-bbatch[ibatch].first;
            boost::shared_ptr<Wavefunction> VabSing = Vab_S[b0];
            boost::shared_ptr<Wavefunction> SigmaSing = SigmaS[b0];
            vector<boost::shared_ptr<Wavefunction> > VabTrip(Vab_T[b0]);
            vector<boost::shared_ptr<Wavefunction> > SigmaTrip(SigmaT[b0]);
            //---------------------------
            //evaluate the overlap of Vab
            //---------------------------
            //get the integral
            boost::shared_ptr<Matrix> Kab = AIBJ->GetMatrix(a,b);
            N = 0.0;
            if (ConventionalOverlap){
              for (t_=0;t_<NActive;t_++){
                for (u_=0;u_<NActive;u_++){
                  x = Kab->element(u_+NInternal,t_+NInternal);
                  for (t=0;t<NActive;t++){
                    for (u=0;u<NActive;u++){
                      y  = Kab->element(u+NInternal,t+NInternal);
                      if (a==b){
                        N += 0.5*x*y*D2(t_,u_,t,u);
                      }
                      else {
                        N += x*y*D2(t_,u_,t,u);
                      }
                    }//u
                  }//t
                }//u_
              }//t_
            }//ConventionalOverlap
            else{
              //check if one of the triplet quanta equals the singlet quantum
              int EqualElement = CheckEquality(VQSing[0],VQTrip);
              //if yes, add the corresponding Vab and Sigma functions together
              if (EqualElement!=-1){
                ScaleAdd(1.0,*(VabTrip[EqualElement]),*VabSing);
                ScaleAdd(1.0,*(SigmaTrip[EqualElement]),*SigmaSing);
                //avoid double counting
                VabTrip[EqualElement]->Clear();
                SigmaTrip[EqualElement]->Clear();
              }
              double NSing,NTrip;
              //evaluate the overlap
              NSing = DotProduct(*VabSing,*VabSing);
              NTrip = DotProduct(*VabTrip[0].get(),*VabTrip[0].get());
              if (VabTrip.size()>1) NTrip += DotProduct(*VabTrip[1].get(),*VabTrip[1].get());
              if (VabTrip.size()>2) NTrip += DotProduct(*VabTrip[2].get(),*VabTrip[2].get());
              N = NSing + NTrip;
              if (a==b) N *= 0.5;
            }//!ConventionalOverlap

            if (N<_NEV_ZERO) continue;
            //---------------------------------------------------
            //evaluate the active contribution to the denominator
            //---------------------------------------------------
            D = 0.0;
            D += DotProduct(*VabSing,*SigmaSing);
            D += DotProduct(*VabTrip[0].get(),*SigmaTrip[0]); 
            if (VabTrip.size()>1) D += DotProduct(*VabTrip[1].get(),*SigmaTrip[1]);
            if (VabTrip.size()>2) D += DotProduct(*VabTrip[2].get(),*SigmaTrip[2]); 
            if (a==b) D *= 0.5;
            D *= 1/N;
            //-----------------------------------------------------
            //evaluate the external contribution to the denominator
            //-----------------------------------------------------
            D += EOrb[a+a0] + EOrb[b+a0] - E0;
            //sprintf(msg,"\n%i %i  %4.12lf %4.12lf  %4.12lf  %4.12lf %4.12lf",a,b,N,D-EOrb[a+a0]-EOrb[b+a0]+E0,EOrb[a+a0]+EOrb[b+a0],E0,-N/D);pout << msg;
            //---------------------------------------
            //calculate the total energy contribution
            //---------------------------------------
            E -= N/D;
          }//b
          //clear memory
          for (int btmp=0;btmp<Vab_S.size();btmp++){
            Vab_S[btmp]->Clear();
            SigmaS[btmp]->Clear();
            for (int iquanta=0;iquanta<Vab_T[btmp].size();iquanta++){
              Vab_T[btmp][iquanta]->Clear();
              SigmaT[btmp][iquanta]->Clear();
            }//btmp
          }//iquanta
          Vab_S.clear();
          Vab_T.clear();
          SigmaS.clear();
          SigmaT.clear();
        }//ibatch
      }//a
      //close the container
      AIBJ->CloseFileRead();
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->CloseFileRead();
        SigmaArray[i]->CloseFileRead();
      }
      //--------
      //clean up
      //--------
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->Clear();
        SigmaArray[i]->Clear();
      }
      AIBJ->Clear();
    }//not conventional
    else{
      //---------------------------------
      // construct the auxiliary K matrix
      //---------------------------------
      IKJL.OpenFileRead();
      int d,e,f;
      array_4d<double> AuxK;
      AuxK.resize(NActive,NActive,NActive,NActive);
      Initialize(AuxK);
      for (t_=0;t_<NActive;t_++){
        for (u_=0;u_<NActive;u_++){
          for (t=0;t<NActive;t++){
            for (u=0;u<NActive;u++){
              for (d=0;d<NActive;d++){
                AuxK(t_,u_,t,u) -= heff.element(t+NInternal,d+NInternal) * D2(t_,u_,d,u);
                AuxK(t_,u_,t,u) -= heff.element(u+NInternal,d+NInternal) * D2(t_,u_,t,d);
              }//d
            }//u
          }//t
        }//u_
      }//t_
      for (d=0;d<NActive;d++){
        for (t=0;t<NActive;t++){
          boost::shared_ptr<Matrix> Kdt = IKJL.GetMatrix(d+NInternal,t+NInternal);
          for (u=0;u<NActive;u++){
            boost::shared_ptr<Matrix> Kdu = IKJL.GetMatrix(d+NInternal,u+NInternal);
            boost::shared_ptr<Matrix> Ktu = IKJL.GetMatrix(t+NInternal,u+NInternal);
            for (t_=0;t_<NActive;t_++){
              for (u_=0;u_<NActive;u_++){
                for (e=0;e<NActive;e++){
                  for (f=0;f<NActive;f++){
                    x = Kdt->element(e+NInternal,f+NInternal);
                    AuxK(t_,u_,t,u) -= x * D3(t_,u_,d,f,u,e);
                    y = Kdu->element(e+NInternal,f+NInternal);
                    AuxK(t_,u_,t,u) -= y * D3(t_,u_,d,t,f,e);
                  }//f
                  AuxK(t_,u_,t,u) -= Ktu->element(d+NInternal,e+NInternal) * D2(t_,u_,d,e);
                }//e
              }//u
            }//t
          }//u_
        }//t_
      }//d
      IKJL.CloseFileRead();
      //---------------------------------
      //Calculate the energy contribution
      //---------------------------------
      //invert the index order of the two-electron integrals
      sprintf(msg,"%s.block.AIBJ.tmp",BaseName);
      boost::shared_ptr<IntegralContainer> AIBJ=IAJB.ReverseOrder(msg);
      //open the new integral container
      AIBJ->OpenFileRead();
      for (a=0;a<NExternal;a++){
        for (b=a;b<NExternal;b++){
          bool ab_equal = (a==b);
          //get the integral
          boost::shared_ptr<Matrix> Kab = AIBJ->GetMatrix(a,b);
          //reset the overlap and denominator 
          N = 0.0;
          D = 0.0;
          for (t_=0;t_<NActive;t_++){
            for (u_=0;u_<NActive;u_++){
              x_ = Kab->element(u_+NInternal,t_+NInternal);
              for (t=0;t<NActive;t++){
                for (u=0;u<NActive;u++){
                  x = Kab->element(u+NInternal,t+NInternal);
                  if (ab_equal){
                    N += x*x_*D2(t_,u_,t,u) * 0.5;
                    D += x*x_*AuxK(t_,u_,t,u) * 0.5;
                  }//a==b
                  else {
                    N += x*x_*D2(t_,u_,t,u);
                    D += x*x_*AuxK(t_,u_,t,u);
                  }
                }//u
              }//t
            }//u_
          }//t_
          if (N<_NEV_ZERO) continue;
          D *= 1/N;
          D += EOrb[a+a0] + EOrb[b+a0]; 
          E -= N/D;
        }//b
      }//a
      //close the integral container and clean up disk
      AIBJ->CloseFileRead();
      AIBJ->Clear();
    }//conventional treatment
    
    return E;
    
  }
  
  //============================================================================
  // Evaluate the energy contribution from class V(2,ij)
  //============================================================================
  double V_ij(int *OrbWin, const vector<double> &EOrb, const SpinBlock &big,
              Wavefunction &WF, IntegralContainer &IKJL, const array_4d<double> &D2_, 
              const array_6d &D3_, const Matrix &heff_, bool Conventional, double E0, 
              bool ConventionalOverlap, int MaxCore){
    
    char msg[512];
    int t,u;
    int i0,i1,t0,t1,a0,a1;
    int a,b;
    int tu;
    int i,j;
    double x,x_,y,y_;
    double N=0.0;
    double D=0.0;
    double E=0.0;
    int t_,u_;

    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    
    //define the orbital spaces
    i0 = OrbWin[0];
    i1 = OrbWin[1];
    t0 = OrbWin[2];
    t1 = OrbWin[3];
    a0 = OrbWin[4];
    a1 = OrbWin[5];
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    
    if (!Conventional){
      //-----------------------------------------------
      //get the quanta of the ground state Wavefunction
      //-----------------------------------------------
      SpinQuantum WFQ = WF.get_deltaQuantum(0);
      IrrepSpace Dummy(0);
      SpinQuantum TripOpQ(2,SpinSpace(2),Dummy);
      SpinQuantum SingOpQ(2,SpinSpace(0),Dummy);
      vector<SpinQuantum> VQTrip = WFQ + TripOpQ;
      vector<SpinQuantum> VQSing = WFQ + SingOpQ;

      //-----------------------------
      //Generate the V(t,u) container
      //-----------------------------
      //the container that hold the V(t,u) operators
      vector<boost::shared_ptr<WavefunctionArray> > VtuArray;
      boost::shared_ptr<WavefunctionArray> VtuSing(new WavefunctionArray);
      sprintf(msg, "%s%s", dmrginp.save_prefix().c_str(), "/VtuSing.tmp");
      VtuSing->Initialize(NActive,msg,mpi_world_size());
      VtuArray.push_back(VtuSing);
      for (int iquanta=0;iquanta<VQTrip.size();iquanta++){
        boost::shared_ptr<WavefunctionArray> VtuTrip(new WavefunctionArray);
        sprintf(msg, "%s%s.%i.tmp", dmrginp.save_prefix().c_str(), "/VtuTrip",iquanta);
        VtuTrip->Initialize(NActive,msg,mpi_world_size());
        VtuArray.push_back(VtuTrip);
      }
      //the matrix that holds the signs of the operators
      Matrix Sign(NActive,NActive);
      for (t=0;t<NActive;t++){
        for (u=0;u<NActive;u++){
          Sign.element(t,u) = 1.0;
        }//u
      }//t
      //evaluate the size of the buffers
      int MaxM = dmrginp.sweep_state_schedule()[dmrginp.sweep_qstate_schedule().size()];
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->CalcBufferSize(MaxCore,MaxM);
      }
      //--------------------------------------------------------------------------
      // generate V(t,u) = a+(t)a+(u)|psi> 
      //--------------------------------------------------------------------------
      //open the container
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->OpenFileWrite();
      }
      //--------------------------------------------------------
      //case 1: both indices are on the left side of the lattice
      //--------------------------------------------------------
      for (tu=0;tu<leftBlock->get_op_array(CRE_CRE).get_size();tu++){
        //------------------------
        //get the triplet operator
        //------------------------
        boost::shared_ptr<SparseMatrix> TripOp = leftBlock->get_op_array(CRE_CRE).get_local_element(tu)[1]->getworkingrepresentation(leftBlock);
        //get the orbital indices
        t = TripOp->get_orbs(0);
        u = TripOp->get_orbs(1);
        //----------------------------
        //generate the possible Quanta
        //----------------------------
        //get the Quanta of the Wavefunction
        SpinQuantum WFQ=WF.get_deltaQuantum(0);
        //the operator quanta
        SpinQuantum opQ=TripOp->get_deltaQuantum(0);
        //the vector that holds all possible combinations of wavefunction plus operator
        std::vector<SpinQuantum> VQ  = opQ  + WFQ;
        //-------------------------------------------
        //Multiply the operator with the wavefunction
        //-------------------------------------------
        for (int iquanta=0;iquanta<VQ.size();iquanta++){
          //generate the wavefunction that holds the result of the multiplication
          boost::shared_ptr<Wavefunction> Vtu(new Wavefunction(VQ[iquanta],&big,true));
          //multiply the wavefunction with the triplet operator
          operatorfunctions::TensorMultiply(leftBlock,*TripOp,&big,WF,*Vtu,opQ,1.0);
          
          //------------------
          //store the operator
          //------------------
          VtuArray[iquanta+1]->AppendOperator(Vtu,t,u);
          if (t!=u){
            VtuArray[iquanta+1]->DuplicateAddress(t,u,u,t);
            Sign.element(u,t)=-1.0;
          }//t!=u
        }//iquanta
        VQ.clear();
        //------------------------
        //get the singlet operator
        //------------------------
        boost::shared_ptr<SparseMatrix> SingOp = leftBlock->get_op_array(CRE_CRE).get_local_element(tu)[0]->getworkingrepresentation(leftBlock);
        //get the orbital indices
        t = SingOp->get_orbs(0);
        u = SingOp->get_orbs(1);
        //the operator quanta
        opQ = SingOp->get_deltaQuantum(0); 
        VQ = opQ + WFQ;
        //generate the wavefunction that holds the result of the multiplication
        boost::shared_ptr<Wavefunction> VtuSing(new Wavefunction(VQ[0],&big,true));
        //multiply the wavefunction with the singlet operator
        operatorfunctions::TensorMultiply(leftBlock,*SingOp,&big,WF,*VtuSing,opQ,1.0);
        //------------------
        //store the operator
        //------------------
        VtuArray[0]->AppendOperator(VtuSing,t,u);
        if (t!=u) VtuArray[0]->DuplicateAddress(t,u,u,t);
      }//tu
      //---------------------------------------------------------
      //case 2: both indices are on the right side of the lattice
      //---------------------------------------------------------
      for (tu=0;tu<rightBlock->get_op_array(CRE_CRE).get_size();tu++){
        //------------------------
        //get the triplet operator
        //------------------------
        boost::shared_ptr<SparseMatrix> TripOp = rightBlock->get_op_array(CRE_CRE).get_local_element(tu)[1]->getworkingrepresentation(rightBlock);
        //get the orbital indices
        t = TripOp->get_orbs(0);
        u = TripOp->get_orbs(1);
        //----------------------------
        //generate the possible Quanta
        //----------------------------
        //get the Quanta of the Wavefunction
        SpinQuantum WFQ=WF.get_deltaQuantum(0);
        //the operator quanta
        SpinQuantum opQ=TripOp->get_deltaQuantum(0);
        //the vector that holds all possible combinations of wavefunction plus operator
        std::vector<SpinQuantum> VQ  = opQ  + WFQ;
        //-------------------------------------------
        //Multiply the operator with the wavefunction
        //-------------------------------------------
        for (int iquanta=0;iquanta<VQ.size();iquanta++){
          //generate the wavefunction that holds the result of the multiplication
          boost::shared_ptr<Wavefunction> Vtu(new Wavefunction(VQ[iquanta],&big,true));
          //multiply the wavefunction with the triplet operator
          operatorfunctions::TensorMultiply(rightBlock,*TripOp,&big,WF,*Vtu,opQ,1.0);
                  
          //------------------
          //store the operator
          //------------------
          VtuArray[iquanta+1]->AppendOperator(Vtu,t,u);
          if (t!=u){
            VtuArray[iquanta+1]->DuplicateAddress(t,u,u,t);
            Sign.element(u,t)=-1.0;
          }//t!=u
        }//iquanta
        VQ.clear();
        //------------------------
        //get the singlet operator
        //------------------------
        boost::shared_ptr<SparseMatrix> SingOp = rightBlock->get_op_array(CRE_CRE).get_local_element(tu)[0]->getworkingrepresentation(rightBlock);
        //get the orbital indices
        t = SingOp->get_orbs(0);
        u = SingOp->get_orbs(1);
        //the operator quanta
        opQ = SingOp->get_deltaQuantum(0); 
        VQ = opQ + WFQ;
        //generate the wavefunction that holds the result of the multiplication
        boost::shared_ptr<Wavefunction> VtuSing(new Wavefunction(VQ[0],&big,true));
        //multiply the wavefunction with the singlet operator
        operatorfunctions::TensorMultiply(rightBlock,*SingOp,&big,WF,*VtuSing,opQ,1.0);
        //------------------
        //store the operator
        //------------------
        VtuArray[0]->AppendOperator(VtuSing,t,u);
        if (t!=u) VtuArray[0]->DuplicateAddress(t,u,u,t);
      }//tu
      //--------------------------------------------------------
      //case 3: the indices are on the both sides of the lattice
      //--------------------------------------------------------
      int tstart,tstop;
      //divide the loop over the parallel processes
      PALDivideLoop(tstart,tstop,0,leftBlock->get_op_array(CRE).get_size());
      //loop over local elements on the left side
      for (int tcount=tstart;tcount<tstop;tcount++){
      //for (int tcount=0;tcount<leftBlock->get_op_array(CRE).size();tcount++){
        //get the operator
        boost::shared_ptr<SparseMatrix> lop = leftBlock->get_op_array(CRE).get_local_element(tcount)[0]->getworkingrepresentation(leftBlock);
        //get the orbital index
        t = lop->get_orbs(0);
        for (int ucount=0;ucount<rightBlock->get_op_array(CRE).size();ucount++){
          //get the operator
          boost::shared_ptr<SparseMatrix> rop = rightBlock->get_op_array(CRE).get_local_element(ucount)[0]->getworkingrepresentation(rightBlock);
          //get the orbital index
          u = rop->get_orbs(0);
          //get the Quanta of the Wavefunction
          SpinQuantum WFQ=WF.get_deltaQuantum(0);
          //the possible operator quanta -singlet and triplet-
          std::vector<SpinQuantum> opQu = lop->get_deltaQuantum(0) + rop->get_deltaQuantum(0); 
          //----------------
          //triplet operator
          //----------------
          SpinQuantum opQ = opQu[1];
          //the vector that holds all possible combinations of wavefunction plus operator
          std::vector<SpinQuantum> VQ = opQ + WFQ;
          for (int iquanta=0;iquanta<VQ.size();iquanta++){
            //generate the wavefunction that holds the result of the multiplication
            boost::shared_ptr<Wavefunction> Vtu(new Wavefunction(VQ[iquanta],&big,true));
            //multiply the wavefunction with the triplet operator
            operatorfunctions::TensorMultiply(leftBlock,*lop,*rop,&big,WF,*Vtu,opQ,1.0);
            //------------------
            //store the operator
            //------------------
            VtuArray[iquanta+1]->AppendOperator(Vtu,t,u);
            if (t!=u){
              VtuArray[iquanta+1]->DuplicateAddress(t,u,u,t);
              Sign.element(u,t)=-1.0;
            }//t!=u
          }//iquanta
          VQ.clear();
          //----------------
          //singlet operator
          //----------------
          //the operator quanta
          opQ = opQu[0];
          VQ = opQ + WFQ;
          //generate the wavefunction that holds the result of the multiplication
          boost::shared_ptr<Wavefunction> VtuSing(new Wavefunction(VQ[0],&big,true));
          //multiply the wavefunction with the triplet operator
          operatorfunctions::TensorMultiply(leftBlock,*lop,*rop,&big,WF,*VtuSing,opQ,1.0);
          //------------------
          //store the operator
          //------------------
          VtuArray[0]->AppendOperator(VtuSing,t,u);
          if (t!=u) VtuArray[0]->DuplicateAddress(t,u,u,t);
        }//u
      }//t
      //close the container
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->CloseFileWrite();
      }
      
      //------------------------------------------------------------------------
      // Construct the action of the Hamiltonian on the V(t,u) functions
      // sigma(t,u) = H * V(t,u)
      //------------------------------------------------------------------------
      //the container that hold the sigma(t,u) functions
      vector<boost::shared_ptr<WavefunctionArray> > SigmaArray;
      boost::shared_ptr<WavefunctionArray> SigmaSing(new WavefunctionArray);
      sprintf(msg, "%s%s", dmrginp.save_prefix().c_str(), "/VtuSing.tmp");
      SigmaSing->Initialize(NActive,msg,mpi_world_size());
      SigmaSing->CalcBufferSize(MaxCore,MaxM);
      SigmaArray.push_back(SigmaSing);
      for (int iquanta=0;iquanta<VQTrip.size();iquanta++){
        boost::shared_ptr<WavefunctionArray> SigmaTrip(new WavefunctionArray);
        sprintf(msg, "%s%s.%i.tmp", dmrginp.save_prefix().c_str(), "/SigmatuTrip",iquanta);
        SigmaTrip->Initialize(NActive,msg,mpi_world_size());
        SigmaTrip->CalcBufferSize(MaxCore,MaxM);
        SigmaArray.push_back(SigmaTrip);
      }
      //generate the sigma(t,u) functions
      GenerateActionOfH(WF,VtuArray,SigmaArray,TripOpQ,SingOpQ,big,OrbWin);
      
      //--------------------------------------------------------------------------
      // construct Vij = gamma(ij){fac * (iu|jt) * Ttu}
      // and calculate the energy contribution
      //--------------------------------------------------------------------------
      //---------------
      // the prefactors
      //---------------
      vector<double> TripFac;
      for (i=0;i<VQTrip.size();i++){
        EvalCompFactors(TripFac, 2, VQTrip[i].get_s().getirrep(), WFQ.get_s().getirrep());
        /*double twoS_ = (double) VQTrip[i].get_s().getirrep();
        double twoS  = (double) WFQ.get_s().getirrep();
        double S_    = twoS_ / 2.0;
        double S     = twoS / 2.0;
        TripFac.push_back(sqrt((2.0*S_+1)/3.0*(2.0*S+1)));
        */ 
      }
      //open the integral and operator container
      IKJL.OpenFileRead();
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->OpenFileRead();VtuArray[i]->ResetBuffer();
        SigmaArray[i]->OpenFileRead();SigmaArray[i]->ResetBuffer();
      }
      for (i=0;i<NInternal;i++){
        //for (j=i;j<NInternal;j++){
        //establish batching over second index j that runs from i to NInternal-1
        vector<pair<int,int> > jbatch;
        EstablishBatching(jbatch,MaxCore,MaxM,i,NInternal,VQTrip.size()+1);
        for (int ibatch=0;ibatch<jbatch.size();ibatch++){
          mpi_barrier();
          //generate the vector that holds Vsing,VTrip,SigmaSing and sigmaTrip
          //for the entire batch
          vector<boost::shared_ptr<Wavefunction> > Vij_S;
          vector<boost::shared_ptr<Wavefunction> > SigmaS;
          vector<vector<boost::shared_ptr<Wavefunction> > > Vij_T;
          vector<vector<boost::shared_ptr<Wavefunction> > > SigmaT;
          for (j=jbatch[ibatch].first;j<jbatch[ibatch].second;j++){
            //------------------------------
            //the singlet perturber function
            //------------------------------
            boost::shared_ptr<Wavefunction> VijSing(new Wavefunction(VQSing[0],&big,true));
            boost::shared_ptr<Wavefunction> SigmaSing(new Wavefunction(VQSing[0],&big,true));
            //-------------------------------
            //the triplet perturber functions
            //-------------------------------
            vector<boost::shared_ptr<Wavefunction> > VijTrip;
            vector<boost::shared_ptr<Wavefunction> > SigmaTrip;
            for (a=0;a<VQTrip.size();a++){
              boost::shared_ptr<Wavefunction> VijT(new Wavefunction(VQTrip[a],&big,true));
              VijTrip.push_back(VijT);
              boost::shared_ptr<Wavefunction> sigmaTrip(new Wavefunction(VQTrip[a],&big,true));
              SigmaTrip.push_back(sigmaTrip);
            }
            //----------------
            //get the integral
            //----------------
            boost::shared_ptr<Matrix> Kij = IKJL.GetMatrix(i,j);

            //-------------
            //construct Vij
            //-------------
            /*
            for (t=0;t<NActive;t++){
              for (u=0;u<NActive;u++){
                //get the operators in active space
                boost::shared_ptr<Wavefunction> Ttusing  = T_int1.GetOperator(t,u);
                boost::shared_ptr<Wavefunction> Ttutrip1 = T_int2.GetOperator(t,u);
                boost::shared_ptr<Wavefunction> Ttutrip2 = T_int3.GetOperator(t,u);
                boost::shared_ptr<Wavefunction> Ttutrip3 = T_int4.GetOperator(t,u);

                //the integral (au|bt) times the sign of the operator
                x = Kij->element(u+NInternal,t+NInternal) * Sign.element(u,t);
                y = Kij->element(u+NInternal,t+NInternal);

                //add them to the according perturber functions
                ScaleAdd(y,*Ttusing,*VijSing);
                ScaleAdd(x*TripFac[0],*Ttutrip1,*(VijTrip[0]));
                if (VijTrip.size()>1) ScaleAdd(x*TripFac[1],*Ttutrip2,*(VijTrip[1]));
                if (VijTrip.size()>2) ScaleAdd(x*TripFac[2],*Ttutrip3,*(VijTrip[2]));
              }//u
            }//t
             */
            for (int iquanta=0;iquanta<VtuArray.size();iquanta++){
              bool EndOfArray=false;
              bool EndOfArray_=false;
              while (!EndOfArray){
                boost::shared_ptr<Wavefunction> Vtu  = VtuArray[iquanta]->GetOpFromBuffer(t,u,EndOfArray);
                boost::shared_ptr<Wavefunction> Sigmatu = SigmaArray[iquanta]->GetOpFromBuffer(t,u,EndOfArray_);
                if (EndOfArray!=EndOfArray_){sprintf(msg,"\n!!!Error NEVPT2 in V(i,j): VtuArray and SigmaArray are not synchronized!!!");pout << msg;}
                if (!EndOfArray){
                  //the integral (au|bt) times the sign of the operator
                  if (t!=u){
                    x = Kij->element(u+NInternal,t+NInternal) * Sign.element(u,t) + Kij->element(t+NInternal,u+NInternal) * Sign.element(t,u); 
                    y = Kij->element(u+NInternal,t+NInternal) + Kij->element(t+NInternal,u+NInternal);
                    }
                  else{
                    x = Kij->element(u+NInternal,t+NInternal) * Sign.element(u,t);
                    y = Kij->element(u+NInternal,t+NInternal);
                  }
                  //the singlet function
                  if (iquanta==0){
                    ScaleAdd(y,*Vtu,*VijSing);
                    ScaleAdd(y,*Sigmatu,*SigmaSing);
                  }//singlet
                  //the triplet functions
                  else{
                    int itmp = iquanta-1;
                    ScaleAdd(x*TripFac[itmp],*Vtu,*(VijTrip[itmp]));
                    ScaleAdd(x*TripFac[itmp],*Sigmatu,*(SigmaTrip[itmp]));
                  }//triplet
                }
              }//tu
            }//iquanta
            //gather all functions
            Vij_S.push_back(VijSing);
            SigmaS.push_back(SigmaSing);
            Vij_T.push_back(VijTrip);
            SigmaT.push_back(SigmaTrip);
          }//j
          //-----------------------------------------------------
          //Add up all contributions from the different processes
          //-----------------------------------------------------
          AddPalWavefunctions(Vij_T,SigmaT,Vij_S,SigmaS);
          for (j=jbatch[ibatch].first;j<jbatch[ibatch].second;j++){
            //get the synchronized vectors
            int j0 = j-jbatch[ibatch].first;
            boost::shared_ptr<Wavefunction> VijSing = Vij_S[j0];
            boost::shared_ptr<Wavefunction> SigmaSing = SigmaS[j0];
            vector<boost::shared_ptr<Wavefunction> > VijTrip(Vij_T[j0]);
            vector<boost::shared_ptr<Wavefunction> > SigmaTrip(SigmaT[j0]);
            //---------------------------
            //evaluate the overlap of Vab
            //---------------------------
            //get the integral
            boost::shared_ptr<Matrix> Kij = IKJL.GetMatrix(i,j);
            N = 0.0;
            if (ConventionalOverlap){
              for (t_=0;t_<NActive;t_++){
                for (u_=0;u_<NActive;u_++){
                  x = Kij->element(u_+NInternal,t_+NInternal);
                  for (t=0;t<NActive;t++){
                    for (u=0;u<NActive;u++){
                      y  = Kij->element(u+NInternal,t+NInternal);
                      if (i==j){
                        N += 0.5*x*y*D2_(t_,u_,t,u);
                      }
                      else {
                        N += x*y*D2_(t_,u_,t,u);
                      }
                    }//u
                  }//t
                }//u_
              }//t_
            }//ConventionalOverlap
            else{
              //check if one of the triplet quanta equals the singlet quantum
              int EqualElement = CheckEquality(VQSing[0],VQTrip);
              //if yes, add the corresponding Vab and Sigma functions together
              if (EqualElement!=-1){
                ScaleAdd(1.0,*(VijTrip[EqualElement]),*VijSing);
                ScaleAdd(1.0,*(SigmaTrip[EqualElement]),*SigmaSing);
                //avoid double counting
                VijTrip[EqualElement]->Clear();
                SigmaTrip[EqualElement]->Clear();
              }
              double NSing,NTrip;
              //evaluate the overlap
              NSing = DotProduct(*VijSing,*VijSing);
              NTrip = DotProduct(*VijTrip[0].get(),*VijTrip[0].get());
              if (VijTrip.size()>1) NTrip+= DotProduct(*VijTrip[1].get(),*VijTrip[1].get());
              if (VijTrip.size()>2) NTrip+= DotProduct(*VijTrip[2].get(),*VijTrip[2].get());
              N = NSing + NTrip;
              if (i==j) N *= 0.5;
            }//!ConventionalOverlap
            if (N<_NEV_ZERO) continue;
            //---------------------------------------------------
            //evaluate the active contribution to the denominator
            //---------------------------------------------------
            D = 0.0;
            D += DotProduct(*VijSing,*SigmaSing);
            D += DotProduct(*VijTrip[0].get(),*SigmaTrip[0]); 
            if (VijTrip.size()>1) D += DotProduct(*VijTrip[1].get(),*SigmaTrip[1]);
            if (VijTrip.size()>2) D += DotProduct(*VijTrip[2].get(),*SigmaTrip[2]); 
            if (i==j) D *= 0.5;
            D *= 1/N;
            //-----------------------------------------------------
            //evaluate the external contribution to the denominator
            //-----------------------------------------------------
            D -= EOrb[i+i0] + EOrb[j+i0] + E0;
            //sprintf(msg,"\n%i %i  %4.12lf %4.12lf  %4.12lf",i,j,N,D+EOrb[i+i0]+EOrb[j+i0]+E0,EOrb[i+i0]+EOrb[j+i0]);pout << msg;
            //---------------------------------------
            //calculate the total energy contribution
            //---------------------------------------
            E -= N/D;
          }//j
          //clear memory
          for (int jtmp=0;jtmp<Vij_S.size();jtmp++){
            Vij_S[jtmp]->Clear();
            SigmaS[jtmp]->Clear();
            for (int iquanta=0;iquanta<Vij_T[jtmp].size();iquanta++){
              Vij_T[jtmp][iquanta]->Clear();
              SigmaT[jtmp][iquanta]->Clear();
            }//btmp
          }//iquanta
          Vij_S.clear();
          Vij_T.clear();
          SigmaS.clear();
          SigmaT.clear();
        }//ibatch
      }//i
      //close the container
      IKJL.CloseFileRead();
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->CloseFileRead();
        SigmaArray[i]->CloseFileRead();
      }
      //--------
      //clean up
      //--------
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->Clear();
        SigmaArray[i]->Clear();
      }
      
      return E;
    }
    else{
      //---------------------------------
      // construct the auxiliary K matrix
      //---------------------------------
      array_4d<double> AuxK;
      AuxK.resize(NActive,NActive,NActive,NActive);
      Initialize(AuxK);
      int c,d,e;
      int t__,u__,c_,d_,e_;
      //open the integral container
      IKJL.OpenFileRead();
      for (t_=0;t_<NActive;t_++){
        for (u_=0;u_<NActive;u_++){
          for (t=0;t<NActive;t++){
            for (u=0;u<NActive;u++){
              //the one-electron part
              for (c=0;c<NActive;c++){
                AuxK(t_,u_,t,u) += heff_.element(c+NInternal,t+NInternal) * D2_(t_,u_,c,u);
                AuxK(t_,u_,t,u) += heff_.element(c+NInternal,u+NInternal) * D2_(t_,u_,t,c);
              }//c
            }//u
          }//t
        }//u_
      }//t_
      for (c=0;c<NActive;c++){
        c_ = c + NInternal;
        for (d=0;d<NActive;d++){
          d_ = d + NInternal;
          //get the integral
          boost::shared_ptr<Matrix> Kcd = IKJL.GetMatrix(c_,d_);
          for (t_=0;t_<NActive;t_++){
            for (u_=0;u_<NActive;u_++){
              for (t=0;t<NActive;t++){
                t__ = t + NInternal;
                for (u=0;u<NActive;u++){
                  u__ = u+ NInternal;
                  for (e=0;e<NActive;e++){
                    e_ = e + NInternal;
                    AuxK(t_,u_,t,u) -= Kcd->element(e_,t__) * D3_(t_,u_,e,d,u,c);
                    AuxK(t_,u_,t,u) -= Kcd->element(e_,u__) * D3_(t_,u_,e,t,d,c);
                  }//e
                  AuxK(t_,u_,t,u) += 2.0 * Kcd->element(c_,t__) * D2_(t_,u_,d,u);
                  AuxK(t_,u_,t,u) -= 0.5 * Kcd->element(d_,t__) * D2_(t_,u_,c,u);
                  AuxK(t_,u_,t,u) -= 0.5 * Kcd->element(u__,t__) * D2_(t_,u_,d,c);
                  
                  AuxK(t_,u_,t,u) += 2.0 * Kcd->element(c_,u__) * D2_(t_,u_,t,d);
                  AuxK(t_,u_,t,u) -= 0.5 * Kcd->element(t__,u__) * D2_(t_,u_,c,d);
                  AuxK(t_,u_,t,u) -= 0.5 * Kcd->element(d_,u__) * D2_(t_,u_,t,c);
                }//d
              }//c
            }//u
          }//t
        }//u_
      }//t_
      //---------------------------------
      //Calculate the energy contribution
      //---------------------------------
      for (i=0;i<NInternal;i++){
        for (j=i;j<NInternal;j++){
          //are i and j equal?
          bool ij_equal = (i==j);
          //get the integral
          boost::shared_ptr<Matrix> Kij = IKJL.GetMatrix(i-i0,j-i0);
          //reset the overlap and denominator
          N= 0.0;
          D = 0.0;
          for (t_=0;t_<NActive;t_++){
            for (u_=0;u_<NActive;u_++){
              x = Kij->element(u_+NInternal,t_+NInternal);
              for (t=0;t<NActive;t++){
                for (u=0;u<NActive;u++){
                  x_ = Kij->element(u+NInternal,t+NInternal);
                  if (ij_equal){
                    N += 0.5 * x * x_ * D2_(t_,u_,t,u);
                    D += 0.5 * x * x_ * AuxK(t_,u_,t,u);
                  }//i==j
                  else{
                    N += x * x_ * D2_(t_,u_,t,u);
                    D += x * x_ * AuxK(t_,u_,t,u);
                  }//i!=j
                }//u
              }//t
            }//u_
          }//t_
          if (N<_NEV_ZERO) continue;
          D *= 1/N;
          D -= EOrb[i+i0] + EOrb[j+i0];
          E -= N/D;
        }//j
      }//i
      //close the integral container
      IKJL.CloseFileRead();
      return E;
    }
  }
  
  //============================================================================
  // Evaluate the energy contribution from class V(0,ia)
  //============================================================================
  double V_ia(int *OrbWin, const vector<double> &EOrb, Wavefunction &WF, const SpinBlock &big,
              IntegralContainer &IKJA, IntegralContainer &IKJL, IntegralContainer &IJKA,
              boost::shared_ptr<IntegralContainer> KIAJ, const Matrix &D1, 
              const array_4d<double> &DC2, const array_6d &DC3, const Matrix &heff,
              const Matrix &heff_, bool Conventional,char *BaseName, double E0, 
              bool ConventionalOverlap,int MaxCore){
    
    char msg[512];
    int t,u;
    int i0,i1,t0,t1,a0,a1;
    int a,b;
    int tu;
    int i,j;
    double x,x_,y,y_,z,z_;
    double N=0.0;
    double D=0.0;
    double E=0.0;
    int t_,u_;

    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    
    //define the orbital spaces
    i0 = OrbWin[0];
    i1 = OrbWin[1];
    t0 = OrbWin[2];
    t1 = OrbWin[3];
    a0 = OrbWin[4];
    a1 = OrbWin[5];
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    
    if (!Conventional){
      
      //-----------------------------------------------
      //get the quanta of the ground state Wavefunction
      //-----------------------------------------------
      SpinQuantum WFQ = WF.get_deltaQuantum(0);
      IrrepSpace Dummy(0);
      SpinQuantum TripOpQ(0,SpinSpace(2),Dummy);
      SpinQuantum SingOpQ(0,SpinSpace(0),Dummy);
      vector<SpinQuantum> VQTrip = WFQ + TripOpQ;
      vector<SpinQuantum> VQSing = WFQ + SingOpQ;
      
      //-----------------------------
      //Generate the T(t,u) container
      //-----------------------------
      //the container that hold the V(t,u) operators
      vector<boost::shared_ptr<WavefunctionArray> > VtuArray;
      boost::shared_ptr<WavefunctionArray> VtuSing(new WavefunctionArray);
      sprintf(msg, "%s%s", dmrginp.save_prefix().c_str(), "/VtuSing.tmp");
      VtuSing->Initialize(NActive,msg,mpi_world_size());
      VtuArray.push_back(VtuSing);
      for (int iquanta=0;iquanta<VQTrip.size();iquanta++){
        boost::shared_ptr<WavefunctionArray> VtuTrip(new WavefunctionArray);
        sprintf(msg, "%s%s.%i.tmp", dmrginp.save_prefix().c_str(), "/VtuTrip",iquanta);
        VtuTrip->Initialize(NActive,msg,mpi_world_size());
        VtuArray.push_back(VtuTrip);
      }
      //evaluate the size of the buffers
      int MaxM = dmrginp.sweep_state_schedule()[dmrginp.sweep_qstate_schedule().size()];
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->CalcBufferSize(MaxCore,MaxM);
      }
      
      //--------------------------------------------------------------------------
      // generate V(t,u) = a+(t)a(u)|psi> 
      //--------------------------------------------------------------------------
      //open the container
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->OpenFileWrite();
      }
      
      //--------------------------------------------------------
      //case 1: both indices are on the left side of the lattice
      //--------------------------------------------------------
      for (tu=0;tu<leftBlock->get_op_array(CRE_DES).get_size();tu++){
        //------------------------
        //get the triplet operator
        //------------------------
        boost::shared_ptr<SparseMatrix> TripOp = leftBlock->get_op_array(CRE_DES).get_local_element(tu)[1]->getworkingrepresentation(leftBlock);
        //get the orbital indices
        t = TripOp->get_orbs(0);
        u = TripOp->get_orbs(1);
        //----------------------------
        //generate the possible Quanta
        //----------------------------
        //get the Quanta of the Wavefunction
        SpinQuantum WFQ=WF.get_deltaQuantum(0);
        //the operator quanta
        SpinQuantum opQ=TripOp->get_deltaQuantum(0);
        //the vector that holds all possible combinations of wavefunction plus operator
        std::vector<SpinQuantum> VQ  = opQ  + WFQ;
        //-------------------------------------------
        //Multiply the operator with the wavefunction
        //-------------------------------------------
        for (int iquanta=0;iquanta<VQ.size();iquanta++){
          //generate the wavefunction that holds the result of the multiplication
          boost::shared_ptr<Wavefunction> Vtu(new Wavefunction(VQ[iquanta],&big,true));
          boost::shared_ptr<Wavefunction> VtuT(new Wavefunction(VQ[iquanta],&big,true));
          //multiply the wavefunction with the triplet operator
          operatorfunctions::TensorMultiply(leftBlock,*TripOp,&big,WF,*Vtu,opQ,1.0);
          if (t!=u){
            operatorfunctions::TensorMultiply(leftBlock,Transposeview(*TripOp),&big,WF,*VtuT,opQ,1.0);
          }//t!=u
          //-------------------
          //store the operators
          //-------------------
          VtuArray[iquanta+1]->AppendOperator(Vtu,t,u);
          if (t!=u) VtuArray[iquanta+1]->AppendOperator(VtuT,u,t);
        }//iquanta
        VQ.clear();
        //------------------------
        //get the singlet operator
        //------------------------
        boost::shared_ptr<SparseMatrix> SingOp = leftBlock->get_op_array(CRE_DES).get_local_element(tu)[0]->getworkingrepresentation(leftBlock);
        //get the orbital indices
        t = SingOp->get_orbs(0);
        u = SingOp->get_orbs(1);
        //the operator quanta
        opQ = SingOp->get_deltaQuantum(0); 
        VQ = opQ + WFQ;
        //generate the wavefunction that holds the result of the multiplication
        boost::shared_ptr<Wavefunction> VtuSing(new Wavefunction(VQ[0],&big,true));
        boost::shared_ptr<Wavefunction> VtuSingT(new Wavefunction(VQ[0],&big,true));
        //multiply the wavefunction with the singlet operator
        operatorfunctions::TensorMultiply(leftBlock,*SingOp,&big,WF,*VtuSing,opQ,1.0);
        if (t!=u){
          operatorfunctions::TensorMultiply(leftBlock,Transposeview(*SingOp),&big,WF,*VtuSingT,opQ,1.0);
        }
        //------------------
        //store the operator
        //------------------
        VtuArray[0]->AppendOperator(VtuSing,t,u);
        if (t!=u) VtuArray[0]->AppendOperator(VtuSingT,u,t);
      }//tu
      //---------------------------------------------------------
      //case 2: both indices are on the right side of the lattice
      //---------------------------------------------------------
      for (tu=0;tu<rightBlock->get_op_array(CRE_DES).get_size();tu++){
        //------------------------
        //get the triplet operator
        //------------------------
        boost::shared_ptr<SparseMatrix> TripOp = rightBlock->get_op_array(CRE_DES).get_local_element(tu)[1]->getworkingrepresentation(rightBlock);
        //get the orbital indices
        t = TripOp->get_orbs(0);
        u = TripOp->get_orbs(1);
        //----------------------------
        //generate the possible Quanta
        //----------------------------
        //get the Quanta of the Wavefunction
        SpinQuantum WFQ=WF.get_deltaQuantum(0);
        //the operator quanta
        SpinQuantum opQ=TripOp->get_deltaQuantum(0);
        //the vector that holds all possible combinations of wavefunction plus operator
        std::vector<SpinQuantum> VQ  = opQ  + WFQ;
        //-------------------------------------------
        //Multiply the operator with the wavefunction
        //-------------------------------------------
        for (int iquanta=0;iquanta<VQ.size();iquanta++){
          //generate the wavefunction that holds the result of the multiplication
          boost::shared_ptr<Wavefunction> Vtu(new Wavefunction(VQ[iquanta],&big,true));
          boost::shared_ptr<Wavefunction> VtuT(new Wavefunction(VQ[iquanta],&big,true));
          //multiply the wavefunction with the triplet operator
          operatorfunctions::TensorMultiply(rightBlock,*TripOp,&big,WF,*Vtu,opQ,1.0);
          if (t!=u){
            operatorfunctions::TensorMultiply(rightBlock,Transposeview(*TripOp),&big,WF,*VtuT,opQ,1.0);
          }//t!=u
          //-------------------
          //store the operators
          //-------------------
          VtuArray[iquanta+1]->AppendOperator(Vtu,t,u);
          if (t!=u) VtuArray[iquanta+1]->AppendOperator(VtuT,u,t);
        }//iquanta
        VQ.clear();
        //------------------------
        //get the singlet operator
        //------------------------
        boost::shared_ptr<SparseMatrix> SingOp = rightBlock->get_op_array(CRE_DES).get_local_element(tu)[0]->getworkingrepresentation(rightBlock);
        //get the orbital indices
        t = SingOp->get_orbs(0);
        u = SingOp->get_orbs(1);
        //the operator quanta
        opQ = SingOp->get_deltaQuantum(0); 
        VQ = opQ + WFQ;
        //generate the wavefunction that holds the result of the multiplication
        boost::shared_ptr<Wavefunction> VtuSing(new Wavefunction(VQ[0],&big,true));
        boost::shared_ptr<Wavefunction> VtuSingT(new Wavefunction(VQ[0],&big,true));
        //multiply the wavefunction with the singlet operator
        operatorfunctions::TensorMultiply(rightBlock,*SingOp,&big,WF,*VtuSing,opQ,1.0);
        if (t!=u){
          operatorfunctions::TensorMultiply(rightBlock,Transposeview(*SingOp),&big,WF,*VtuSingT,opQ,1.0);
        }
        //------------------
        //store the operator
        //------------------
        VtuArray[0]->AppendOperator(VtuSing,t,u);
        if (t!=u) VtuArray[0]->AppendOperator(VtuSingT,u,t);
      }//tu
      //--------------------------------------------------------
      //case 3: the indices are on the both sides of the lattice
      //--------------------------------------------------------
      int tstart,tstop;
      //divide the loop over the parallel processes
      PALDivideLoop(tstart,tstop,0,leftBlock->get_op_array(CRE).get_size());
      //loop over local elements on the left side
      for (int tcount=tstart;tcount<tstop;tcount++){
      //for (int tcount=0;tcount<leftBlock->get_op_array(CRE).size();tcount++){
        //get the operator
        boost::shared_ptr<SparseMatrix> lop = leftBlock->get_op_array(CRE).get_local_element(tcount)[0]->getworkingrepresentation(leftBlock);
        //get the orbital index
        t = lop->get_orbs(0);
        for (int ucount=0;ucount<rightBlock->get_op_array(CRE).size();ucount++){
          //get the operator
          boost::shared_ptr<SparseMatrix> rop = rightBlock->get_op_array(CRE).get_local_element(ucount)[0]->getworkingrepresentation(rightBlock);
          //get the orbital index
          u = rop->get_orbs(0);
          //get the Quanta of the Wavefunction
          SpinQuantum WFQ=WF.get_deltaQuantum(0);
          //the possible operator quanta -singlet and triplet-
          std::vector<SpinQuantum> opQu = lop->get_deltaQuantum(0) + Transposeview(*rop).get_deltaQuantum(0);
          //----------------
          //triplet operator
          //----------------
          SpinQuantum opQ = opQu[1];
          //the vector that holds all possible combinations of wavefunction plus operator
          std::vector<SpinQuantum> VQ = opQ + WFQ;
          for (int iquanta=0;iquanta<VQ.size();iquanta++){
            //generate the wavefunction that holds the result of the multiplication
            boost::shared_ptr<Wavefunction> Vtu(new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> VtuT(new Wavefunction(VQ[iquanta],&big,true));
            //multiply the wavefunction with the triplet operator
            operatorfunctions::TensorMultiply(leftBlock,*lop,Transposeview(*rop),&big,WF,*Vtu,opQ,1.0);
            operatorfunctions::TensorMultiply(leftBlock,Transposeview(*lop),*rop,&big,WF,*VtuT,opQ,-1.0);
            //------------------
            //store the operator
            //------------------
            VtuArray[iquanta+1]->AppendOperator(Vtu,t,u);
            VtuArray[iquanta+1]->AppendOperator(VtuT,u,t);
          }//iquanta
          VQ.clear();
          //----------------
          //singlet operator
          //----------------
          //the operator quanta
          opQ = opQu[0];
          VQ = opQ + WFQ;
          //generate the wavefunction that holds the result of the multiplication
          boost::shared_ptr<Wavefunction> VtuSing(new Wavefunction(VQ[0],&big,true));
          boost::shared_ptr<Wavefunction> VtuSingT(new Wavefunction(VQ[0],&big,true));
          //multiply the wavefunction with the triplet operator
          operatorfunctions::TensorMultiply(leftBlock,*lop,Transposeview(*rop),&big,WF,*VtuSing,opQ,1.0);
          operatorfunctions::TensorMultiply(leftBlock,Transposeview(*lop),*rop,&big,WF,*VtuSingT,opQ,1.0);
          //------------------
          //store the operator
          //------------------
          VtuArray[0]->AppendOperator(VtuSing,t,u);
          VtuArray[0]->AppendOperator(VtuSingT,u,t);
        }//u
      }//t
      //close the container
      //open the container
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->CloseFileWrite();
      }
            
      //------------------------------------------------------------------------
      // Construct the action of the Hamiltonian on the V(t,u) functions
      // sigma(t,u) = H * V(t,u)
      //------------------------------------------------------------------------
      //the container that hold the sigma(t,u) functions
      vector<boost::shared_ptr<WavefunctionArray> > SigmaArray;
      boost::shared_ptr<WavefunctionArray> SigmaSing(new WavefunctionArray);
      sprintf(msg, "%s%s", dmrginp.save_prefix().c_str(), "/VtuSing.tmp");
      SigmaSing->Initialize(NActive,msg,mpi_world_size());
      SigmaSing->CalcBufferSize(MaxCore,MaxM);
      SigmaArray.push_back(SigmaSing);
      for (int iquanta=0;iquanta<VQTrip.size();iquanta++){
        boost::shared_ptr<WavefunctionArray> SigmaTrip(new WavefunctionArray);
        sprintf(msg, "%s%s.%i.tmp", dmrginp.save_prefix().c_str(), "/SigmatuTrip",iquanta);
        SigmaTrip->Initialize(NActive,msg,mpi_world_size());
        SigmaTrip->CalcBufferSize(MaxCore,MaxM);
        SigmaArray.push_back(SigmaTrip);
      }
      //generate the sigma(t,u) functions
      GenerateActionOfH(WF,VtuArray,SigmaArray,TripOpQ,SingOpQ,big,OrbWin);
      //generate the regular sigma function
      boost::shared_ptr<Wavefunction> RegSigma (new Wavefunction (WFQ,&big,true));
      big.multiplyH(WF,RegSigma.get(),MAX_THRD);
      //------------------------------------------------------------------------
      // construct Vij = gamma(ia){fac * (iu|jt) * Ttu}
      // and calculate the energy contribution
      //------------------------------------------------------------------------
      double tmp=0.0;
      double h=0.0;
      //---------------
      // the prefactors
      //---------------
      vector<double> TripFac;
      for (i=0;i<VQTrip.size();i++){
        EvalCompFactors(TripFac, 2, VQTrip[i].get_s().getirrep(), WFQ.get_s().getirrep());
      }
      //open the integral and operator container
      sprintf(msg,"%s.block.KIAJ.tmp",BaseName);
      boost::shared_ptr<IntegralContainer> KIAJ = IKJA.ReverseOrder(msg);
      sprintf(msg,"%s.block.KAIJ.tmp",BaseName);
      boost::shared_ptr<IntegralContainer> KAIJ = IJKA.ReverseOrder(msg);
      KIAJ->OpenFileRead();
      KAIJ->OpenFileRead();
      //open the container
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->OpenFileRead();VtuArray[i]->ResetBuffer();
        SigmaArray[i]->OpenFileRead();SigmaArray[i]->ResetBuffer();
      }
      for (i=0;i<NInternal;i++){
        int i_ = i+i0;
        //for (a=0;a<NExternal;a++){
        //establish batching over second index j that runs from 0 to NExternal-1
        vector<pair<int,int> > abatch;
        EstablishBatching(abatch,MaxCore,MaxM,0,NExternal,VQTrip.size()+1);
        for (int ibatch=0;ibatch<abatch.size();ibatch++){
          mpi_barrier();
          //generate the vector that holds Vsing,VTrip,SigmaSing and sigmaTrip
          //for the entire batch
          vector<boost::shared_ptr<Wavefunction> > Via_S;
          vector<boost::shared_ptr<Wavefunction> > SigmaS;
          vector<vector<boost::shared_ptr<Wavefunction> > > Via_T;
          vector<vector<boost::shared_ptr<Wavefunction> > > SigmaT;
          for (a=abatch[ibatch].first;a<abatch[ibatch].second;a++){
            int a_ = a+a0;
            h = heff.element(a_,i_);

            //------------------------------
            //the singlet perturber function
            //------------------------------
            boost::shared_ptr<Wavefunction> ViaSing(new Wavefunction(VQSing[0],&big,true));
            boost::shared_ptr<Wavefunction> SigmaSing(new Wavefunction(VQSing[0],&big,true));
            //-------------------------------
            //the triplet perturber functions
            //-------------------------------
            vector<boost::shared_ptr<Wavefunction> > ViaTrip;
            vector<boost::shared_ptr<Wavefunction> > SigmaTrip;
            for (j=0;j<VQTrip.size();j++){
              boost::shared_ptr<Wavefunction> ViaT(new Wavefunction(VQTrip[j],&big,true));
              ViaTrip.push_back(ViaT);
              boost::shared_ptr<Wavefunction> SigmaiaT(new Wavefunction(VQTrip[j],&big,true));
              SigmaTrip.push_back(SigmaiaT);
            }
            //----------------
            //get the integral
            //----------------
            boost::shared_ptr<Matrix> Kia = KIAJ->GetMatrix(i,a);
            boost::shared_ptr<Matrix> Jia = KAIJ->GetMatrix(i,a);
            //-------------
            //construct Via
            //-------------
            /*
            for (t=0;t<NActive;t++){
              for (u=0;u<NActive;u++){
                //get the operators in active space
                boost::shared_ptr<Wavefunction> Ttusing  = T_int1.GetOperator(t,u);
                boost::shared_ptr<Wavefunction> Ttutrip1 = T_int2.GetOperator(t,u);
                boost::shared_ptr<Wavefunction> Ttutrip2 = T_int3.GetOperator(t,u);
                boost::shared_ptr<Wavefunction> Ttutrip3 = T_int4.GetOperator(t,u);

                //the integral (it|au)
                x = Kia->element(t+NInternal,u+NInternal);
                y = Jia->element(t+NInternal,u+NInternal);

                //add them to the according perturber functions
                ScaleAdd(2.0*y-x,*Ttusing,*ViaSing);
                ScaleAdd(x*TripFac[0],*Ttutrip1,*(ViaTrip[0]));
                if (ViaTrip.size()>1) ScaleAdd(x*TripFac[1],*Ttutrip2,*(ViaTrip[1]));
                if (ViaTrip.size()>2) ScaleAdd(x*TripFac[2],*Ttutrip3,*(ViaTrip[2]));
              }//u
            }//t
            //add the one-particle term
            ScaleAdd(sqrt(2.0)*h,WF,*ViaSing);
             */
            for (int iquanta=0;iquanta<VtuArray.size();iquanta++){
              bool EndOfArray=false;
              bool EndOfArray_=false;
              while (!EndOfArray){
                boost::shared_ptr<Wavefunction> Vtu = VtuArray[iquanta]->GetOpFromBuffer(t,u,EndOfArray);
                boost::shared_ptr<Wavefunction> Sigmatu = SigmaArray[iquanta]->GetOpFromBuffer(t,u,EndOfArray_);
                if (EndOfArray!=EndOfArray_){sprintf(msg,"\n!!!Error NEVPT2 in V(i,j): VtuArray and SigmaArray are not synchronized!!!");pout << msg;}
                if (!EndOfArray){
                  //the integral (au|bt) times the sign of the operator
                  x = Kia->element(t+NInternal,u+NInternal);
                  y = Jia->element(t+NInternal,u+NInternal);
                  //singlet functions
                  if (iquanta==0){
                    ScaleAdd(2.0*y-x,*Vtu,*ViaSing);
                    ScaleAdd(2.0*y-x,*Sigmatu,*SigmaSing);
                  }//singlet
                  //triplet functions
                  else{
                    int itmp = iquanta-1;
                    ScaleAdd(x*TripFac[itmp],*Vtu,*(ViaTrip[itmp]));
                    ScaleAdd(x*TripFac[itmp],*Sigmatu,*(SigmaTrip[itmp]));
                  }//triplet
                }
              }//tu
            }//iquanta
            //add the one-particle term
            if (mpi_rank()==0) ScaleAdd(sqrt(2.0)*h,WF,*ViaSing);
            if (mpi_rank()==0) ScaleAdd(sqrt(2.0)*h,*RegSigma,*SigmaSing);
            //gather all functions
            Via_S.push_back(ViaSing);
            SigmaS.push_back(SigmaSing);
            Via_T.push_back(ViaTrip);
            SigmaT.push_back(SigmaTrip);
          }//a
          //-----------------------------------------------------
          //Add up all contributions from the different processes
          //-----------------------------------------------------
          AddPalWavefunctions(Via_T,SigmaT,Via_S,SigmaS);
          for (a=abatch[ibatch].first;a<abatch[ibatch].second;a++){
            //get the synchronized vectors
            int a0_ = a-abatch[ibatch].first;
            boost::shared_ptr<Wavefunction> ViaSing = Via_S[a0_];
            boost::shared_ptr<Wavefunction> SigmaSing = SigmaS[a0_];
            vector<boost::shared_ptr<Wavefunction> > ViaTrip(Via_T[a0_]);
            vector<boost::shared_ptr<Wavefunction> > SigmaTrip(SigmaT[a0_]);
            //---------------------------
            //evaluate the overlap of Vab
            //---------------------------
            N = 0.0;
            if (ConventionalOverlap){
              //get the integrals
              int a_ = a+a0;
              h = heff.element(a_,i_);
              boost::shared_ptr<Matrix> Kia = KIAJ->GetMatrix(i,a);
              boost::shared_ptr<Matrix> Jia = KAIJ->GetMatrix(i,a);
              for (t_=0;t_<NActive;t_++){
                int at_ = t_ + NInternal;
                for (u_=0;u_<NActive;u_++){
                  int au_ = u_ + NInternal;
                  x_ = Jia->element(at_,au_);
                  y_ = Kia->element(at_,au_);
                  for (t=0;t<NActive;t++){
                    int at = t + NInternal;
                    for (u=0;u<NActive;u++){
                      int au = u + NInternal;
                      x = Jia->element(at,au);
                      y = Kia->element(at,au);
                      tmp = 2.0*x_*x - x_*y - y_*x;
                      N += tmp * DC2(u_,t,t_,u); 
                      tmp = y_*y;
                      N -= tmp * DC2(u_,t,u,t_);
                    }//u
                    z_ = Kia->element(at_,au_);
                    z = Kia->element(at_,at);
                    y = Kia->element(at,at);
                    N += 2.0 * z_*z * D1.element(u_,t);
                    N += y_*y * D1.element(u_,t_);
                  }//t
                  tmp = 4.0 * x_ - 2.0 * y_;
                  N += tmp * h * D1.element(u_,t_);
                }//u_
              }//t_
              N += 2.0 * h * h;
            }//ConventionalOverlap
            else{
              //check if one of the triplet quanta equals the singlet quantum
              int EqualElement = CheckEquality(VQSing[0],VQTrip);
              //if yes, add the corresponding Vab and Sigma functions together
              if (EqualElement!=-1){
                ScaleAdd(1.0,*(ViaTrip[EqualElement]),*ViaSing);
                ScaleAdd(1.0,*(SigmaTrip[EqualElement]),*SigmaSing);
                //avoid double counting
                ViaTrip[EqualElement]->Clear();
                SigmaTrip[EqualElement]->Clear();
              }
              double NSing,NTrip;
              //evaluate the overlap
              NSing = DotProduct(*ViaSing,*ViaSing);
              NTrip = DotProduct(*ViaTrip[0].get(),*ViaTrip[0].get());
              if (ViaTrip.size()>1) NTrip = DotProduct(*ViaTrip[1].get(),*ViaTrip[1].get());
              if (ViaTrip.size()>2) NTrip = DotProduct(*ViaTrip[2].get(),*ViaTrip[2].get());
              N = NSing + NTrip;
            }//!ConventionalOverlap
            //sprintf(msg,"\nmy i=%i a=%i N(conv)=%4.6lf N=%4.6lf NSing=%4.6lf NTrip=%4.6lf NH=%4.6lf quot=%4.6lf",i,a,N,Nij,NSing,NTrip,NH,(N-Nij)/NTrip);pout << msg;
            if (N<_NEV_ZERO) continue;
            //---------------------------------------------------
            //evaluate the active contribution to the denominator
            //---------------------------------------------------
            D = 0.0;
            D += DotProduct(*ViaSing,*SigmaSing);
            D += DotProduct(*ViaTrip[0].get(),*SigmaTrip[0]); 
            if (ViaTrip.size()>1) D += DotProduct(*ViaTrip[1].get(),*SigmaTrip[1]);
            if (ViaTrip.size()>2) D += DotProduct(*ViaTrip[2].get(),*SigmaTrip[2]); 
            D *= 1/N;
            //-----------------------------------------------------
            //evaluate the external contribution to the denominator
            //-----------------------------------------------------
            D += EOrb[a+a0] - EOrb[i+i0] - E0;
            //sprintf(msg,"\n%i %i  %4.12lf %4.12lf  %4.12lf  %4.12lf",i,a,N,D+EOrb[i+i0]-EOrb[a+a0]+E0,-EOrb[i+i0]+EOrb[j+i0],-N/D);pout << msg;
            //---------------------------------------
            //calculate the total energy contribution
            //---------------------------------------
            E -= N/D;
          }//a
          //clear memory
          for (int btmp=0;btmp<Via_S.size();btmp++){
            Via_S[btmp]->Clear();
            SigmaS[btmp]->Clear();
            for (int iquanta=0;iquanta<Via_T[btmp].size();iquanta++){
              Via_T[btmp][iquanta]->Clear();
              SigmaT[btmp][iquanta]->Clear();
            }//btmp
          }//iquanta
          Via_S.clear();
          Via_T.clear();
          SigmaS.clear();
          SigmaT.clear();
        }//ibatch
      }//i
      
      //close the container
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->CloseFileRead();
        SigmaArray[i]->CloseFileRead();
      }
      //--------
      //clean up
      //--------
      for (i=0;i<VtuArray.size();i++){
        VtuArray[i]->Clear();
        SigmaArray[i]->Clear();
      }
      return E;
    }
    else{
      int i_,a_,c,d,e;
      double z,z_;
      double tmp;
      double h;
      int at,at_,au,au_,ac,ad,ae;
      //generate the auxiliary matrices A and D
      array_4d<double> AuxA;
      array_4d<double> AuxD;
      AuxA.resize(NActive,NActive,NActive,NActive);
      AuxD.resize(NActive,NActive,NActive,NActive);
      Initialize(AuxA);
      Initialize(AuxD);

      //-------------------------------
      //Generate the auxiliary A matrix
      //-------------------------------
      /*//start with the auxiliary h_eff_ matrix
      Matrix heff_(heff);
      for (t=0;t<NActive;t++){
        for (u=0;u<NActive;u++){
          //copy the old heff matrix
          heff_.element(t,u) = heff.element(t,u);
          
        }//u
      }//t
      //open the integral container
      IKJL.OpenFileRead();
      for (c=0;c<NActive;c++){
        //get the K integral
        boost::shared_ptr<Matrix> Kcc = IKJL.GetMatrix(c+NInternal,c+NInternal);
        for (t=0;t<NActive;t++){
          t_ = t + NInternal;
          for (u=0;u<NActive;u++){
            u_ = u + NInternal;
            heff_.element(t,u) -= 0.5 * Kcc->element(t_,u_);
          }//u
        }//t
      }//c
      //close the integral container
      IKJL.CloseFileRead();
      PrintMatrix(heff_,"heff__.tmp");*/
      
      //open the integral container
      IKJL.OpenFileRead();
      for (c=0;c<NActive;c++){
        ac = c + NInternal;
        for (t=0;t<NActive;t++){
          at = t + NInternal;
          for (u=0;u<NActive;u++){
            au = u + NInternal;
            //get the integrals
            boost::shared_ptr<Matrix> Ktc = IKJL.GetMatrix(t+t0,c+t0);
            boost::shared_ptr<Matrix> Kcu = IKJL.GetMatrix(c+t0,u+t0);
            for (t_=0;t_<NActive;t_++){
              at_ = t_ + NInternal;
              for (u_=0;u_<NActive;u_++){
                au_ = u_ + NInternal;
                //the one-electron part
                AuxA(t_,u_,t,u) += heff_.element(c+t0,t+t0)*DC2(u_,c,t_,u);
                AuxA(t_,u_,t,u) -= heff_.element(u+t0,c+t0)*DC2(u_,t,t_,c);
                for (d=0;d<NActive;d++){
                  ad = d+NInternal;
                  for (e=0;e<NActive;e++){
                    ae = e + NInternal;
                    //the two-electron part
                    AuxA(t_,u_,t,u) += 0.5*Ktc->element(ad,ae) * DC3(u_,c,d,t_,e,u);
                    AuxA(t_,u_,t,u) += 0.5*Ktc->element(ad,ae) * DC3(u_,d,c,t_,u,e);
                    AuxA(t_,u_,t,u) -= 0.5*Kcu->element(ad,ae) * DC3(u_,t,c,t_,e,d);
                    AuxA(t_,u_,t,u) -= 0.5*Kcu->element(ad,ae) * DC3(u_,c,t,t_,d,e);
                  }//e
                }//d
              }//c
            }//u
          }//t
        }//u_
      }//u_
      //-------------------------------
      //Generate the auxiliary D matrix
      //-------------------------------
      for (c=0;c<NActive;c++){
        ac = c + NInternal;
        for (t=0;t<NActive;t++){
          at = t + NInternal;
          for (u=0;u<NActive;u++){
            au = u + NInternal;
            //get the integrals
            boost::shared_ptr<Matrix> Ktc = IKJL.GetMatrix(t+t0,c+t0);
            boost::shared_ptr<Matrix> Kuc = IKJL.GetMatrix(u+t0,c+t0);
            for (t_=0;t_<NActive;t_++){
              at_ = t_ + NInternal;
              for (u_=0;u_<NActive;u_++){
                au_ = u_ + NInternal;
                //the one-electron part
                AuxD(t_,u_,t,u) +=       heff_.element(u+t0,c+t0) *  DC2(t,u_,t_,c);
                AuxD(t_,u_,t,u) -=       heff_.element(c+t0,t+t0) *  DC2(c,u_,t_,u);
                
                for (d=0;d<NActive;d++){
                  ad = d+NInternal;
                  for (e=0;e<NActive;e++){
                    ae = e + NInternal;
                    //get the integral values
                    x = -0.5 * Kuc->element(ad,ae);
                    y =  0.5 * Ktc->element(ad,ae);
                    AuxD(t_,u_,t,u) -= x * DC3(c,t,u_,e,t_,d);
                    AuxD(t_,u_,t,u) -= x * DC3(t,u_,c,t_,d,e);
                    AuxD(t_,u_,t,u) -= y * DC3(c,d,u_,e,t_,u);
                    AuxD(t_,u_,t,u) -= y * DC3(d,u_,c,t_,u,e);
                  }//e
                  x = -0.5 * Kuc->element(ad,at_);
                  y =  0.5 * Ktc->element(at_,ad);
                  z =  0.5 * Ktc->element(ad,at_);
                  
                  AuxD(t,u_,t,u) += 2.0 * x * DC2(c,u_,t_,d);
                  AuxD(t,u_,t,u) += 2.0 * x * DC2(u_,c,d,t_);
                  AuxD(t_,u_,t,u)+= 2.0 * y * DC2(c,u_,d,u);
                  AuxD(t_,u_,t,u)+= 2.0 * y * DC2(u_,c,u,d);
                  AuxD(u_,u_,t,u)+=       x * DC2(c,t,t_,d);
                  AuxD(u_,u_,t,u)+=       x * DC2(t,c,d,t_);
                  AuxD(u_,u_,t,u)+=       z * DC2(c,d,t_,u);
                  AuxD(u_,u_,t,u)+=       z * DC2(d,c,u,t_);
                  
                  AuxD(c,u_,t,u) -=       x * DC2(t,u_,t_,d);
                  if (t_==t) AuxD(c,u_,t,u) += 2.0 * x * D1.element(u_,d);
                  AuxD(u_,t_,t,u)+=       x * DC2(t,c,u_,d);
                  if (t==u_) AuxD(u_,t_,t,u) -= 2.0 * x * D1.element(c,d);
                  AuxD(c,u_,t,u) -=       z * DC2(d,u_,t_,u);
                  if (d==t_) AuxD(c,u_,t,u) += 2.0 * z * D1.element(u_,u);
                  AuxD(u_,t_,t,u)+=       z *  DC2(d,c,u_,u);
                  if (d==u_) AuxD(u_,t_,t,u)-= 2.0 * z *  D1.element(c,u);
                }//d
              }//u_
              //the one-electron part (cont'd)
              AuxD(t,t_,t,u) -= 2.0 * heff_.element(u+t0,c+t0) *  D1.element(t_,c);
              AuxD(t_,t_,t,u)-=       heff_.element(u+t0,c+t0) *  D1.element(c,t);
              AuxD(c,t_,t,u) += 2.0 * heff_.element(c+t0,t+t0) *  D1.element(t_,u);
              AuxD(t_,t_,t,u)+=       heff_.element(c+t0,t+t0) *  D1.element(u,c);
            }//t_
          }//u
        }//t
      }//c
      //close the integral container
      IKJL.CloseFileRead();
      
      //---------------------------------
      //Calculate the energy contribution
      //---------------------------------
      //invert the order of the half external integrals
      sprintf(msg,"%s.block.KAIJ.tmp",BaseName);
      boost::shared_ptr<IntegralContainer> KAIJ = IJKA.ReverseOrder(msg);
      //open the integral container
      KAIJ->OpenFileRead();
      KIAJ->OpenFileRead();
      for (i=0;i<NInternal;i++){
        i_ = i + i0;
        for (a=0;a<NExternal;a++){
          a_ = a+a0;
          h = heff.element(a_,i_);
          //get the integral
          boost::shared_ptr<Matrix> Jia = KAIJ->GetMatrix(i,a);
          boost::shared_ptr<Matrix> Kia = KIAJ->GetMatrix(i,a);
          //reset the overlap and denominator
          N= 0.0;
          D = 0.0;
          for (t_=0;t_<NActive;t_++){
            at_ = t_ + NInternal;
            for (u_=0;u_<NActive;u_++){
              au_ = u_ + NInternal;
              x_ = Jia->element(at_,au_);
              y_ = Kia->element(at_,au_);
              for (t=0;t<NActive;t++){
                at = t + NInternal;
                for (u=0;u<NActive;u++){
                  au = u + NInternal;
                  x = Jia->element(at,au);
                  y = Kia->element(at,au);
                  tmp = 2.0*x_*x - x_*y - y_*x;
                  N += tmp * DC2(u_,t,t_,u); 
                  D += tmp * AuxA(t_,u_,t,u);
                  tmp = y_*y;
                  N -= tmp * DC2(u_,t,u,t_);
                  D += tmp * AuxD(t_,u_,t,u);
                }//u
                z_ = Kia->element(at_,au_);
                z = Kia->element(at_,at);
                y = Kia->element(at,at);
                N += 2.0 * z_*z * D1.element(u_,t);
                N += y_*y * D1.element(u_,t_);
              }//t
              tmp = 4.0 * x_ - 2.0 * y_;
              N += tmp * h * D1.element(u_,t_);
            }//u_
          }//t_
          N += 2.0 * h * h;
          if (N<_NEV_ZERO) continue;
          D *= 1/N;
          D += EOrb[a+a0] - EOrb[i+i0];
          E -= N/D;
          
        }//a
      }//i
      //close the integral container
      KAIJ->CloseFileRead();
      KIAJ->CloseFileRead();
      //clean up disk
      KAIJ->Clear();
    }
    
    return E;
  }
  
  
  //============================================================================
  // The main driver for BLOCK-NEVPT2 calculations
  // INPUT      big     the spinblock of the whole lattice
  //            WF      the vector wavefunctions for which you want to evaluate
  //                    the NEVPT2 energy
  // also you need to supply an input file with name "dmrg.nevpt2.inp" that 
  // specifies (in that order) the name of your calculation, first and last
  // indices for internal, active and external orbital spaces, the static energy
  // of the nuclei and if the Overlap should be calculated conventionally or not. 
  //============================================================================
  void NEVPT2_Driver(SpinBlock &big, std::vector<Wavefunction> &WF, NEVPT2Info &Info){
    
    char msg[512];
    char BaseName[512];
    double E,E0,ENuc,EEl,Eval;
    int nroots = WF.size();
    vector<double> EOrb;
    vector< vector<double> > NevE;
    NevE.resize(nroots);
    bool Conventional = dmrginp.read_higherpdm();
    bool ConventionalOverlap = false;
    double densTime=0.0;
    double CoreEnTime=0.0;
    double NEVTime=0.0;
    double totalTimeStart = GetTime();
    double prepTimeStart = GetTime();
    int MaxCore = Info.GetMaxCore();
    
    //print a header
    if (Conventional){
      sprintf(msg,"\n");pout << msg;
      sprintf(msg,"\n=============================================================");pout << msg;
      sprintf(msg,"\n                   BLOCK-NEVPT2 Calculation");pout << msg;
      sprintf(msg,"\n=============================================================");pout << msg;
      sprintf(msg,"\n");pout << msg;
    }
    
    //define the orbital spaces
    int i0,i1,t0,t1,a0,a1;
    int OrbWin[6];
    Info.getBaseName(BaseName);
    Info.GetOrbWin(OrbWin);
    ConventionalOverlap = Info.ConventionalOverlap();
    ENuc = Info.GetENuc();
    i0 = OrbWin[0];
    i1 = OrbWin[1];
    t0 = OrbWin[2];
    t1 = OrbWin[3];
    a0 = OrbWin[4];
    a1 = OrbWin[5];
    int NInternal = i1-i0+1;
    int NActive = big.get_sites().size();
    int NExternal = a1-a0+1;
    int OrbDim = NInternal+NActive+NExternal;
    
    //take into account that active orbitals might have been reordered
    vector<int> ReOrder;
    Info.GetOrbOrder(ReOrder);
    
    //continue header
    if (Conventional){
      sprintf(msg,"\nInternal Orbitals    [%4i ; %4i]",i0,i1);pout << msg;
      sprintf(msg,"\nActive Orbitals      [%4i ; %4i]",t0,t1);pout << msg;
      sprintf(msg,"\nExternal Orbitals    [%4i ; %4i]",a0,a1);pout << msg;
      sprintf(msg,"\nOverlap is calculated in ");pout << msg;
      if (ConventionalOverlap) {sprintf(msg,"conventional ");pout << msg;}
      else {sprintf(msg,"non-conventional ");pout << msg;}
      sprintf(msg,"fashion!");pout << msg;
    }
    
    //check if the number of sites in DMRG coincides with the given orbital space size
    if (NActive!=t1-t0+1){
      sprintf(msg,"\n\nERROR NEVPT2: The number of sites in DMRG does not match the number of active orbitals!");pout << msg;
      sprintf(msg,"\nLeaving NEVPT2 section.......\n\n");pout << msg;
      exit(1);
    }
    
    //the one-electron matrix
    Matrix h,h_eff,h_eff_;
    h.ReSize(OrbDim,OrbDim);
    h_eff.ReSize(OrbDim,OrbDim);
    h_eff_.ReSize(OrbDim,OrbDim);
    
    //the one-electron and two-electron density matrices
    Matrix D,D_;
    D.ReSize(NActive,NActive);
    D_.ReSize(NActive,NActive);
    array_4d<double> D2,D2_,DC2;
    D2.resize(NActive,NActive,NActive,NActive);
    D2_.resize(NActive,NActive,NActive,NActive);
    DC2.resize(NActive,NActive,NActive,NActive);
    
    //the 3pdm
    array_6d D3,DC3,D3_,AuxA,AuxA_;
    if (Conventional){
      D3.resize(NActive);
      DC3.resize(NActive);
      D3_.resize(NActive);
      AuxA.resize(NActive);
      AuxA_.resize(NActive);
      D3.initialize();
      DC3.initialize();
      D3.initialize();
      AuxA.initialize();
      AuxA_.initialize();
    }
    //if not requested, free the memory immediately
    if (!Conventional){
      D3.clear();
      DC3.clear();
      D3_.clear();
      AuxA.clear();
      AuxA_.clear();
    }
    
    //----------------------------
    //Read the Integrals from disk
    //----------------------------
    IntegralContainer IJKL(t1-i0+1,t1-i0+1,t1-i0+1,t1-i0+1,_COULOMB_);
    IntegralContainer IKJL(t1-i0+1,t1-i0+1,t1-i0+1,t1-i0+1,_EXCHANGE_);
    IntegralContainer IAJB(t1-i0+1,t1-i0+1,a1-a0+1,a1-a0+1,_EXCHANGE_);
    IntegralContainer IJAB(t1-i0+1,t1-i0+1,a1-a0+1,a1-a0+1,_COULOMB_);
    IntegralContainer IJKA(t1-i0+1,t1-i0+1,t1-i0+1,a1-a0+1,_NO_SYMM_);
    IntegralContainer IKJA(t1-i0+1,t1-i0+1,t1-i0+1,a1-a0+1,_NO_SYMM_);
    sprintf(msg,"%s.block.IJKL.tmp",BaseName);IJKL.SetFileName(msg);
    sprintf(msg,"%s.block.IKJL.tmp",BaseName);IKJL.SetFileName(msg);
    sprintf(msg,"%s.block.IAJB.tmp",BaseName);IAJB.SetFileName(msg);
    sprintf(msg,"%s.block.IJAB.tmp",BaseName);IJAB.SetFileName(msg);
    sprintf(msg,"%s.block.IJKA.tmp",BaseName);IJKA.SetFileName(msg);
    sprintf(msg,"%s.block.IKJA.tmp",BaseName);IKJA.SetFileName(msg);
    ReadIntegrals(IJKL,IKJL,IAJB,IJAB,IJKA,IKJA,h,OrbWin,BaseName,ReOrder);
    sprintf(msg,"%s.block.KIAJ.tmp",BaseName);
    boost::shared_ptr<IntegralContainer> KIAJ = IKJA.ReverseOrder(msg);
    
    //-----------------------------------------------------
    //construct or read the effective one-electron matrices
    //-----------------------------------------------------
    //GenerateHeff(OrbWin,h,h_eff,h_eff_,IJKL,IKJL,IJAB,IAJB,IJKA,IKJA);
    //ReadHeff(OrbWin,h_eff,h_eff_,BaseName,ReOrder);
    //Info.AddH(h);Info.AddH(h_eff);Info.AddH(h_eff_);
    Info.GetH(1,h_eff);
    Info.GetH(2,h_eff_);

    //-------------------------
    //Read the orbital energies
    //-------------------------
    Info.GetOrbEnergies(EOrb);
    double prepTimeEnd = GetTime();
    
    for (int iroot=0;iroot<nroots;iroot++){
      
      //----------------------------------------------
      //read in the 1-body and 2-body density matrices
      //----------------------------------------------
      double densTimeStart = GetTime();
      ReadDensityMatrices(D,D2,iroot,ReOrder);
      
      //------------------------
      //build the hole densities
      //------------------------
      BuildHoleDensities(OrbWin,D,D_,D2,D2_);
      
      //-------------------------------------------
      //if conventional, read 3-body density matrix 
      //-------------------------------------------
      if (Conventional){
        sprintf(msg, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/spatial_threepdm_RI_bin.", iroot,".rsdm");
        Read3PDM(D3,msg);
        sprintf(msg, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/A_bin.", iroot,".rsdm");
        Read3PDM(AuxA,msg);
        sprintf(msg, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/A_prime_bin.", iroot,".rsdm");
        Read3PDM(AuxA_,msg);
      }
      ConstructAuxPDM(DC3,D3_,D3,D2,DC2,D2_,D,D_,Conventional);
      double densTimeEnd = GetTime();
      densTime += densTimeEnd - densTimeStart;
      //------------------------------
      //calculate the orbital energies
      //------------------------------
      //CalcOrbEnergies(OrbWin,h,IJKL,IKJL,IJAB,IAJB,EOrb,D);
      
      //-------------------------------
      //Evaluate the unperturbed energy
      //-------------------------------
      double CoreEnTimeStart = GetTime();
      CalcCoreEnergy(EEl,h,D,i0,i1,t0,t1,IJKL,IKJL,D2,Eval);
      E0 = EEl + ENuc;
      Info.SetE0(iroot,E0);
      Info.SetEVal(iroot,Eval);
      if (Conventional){
        sprintf(msg,"\n");pout << msg;
        sprintf(msg,"\n");pout << msg;
        sprintf(msg,"    **********\n");pout << msg;
        sprintf(msg,"      ROOT %i\n",iroot);pout << msg;
        sprintf(msg,"    **********\n");pout << msg;
        sprintf(msg,"\n");pout << msg;
        sprintf(msg,"-------------------------\n");pout << msg;
        sprintf(msg,"E0 = %6.8lf\n",E0);pout << msg;
        sprintf(msg,"-------------------------\n");pout << msg;
      }
      double CoreEnTimeEnd = GetTime();
      CoreEnTime += CoreEnTimeEnd - CoreEnTimeStart;
      //--------------------------------------------------------
      //evaluate the eight contributions to the perturbed energy
      //--------------------------------------------------------
      double NEVTimeStart = GetTime();
      NevE[iroot].resize(8);
      //V(ijab)
      double VijabTimeStart = GetTime();
      NevE[iroot][0] = V_ijab(OrbWin,EOrb,big,WF[iroot],IAJB);
      double VijabTimeEnd = GetTime();
      Info.SetTimePerClass(0,VijabTimeEnd-VijabTimeStart);
      //V(iab)
      double ViabTimeStart = GetTime();
      NevE[iroot][1] = V_iab(OrbWin,EOrb,big,WF[iroot],IAJB,IKJL,D,D2,h_eff);
      double ViabTimeEnd = GetTime();
      Info.SetTimePerClass(1,ViabTimeEnd-ViabTimeStart);
      //V(ija)
      double VijaTimeStart = GetTime();
      NevE[iroot][2] = V_ija(OrbWin,EOrb,big,WF[iroot],IKJA,IKJL,IJKL,D,D_,D2,h_eff);
      double VijaTimeEnd = GetTime();
      Info.SetTimePerClass(2,VijaTimeEnd-VijaTimeStart);
      //V(ab)
      double VabTimeStart = GetTime();
      if (mpi_rank()==0) {sprintf(msg,"\nStarting Vab...");pout << msg;}
      NevE[iroot][3] = V_ab(OrbWin,EOrb,big,WF[iroot],IAJB,IKJL,D2,D3,h_eff,Conventional,BaseName,Eval,ConventionalOverlap,MaxCore);
      double VabTimeEnd = GetTime();
      Info.SetTimePerClass(3,VabTimeEnd-VabTimeStart);
      //V(ij)
      double VijTimeStart = GetTime();
      if (mpi_rank()==0) {sprintf(msg,"Done\nStarting Vij...");pout << msg;}
      NevE[iroot][4] = V_ij(OrbWin,EOrb,big,WF[iroot],IKJL,D2_,D3_,h_eff_,Conventional,Eval,ConventionalOverlap,MaxCore);
      double VijTimeEnd = GetTime();
      Info.SetTimePerClass(4,VijTimeEnd-VijTimeStart);
      //V(ia)
      double ViaTimeStart = GetTime();
      if (mpi_rank()==0) {sprintf(msg,"Done\nStarting Via...");pout << msg;}
      NevE[iroot][5] = V_ia(OrbWin,EOrb,WF[iroot],big,IKJA,IKJL,IJKA,KIAJ,D,DC2,DC3,h_eff,h_eff_, Conventional,BaseName,Eval,ConventionalOverlap,MaxCore);
      if (mpi_rank()==0) {sprintf(msg,"Done");pout << msg;}
      double ViaTimeEnd = GetTime();
      Info.SetTimePerClass(5,ViaTimeEnd-ViaTimeStart);
      //V(a)
      double VaTimeStart = GetTime();
      NevE[iroot][6] = V_a(OrbWin,EOrb,big,WF[iroot],IKJL,IKJA,KIAJ,DC3,DC2,D,AuxA,h_eff_,h_eff,Conventional,ConventionalOverlap);
      double VaTimeEnd = GetTime();
      if (Conventional) Info.SetTimePerClass(6,VaTimeEnd-VaTimeStart);
      //V(i)
      double ViTimeStart = GetTime();
      NevE[iroot][7] = V_i(OrbWin,EOrb,big,WF[iroot],IKJL,DC3,DC2,D,D_,AuxA_,h_eff_,h_eff,Conventional,ConventionalOverlap);
      double ViTimeEnd = GetTime();
      if (Conventional) Info.SetTimePerClass(7,ViTimeEnd-ViTimeStart);
      
      double NEVTimeEnd = GetTime();
      NEVTime += NEVTimeEnd - NEVTimeStart;
      Info.SetE(iroot,NevE[iroot]);
      //-----------------
      //print the results
      //-----------------
      double TotEnergy = E0;
      if (Conventional&&mpi_rank()==0){
        sprintf(msg,"\nE(0,ijab) = %5.12lf",NevE[iroot][0]);pout << msg;TotEnergy += NevE[iroot][0];
        sprintf(msg,"\nE(-1,iab) = %5.12lf",NevE[iroot][1]);pout << msg;TotEnergy += NevE[iroot][1];
        sprintf(msg,"\nE(1,ija)  = %5.12lf",NevE[iroot][2]);pout << msg;TotEnergy += NevE[iroot][2];
        sprintf(msg,"\nE(-2,ab)  = %5.12lf",NevE[iroot][3]);pout << msg;TotEnergy += NevE[iroot][3];
        sprintf(msg,"\nE(2,ij)   = %5.12lf",NevE[iroot][4]);pout << msg;TotEnergy += NevE[iroot][4];
        sprintf(msg,"\nE(0,ia)   = %5.12lf",NevE[iroot][5]);pout << msg;TotEnergy += NevE[iroot][5];
        sprintf(msg,"\nE(-1,a)   = %5.12lf",NevE[iroot][6]);pout << msg;TotEnergy += NevE[iroot][6];
        sprintf(msg,"\nE(1,i)    = %5.12lf",NevE[iroot][7]);pout << msg;TotEnergy += NevE[iroot][7];
        sprintf(msg,"\n");pout << msg;
        sprintf(msg,"\n-------------------------------");pout << msg;
        sprintf(msg,"\nTOTAL ENERGY = %6.12lf",TotEnergy);pout << msg;
        sprintf(msg,"\n-------------------------------");pout << msg;
        sprintf(msg,"\n");pout << msg;
        sprintf(msg,"\n");pout << msg;
      }//conventional
      
    }//iroot
    double totalTimeEnd = GetTime();
    if (Conventional&&mpi_rank()==0){
      sprintf(msg,"\n");pout << msg;
      sprintf(msg,"\nTimings within BLOCK-NEVPT2:");pout << msg;
      sprintf(msg,"\nPreparation Time          = %7.6lf sec",prepTimeEnd-prepTimeStart);pout << msg;
      sprintf(msg,"\nDensity Formation Time    = %7.6lf sec",densTime);pout << msg;
      sprintf(msg,"\nCoreEnergy Formation Time = %7.6lf sec",CoreEnTime);pout << msg;
      sprintf(msg,"\nNEVPT2 Time               = %7.6lf sec",NEVTime);pout << msg;
      sprintf(msg,"\n----------------------------------------");pout << msg;
      sprintf(msg,"\nTotal Time                = %lf sec",totalTimeEnd-totalTimeStart);pout << msg;
      sprintf(msg,"\n----------------------------------------");pout << msg;
    }
    else{
      //store the timings
      Info.SetTime(0,prepTimeEnd-prepTimeStart);
      Info.SetTime(1,densTime);
      Info.SetTime(2,CoreEnTime);
      Info.SetTime(4,NEVTime);
    }
    
    //---------------------
    //clear memory and disk
    //---------------------
    if (Conventional){
      D3.clear();
      DC3.clear();
      D3_.clear();
      AuxA.clear();
      AuxA_.clear();
    }
    IJKL.Clear();
    IKJL.Clear();
    IJKA.Clear();
    IKJA.Clear();
    IJAB.Clear();
    IAJB.Clear();
    KIAJ->Clear();
    D.CleanUp();
    D2.Clear();
    DC2.Clear();
    D2_.Clear();
    h.CleanUp();
    h_eff.CleanUp();
    h_eff_.CleanUp();
    
  }
  
  //============================================================================
  // Evaluate the V(i) contribution to the NEVPT2 energy in unconventional mode
  //============================================================================
  double V_i(vector<WavefunctionArray> &T, SpinBlock &big, Wavefunction &WF, NEVPT2Info &Info, int iroot){
    
    char msg[512];
    int OrbWin[6];
    Info.GetOrbWin(OrbWin);
    int t,u,v,i,a,tu;
    int dummy;
    
    //the orbital spaces
    int i0 = OrbWin[0];
    int i1 = OrbWin[1];
    int t0 = OrbWin[2];
    int t1 = OrbWin[3];
    int a0 = OrbWin[4];
    int a1 = OrbWin[5];
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    int OrbDim = NInternal+NActive+NExternal;
    
    //the NEVPT2 quantities
    double E=0;
    vector<double> N,D;
    N.resize(NInternal,0.0);
    D.resize(NInternal,0.0);
    double Norm,Denominator;
    
    //take time
    double Start=GetTime();
    
    //get the orbital energies
    vector<double> EOrb;
    Info.GetOrbEnergies(EOrb);
    
    //the possible spin quanta
    SpinQuantum opQ(1,SpinSpace(1),IrrepSpace(0));
    SpinQuantum WFQ = WF.get_deltaQuantum(0);
    vector<SpinQuantum> VQ = opQ + WFQ;
    
    //-------------------------------------------
    //get the wavefunction V(i)|psi> and evaluate 
    //  a) the norm 
    //  b) the expectation value
    //-------------------------------------------
    int NumProcs = mpi_world_size();
    for (int iproc=0;iproc<NumProcs;iproc++){
      bool Sender = (mpi_rank()==iproc);
      for (int iquanta=0;iquanta<T.size();iquanta++){
        T[iquanta].OpenFileRead();
        bool EndOfArray = false;
        T[iquanta].ResetBuffer();
        while (!EndOfArray){
          boost::shared_ptr<Wavefunction> Vi;
          if (Sender){
            Vi = T[iquanta].GetOpFromBuffer(dummy,i,EndOfArray);
          }//sender
          else {
            Vi = boost::make_shared<Wavefunction>(VQ[iquanta],&big,true);
          }//receiver
#ifndef SERIAL
          mpi::communicator world;
          mpi::broadcast(world,EndOfArray,iproc);
          mpi::broadcast(world,i,iproc);
#endif
          if (!EndOfArray){
#ifndef SERIAL
            //broadcast the V(a) function
            mpi::broadcast(world,*Vi,iproc);
#endif
  //Evaluate the overlap
            N[i] += DotProduct(*Vi,*Vi);
            //Evaluate the energy expectation value <psi|V(i)+ H V(i)|psi>
            boost::shared_ptr<Wavefunction> sigma (new Wavefunction(VQ[iquanta],&big,true));
            big.multiplyH_Q(*Vi,sigma.get(),MAX_THRD,VQ[iquanta]);
            D[i] += DotProduct(*Vi,*sigma);
          }//!EndOfArray
        }//i
        T[iquanta].CloseFileRead();
      }//iquanta
    }//iproc
    //----------------------------------------------------
    //Evaluate the contribution to the perturbation energy
    //----------------------------------------------------
    double E0 = Info.getEVal(iroot);
    for (i=0;i<NInternal;i++){
      Norm = N[i];
      Denominator = D[i];
      if (Norm<_NEV_ZERO) continue;
      Denominator *= 1/Norm;
      Denominator += -EOrb[i] - E0;
      E -= Norm/Denominator;
      //sprintf(msg,"\n%i %4.12lf %4.12lf  %4.12lf  %4.12lf",i,Norm,Denominator+EOrb[i]+E0,-EOrb[i],-Norm/Denominator);pout << msg;
    }

    //take time
    double Finish = GetTime();
    double Time = Finish-Start;
    Info.AddTime(4,Time);
    
    //return the energy
    return E;
  }
  
  //============================================================================
  // Evaluate the V(a) contribution to the NEVPT2 energy in unconventional mode
  //============================================================================
  double V_a(vector<WavefunctionArray> &T, SpinBlock &big, Wavefunction &WF, NEVPT2Info &Info, int iroot){
    char msg[512];
    int OrbWin[6];
    Info.GetOrbWin(OrbWin);
    int t,u,v,i,a,tu;
    int dummy;
    
    //the orbital spaces
    int i0 = OrbWin[0];
    int i1 = OrbWin[1];
    int t0 = OrbWin[2];
    int t1 = OrbWin[3];
    int a0 = OrbWin[4];
    int a1 = OrbWin[5];
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    int OrbDim = NInternal+NActive+NExternal;
    
    //the NEVPT2 quantities
    double E=0;
    vector<double> N,D;
    N.resize(NExternal,0.0);
    D.resize(NExternal,0.0);
    double Norm,Denominator;
    for (a=0;a<NExternal;a++){
      N[a]=0.0;
      D[a]=0.0;
    }//a
    
    double Start = GetTime();
    
    //get the orbital energies
    vector<double> EOrb;
    Info.GetOrbEnergies(EOrb);
       
    //the possible spin quanta
    SpinQuantum opQ(-1,SpinSpace(1),IrrepSpace(0));
    SpinQuantum WFQ = WF.get_deltaQuantum(0);
    vector<SpinQuantum> VQ = opQ + WFQ;
    
    //-------------------------------------------
    //get the wavefunction V(a)|psi> and evaluate 
    //  a) the norm 
    //  b) the expectation value
    //-------------------------------------------
    int NumProcs = mpi_world_size();
    for (int iproc=0;iproc<NumProcs;iproc++){
      for (int iquanta=0;iquanta<T.size();iquanta++){
        T[iquanta].OpenFileRead();
        bool Sender = (mpi_rank()==iproc);
        bool EndOfArray = false;
        T[iquanta].ResetBuffer();
        while (!EndOfArray){
          boost::shared_ptr<Wavefunction> Va;
          if (Sender){
            Va = T[iquanta].GetOpFromBuffer(dummy,a,EndOfArray);
          }//sender
          else {
            Va = boost::make_shared<Wavefunction>(VQ[iquanta],&big,true);
          }//receiver
#ifndef SERIAL
          mpi::communicator world;
          mpi::broadcast(world,EndOfArray,iproc);
          mpi::broadcast(world,a,iproc);
#endif
          if (!EndOfArray){
#ifndef SERIAL
            //broadcast the V(a) function
            mpi::broadcast(world,*Va,iproc);
#endif
            //Evaluate the overlap
            N[a] += DotProduct(*Va,*Va);
            //Evaluate the energy expectation value <psi|V(i)+ H V(i)|psi>
            boost::shared_ptr<Wavefunction> sigma (new Wavefunction(VQ[iquanta],&big,true));
            big.multiplyH_Q(*Va,sigma.get(),MAX_THRD,VQ[iquanta]);
            D[a] += DotProduct(*Va,*sigma);
          }//!EndOfArray
        }//a
        T[iquanta].CloseFileRead();
      }//iquanta
    }//iproc
    
    //----------------------------------------------------
    //Evaluate the contribution to the perturbation energy
    //----------------------------------------------------
    double E0 = Info.getEVal(iroot);
    for (a=0;a<NExternal;a++){
      Norm = N[a];
      Denominator = D[a];
      //if (Norm<_NEV_ZERO) continue;
      Denominator *= 1/Norm;
      Denominator += EOrb[a+a0] - E0;
      if (Norm<_NEV_ZERO) continue;
      E -= Norm/Denominator;
      //sprintf(msg,"\n%i %4.12lf %4.12lf  %4.12lf  %4.12lf",a,Norm,(Denominator-EOrb[a+a0]+E0)*Norm,EOrb[a+a0],-Norm/Denominator);pout << msg;
    }

    //take time
    double Finish = GetTime();
    double Time = Finish-Start;
    Info.AddTime(4,Time);

    //return the energy
    return E;
  }


}




