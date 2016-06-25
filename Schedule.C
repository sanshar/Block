#include <iostream>
#include <fstream>
#include <string.h>
#include <boost/algorithm/string.hpp>
#include <ctype.h>
#include "global.h"
#include "input.h"
#include "pario.h"


using namespace std;

void SpinAdapted::Input::generateDefaultSchedule(){
#ifndef SERIAL
  if (mpigetrank() == 0) {
#endif

   if (m_sweep_tol <= 0.0) {
      pout << "Using the default tolerance sweep tolerance of 1.0e-5."<<endl;
      m_sweep_tol = 1.0e-5;
   }

   if (m_schedule_type_default && !m_schedule_type_backward) {
      int nentry = 0;
      m_sweep_iter_schedule.resize(nentry);
      m_sweep_state_schedule.resize(nentry);
      m_sweep_qstate_schedule.resize(nentry,0);
      m_sweep_tol_schedule.resize(nentry);
      m_sweep_noise_schedule.resize(nentry);
      m_sweep_additional_noise_schedule.resize(nentry,0);

      double sweeptol = m_sweep_tol;
      int lastiter = 0;
      int firstSched = 0;
      int sweepCount = 0;
      int nSched=14;
      int defM [] = {50, 100, 250, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000};
      int defIter [] = {8, 8, 8, 8, 8, 4, 4, 4, 4, 4, 4, 4, 4, 4};
      double defNoise [] = {1.0e-4, 1.0e-4, 1e-4, 1e-4, 5e-5, 5e-5, 5e-5, 5e-5, 5e-5, 5e-5, 5e-5, 5e-5, 5e-5, 5e-5};
      double defTol [] = {1.0e-5, 1.0e-5, 1e-5, 1e-5, 5e-5, 5e-6, 5e-6, 5e-6, 5e-6, 5e-6, 5e-6, 5e-6, 5e-6, 5e-6};

      for (int i=0; i<nSched;++i){
         if (firstSched==0){
            if (m_startM == m_maxM){
               //pout << sweepCount << " " << m_startM << " " << defTol[i] << " " << defNoise[i] << endl;
               m_sweep_iter_schedule.push_back(sweepCount); m_sweep_state_schedule.push_back(m_startM); m_sweep_tol_schedule.push_back(1E-5);  m_sweep_noise_schedule.push_back(1E-4);
               sweepCount += 8;
               m_sweep_iter_schedule.push_back(sweepCount); m_sweep_state_schedule.push_back(m_startM); m_sweep_tol_schedule.push_back(5E-6);  m_sweep_noise_schedule.push_back(5E-5);
               break;
            }

            else if (m_startM == defM[i]){
               //pout << sweepCount << " " << defM[i] << " " << defTol[i] << " " << defNoise[i] << endl;
               firstSched = 1;
               m_sweep_iter_schedule.push_back(sweepCount); m_sweep_state_schedule.push_back(defM[i]); m_sweep_tol_schedule.push_back(defTol[i]);  m_sweep_noise_schedule.push_back(defNoise[i]);
               sweepCount += defIter[i];
            }
            else if (m_startM < defM[i]){
               firstSched = 1;
               //pout << sweepCount << " " << m_startM << " " << defTol[i-1] << " " << defNoise[i-1] << endl;
               m_sweep_iter_schedule.push_back(sweepCount); m_sweep_state_schedule.push_back(m_startM); m_sweep_tol_schedule.push_back(defTol[i-1]);  m_sweep_noise_schedule.push_back(defNoise[i-1]);
               sweepCount += defIter[i-1];

               //pout << sweepCount << " " << defM[i] << " " << defTol[i] << " " << defNoise[i] << endl;
               m_sweep_iter_schedule.push_back(sweepCount); m_sweep_state_schedule.push_back(defM[i]); m_sweep_tol_schedule.push_back(defTol[i]);  m_sweep_noise_schedule.push_back(defNoise[i]);
               sweepCount += defIter[i];
            }
         }
         else{//After first iteration
            if (defM[i]>=m_maxM){
               //pout << sweepCount << " " << m_maxM << " " << defTol[i] << " " << defNoise[i] << endl;
               m_sweep_iter_schedule.push_back(sweepCount); m_sweep_state_schedule.push_back(m_maxM); m_sweep_tol_schedule.push_back(defTol[i]);  m_sweep_noise_schedule.push_back(defNoise[i]);
               sweepCount += defIter[i];
               break;
            }
            else{
               //pout << sweepCount << " " << defM[i] << " " << defTol[i] << " " << defNoise[i] << endl;
               m_sweep_iter_schedule.push_back(sweepCount); m_sweep_state_schedule.push_back(defM[i]); m_sweep_tol_schedule.push_back(defTol[i]);  m_sweep_noise_schedule.push_back(defNoise[i]);
               sweepCount += defIter[i];
            }
         }
      }
         lastiter = m_sweep_iter_schedule.back();
         //lastiter = sweepCount;
         m_sweep_iter_schedule.push_back(lastiter+2); m_sweep_state_schedule.push_back(m_maxM); m_sweep_tol_schedule.push_back(sweeptol/10.0);  m_sweep_noise_schedule.push_back(0.0e-5);


      if (m_twodot_to_onedot_iter < 18 && m_algorithm_type == TWODOT_TO_ONEDOT) {
        if (m_twodot_to_onedot_iter <= 0)
           pout << "Sweep at which the switch from twodot to onedot will happen -> "<<lastiter+4<<endl;
        m_twodot_to_onedot_iter = lastiter+4;
      }
      if (m_maxiter <= m_sweep_iter_schedule.back()) {
        //pout << "With the default schedule and maxM specified, maxiter has to be at least "<<lastiter+6<<endl;
        //pout << "changing maxiter to "<<lastiter+6<<endl;
        m_maxiter = lastiter+6;
      }
    }
 else if (m_schedule_type_backward) {

   int nentry = 0;
   m_sweep_iter_schedule.resize(nentry);
   m_sweep_state_schedule.resize(nentry);
   m_sweep_qstate_schedule.resize(nentry,0);
   m_sweep_tol_schedule.resize(nentry);
   m_sweep_noise_schedule.resize(nentry);
   m_sweep_additional_noise_schedule.resize(nentry,0);

   double sweeptol = m_sweep_tol;
   int lastiter = 0;
   int firstSched = 0;
   int sweepCount = 0;
   int nSched=17;
   int lastM=0; //To be changed to a variable
   lastM = m_lastM;
   if (lastM==0) lastM=50;
   if (lastM > m_startM) {
      pout << "lastM is larger than startM" << endl;
      pout << "Make sure you specify a lastM larger than " << lastM << endl;
      pout << "or specify a larger startM " << endl;
      abort();
   }
   double bNoise = 0.0;
   double bTol = 5e-6;

   int defIter [] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
   int defM [] = {10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, 500, 250, 100, 50, 25, 10, 1};

   for (int i=0; i<nSched;++i){
      if (firstSched==0){
         if (lastM == m_startM){
            m_sweep_iter_schedule.push_back(sweepCount); m_sweep_state_schedule.push_back(m_startM); m_sweep_tol_schedule.push_back(bTol);  m_sweep_noise_schedule.push_back(bNoise);
            sweepCount += 3;
            break;
         }
         else if(m_startM==defM[i]){ 
            firstSched = 1;
            m_sweep_iter_schedule.push_back(sweepCount); m_sweep_state_schedule.push_back(defM[i]); m_sweep_tol_schedule.push_back(bTol);  m_sweep_noise_schedule.push_back(bNoise);
            sweepCount += defIter[i];
      }
         else if (m_startM > defM[i]){
            firstSched = 1;
            m_sweep_iter_schedule.push_back(sweepCount); m_sweep_state_schedule.push_back(m_startM); m_sweep_tol_schedule.push_back(bTol);  m_sweep_noise_schedule.push_back(bNoise);
            sweepCount += defIter[i-1];
            m_sweep_iter_schedule.push_back(sweepCount); m_sweep_state_schedule.push_back(defM[i]); m_sweep_tol_schedule.push_back(bTol);  m_sweep_noise_schedule.push_back(bNoise);
            sweepCount += defIter[i];
         }
      }
      else{//After first iteration
         if (defM[i]<=lastM){
            m_sweep_iter_schedule.push_back(sweepCount); m_sweep_state_schedule.push_back(lastM); m_sweep_tol_schedule.push_back(bTol);  m_sweep_noise_schedule.push_back(bNoise);
            sweepCount += defIter[i];
            break;
         }
         else{
            m_sweep_iter_schedule.push_back(sweepCount); m_sweep_state_schedule.push_back(defM[i]); m_sweep_tol_schedule.push_back(bTol);  m_sweep_noise_schedule.push_back(bNoise);
            sweepCount += defIter[i];
         }

       }
   }
      lastiter = m_sweep_iter_schedule.back();
      m_maxiter = lastiter+3;
      pout << "maxiter " << m_maxiter << endl;
 }

 // Questions... What happens to twodottoonedot?
 if(m_algorithm_type == TWODOT_TO_ONEDOT && m_twodot_to_onedot_iter == 0)
   m_twodot_to_onedot_iter = min(m_sweep_iter_schedule.back()+2, m_maxiter-1);

  if (m_maxiter < m_sweep_iter_schedule.back()) {
    pout << "maximum iterations allowed is less than the last sweep iteration in your schedule."<<endl;
    pout << m_maxiter <<" < "<< (m_sweep_iter_schedule.back())<<endl;
    pout << "either increase the max_iter or reduce the number of sweeps"<<endl;
    abort();
  }
#ifndef SERIAL
  }
#endif
}



