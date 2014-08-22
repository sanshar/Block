#!/usr/bin/python

import os
from math import sqrt
import numpy as N
import struct

def run():

  nelec = 8.0 #the number of electrons
  
  fileint1 = open("qcdmrg.int1", "r")
  fileint2 = open("qcdmrg.int2", "r")
  norb = int(fileint1.readline().split()[0])
  int1 = N.zeros((norb, norb))
  int2 = N.zeros((norb, norb, norb, norb))

  #read onelectron integral
  for line in fileint1.readlines():
    linesp = line.split()
    int1[int(linesp[0]),int(linesp[1])] = float(linesp[2])
    int1[int(linesp[1]),int(linesp[0])] = float(linesp[2])

  #read twoelectron integrals
  fileint2.readline()
  for line in fileint2.readlines():
    linesp = line.split()
    int2[int(linesp[0]),int(linesp[1]),int(linesp[2]),int(linesp[3])] = float(linesp[4])
    int2[int(linesp[0]),int(linesp[3]),int(linesp[2]),int(linesp[1])] = float(linesp[4])
    int2[int(linesp[2]),int(linesp[1]),int(linesp[0]),int(linesp[3])] = float(linesp[4])
    int2[int(linesp[2]),int(linesp[3]),int(linesp[0]),int(linesp[1])] = float(linesp[4])
    int2[int(linesp[3]),int(linesp[0]),int(linesp[1]),int(linesp[2])] = float(linesp[4])
    int2[int(linesp[3]),int(linesp[2]),int(linesp[1]),int(linesp[0])] = float(linesp[4])
    int2[int(linesp[1]),int(linesp[0]),int(linesp[3]),int(linesp[2])] = float(linesp[4])
    int2[int(linesp[1]),int(linesp[2]),int(linesp[3]),int(linesp[0])] = float(linesp[4])


  #test energy for two lowest states
  energyfile = open("dmrg.e","rb")
  for state in range(2):
    file1 = open("spatial_twopdm.%i.%i.txt"%(state,state),"r")
    rdm2 = N.zeros((norb,norb,norb,norb))

    file1.readline()
    for line in file1.readlines():
      linesp = line.split()
      rdm2[int(linesp[0]),int(linesp[1]),int(linesp[2]),int(linesp[3])] = float(linesp[4])*2.0

    rdm1 = N.zeros((norb, norb)) 
    for i in range(norb):
      for j in range(norb):
        for k in range(norb):
            rdm1[i,j] += rdm2[i,k,k,j]/(nelec-1.0)

    #trace over int1 and int2 to get energy
    energy1e = 0.0;
    energy2e = 0.0;
    electron = 0.0

    for i in range(norb):
      for j in range(norb):
        energy1e += int1[i,j]*rdm1[i,j]

    for i in range(norb):
      for j in range(norb):
        for k in range(norb):
          for l in range(norb):
            energy2e += 0.5*int2[i,j,k,l] * rdm2[i,j,l,k]
            
    calc_e = struct.unpack('d', energyfile.read(8))[0]

    if (abs(calc_e-energy1e-energy2e+57.9042346988) > 1.0e-8):
      print "FAILED ...."
    else:
      print "PASSED ...."

if __name__=="__main__":
    run()
