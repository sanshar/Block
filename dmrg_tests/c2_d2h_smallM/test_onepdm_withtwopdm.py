#!/usr/bin/python

import os
from math import sqrt
import numpy as N
import struct

def run():

  nelec = 8.0 #the number of electrons
  fileint1 = open("qcdmrg.int1", "r")
  norb = int(fileint1.readline().split()[0])
  
  for state in range(2):
    file1 = open("node0/spatial_twopdm.%i.%i.txt"%(state,state),"r")
    rdm2 = N.zeros((norb,norb,norb,norb))

    file1.readline()
    for line in file1.readlines():
      linesp = line.split()
      rdm2[int(linesp[0]),int(linesp[1]),int(linesp[2]),int(linesp[3])] = float(linesp[4])*2.0

    rdm1 = N.zeros((norb, norb)) 

    rdm1r = N.zeros((norb, norb)) 
    file2 = open("node0/spatial_onepdm.%i.%i.txt"%(state,state),"r")
    file2.readline()
    for line in file2.readlines():
      linesp = line.split()
      rdm1r[int(linesp[0]),int(linesp[1])] = float(linesp[2])


    for i in range(norb):
      for j in range(norb):
        for k in range(norb):
            rdm1[i,j] += rdm2[i,k,k,j]/(nelec-1.0)

    error = 0.0
    for i in range(norb):
      for j in range(norb):
        error += (rdm1[i,j] - rdm1r[i,j])**2

    if (abs(error) > 1.0e-4):
      print "FAILED ...."
    else:
      print "PASSED ...."

if __name__=="__main__":
    run()
