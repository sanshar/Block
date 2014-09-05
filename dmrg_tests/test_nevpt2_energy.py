#!/usr/bin/python

import os
from math import sqrt
import numpy as N
import struct

def run(args):
  #the tolerance
  tol = float(args[2])
  
  #read the reference values
  RefValues = []
  RefFile = open(args[1],"r")
  for line in RefFile:
    linesp = line.split()
    if len(linesp) >=3:
      RefValues.append(float(linesp[2]))

  #read and compare the calculated values
  OutFile = open("dmrg.out","r")
  for line in OutFile:
    linesp = line.split()
    #E(ijab)
    if len(linesp)==3:
      if linesp[0]=="E(0,ijab)":
        calc_e = float(linesp[2])
        if abs(calc_e-RefValues[0]) > tol:
          print calc_e,"-", RefValues[0], " > ", tol
          print "FAILED ...."
        else:
          print "PASSED ...."  
    #E(iab)
    if len(linesp)==3:
      if linesp[0]=="E(-1,iab)":
        calc_e = float(linesp[2])
        if abs(calc_e-RefValues[1]) > tol:
          print calc_e,"-", RefValues[1], " > ", tol
          print "FAILED ...."
        else:
          print "PASSED ...."  
    #E(ija)
    if len(linesp)==3:
      if linesp[0]=="E(1,ija)":
        calc_e = float(linesp[2])
        if abs(calc_e-RefValues[2]) > tol:
          print calc_e,"-", RefValues[2], " > ", tol
          print "FAILED ...."
        else:
          print "PASSED ...."  
    #E(ab)
    if len(linesp)==3:
      if linesp[0]=="E(-2,ab)":
        calc_e = float(linesp[2])
        if abs(calc_e-RefValues[3]) > tol:
          print calc_e,"-", RefValues[3], " > ", tol
          print "FAILED ...."
        else:
          print "PASSED ...."  
    #E(ij)
    if len(linesp)==3:
      if linesp[0]=="E(2,ij)":
        calc_e = float(linesp[2])
        if abs(calc_e-RefValues[4]) > tol:
          print calc_e,"-", RefValues[4], " > ", tol
          print "FAILED ...."
        else:
          print "PASSED ...."  
    #E(ia)
    if len(linesp)==3:
      if linesp[0]=="E(0,ia)":
        calc_e = float(linesp[2])
        if abs(calc_e-RefValues[5]) > tol:
          print calc_e,"-", RefValues[5], " > ", tol
          print "FAILED ...."
        else:
          print "PASSED ...."  
    #E(a)
    if len(linesp)==3:    
      if linesp[0]=="E(-1,a)":
        calc_e = float(linesp[2])
        if abs(calc_e-RefValues[6]) > tol:
          print calc_e,"-", RefValues[6], " > ", tol
          print "FAILED ...."
        else:
          print "PASSED ...."  
    #E(i)
    if len(linesp)==3:
      if linesp[0]=="E(0,i)":
        calc_e = float(linesp[2])
        if abs(calc_e-RefValues[7]) > tol:
          print calc_e,"-", RefValues[7], " > ", tol
          print "FAILED ...."
        else:
          print "PASSED ...."  

if __name__=="__main__":
    import sys
    run(sys.argv)
