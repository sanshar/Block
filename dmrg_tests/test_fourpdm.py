#!/usr/bin/python

import os
from math import sqrt
import numpy as N

factor = 1.0

def run(file1, file2, tol):

  #int1
  filer1 = open(file1,"r")
  filer2 = open(file2,"r")
  sz = int(filer1.readline().split()[0])
  sz2 = int(filer2.readline().split()[0])
  print "npdm dimensions = ", sz, sz2

  mat1 = N.zeros((sz,sz,sz,sz,sz,sz,sz,sz))
  mat2 = N.zeros((sz2,sz2,sz2,sz2,sz2,sz2,sz2,sz2))
  for line in filer1.readlines():
    linesp = line.split()
    mat1[int(linesp[0]),int(linesp[1]),int(linesp[2]),int(linesp[3]),int(linesp[4]),int(linesp[5]),int(linesp[6]),int(linesp[7])] = float(linesp[8])

  for line in filer2.readlines():
    linesp = line.split()
    mat2[int(linesp[0]),int(linesp[1]),int(linesp[2]),int(linesp[3]),int(linesp[4]),int(linesp[5]),int(linesp[6]),int(linesp[7])] = float(factor)*float(linesp[8])

  filer1.close()
  filer2.close()

  # Loop over mats
  val = 0.0
  for i in xrange(0,sz):
    for j in xrange(0,sz):
      for k in xrange(0,sz):
        for l in xrange(0,sz):
          for m in xrange(0,sz):
            for n in xrange(0,sz):
              for p in xrange(0,sz):
                for q in xrange(0,sz):
                  res = (mat1[i,j,k,l,m,n,p,q] - mat2[i,j,k,l,m,n,p,q]) * (mat1[i,j,k,l,m,n,p,q] - mat2[i,j,k,l,m,n,p,q])
                  val = val + sqrt(res)

  
  dev = val/sz**8
  print "mean deviation",dev
  if (dev > float(tol)):
    print "FAILED ...."
  else:
    print "PASSED ...."

if __name__=="__main__":
    import sys
    run(sys.argv[1], sys.argv[2], sys.argv[3])
