import numpy
import scipy
import pdb

''' this function takes one electron and two electron integral files and makes new integral files in which the basis functions are eigenfunctions of 
D_{inf}h point group. '''
'''
Lz0 = [0, 1, 4, 5, 6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 2, 3, 38, 39, 40, 41, 42, 44, 45, 47, 48, 49, 50, 51]
Lz1a = [[20, 21, 22, 23, 24, 25, 26, 27], [28, 29, 30, 31, 32, 33, 34, 35]]
Lz1b = [[52, 53, 54, 55, 56, 57, 58, 59], [60, 61, 62, 63, 64, 65, 66, 67]]
Lz2 = [[9, 13, 36, 37], [18, 19, 43, 46]]
Larray = [Lz0, Lz1a, Lz1b, Lz2];
mLz = [0,1,1,2]
'''
''' pvdz'''
norbs = 11
Lz0a = [0, 2]
Lz0b = [3]
Lz1b = [[5], [9]]
Lz1 = [[1], [7]]
Lz2 = [[4],[8]]
Lz2b = [[6],[10]]
Larray = [Lz1, Lz2, Lz1b, Lz2b, Lz0a, Lz0b];
mLz = [1,2,1,2,0,0]

#negatives = [50] + [95] + [28] + [68] + [53] + [103]
#negatives = [51,62,63,  106,107,108,    28,29,30,32,33,  64,65,66,68       ,53, 103]
negatives = [9]

screen_error = 1e-7

int2 = numpy.zeros(shape=(norbs, norbs, norbs, norbs));
int1 = numpy.zeros(shape=(norbs, norbs));

coeffs = numpy.zeros(shape=(norbs, norbs)).astype(complex);

index = 0;
m = 0;
lindex = 0
for Li in Larray:
   m = mLz[lindex]
   if type(Li[0]) == type(1):
       print "transforming orbs with L = 0"
       for i in Li:
           coeffs[index, i] = 1.0
           index =index+1
   else :
       for j in range(len(Li[0])):
           scale = 1.0
           if (Li[1][j] in negatives):
              scale = -1.0
              print "index with negatives ", index, index+1

           coeffs[index, Li[0][j] ] =  (-1)**m*1.0/(2.0**0.5)
           #coeffs[index, Li[1][j] ] = -1.0/(2.0**0.5)
           coeffs[index, Li[1][j] ] =  (-1)**m*scale*1.0j/(2.0**0.5)
           index = index+1

           coeffs[index, Li[0][j] ] = 1.0/(2.0**0.5)
           #coeffs[index, Li[1][j] ] = ((-1)**m)*1.0/(2.0**0.5)
           coeffs[index, Li[1][j] ] = scale*-1.0j/(2.0**0.5)
           index = index+1
   lindex=lindex+1

print "read two ints"
#print numpy.dot(coeffs, coeffs.conj().transpose())
#now read integrals
f = open("qcdmrg.int2","r")
lines = f.readlines()
n = int(lines[0].split()[0])
if (n != norbs) :
    print "number of orbitals don't match", n, norbs
    exit(-1)
for line in lines[1:]:
    tokens = line.split()
    int2[int(tokens[0]),int(tokens[1]),int(tokens[2]),int(tokens[3])] = float(tokens[4])
    int2[int(tokens[2]),int(tokens[1]),int(tokens[0]),int(tokens[3])] = float(tokens[4])
    int2[int(tokens[0]),int(tokens[3]),int(tokens[2]),int(tokens[1])] = float(tokens[4])
    int2[int(tokens[2]),int(tokens[3]),int(tokens[0]),int(tokens[1])] = float(tokens[4])

    int2[int(tokens[1]),int(tokens[0]),int(tokens[3]),int(tokens[2])] = float(tokens[4])
    int2[int(tokens[3]),int(tokens[0]),int(tokens[1]),int(tokens[2])] = float(tokens[4])
    int2[int(tokens[1]),int(tokens[2]),int(tokens[3]),int(tokens[0])] = float(tokens[4])
    int2[int(tokens[3]),int(tokens[2]),int(tokens[1]),int(tokens[0])] = float(tokens[4])

f.close()

#pdb.set_trace()

print "reading one ints"
f = open("qcdmrg.int1","r")
lines = f.readlines()
n = int(lines[0].split()[0])
if (n != norbs) :
    exit(-1)
for line in lines[1:]:
    tokens = line.split()
    int1[int(tokens[0]),int(tokens[1])] = float(tokens[2])
    int1[int(tokens[1]),int(tokens[0])] = float(tokens[2])

print "one int transform"
#newint1 = numpy.dot(numpy.dot(coeffs,int1),coeffs.conj().transpose())
newintt = numpy.tensordot(coeffs.conj(), int1, axes=([1],[0]))
newint1 = numpy.tensordot(newintt, coeffs, axes=([1],[1]))
f.close()

#norbs = 6
f = open("new.int1","w")
f.write(str(norbs)+'\n')
for i in range(norbs):
  for j in range(norbs):
     if (abs(newint1[i,j]) >=screen_error):
        if (newint1[i,j].imag >= screen_error):
           print 'problem in ', i, j, newint1[i,j]
        s = str(i)+'  '+str( j)+"  "+str( newint1[i,j].real)+'\n'
        f.write(s)
f.close()

print "step1"
newintt = numpy.tensordot(coeffs.conj(), int2, axes=([1],[1]))
print "step2"
newint2 = numpy.tensordot(coeffs.conj(), newintt, axes=([1],[1]))
print "step3"
newintt = numpy.tensordot(newint2, coeffs, axes=([2],[1]))
print "step4"
newint2 = numpy.tensordot(newintt, coeffs, axes=([2],[1]))

print "writing ints"

f = open("new.int2","w")
f.write(str(norbs)+'\n')
for i in range(norbs):
  for j in range(norbs):
      for k in range(norbs):
          for l in range(norbs):
             if (abs(newint2[i,j,k,l]) >=screen_error ):
                if (newint2[i,j,k,l].imag >= screen_error):
                   print 'problem in ', i, j, k,l,newint2[i,j,k,l]

                s = str(i)+'  '+str( j)+"  "+str( k)+"  "+str( l)+"  "+str( newint2[i,j, k, l].real)+'\n'
                f.write(s)
f.close()

'''
err = 1.0e-7
for i in range(norbs):
  for j in range(norbs):
      for k in range(norbs):
          for l in range(norbs):
             if(abs(newint2[i,j,k,l]) < 1.0e-14):
                continue;
             if (abs(newint2[i,j,k,l] - newint2[i,l,k,j]) >err 
                 or abs(newint2[i,j,k,l] - newint2[k,j,i,l]) >err
                 or abs(newint2[i,j,k,l] - newint2[k,l,i,j]) >err
                 or abs(newint2[i,j,k,l] - newint2[l,i,j,k]) >err
                 or abs(newint2[i,j,k,l] - newint2[l,k,j,i]) >err
                 or abs(newint2[i,j,k,l] - newint2[j,i,l,k]) >err
                 or abs(newint2[i,j,k,l] - newint2[j,k,l,i]) >err) :
                print i, j, k, l, newint2[i,j,k,l]
'''



