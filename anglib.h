/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

/*This has been converted from original anglib.f90 code to a c code. The 
the license agreement is reproduced below

!    anglib.f90: angular momentum coupling coefficients in Fortran 90
!    Copyright (C) 1998  Paul Stevenson
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  
!    02110-1301,USA
*/
#ifndef ANGLIB_H
#define ANGLIB_H

double factorial(int n) ;

double binom(int n, int r);

double cleb_(int j1, int m1, int j2, int m2, int j, int m);

double sixj_(int a, int b, int c, int d, int e, int f);

double angdelta(int a,int b,int c);

double ninej_(int a, int b, int c, int d, int e, int f, int g, int h, int i);

#endif
