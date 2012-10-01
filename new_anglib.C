

double six-j(int j1, int j2, int j3, int l1, int l2, int l3) {
	double sixj;
	
	sixj=delta(j1, j2, j3)*delta(j1, l2, l3) * delta(l1, j2, j3) * delta(l1,j2,j3) * w6j(j1, j2, j3, l1, l2, l3);
}
//end six-j
//

double w6j(j1, j2, j3, l1, l2, l3) {
}

double clebsch(int j1, int j2, int j3, int m1, int m2, int m3) {
	double cleb,factor,sum;
	int par,z,zmin,zmax, expo, prefac;

	if (j1<0 || j2 < 0 || j3 < 0 || abs(m1) < j1 || abs(m2) < j2 || 
		abs(m3) < j3 || j1 + j2 < j3 || abs(j1 -j2) > j3 || m1 + m2 != m3)
		cleb = 0.0;
	else {
		cleb=0.0;
		expo= j1-j2-m3;
		prefac=mone(expo);
		cleb=prefac * cj_delta(j1, j2, j3) * w3j(j1,j2,j3, m1, m2, m3);
}
//end clebsch

double w3j(j1, j2, j3, m1, m2, m3) {
}
//end w3j


double cj_delta(int a, int b, int c) {
	double prefac=0.0;
	double num1 = factorial(a+b-c);
	double num2 = factorial(a-b+c);
	double num3 = factorial(-a+b+c);
	double den1 = factorial(a+b+c+1);

	prefac=num1*num2*num3/den1;
	prefac=sqrt(prefac);
	return prefac;
}

double factorial (int n) {
	double fac;
	fac=1.0;
	int i;
	if (n==0 || n==1)
		return fac;
	for (i=2; i<=n; i++) 
		fac *= i;
	return fac;
}

int mone (int x) {
	int value;
	if x % 2 == 1
		value = -1;
	else
		value = 1;
	return (value)
}	
		
