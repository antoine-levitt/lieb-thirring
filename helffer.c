//Local Variables:
//compile-command: "gcc helffer.c -o helffer -lm -Wextra -Wall -O3 -ffast-math -std=c99 -g"
//End:

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double gam = .9;
double h = 0.01;

long long C(int k, int n){
	if(k== 0)
		return 1;
	else if(k == 1)
		return n;
	else if(k == 2)
		return (n*(n-1))/2;
	else if(k == 3)
		return (n*(n-1)*(n-2))/6;
	return 0;
}

// sum_{x_1 + x_2 ... + x_d = n} 1
int coeff(int n, int d){
	if(d == 1)
		return 1;
	if(d == 2)
		return n+1;
	if(d == 3)
		return ((n+1)*(n+2))/2;
	return 0;
}

double compute_r(int d){
	double sum = 0, num = 0;
	for(int i = 0;;i++){
		num = (1-(2*i+d)*h);
		if(num <= 0){
			break;
		}
		else{
			//sum += pow(num,gam)*C(d-1,i+d-1);
			/* printf("%g\n", h*h*h*sum); */
			sum += pow(num,gam)*coeff(i,d);
		}
	}
	return sum;
}

// OK
double lsc(int d){
	return tgamma(gam+1.0)/(pow(2*sqrt(M_PI),d)*tgamma(gam+1.0+((double) d)/2));
}

int main(){
	double lsc1D = lsc(1);
	double rhs1D = 0.5*tgamma(gam+1.5)*sqrt(M_PI)/tgamma(gam+2.0)*2;
	double limit1D = lsc1D*rhs1D;
	double r1 = compute_r(1);
	printf("R 1D %g\n",h*r1/limit1D - 1);

	double lsc2D = lsc(2);
	double rhs2D = M_PI/(gam+2.0);
	double limit2D = rhs2D * lsc2D;
	double r2 = compute_r(2);
	printf("R 2D %g\n",pow(h,2)*r2/limit2D - 1);

	double lsc3D = lsc(3); // OK
	double rhs3D = 2*M_PI*M_PI*tgamma(1.5)*tgamma(gam+2.5)/tgamma(gam+2.5 + 1.5);
	double limit3D = rhs3D * lsc3D;
	double r3 = compute_r(3);
	/* printf("%g\n",rhs3D); */
	/* printf("%g\n",lsc3D); */
	printf("R 3D %g\n",h*h*h*r3/limit3D);

	return 0;
}
