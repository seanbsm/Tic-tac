
#include "gauss_legendre.h"

void gauss(float* x, float* w, int N){
	
	float p3	= 0;								// Legendre polynomial, P(xi), for i=3
	float eps = 3.e-7;
	int m	= (N+1)/2;
	
	for (int i=1; i<m+1; i++){
		float t  = cos(M_PI*(i-0.25)/(N+0.5));		// xi
		float t1 = 1;								// old xi
		float pp = 0;
		
		while (std::abs(t-t1) >= eps){
			float p1 = 1;							// Legendre polynomial, P(xi), for i=1
			float p2 = 0;							// Legendre polynomial, P(xi), for i=2
			
			for (int j=1; j<N+1; j++){
				p3 = p2;
				p2 = p1;
				p1 = ((2*j-1)*t*p2 - (j-1)*p3)/j;	// recurrence relation for Legendre polynomials
			}
			pp = N*(t*p1-p2)/(t*t-1);				// identity for P'(xi)
			t1 = t;
			t  = t1 - p1/pp; 						// Newton's method for finding roots
		}
		
		x[i-1] = -t;
		x[N-i] = t;
		
		w[i-1] = 2./((1-t*t)*pp*pp);
		w[N-i] = w[i-1];
	}
}
void gauss(double* x, double* w, int N){
	
	double p3	= 0;								// Legendre polynomial, P(xi), for i=3
	double eps = 3.e-16;
	int m	= (N+1)/2;
	
	for (int i=1; i<m+1; i++){
		double t  = cos(M_PI*(i-0.25)/(N+0.5));		// xi
		double t1 = 1;								// old xi
		double pp = 0;
		
		while (std::abs(t-t1) >= eps){
			double p1 = 1;							// Legendre polynomial, P(xi), for i=1
			double p2 = 0;							// Legendre polynomial, P(xi), for i=2
			
			for (int j=1; j<N+1; j++){
				p3 = p2;
				p2 = p1;
				p1 = ((2*j-1)*t*p2 - (j-1)*p3)/j;	// recurrence relation for Legendre polynomials
			}
			pp = N*(t*p1-p2)/(t*t-1);				// identity for P'(xi)
			t1 = t;
			t  = t1 - p1/pp; 						// Newton's method for finding roots
		}
		
		x[i-1] = -t;
		x[N-i] = t;
		
		w[i-1] = 2./((1-t*t)*pp*pp);
		w[N-i] = w[i-1];
	}
}

// changes range from (-1,1) to (0,inf)
void rangeChange_0_inf(float* x, float* w, float scale, int N){
	for (int i=0; i<N; i++){
		w[i] = scale*w[i]*M_PI/(4*(cos(M_PI*(x[i]+1)/4))*(cos(M_PI*(x[i]+1)/4)));
		x[i] = scale*tan(M_PI*(x[i]+1)/4);
	}
}
void rangeChange_0_inf(double* x, double* w, double scale, int N){
	for (int i=0; i<N; i++){
		w[i] = scale*w[i]*M_PI/(4*(cos(M_PI*(x[i]+1)/4))*(cos(M_PI*(x[i]+1)/4)));
		x[i] = scale*tan(M_PI*(x[i]+1)/4);
	}
}

// changes range from (-1,1) to (a,b)
void updateRange_a_b(float* x, float* w, float a, float b, int N){
	for (int i=0; i<N; i++){
		x[i] = 0.5*(b-a)*x[i] + 0.5*(b+a);
		w[i] = 0.5*(b-a)*w[i];
	}
}
void updateRange_a_b(double* x, double* w, double a, double b, int N){
	for (int i=0; i<N; i++){
		x[i] = 0.5*(b-a)*x[i] + 0.5*(b+a);
		w[i] = 0.5*(b-a)*w[i];
	}
}
