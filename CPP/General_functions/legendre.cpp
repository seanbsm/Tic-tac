
#include "legendre.h"

LegendrePolynomial::LegendrePolynomial(){
}

void LegendrePolynomial::findRoot(int d, float x, float &root){
	
	if (d==0){
		root = 1.;
	}
	else if (d==1){
		root = x;
	}
	else if (d>=2){
		/* Using recurrance relation for Legendre polynomials */
		float Pnext = 0;
		float Pprev = 1;
		float Pcurr = x;
		int c = 1;
		while (c<d){
			Pnext = ((2*c+1)*x*Pcurr - c*Pprev)/(c+1);	//recurrence relation
			Pprev = Pcurr;
			Pcurr = Pnext;
			c += 1;
			root = Pnext;
		}
	}
	else{
		throw std::invalid_argument("Function findRoot() in Legendre.cpp recieved negative integer.");
	}
}

void LegendrePolynomial::findRoot(int d, double x, double &root){
	
	if (d==0){
		root = 1.;
	}
	else if (d==1){
		root = x;
	}
	else if (d>=2){
		/* Using recurrance relation for Legendre polynomials */
		double Pnext = 0;
		double Pprev = 1;
		double Pcurr = x;
		int c = 1;
		while (c<d){
			Pnext = ((2*c+1)*x*Pcurr - c*Pprev)/(c+1);	//recurrence relation
			Pprev = Pcurr;
			Pcurr = Pnext;
			c += 1;
			root = Pnext;
		}
	}
	else{
		throw std::invalid_argument("Function findRoot() in Legendre.cpp recieved negative integer.");
	}
}

void LegendrePolynomial::findRoots(int L, std::vector<float> &z){
	
	/* List with roots per l, for each z[i] */
	std::vector<float> currRootList; (z.size());
	
	/* Set all roots to 1 for l=0 */
	currRootList.assign(z.size(), 1);
	roots_f.push_back(currRootList);
	
	/* Set all roots to z[i] for l=1 */
	roots_f.push_back(z);
	
	/* Use recurrance relation for Legendre polynomials
	 * to set remainder of roots */
	for (int l=1; l<L; l++){
		for (size_t i=0; i<z.size(); i++){
			float Pcurr = roots_f[l][i];
			float Pprev = roots_f[l-1][i];
			float root = ( (2*l+1)*z[i]*Pcurr - l*Pprev ) / (l+1);
			currRootList[i] = root;
		}
		roots_f.push_back(currRootList);
	}
}

void LegendrePolynomial::findRoots(int L, std::vector<double> &z){
	
	/* List with roots per l, for each z[i] */
	std::vector<double> currRootList; (z.size());
	
	/* Set all roots to 1 for l=0 */
	currRootList.assign(z.size(), 1);
	roots_d.push_back(currRootList);
	
	/* Set all roots to z[i] for l=1 */
	roots_d.push_back(z);
	
	/* Use recurrance relation for Legendre polynomials
	 * to set remainder of roots */
	for (int l=1; l<L; l++){
		for (size_t i=0; i<z.size(); i++){
			double Pcurr = roots_d[l][i];
			double Pprev = roots_d[l-1][i];
			double root = ( (2*l+1)*z[i]*Pcurr - l*Pprev ) / (l+1);
			currRootList[i] = root;
		}
		roots_d.push_back(currRootList);
	}
}

void LegendrePolynomial::fetchRoot(int d, int i, float &root){
	root = roots_f[d][i];
}

void LegendrePolynomial::fetchRoot(int d, int i, double &root){
	root = roots_d[d][i];
}
