#include <Rcpp.h>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace Rcpp;

//' @title A Gaussian distribution random number generater
//' @description A Gaussian distribution random number generater using Rcpp
//' @param mu the mean of Gaussian distribution
//' @param sigma the standard deviation of Gaussian distribution
//' @return a random number form N(mu,sigma^2)
//' @useDynLib StatComp20056
//' @examples
//' \dontrun{
//' Guassian(0,1)
//' }
//' @export
// [[Rcpp::export]]
double Gaussian(double mu, double sigma)
 {
    const double epsilon = std::numeric_limits<double>::min();
    const double two_pi = 2.0*3.14159265358979323846;

    static double z0, z1;
    static bool generate;
    generate = !generate;

    if (!generate)
       return z1 * sigma + mu;

    double u1, u2;
    do
     {
       u1 = rand() * (1.0 / RAND_MAX);
       u2 = rand() * (1.0 / RAND_MAX);
     }
    while ( u1 <= epsilon );

    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
    return z0 * sigma + mu;
}

#include <Rcpp.h>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace Rcpp;

//' @title Laplace density function
//' @description Standard Laplace distribution density  using Rcpp
//' @param x a real number
//' @return the density at number x
//' @useDynLib StatComp20056
//' @examples
//' \dontrun{
//' lap(0)
//' }
//' @export
// [[Rcpp::export]]
double lap(double x){
  return exp(-fabs(x))/2.0;
}

#include <Rcpp.h>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace Rcpp;

//' @title Random Walk Metropolis Sampler
//' @description  Generate  a sample that obeys standard Laplace distribution by Random walk Metropolis, which uses proposal distribution normal (X_t,sigma^2). The data generation process records the acceptance rate.
//' @param n the size of sample
//' @param sigma the standard deviation of proposal Gaussian distribution
//' @param x0 the initial value of random walk
//' @return A sample that obeys standard Laplace distribution and acceptance rates.
//' @useDynLib StatComp20056
//' @examples
//' \dontrun{
//' a1<-as.numeric(metro(4000,1,10))
//' a1[1]
//' }
//' @export
// [[Rcpp::export]]
List metro( int n, double sigma, double x0){
	List out(n+1);
	double x[n+2];
    int k=0;
    double u;
    double y;
    x[0]=x0;
    for(int i=1;i<n+1;++i){
    	u = rand()*(1.0/RAND_MAX);
    	y = Gaussian(x[i-1],sigma);
    	if(u <= (lap(y) / lap(x[i-1]))){
    		x[i]=y;
    		k=k+1;
		}
		else{
			x[i]=x[i-1];
		}
		out[i]=x[i];	
	}
	out[0]=k*(1.0/n);
	return out;
}