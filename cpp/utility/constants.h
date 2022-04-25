

#ifndef CPP_CONSTANTS_H
#define CPP_CONSTANTS_H

#include <math.h>
#include <vector>
#include <complex>
#include <string>
#include <utility>
#include "../matricies/matrix4d.h"
#include "../matricies/matrix3d.h"
#include "../matricies/matrix2d.h"


namespace constants{
    const double PI = 3.1415926535897932385;
    const double eps0 = 8.85419e-12;
    const double mu0 = 4 * PI * 1e-7;
    const double z0 = sqrt(mu0 / eps0);
    //const std::pair<double, double> cone = std::make_pair(1, 0);
	//std::complex<double> cone = std::complex<double>(1,0); //cam did this yell at him if you don't like it; I fixed it
    //const std::pair<double, double> cj = std::make_pair(0, 1);
	//const std::complex<double> cj = std::complex<double>(0,1); //cam did this too

}



#endif //CPP_CONSTANTS_H
