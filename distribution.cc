/* 
 * File:   distribution.cc
 * Author: irina
 * 
 * Created on June 2, 2014, 2:09 PM
 */

#include "distribution.h"

#include <vector>

#include "errors.h"
#include "precision.h"

using std::list;
using std::vector;

NormalDistribution::NormalDistribution(double center_x, double center_y, double sigma)
    : center_x_(center_x), center_y_(center_y), sigma_(sigma) {
}

NormalDistribution::~NormalDistribution() {
}

//void NormalDistribution::approximate(double error, std::list<WeightedAnnulus>* annuli) {
//    double eps = error/2;
//
//    vector<double> rs;
//    for (int i = -1; i < log(1/eps)/log(1+eps); ++i) {
//        double r = sigma_*sqrt(2*(log(1/eps)-log(1+eps)*i));
//        rs.push_back(r);
//    }
//
//    for (int i = 0; i < rs.size()-1; ++i) {
//        WeightedAnnulus annulus;
//        annulus.r_min = rs[i+1];
//        annulus.r_max = rs[i];
//        annulus.weight = (eps*pow(1.0+eps,i+1))/(M_PI*2.0*sigma_*sigma_); // == F(r_min)
//        annuli->push_back(annulus);
//    }
//}

void NormalDistribution::approximate(double error, std::list<WeightedAnnulus>* annuli) {
    double k = ceil(1.0/(exp(error)-1.0));

    for (int i = 0; i < k; ++i) {
        double r_odd = 2.0*sigma_*sqrt(0.5*log((k*(k+1))/pow(k-i,2)));  // r_{2i+1}
        double rho_minus = sigma_*sqrt(2.0*log((2.0*k*(k+1))/((2*k+1-2*i)*(k-i))));
        double rho_plus = sigma_*sqrt(2.0*log((2.0*k*(k+1))/((2*k-1-2*i)*(k-i))));

        WeightedAnnulus annulus1;
        annulus1.r_min = rho_minus;
        annulus1.r_max = r_odd;
        annulus1.weight = ((k-i)*(k+1-i))/(2.0*k*(k+1)*M_PI*sigma_*sigma_); // == F(r_{2i})
        annuli->push_back(annulus1);

        WeightedAnnulus annulus2;
        annulus2.r_min = r_odd;
        annulus2.r_max = rho_plus;
        annulus2.weight = (k-i)*(k-i)/(2.0*k*(k+1)*M_PI*sigma_*sigma_); // == F(r_{2i+1})
        annuli->push_back(annulus2);
    }
}

