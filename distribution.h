/* 
 * File:   distribution.h
 * Author: irina
 *
 * Created on June 2, 2014, 2:09 PM
 */

#ifndef DISTRIBUTION_H
#define	DISTRIBUTION_H

#include <list>

#include "geometry.h"
#include "polygon.h"

struct WeightedAnnulus {
    double r_min;
    double r_max;
    double weight;
};

class NormalDistribution {
public:
    NormalDistribution(double center_x, double center_y, double sigma);
    virtual ~NormalDistribution();

    void approximate(double error, std::list<WeightedAnnulus>* annuli);
private:
    double center_x_;
    double center_y_;
    double sigma_;
};

#endif	/* DISTRIBUTION_H */

