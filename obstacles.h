/* 
 * File:   obstacles.h
 * Author: irina
 *
 * Created on May 7, 2014, 3:47 PM
 */

#ifndef OBSTACLES_H
#define	OBSTACLES_H

#include <list>

#include "dual_representation.h"
#include "polygon.h"

struct Obstacle {
    SimplePolygon obstacle;
    Cell dual_cell;
};

class Obstacles {
public:
    Obstacles();
    virtual ~Obstacles();

    const std::list<Obstacle>& obstacles() const;

    /* assuming no vertical segments */
    void addObstacle(const SimplePolygon& ob);
    
    void intersect(const SimplePolygon& ch);
private:
    std::list<Obstacle> obstacles_;
};

#endif	/* OBSTACLES_H */

