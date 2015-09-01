/* 
 * File:   obstacles.cc
 * Author: irina
 * 
 * Created on May 7, 2014, 3:47 PM
 */

#include "obstacles.h"

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/convex_hull_2.h>

using std::list;

Obstacles::Obstacles() {
}

Obstacles::~Obstacles() {
}

const list<Obstacle>& Obstacles::obstacles() const {
    return obstacles_;
}

void Obstacles::addObstacle(const SimplePolygon& ob) {
    Obstacle obstacle;
    CGAL::convex_hull_2(ob.vertices_begin(), ob.vertices_end(), std::back_inserter(obstacle.obstacle));

    obstacle.dual_cell.Initialize(obstacle.obstacle);

    obstacles_.push_back(obstacle);
}

void Obstacles::intersect(const SimplePolygon& ch) {
    std::list<Obstacle> trimmed_obstacles;
    for (list<Obstacle>::iterator it = obstacles_.begin(); it != obstacles_.end(); ++it)
        if (CGAL::do_intersect(SimplePolygon(it->obstacle), ch))
            trimmed_obstacles.push_back(*it);

    std::swap(obstacles_, trimmed_obstacles);
}
