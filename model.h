//
//  model.h
//  monkeys
//
//  Created by Irina Kostitsyna on 6/9/14.
//
//

#ifndef MODEL_H
#define MODEL_H


#include <string>
#include <list>
#include <map>

#include <CGAL/Random.h>

#include "geometry.h"
#include "obstacles.h"
#include "polygon.h"

class Model {
public:
    explicit Model(double alpha);

    void set_obstacles(const std::list<SimplePolygon>& obs);
    void read_obstacles_from_file(const std::string& file_name);
    void write_obstacles_to_file(const std::string& file_name);
    void set_first_monkey(double x, double y, double sigma);
    void set_second_monkey(double x, double y, double sigma);

    void trim_obstacles();

    const std::list<ProbPolygon>& monkey(int index) const;
    const Obstacles& obstacles() const;
private:
    void set_monkey(int index, double x, double y, double sigma);
    bool check_overlap_with_x_set(const SimplePolygon& p, std::map<double, Point>* new_xs, std::list<Segment>* new_segs);
    bool check_overlap_with_x_set(const std::list<SimplePolygon>& polys, std::map<double, Point>* new_xs, std::list<Segment>* new_segs);
    void update_x_set(const std::map<double, Point>& new_xs, const std::list<Segment>& new_segs);
    bool update_x_set_if_checks(const SimplePolygon& p);
private:
    Obstacles obstacles_;
    std::list<ProbPolygon> monkeys_[2];

    double error_;

    std::map<double, Point> xs_;
    std::list<Segment> segs_;
};

#endif // MODEL_H
