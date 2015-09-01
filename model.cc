//
//  model.cc
//  monkeys
//
//  Created by Irina Kostitsyna on 6/9/14.
//
//

#include "model.h"

#include <fstream>
#include <vector>

#include <CGAL/convex_hull_2.h>

#include "distribution.h"
//#include "precision.h"

using std::list;
using std::vector;
using std::map;

Model::Model(double alpha) :  error_(alpha) {
}

void Model::set_obstacles(const list<SimplePolygon>& obs) {
    for (list<SimplePolygon>::const_iterator it = obs.begin(); it != obs.end(); ++it) {
//        if (update_x_set_if_checks(*it))
            obstacles_.addObstacle(*it);
    }
}

void Model::read_obstacles_from_file(const std::string& file_name) {
    std::ifstream f(file_name);
    SimplePolygon ob;
    while (f >> ob) {
//        if (update_x_set_if_checks(ob))
            obstacles_.addObstacle(ob);
    }
    f.close();
}

void Model::write_obstacles_to_file(const std::string& file_name) {
    std::ofstream f(file_name);
    for (list<Obstacle>::const_iterator it = obstacles_.obstacles().begin(); it != obstacles_.obstacles().end(); ++it)
        f << it->obstacle;
    f.close();
}

void Model::set_first_monkey(double x, double y, double sigma) {
    set_monkey(0, x, y, sigma);
}

void Model::set_second_monkey(double x, double y, double sigma) {
    set_monkey(1, x, y, sigma);
}

void Model::trim_obstacles() {
    vector<Point> pts;
    for (list<ProbPolygon>::iterator it = monkeys_[0].begin(); it != monkeys_[0].end(); ++it) {
        pts.insert(pts.end(), it->poly.vertices_begin(), it->poly.vertices_end());
    }
    for (list<ProbPolygon>::iterator it = monkeys_[1].begin(); it != monkeys_[1].end(); ++it) {
        pts.insert(pts.end(), it->poly.vertices_begin(), it->poly.vertices_end());
    }

    SimplePolygon ch;
    CGAL::convex_hull_2(pts.begin(), pts.end(), std::back_inserter(ch));

    obstacles_.intersect(ch);
}

const list<ProbPolygon>& Model::monkey(int index) const {
    return monkeys_[index];
}

const Obstacles& Model::obstacles() const {
    return obstacles_;
}

void Model::set_monkey(int index, double x, double y, double sigma) {
    NormalDistribution distr(x, y, sigma);
    list<WeightedAnnulus> annuli;
    distr.approximate(error_, &annuli);

    double prev_weight(0);
    for (list<WeightedAnnulus>::reverse_iterator it = annuli.rbegin(); it != annuli.rend(); ++it) {
        list<ProbPolygon> prob_polys;
        bool generate_success;
        do {
            generate_success = true;
            prob_polys.clear();

            SimplePolygon p;
            generate_polygon(x, y, it->r_min, it->r_max, &p);
            PolygonSet p_set(p);
            for (list<Obstacle>::const_iterator ob_it = obstacles_.obstacles().begin();
                 ob_it != obstacles_.obstacles().end(); ++ob_it)
                p_set.difference(ob_it->obstacle);

            list<Polygon> poly_with_holes;
            p_set.polygons_with_holes(std::back_inserter(poly_with_holes));

            for (list<Polygon>::iterator p_it = poly_with_holes.begin();p_it != poly_with_holes.end(); ++p_it) {
                list<SimplePolygon> convex_polys;
                partition(*p_it, &convex_polys);
                map<double, Point> temp_xs;
                list<Segment> temp_segs;
//                if (!check_overlap_with_x_set(convex_polys, &temp_xs, &temp_segs)) {
//                    generate_success = false;
//                    break;
//                }

                for (list<SimplePolygon>::iterator conv_it = convex_polys.begin(); conv_it != convex_polys.end(); ++conv_it) {
                    ProbPolygon prob_poly;
                    prob_poly.weight = it->weight - prev_weight;
                    prob_poly.poly = *conv_it;
                    prob_polys.push_back(prob_poly);
                }
            }
        } while (!generate_success);

        prev_weight = it->weight;
        monkeys_[index].insert(monkeys_[index].end(), prob_polys.begin(), prob_polys.end());
    }
}

bool Model::check_overlap_with_x_set(const SimplePolygon& poly, std::map<double, Point>* new_xs, std::list<Segment>* new_segs) {
    for (SimplePolygon::Vertex_const_iterator it = poly.vertices_begin(); it != poly.vertices_end(); ++it) {
        double x = CGAL::to_double(it->x());
//        double x = it->x().get_d();
        if ((xs_.find(x) != xs_.end() && xs_[x] != *it) ||
            (new_xs->find(x) != new_xs->end() && (*new_xs)[x] != *it))
            return false;

        (*new_xs)[x] = *it;
    }

    for (SimplePolygon::Edge_const_iterator edge_it = poly.edges_begin(); edge_it != poly.edges_end(); ++edge_it) {
        for (list<Segment>::iterator seg_it = new_segs->begin(); seg_it != new_segs->end(); ++seg_it) {
            CGAL::Object obj = CGAL::intersection(*edge_it, *seg_it);
            Point p;
            if (CGAL::assign(p, obj)) {
                double x = CGAL::to_double(p.x());
//                double x = p.x().get_d();
                if ((xs_.find(x) != xs_.end() && xs_[x] != p) ||
                    (new_xs->find(x) != new_xs->end() && (*new_xs)[x] != p))
                    return false;

                (*new_xs)[x] = p;
            }
        }
        for (list<Segment>::iterator seg_it = segs_.begin(); seg_it != segs_.end(); ++seg_it) {
            CGAL::Object obj = CGAL::intersection(*edge_it, *seg_it);
            Point p;
            if (CGAL::assign(p, obj)) {
                double x = CGAL::to_double(p.x());
//                double x = p.x().get_d();
                if ((xs_.find(x) != xs_.end() && xs_[x] != p) ||
                    (new_xs->find(x) != new_xs->end() && (*new_xs)[x] != p))
                    return false;

                (*new_xs)[x] = p;
            }
        }

        new_segs->push_back(*edge_it);
    }

    return true;
}

bool Model::check_overlap_with_x_set(const list<SimplePolygon>& polys, map<double, Point>* new_xs, list<Segment>* new_segs) {
    for (list<SimplePolygon>::const_iterator poly_it = polys.begin(); poly_it != polys.end(); ++poly_it) {
        const SimplePolygon &p = *poly_it;
        for (SimplePolygon::Vertex_const_iterator it = p.vertices_begin(); it != p.vertices_end(); ++it) {
            double x = CGAL::to_double(it->x());
//            double x = it->x().get_d();
            if ((xs_.find(x) != xs_.end() && xs_[x] != *it) ||
                (new_xs->find(x) != new_xs->end() && (*new_xs)[x] != *it))
                return false;

            (*new_xs)[x] = *it;
        }

        for (SimplePolygon::Edge_const_iterator edge_it = p.edges_begin(); edge_it != p.edges_end(); ++edge_it) {
            for (list<Segment>::iterator seg_it = new_segs->begin(); seg_it != new_segs->end(); ++seg_it) {
                CGAL::Object obj = CGAL::intersection(*edge_it, *seg_it);
                Point p;
                if (CGAL::assign(p, obj)) {
                    double x = CGAL::to_double(p.x());
//                    double x = p.x().get_d();
                    if ((xs_.find(x) != xs_.end() && xs_[x] != p) ||
                        (new_xs->find(x) != new_xs->end() && (*new_xs)[x] != p))
                        return false;

                    (*new_xs)[x] = p;
                }
            }
            for (list<Segment>::iterator seg_it = segs_.begin(); seg_it != segs_.end(); ++seg_it) {
                CGAL::Object obj = CGAL::intersection(*edge_it, *seg_it);
                Point p;
                if (CGAL::assign(p, obj)) {
                    double x = CGAL::to_double(p.x());
//                    double x = p.x().get_d();
                    if ((xs_.find(x) != xs_.end() && xs_[x] != p) ||
                        (new_xs->find(x) != new_xs->end() && (*new_xs)[x] != p))
                        return false;
                    
                    (*new_xs)[x] = p;
                }
            }
            
            new_segs->push_back(*edge_it);
        }
    }

    return true;
}

void Model::update_x_set(const map<double, Point>& new_xs, const list<Segment>& new_segs) {
    xs_.insert(new_xs.begin(), new_xs.end());
    segs_.insert(segs_.end(), new_segs.begin(), new_segs.end());
}

bool Model::update_x_set_if_checks(const SimplePolygon& p) {
    map<double, Point> temp_xs;
    list<Segment> temp_segs;
    if (!check_overlap_with_x_set(p, &temp_xs, &temp_segs))
        return false;

    update_x_set(temp_xs, temp_segs);
    return true;
}
