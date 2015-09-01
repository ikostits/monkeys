/* 
 * File:   polygon.h
 * Author: irina
 *
 * Created on May 2, 2014, 1:59 PM
 */

#ifndef POLYGON_H
#define	POLYGON_H

#include <list>
#include <vector>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_set_2.h>

#include "geometry.h"

typedef CGAL::Triangle_2<Kernel> Triangle;

typedef CGAL::Polygon_2<Kernel> SimplePolygon;

typedef CGAL::Polygon_with_holes_2<Kernel> Polygon;

typedef CGAL::Polygon_set_2<Kernel> PolygonSet;
/*
class SimplePolygon : public CGAL::Polygon_2<Kernel> {
 public:
  SimplePolygon();
  template <typename Iterator>
  SimplePolygon(Iterator begin, Iterator end) : CGAL::Polygon_2<Kernel>(begin, end) {}
  virtual ~SimplePolygon();
};
*/
/*
class Polygon : public CGAL::Polygon_with_holes_2<Kernel> {
 public:
  Polygon();
  Polygon(const SimplePolygon& p);
  virtual ~Polygon();

  void triangulate(std::list<Triangle>* triangles) const;
};
*/
/*
class ConvexPolygon : public SimplePolygon {
    //friend void swap(ConvexPolygon& c1, ConvexPolygon& c2);
public:
    ConvexPolygon();
    ConvexPolygon(const SimplePolygon& poly);
    
    template <typename Iterator>
    ConvexPolygon(Iterator begin, Iterator end) : SimplePolygon(begin, end) {
        correct();
    }
    
    virtual ~ConvexPolygon();
    
    void correct();
    
    Vertex_const_iterator left_top_vertex() const;
    Vertex_const_iterator right_bottom_vertex() const;
    
    void intersect(const ConvexPolygon& cell, ConvexPolygon* out) const;
    void get_lower_envelope(std::vector<Point>* out) const;
    void get_upper_envelope(std::vector<Point>* out) const;
    
    bool is_degenerate() const;
};
*/

SimplePolygon::Vertex_const_iterator left_top_vertex(const SimplePolygon& convex_poly);
SimplePolygon::Vertex_const_iterator right_bottom_vertex(const SimplePolygon& convex_poly);

void intersect_convex_polygons(const SimplePolygon& poly1, const SimplePolygon& poly2, SimplePolygon* out);
void get_lower_envelope(const SimplePolygon& convex_poly, std::vector<Point>* out);
void get_upper_envelope(const SimplePolygon& convex_poly, std::vector<Point>* out);
bool is_degenerate(const SimplePolygon& convex_poly);
bool do_intersect_interiors(const SimplePolygon& poly1, const SimplePolygon& poly2);

struct ProbPolygon {
    SimplePolygon poly;
    double weight;
};

struct ProbPolygonInexact {
    CGAL::Polygon_2<InexactKernel> poly;
    double weight;
};

/*
namespace CGAL {
    bool do_intersect(const ConvexPolygon& p1, const ConvexPolygon& p2);
    template<typename OutputIterator >
    OutputIterator intersection(const ConvexPolygon &p1, const ConvexPolygon &p2, OutputIterator intersections) {
        return intersection(dynamic_cast<const SimplePolygon&>(p1), dynamic_cast<const SimplePolygon&>(p2), intersections);
    }
}
*/

void triangulate(const Polygon& poly, std::list<Triangle>* triangles);
void partition(const Polygon& poly, std::list<SimplePolygon>* convex_polys);

void generate_polygon(double center_x, double center_y, double r_min, double r_max, SimplePolygon* out);

NT area(const Polygon& poly);

void print_to_file(const SimplePolygon& poly, std::string file_name);
void print_to_file(const Polygon& poly, std::string file_name);
void print_to_file(const SimplePolygon& poly1, const SimplePolygon& poly2, const std::list<SimplePolygon>& obs, std::string file_name);
void print_to_file(const std::list<SimplePolygon>& poly1, const std::list<SimplePolygon>& poly2, const std::list<SimplePolygon>& obs, std::string file_name);
void print_to_file(double x1, double y1,
                   double x2, double y2,
                   double sigma,
                   const std::list<ProbPolygon>& poly1,
                   const std::list<ProbPolygon>& poly2,
                   double opacity,
                   const std::list<SimplePolygon>& obs,
                   std::string file_name);

#endif	/* POLYGON_H */

