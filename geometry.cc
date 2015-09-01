/* 
 * File:   geometry.cc
 * Author: irina
 * 
 * Created on May 7, 2014, 11:01 AM
 */

#include "geometry.h"

#include "errors.h"

CGAL::Random rand___;

/*
Segment::Segment() : CGAL::Segment_2<Kernel>() {
}

Segment::Segment(const Point& p1, const Point& p2) : CGAL::Segment_2<Kernel>(p1, p2) {
}

bool Segment::operator==(const Segment& s) const {
    return (this->source() == s.source() && this->target() == s.target()) ||
           (this->source() == s.target() && this->target() == s.source());

}

bool Segment::operator!=(const Segment& s) const {
    return !(this->operator==(s));
}
*/
/*
NonVerticalLine::NonVerticalLine() : Line() {
}

NonVerticalLine::NonVerticalLine(const Point& p1, const Point& p2) : Line(p1, p2) {
  if (is_vertical())
    throw integral_error() << err_description(__FILE__) << code_line(__LINE__);
}

NonVerticalLine::NonVerticalLine(const Point& p, const Direction& dir) : Line(p, dir) {
  if (is_vertical())
    throw integral_error() << err_description(__FILE__) << code_line(__LINE__);
}

NonVerticalLine::NonVerticalLine(const Point& dual) : Line(-dual.x(), 1, dual.y()) {
}

NonVerticalLine::~NonVerticalLine() {
}

NT NonVerticalLine::slope() const {
  return -a()/b();
}

NT NonVerticalLine::y_intersept() const {
  return this->y_at_x(0);
}

bool NonVerticalLine::operator<(const NonVerticalLine& l) const {
  return this->slope() < l.slope() || (this->slope() == l.slope() && this->y_intersept() < l.y_intersept());
}

void NonVerticalLine::dual(Point* p) const {
  *p = Point(slope(), -y_intersept());
}
*/

//bool operator==(const Segment& s1, const Segment& s2) {
//    return (s1.source() == s2.source() && s1.target() == s2.target()) ||
//           (s1.source() == s2.target() && s1.target() == s2.source());
//}
//
//bool operator!=(const Segment& s1, const Segment& s2) {
//    return !(s1 == s2);
//}

bool operator<(const NonVerticalLine<Kernel>& l1, const NonVerticalLine<Kernel>& l2) {
    return l1.slope() < l2.slope() || (l1.slope() == l2.slope() && l1.y_intersept() < l2.y_intersept());
}

bool operator<(const NonVerticalLine<InexactKernel>& l1, const NonVerticalLine<InexactKernel>& l2) {
    return l1.slope() < l2.slope() - PRECISION || (CGAL::abs(l1.slope() - l2.slope()) < PRECISION && l1.y_intersept() < l2.y_intersept() - PRECISION);
}

bool operator<(const Line& l1, const Line& l2) {
  if (l1.is_vertical() & l2.is_vertical())
    return -l1.c()/l1.a() < -l2.c()/l2.a();
  else if (l1.is_vertical())
    return false;
  else if (l2.is_vertical())
    return true;

  NonVerticalLine<Kernel> nv_l1(l1.point(0), l1.point(1));
  NonVerticalLine<Kernel> nv_l2(l2.point(0), l2.point(1));
  return nv_l1 < nv_l2;
}

bool operator<(const InexactLine& l1, const InexactLine& l2) {
    if (l1.is_vertical() & l2.is_vertical())
        return -l1.c()/l1.a() < -l2.c()/l2.a() - PRECISION;
    else if (l1.is_vertical())
        return false;
    else if (l2.is_vertical())
        return true;
    
    NonVerticalLine<InexactKernel> nv_l1(l1.point(0), l1.point(1));
    NonVerticalLine<InexactKernel> nv_l2(l2.point(0), l2.point(1));
    return nv_l1 < nv_l2;
}
