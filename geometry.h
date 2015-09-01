/* 
 * File:   geometry.h
 * Author: irina
 *
 * Created on May 7, 2014, 11:01 AM
 */

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* pi             */
#define M_PI_2      1.57079632679489661923132169163975144   /* pi/2           */
#define M_PI_4      0.785398163397448309615660845819875721  /* pi/4           */
#endif

#ifndef GEOMETRY_H
#define	GEOMETRY_H

//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Cartesian.h>
//#include <CGAL/Lazy_exact_nt.h>
//#include <CGAL/Quotient.h>
////#include <CGAL/Gmpzf.h>
//#include <CGAL/Gmpq.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Bbox_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Vector_2.h>

#include <CGAL/Random.h>

//#include <gmpxx.h>
//#include <CGAL/mpq_class.h>

#include "errors.h"
#include "precision.h"

extern CGAL::Random rand___;

//typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > NT;
//typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::Gmpzf> > NT;
//typedef CGAL::Quotient<CGAL::MP_Float> NT;
//typedef CGAL::Gmpzf NT;
//typedef CGAL::Gmpq NT;
//typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::Gmpq> > NT;
//typedef CGAL::Cartesian<NT> Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::FT NT;
//typedef mpq_class NT;
//typedef CGAL::Simple_cartesian<NT> Kernel;

typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel;
typedef CGAL::Point_2<Kernel> Point;
typedef CGAL::Point_2<InexactKernel> InexactPoint;
//typedef CGAL::Segment_2<Kernel> Segment;
//typedef CGAL::Segment_2<InexactKernel> InexactSegment;
typedef CGAL::Line_2<Kernel> Line;
typedef CGAL::Line_2<InexactKernel> InexactLine;
typedef CGAL::Direction_2<Kernel> Direction;
typedef CGAL::Vector_2<Kernel> Vector;
typedef CGAL::Bbox_2 Bbox;
typedef CGAL::Iso_rectangle_2<Kernel> IsoRectangle;
typedef CGAL::Ray_2<Kernel> Ray;

template <class K>
class BaseSegment : public CGAL::Segment_2<K> {
public:
    BaseSegment() : CGAL::Segment_2<K>() {};
    BaseSegment(const CGAL::Segment_2<K>& s) : CGAL::Segment_2<K>(s) {};
    BaseSegment(const CGAL::Point_2<K>& p1, const CGAL::Point_2<K>& p2) : CGAL::Segment_2<K>(p1, p2) {};
    bool operator==(const BaseSegment<K>& s) const {
        return (this->source() == s.source() && this->target() == s.target()) ||
               (this->source() == s.target() && this->target() == s.source());
    }
    bool operator!=(const BaseSegment<K>& s) const {
        return !(this->operator==(s));
    };
};

typedef BaseSegment<Kernel> Segment;
typedef BaseSegment<InexactKernel> InexactSegment;

template <class K>
class NonVerticalLine : public CGAL::Line_2<K> {
 public:
    NonVerticalLine() : CGAL::Line_2<K>() {};
    NonVerticalLine(const CGAL::Point_2<K>& p1, const CGAL::Point_2<K>& p2) : CGAL::Line_2<K>(p1, p2) {
        if (CGAL::Line_2<K>::is_vertical())
            throw boost::enable_current_exception(integral_error()) << err_description(__FILE__) << code_line(__LINE__);
    };
    NonVerticalLine(const CGAL::Point_2<K>& p, const CGAL::Direction_2<K>& dir) : CGAL::Line_2<K>(p, dir) {
        if (CGAL::Line_2<K>::is_vertical())
            throw boost::enable_current_exception(integral_error()) << err_description(__FILE__) << code_line(__LINE__);
    };
    NonVerticalLine(const CGAL::Point_2<K>& dual) : CGAL::Line_2<K>(-dual.x(), 1, dual.y()) {};
    virtual ~NonVerticalLine() {};

    NT slope() const {
        return -a()/b();
    };
    NT y_intersept() const{
        return this->y_at_x(0);
    };

    void dual(CGAL::Point_2<K>* p) const {
        *p = CGAL::Point_2<K>(slope(), -y_intersept());
    };
 private:
    typename K::FT a() const { return CGAL::Line_2<K>::a(); };
    typename K::FT b() const { return CGAL::Line_2<K>::b(); };
};

//bool operator==(const Segment& s1, const Segment& s2);
//bool operator!=(const Segment& s1, const Segment& s2);
bool operator<(const NonVerticalLine<Kernel>& l1, const NonVerticalLine<Kernel>& l2);
bool operator<(const NonVerticalLine<InexactKernel>& l1, const NonVerticalLine<InexactKernel>& l2);
bool operator<(const Line& l1, const Line& l2);
bool operator<(const InexactLine& l1, const InexactLine& l2);

namespace CGAL {
    template<class K>
    bool do_intersect(const BaseSegment<K>& s1, const BaseSegment<K>& s2) {
        return do_intersect(dynamic_cast<const CGAL::Segment_2<K>&>(s1), dynamic_cast<const CGAL::Segment_2<K>&>(s2));
    }
}

#endif	/* GEOMETRY_H */

