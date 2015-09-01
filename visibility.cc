/*
 * File:   visibility.cc
 * Author: irina
 *
 * Created on May 7, 2014, 1:15 PM
 */

#include "visibility.h"

#include <CGAL/IO/io.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>
#include <CGAL/number_utils.h>

#include <set>
#include <list>
#include <iostream>

#include "geometry.h"
#include "errors.h"
#include "obstacles.h"
#include "precision.h"

using std::list;
using std::set;
using std::swap;
using std::vector;

namespace visibility_intl {
    
    bool intersects_in_order(const Line& l, const Segment (&segs)[4]) {
        struct LineSegmentIntersection {
            bool is_point;  // if false then segment
            Point p;
            Segment s;
        } intersections[4];
        
        for (int i = 0; i < 4; ++i) {
            if (!CGAL::do_intersect(l, segs[i]))
                return false;
        }
        
        for (int i = 0; i < 4; ++i) {
            CGAL::Object tmp = CGAL::intersection(l, segs[i]);
            if (CGAL::assign(intersections[i].s, tmp)) {
                intersections[i].is_point = false;
            } else if (CGAL::assign(intersections[i].p, tmp)) {
                intersections[i].is_point = true;
            } else {
                return false;
            }
        }
        
        set<Direction> directions;
        for (int i = 0; i < 4; ++i)
            for (int j = i+1; j < 4; ++j) {
                if (intersections[i].is_point && intersections[j].is_point) {
                    if (intersections[i].p == intersections[j].p) {
                        vector<Point> neg_side;
                        vector<Point> pos_side;
                        for (int k = 0; k < 2; ++k) {
                            if (l.has_on_negative_side(segs[i].vertex(k))) {
                                neg_side.insert(neg_side.begin(), segs[i].vertex(k));
                            } else if (l.has_on_positive_side(segs[i].vertex(k))) {
                                pos_side.insert(pos_side.begin(), segs[i].vertex(k));
                            }
                            
                            if (l.has_on_negative_side(segs[j].vertex(k))) {
                                neg_side.push_back(segs[j].vertex(k));
                            } else if (l.has_on_positive_side(segs[j].vertex(k))) {
                                pos_side.push_back(segs[j].vertex(k));
                            }
                        }
                        
                        if (neg_side.size() < 2 && pos_side.size() < 2)
                            return false;
                        
                        Point p = l.point(0);
                        if (intersections[i].p == p)
                            p = l.point(1);
                        
                        vector<Point>* cur_side;
                        if (neg_side.size() >= 2)
                            cur_side = &neg_side;
                        else
                            cur_side = &pos_side;
                        
                        if ((CGAL::left_turn(p, intersections[i].p, (*cur_side)[0]) &&
                             CGAL::left_turn((*cur_side)[0], intersections[i].p, (*cur_side)[1])) ||
                            (CGAL::right_turn(p, intersections[i].p, (*cur_side)[0]) &&
                             CGAL::right_turn((*cur_side)[0], intersections[i].p, (*cur_side)[1])))
                            directions.insert(Direction(Segment(p, intersections[i].p)));
                        else if (!CGAL::collinear((*cur_side)[0], intersections[i].p, (*cur_side)[1]))
                            directions.insert(Direction(Segment(intersections[i].p, p)));
                    } else {
                        directions.insert(Direction(Segment(intersections[i].p, intersections[j].p)));
                    }
                } else if (intersections[i].is_point != intersections[j].is_point) {
                    Point p = intersections[i].is_point ? intersections[i].p : intersections[j].p;
                    Segment s = !intersections[i].is_point ? intersections[i].s : intersections[j].s;
                    if (!s.has_on(p) || s.vertex(1) == p || s.vertex(0) == p) {
                        Point seg_middle = s.vertex(0) + (s.vertex(1) - s.vertex(0))/2;
                        if (intersections[i].is_point) {
                            directions.insert(Direction(Segment(p, seg_middle)));
                        } else {
                            directions.insert(Direction(Segment(seg_middle, p)));
                        }
                    }
                } else {
                    if (!CGAL::do_intersect(intersections[i].s, intersections[j].s))
                        directions.insert(Direction(CGAL::Segment_2<Kernel>(intersections[i].s.vertex(0), intersections[j].s.vertex(0))));
                }
            }
        
        return directions.size() == 1;
    }
    
    bool intersects_in_order(const InexactLine& l, const InexactSegment (&segs)[4]) {
        struct LineSegmentIntersection {
            bool is_point;  // if false then segment
            InexactPoint p;
            bool end_point;
            InexactPoint p_other;
            InexactSegment s;
        } intersections[4];
        
        for (int i = 0; i < 4; ++i) {
            InexactKernel::FT dist1 = CGAL::squared_distance(segs[i].vertex(0), l);
            InexactKernel::FT dist2 = CGAL::squared_distance(segs[i].vertex(1), l);
            if (dist1 < SQUARED_PRECISION && dist2 < SQUARED_PRECISION) {
                intersections[i].is_point = false;
                intersections[i].s = segs[i];
            } else if (dist1 < SQUARED_PRECISION) {
                intersections[i].is_point = true;
                intersections[i].end_point = true;
                intersections[i].p = segs[i].vertex(0);
                intersections[i].p_other = segs[i].vertex(1);
            } else if (dist2 < SQUARED_PRECISION) {
                intersections[i].is_point = true;
                intersections[i].end_point = true;
                intersections[i].p = segs[i].vertex(1);
                intersections[i].p_other = segs[i].vertex(0);
            } else if (!CGAL::do_intersect(l, segs[i])) {
                return false;
            } else {
                intersections[i].is_point = true;
                intersections[i].end_point = false;
                CGAL::Object tmp = CGAL::intersection(l, segs[i]);
                if (!CGAL::assign(intersections[i].p, tmp))
                    return false;
            }
        }
        
        Bbox box = CGAL::bbox_2(segs, segs + 4);
        InexactPoint pt_ref;
        if (l.is_vertical())
            pt_ref = InexactPoint(l.x_at_y(box.ymin()-1), box.ymin()-1);
        else
            pt_ref = InexactPoint(box.xmin()-1, l.y_at_x(box.xmin()-1));
        
        set<bool> directions;
        for (int i = 0; i < 3; ++i)
            for (int j = i+1; j < 4; ++j) {
                if (i == 1 && j == 2 && segs[i] == segs[j])
                    continue;

                if (intersections[i].is_point && intersections[j].is_point) {
                    if (intersections[i].p == intersections[j].p) {
                        if (!intersections[i].end_point || !intersections[j].end_point) {
//                            std::cout << std::setprecision(10) << "segs[0]" << segs[0] << std::endl
//                                                               << "segs[1]" << segs[1] << std::endl
//                                                               << "segs[2]" << segs[2] << std::endl
//                                                               << "segs[3]" << segs[3] << std::endl << std::endl;
//                            std::cout << "i=" << i << "\tj=" << j << std::endl << std::endl;
                            throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);
                        }
                        
                        if ((l.has_on_positive_side(intersections[i].p_other) && l.has_on_negative_side(intersections[j].p_other)) ||
                            (l.has_on_negative_side(intersections[i].p_other) && l.has_on_positive_side(intersections[j].p_other)))
                            return false;
                        
                        if ((CGAL::left_turn(pt_ref, intersections[i].p, intersections[i].p_other) &&
                             CGAL::left_turn(intersections[i].p_other, intersections[i].p, intersections[j].p_other)) ||
                            (CGAL::right_turn(pt_ref, intersections[i].p, intersections[i].p_other) &&
                             CGAL::right_turn(intersections[i].p_other, intersections[i].p, intersections[j].p_other)))
                            directions.insert(true);
                        else if (!CGAL::collinear(intersections[i].p_other, intersections[i].p, intersections[j].p_other))
                            directions.insert(false);
                    } else {
                        directions.insert(CGAL::collinear_are_strictly_ordered_along_line(pt_ref, intersections[i].p, intersections[j].p));
                    }
                } else if (intersections[i].is_point != intersections[j].is_point) {
                    InexactPoint &p = intersections[i].is_point ? intersections[i].p : intersections[j].p;
                    InexactSegment &s = !intersections[i].is_point ? intersections[i].s : intersections[j].s;
                    if (!CGAL::collinear_are_strictly_ordered_along_line(s.vertex(0), p, s.vertex(1))) {
                        if (intersections[i].is_point) {
                            directions.insert(CGAL::collinear_are_ordered_along_line(pt_ref, p, s.vertex(0)) &&
                                              CGAL::collinear_are_ordered_along_line(pt_ref, p, s.vertex(1)));
                        } else {
                            directions.insert(CGAL::collinear_are_ordered_along_line(pt_ref, s.vertex(0), p) &&
                                              CGAL::collinear_are_ordered_along_line(pt_ref, s.vertex(1), p));
                        }
                    }
                } else {
                    if (!CGAL::do_intersect(intersections[i].s, intersections[j].s))
                        directions.insert(CGAL::collinear_are_ordered_along_line(pt_ref,intersections[i].s.vertex(0), intersections[j].s.vertex(0)));
                }
            }
        
        return directions.size() == 1;
    }
/*
    mpf_class power(const mpf_class& vt, int i) {
        mpf_class result(1.0);
        if (i == 0) {
            return result;
        } else if (i > 0) {
            for (int j = 0; j < i; ++j)
                result *= vt;
        } else {
            for (int j = 0; j < -i; ++j)
                result /= vt;
        }
        
        return result;
    }
    
    mpf_class logarithm(const mpf_class& vt) {
        mpfr_t mpfr_result, source;
        mpfr_init(mpfr_result);
        mpfr_init_set_f(source, vt.get_mpf_t(), GMP_RNDN);
        mpfr_log(mpfr_result, source, GMP_RNDN);

        mpf_t result;
        mpf_init(result);
        mpfr_get_f(result, mpfr_result, GMP_RNDN);

        return mpf_class(result);
//        return mpf_class(log(vt.get_d()));
//        return mpf_class(log(CGAL::to_double(vt)));
    }

    mpf_class to_mpf(const NT& x) {
//          mpf_t num_mpf, den_mpf;
//          mpf_init(num_mpf);
//          mpf_init(den_mpf);
//          mpf_set_q(num_mpf, x.exact().numerator().mpq());
//          mpf_set_q(den_mpf, x.exact().denominator().mpq());
//          mpf_class num(num_mpf);
//          mpf_class den(den_mpf);
        mpf_t result;
        mpf_init(result);
        mpf_set_q(result, x.exact().mpq());
//        mpf_set_q(result, x.mpq());
        return mpf_class(result);
    }
*/
    mpfr_class XiXjXj(const mpfr_class& alpha,
                     const mpfr_class (&A)[2], const mpfr_class (&B)[2],
                     const mpfr_class& ai, const mpfr_class& bi, const mpfr_class& ci,
                     const mpfr_class& aj, const mpfr_class& bj, const mpfr_class& cj) {
        mpfr_class pr(0);
        const mpfr_class& A1 = A[0];
        const mpfr_class& A2 = A[1];
        const mpfr_class& B1 = B[0];
        const mpfr_class& B2 = B[1];
        if (bi.is_zero() && bj.is_zero()) {
            pr = (alpha*ci*cj*cj*(A1*alpha-A2*alpha+2*B1-2*B2))/(2*ai*aj*aj);
        } else if (bj.is_zero()) {  // bi != 0
            if ((A1*ai-B1*bi+ci).is_zero() && (A2*ai-B2*bi+ci).is_zero()) {
                if ((ai+bi*alpha).is_zero()) {
//                    pr = ((B1-B2)*cj*cj*((B1+B2)*bi-2*ci))/(4*aj*aj*bi);
                    pr = (A1*A1-A2*A2)*ai*ai*cj*cj/(4*aj*aj*bi*bi);
                } else {
//                    pr = (alpha*cj*cj*(bi*(A2*(A2*alpha+2*B2)-2*A1*B1+A1*A1*(-alpha))+2*(A1-A2)*ci))/(4*aj*aj*bi);
                    pr = -(A1*A1-A2*A2)*alpha*cj*cj*(2*ai+alpha*bi)/(4*aj*aj*bi);
                }
            } else {
                mpfr_class tmp = abs(ai+alpha*bi);
                pr = -logarithm(tmp)*(cj*cj*((A1-A2)*ai+(B2-B1)*bi)*((A1+A2)*ai-(B1+B2)*bi+2*ci))/(2*aj*aj*bi*bi)+
                (alpha*cj*cj*(2*(A1*A1-A2*A2)*ai+bi*(alpha*(A2*A2-A1*A1)+4*A2*B2-4*A1*B1)+4*(A1-A2)*ci))/(4*aj*aj*bi);
            }
        } else if (bi.is_zero()) {  // bj != 0
            if ((A1*aj-B1*bj+cj).is_zero() && (A2*aj-B2*bj+cj).is_zero()) {
                if ((aj+bj*alpha).is_zero()) {
                    pr = (power(A2,3)-power(A1,3))*aj*aj*ci/(6*ai*bj*bj);
                } else {
                    pr = (power(A1,3)-power(A2,3))*alpha*ci*(2*aj+alpha*bj)/(6*ai*bj);
                }
            } else {
                mpfr_class tmp = abs(aj+alpha*bj);
                pr = (ci*logarithm(tmp)*(2*aj*(A1*A1*(cj-B1*bj)+A2*A2*(B2*bj-cj))+(power(A1,3)-power(A2,3))*aj*aj+A1*power(cj-B1*bj,2)-A2*power(cj-B2*bj,2)))/(ai*bj*bj);
                pr += (ci*(3*aj*(A2*(2*A2*alpha*bj*(cj-B2*bj)+A2*A2*alpha*alpha*bj*bj-2*power(cj-B2*bj,2))+2*A1*A1*alpha*bj*(B1*bj-cj)+2*A1*power(cj-B1*bj,2)+power(A1,3)*(-alpha*alpha)*bj*bj)+2*aj*aj*(A2*A2*(bj*(2*A2*alpha+3*B2)-3*cj)+3*A1*A1*(cj-B1*bj)-2*power(A1,3)*alpha*bj)+2*(power(A1,3)-power(A2,3))*power(aj,3)+bj*(6*bj*cj*((A2*A2-A1*A1)*alpha*alpha+B1*B1-B2*B2)+bj*bj*(6*A1*A1*B1*alpha*alpha-6*A2*A2*B2*alpha*alpha+power(A1,3)*power(alpha,3)-power(A2,3)*power(alpha,3)-2*power(B1,3)+2*power(B2,3))+6*(B2-B1)*cj*cj)))/(6*ai*bj*bj*(aj+alpha*bj));
            }
        } else {  // bi != 0 and bj != 0
            if (ai/bi == aj/bj) {
                if ((ai+alpha*bi).is_zero())
                    throw boost::enable_current_exception(integral_error()) << err_description(__FILE__) << errno_code(ERR_IS_ZERO) << code_line(__LINE__);
                mpfr_class tmp = abs(ai+alpha*bi);
                pr = logarithm(tmp)*(2*ai*bj*(bj*(power(A1,3)*(3*B1*bi-ci)+power(A2,3)*(ci-3*B2*bi))+2*(power(A2,3)-power(A1,3))*bi*cj)+3*(power(A2,4)-power(A1,4))*ai*ai*bj*bj+bi*(power(A1,2)*(B1*bj-cj)*(bi*(cj-3*B1*bj)+2*bj*ci)+power(A2,2)*(B2*bj-cj)*(bi*(3*B2*bj-cj)-2*bj*ci)))/(2*bi*bi*bj*bj);
                pr += (20*power(ai,3)*bj*(bj*(power(A1,3)*(3*B1*bi-ci)+power(A2,3)*(ci-3*B2*bi))+2*(power(A2,3)-power(A1,3))*bi*cj)+18*ai*ai*bi*(A1*A1*(B1*bj-cj)*(bi*(cj-3*B1*bj)+2*bj*ci)+A2*A2*(B2*bj-cj)*(bi*(3*B2*bj-cj)-2*bj*ci))+12*ai*bi*bi*(A1*(B1*bi-ci)*power(cj-B1*bj,2)-A2*(B2*bi-ci)*power(cj-B2*bj,2))+21*(power(A2,4)-power(A1,4))*power(ai,4)*bj*bj+power(bi,3)*(4*ci*(3*(B1*B1-B2*B2)*bj*cj+(power(B2,3)-power(B1,3))*bj*bj+3*(B2-B1)*cj*cj)+bi*(8*(power(B2,3)-power(B1,3))*bj*cj+3*(power(B1,4)-power(B2,4))*bj*bj+6*(B1*B1-B2*B2)*cj*cj)))/(24*bi*bi*bj*bj*power(ai+alpha*bi,2));
                pr += alpha*(8*ai*ai*bj*(bj*(power(A1,3)*(3*B1*bi-ci)+power(A2,3)*(ci-3*B2*bi))+2*(power(A2,3)-power(A1,3))*bi*cj)+12*ai*bi*(A1*A1*(B1*bj-cj)*(bi*(cj-3*B1*bj)+2*bj*ci)+A2*A2*(B2*bj-cj)*(bi*(3*B2*bj-cj)-2*bj*ci))+3*(power(A2,4)-power(A1,4))*power(ai,3)*bj*bj+12*bi*bi*(A1*(B1*bi-ci)*power(cj-B1*bj,2)-A2*(B2*bi-ci)*power(cj-B2*bj,2)))/(12*bi*bj*bj*power(ai+alpha*bi,2));
                pr += alpha*alpha*(ai*(bj*(33*(power(A1,4)-power(A2,4))*ai+16*(power(A1,3)*(ci-3*B1*bi)+power(A2,3)*(3*B2*bi-ci)))+32*(power(A1,3)-power(A2,3))*bi*cj))/(24*bj*power(ai+alpha*bi,2));
                pr += power(alpha,3)*(bi*(bj*(3*(power(A1,4)-power(A2,4))*ai+2*power(A1,3)*(ci-3*B1*bi)+2*power(A2,3)*(3*B2*bi-ci))+4*(power(A1,3)-power(A2,3))*bi*cj))/(6*bj*power(ai+alpha*bi,2));
                pr += power(alpha,4)*((power(A2,4)-power(A1,4))*bi*bi)/(8*power(ai+alpha*bi,2));
            } else {
//                std::cout << A1*aj-B1*bj+cj << std::endl;
//                std::cout << A2*aj-B2*bj+cj << std::endl;
//                std::cout << A1*ai-B1*bi+ci << std::endl;
//                std::cout << A2*ai-B2*bi+ci << std::endl;

                if ((A1*aj-B1*bj+cj).is_zero() && (A2*aj-B2*bj+cj).is_zero()) {
//                    std::cout << aj+alpha*bj << std::endl;
                    if ((aj+alpha*bj).is_zero()) {
                        mpfr_class tmp = abs((ai*bj-aj*bi)/bj);
                        pr = -((aj*bi-ai*bj)*logarithm(tmp)*(3*(power(A1,4)-power(A2,4))*aj*bi+3*(power(A2,4)-power(A1,4))*ai*bj-4*(power(A1,3)-power(A2,3))*(bj*ci-bi*cj)))/(12*bi*bi*bj*bj);
                        pr += (aj*(9*(power(A1,4)-power(A2,4))*aj*bi+6*(power(A2,4)-power(A1,4))*ai*bj-8*(power(A1,3)-power(A2,3))*(bj*ci-bi*cj)))/(24*bi*bj*bj);
//                        pr = alpha*(3*ai*ai*bj*bj*(power(A1,3)*(cj-B1*bj)+power(A2,3)*(B2*bj-cj))+ai*bj*(A1*A1*(B1*bj-cj)*(bi*(9*B1*bj-5*cj)-4*bj*ci)+A2*A2*(B2*bj-cj)*(bi*(5*cj-9*B2*bj)+4*bj*ci))+2*bi*(A1*power(cj-B1*bj,2)*(bi*(cj-3*B1*bj)+2*bj*ci)+A2*power(cj-B2*bj,2)*(bi*(3*B2*bj-cj)-2*bj*ci)))/(12*bi*bj*(aj*bi-ai*bj)*(aj+alpha*bj));
//                        pr += alpha*alpha*(ai*bj*(bj*(power(A1,3)*(21*B1*bi-8*ci)+power(A2,3)*(8*ci-21*B2*bi))+13*(power(A2,3)-power(A1,3))*bi*cj)+6*(power(A2,4)-power(A1,4))*ai*ai*bj*bj+bi*(A1*A1*(B1*bj-cj)*(bi*(7*cj-15*B1*bj)+8*bj*ci)+A2*A2*(B2*bj-cj)*(bi*(15*B2*bj-7*cj)-8*bj*ci)))/(24*bi*(aj*bi-ai*bj)*(aj+alpha*bj));
//                        pr += power(alpha,3)*(bj*((power(A1,4)-power(A2,4))*ai*bj+bi*(power(A1,3)*(cj-B1*bj)+power(A2,3)*(B2*bj-cj))))/(8*(aj*bi-ai*bj)*(aj+alpha*bj));
                    } else {
                        mpfr_class tmp = abs(ai+alpha*bi);
                        pr = -((aj*bi-ai*bj)*logarithm(tmp)*(3*(power(A1,4)-power(A2,4))*aj*bi+3*(power(A2,4)-power(A1,4))*ai*bj-4*(power(A1,3)-power(A2,3))*(bj*ci-bi*cj)))/(12*bi*bi*bj*bj);
                        pr += (alpha*(12*(power(A2,4)-power(A1,4))*aj*bi+6*(power(A1,4)-power(A2,4))*ai*bj+8*power(A1,3)*bj*ci-8*power(A1,3)*bi*cj-8*power(A2,3)*bj*ci+8*power(A2,3)*bi*cj-3*power(A1,4)*alpha*bi*bj+3*power(A2,4)*alpha*bi*bj))/(24*bi*bj);
                    }
                } else if ((A1*ai-B1*bi+ci).is_zero() && (A2*ai-B2*bi+ci).is_zero()) {
//                    pr = -logarithm(abs(tmp))*(((A1-A2)*aj+(B2-B1)*bj)*(aj*aj*bi*(bj*(27*(-A1*A1*B1-(A1*A1+A2*A1+A2*A2)*B2)*bi+(38*A1*A1+11*A2*A1+11*A2*A2)*ci)+16*(A1*A1+A2*A1+A2*A2)*bi*cj)+aj*(bj*bj*(-(22*A2*B2+49*A1*(B1+B2))*bi*ci+27*(A2*B2*B2+A1*(B1*B1+B2*B1+B2*B2))*bi*bi+(23*A1+A2)*ci*ci)-4*bi*bj*cj*(8*(A2*B2+A1*(B1+B2))*bi-(13*A1+5*A2)*ci)+6*(A1+A2)*bi*bi*cj*cj)+9*(A1+A2)*(A1*A1+A2*A2)*power(aj,3)*bi*bi+bj*(bi*bi*(16*(B1*B1+B2*B1+B2*B2)*bj*cj-9*(B1+B2)*(B1*B1+B2*B2)*bj*bj-6*(B1+B2)*cj*cj)+4*bi*ci*(-9*(B1+B2)*bj*cj+5*(B1*B1+B2*B1+B2*B2)*bj*bj+3*cj*cj)-12*bj*ci*ci*((B1+B2)*bj-2*cj))))/(12*bj*bj*power(aj*bi-ai*bj,2));
                    if ((ai+alpha*bi).is_zero()) {
                        mpfr_class tmp = abs((aj*bi-ai*bj)/bi);
                        pr = -(logarithm(tmp)*(2*aj*bi*(9*(power(A2,4)-power(A1,4))*ai*bj-8*(power(A1,3)-power(A2,3))*(bj*ci-bi*cj))+16*(power(A1,3)-power(A2,3))*ai*bj*(bj*ci-bi*cj)+9*(power(A1,4)-power(A2,4))*aj*aj*bi*bi+9*(power(A1,4)-power(A2,4))*ai*ai*bj*bj+6*(A1*A1-A2*A2)*power(bj*ci-bi*cj,2)))/(12*bi*bi*bj*bj);
                        pr += (-16*(power(A1,3)-power(A2,3))*aj*bi*(bi*cj-bj*ci)+3*(A1*A1-A2*A2)*(3*(A1*A1+A2*A2)*ai*ai*bj*bj-4*power(bj*ci-bi*cj,2))+6*(power(A2,4)-power(A1,4))*aj*aj*bi*bi)/(24*bi*bi*bj*bj);
//                        pr += (4*power(aj,3)*(bj*(power(A1,3)*(ci-3*B1*bi)+power(A2,3)*(3*B2*bi-ci))+2*(power(A1,3)-power(A2,3))*bi*cj)+6*aj*aj*(A1*A1*(B1*bj-cj)*(bi*(3*B1*bj-cj)-2*bj*ci)+A2*A2*(B2*bj-cj)*(bi*(cj-3*B2*bj)+2*bj*ci))+12*aj*bj*(A2*(B2*bi-ci)*power(cj-B2*bj,2)-A1*(B1*bi-ci)*power(cj-B1*bj,2))+3*(power(A1,4)-power(A2,4))*power(aj,4)*bi+bj*bj*(4*ci*(3*(B1*B1-B2*B2)*bj*cj+(power(B2,3)-power(B1,3))*bj*bj+3*(B2-B1)*cj*cj)+bi*(8*(power(B2,3)-power(B1,3))*bj*cj+3*(power(B1,4)-power(B2,4))*bj*bj+6*(B1*B1-B2*B2)*cj*cj)))/(12*bj*bj*(ai*bj-aj*bi)*(aj+alpha*bj));
//                        pr += alpha*(aj*(aj*bi*(power(A1,3)*(bi*(8*cj-15*B1*bj)+7*bj*ci)+power(A2,3)*(bi*(15*B2*bj-8*cj)-7*bj*ci))+ai*bj*(bj*(3*(power(A2,4)-power(A1,4))*ai+4*power(A1,3)*(3*B1*bi-ci)+4*power(A2,3)*(ci-3*B2*bi))+8*(power(A2,3)-power(A1,3))*bi*cj)+6*(power(A1,4)-power(A2,4))*aj*aj*bi*bi))/(12*bi*bj*(aj*bi-ai*bj)*(aj+alpha*bj));
//                        pr += alpha*alpha*(2*bj*(3*(power(A1,4)-power(A2,4))*ai+4*power(A1,3)*(ci-3*B1*bi)+4*power(A2,3)*(3*B2*bi-ci))+9*(power(A1,4)-power(A2,4))*aj*bi+16*(power(A1,3)-power(A2,3))*bi*cj)/(24*bi*(aj+alpha*bj));
//                        pr += power(alpha,3)*((power(A2,4)-power(A1,4))*bj)/(8*(aj+alpha*bj));
                    } else {
                        mpfr_class tmp = abs(aj+alpha*bj);
                        pr = -(logarithm(tmp)*(2*aj*bi*(9*(power(A2,4)-power(A1,4))*ai*bj-8*(power(A1,3)-power(A2,3))*(bj*ci-bi*cj))+16*(power(A1,3)-power(A2,3))*ai*bj*(bj*ci-bi*cj)+9*(power(A1,4)-power(A2,4))*aj*aj*bi*bi+9*(power(A1,4)-power(A2,4))*ai*ai*bj*bj+6*(A1*A1-A2*A2)*power(bj*ci-bi*cj,2)))/(12*bi*bi*bj*bj);
                        pr += (aj*bi*(2*ai*bj*(power(A2,3)*(bj*(9*A2*alpha*bi+16*ci)-16*bi*cj)+16*power(A1,3)*(bi*cj-bj*ci)-9*power(A1,4)*alpha*bi*bj)+18*(power(A2,4)-power(A1,4))*ai*ai*bj*bj+A2*A2*(16*A2*alpha*bi*bj*(bj*ci-bi*cj)-9*A2*A2*alpha*alpha*bi*bi*bj*bj+12*power(bj*ci-bi*cj,2))+16*power(A1,3)*alpha*bi*bj*(bi*cj-bj*ci)-12*A1*A1*power(bj*ci-bi*cj,2)+9*power(A1,4)*alpha*alpha*bi*bi*bj*bj)+bj*(-6*(A1*A1-A2*A2)*ai*(3*A1*A1*alpha*alpha*bi*bi*bj*bj+3*A2*A2*alpha*alpha*bi*bi*bj*bj-2*power(bj*ci-bi*cj,2))+16*(power(A1,3)-power(A2,3))*ai*ai*bj*(bj*ci-bi*cj)+6*(power(A1,4)-power(A2,4))*power(ai,3)*bj*bj+alpha*alpha*bi*bi*bj*(power(A2,3)*(bj*(3*A2*alpha*bi+16*ci)-16*bi*cj)+16*power(A1,3)*(bi*cj-bj*ci)-3*power(A1,4)*alpha*bi*bj))+2*aj*aj*bi*bi*(9*(power(A1,4)-power(A2,4))*ai*bj-2*power(A2,3)*(bj*(3*A2*alpha*bi+4*ci)-4*bi*cj)+8*power(A1,3)*(bj*ci-bi*cj)+6*power(A1,4)*alpha*bi*bj)+6*(power(A2,4)-power(A1,4))*power(aj,3)*power(bi,3))/(24*power(bi,3)*bj*bj*(aj+alpha*bj));
//                        pr += (-16*(power(A1,3)-power(A2,3))*aj*bi*(A2*ai*bj+bi*(cj-B2*bj))+3*(A1*A1-A2*A2)*(8*A2*ai*bi*bj*(B2*bj-cj)+(3*A1*A1-A2*A2)*ai*ai*bj*bj-4*bi*bi*power(cj-B2*bj,2))+6*(power(A2,4)-power(A1,4))*aj*aj*bi*bi)/(24*bi*bi*bj*bj);
                    }
                } else {
                    if ((aj+alpha*bj).is_zero()) {
//                        std::cout << A1*aj-B1*bj+cj << std::endl;
//                        std::cout << A2*aj-B2*bj+cj << std::endl;
//                        std::cout << A1*ai-B1*bi+ci << std::endl;
//                        std::cout << A2*ai-B2*bi+ci << std::endl;
//                        std::cout << "A1=" << A1 << std::endl;
//                        std::cout << "A2=" << A2 << std::endl;
//                        std::cout << "B1=" << B1 << std::endl;
//                        std::cout << "B2=" << B2 << std::endl;
//                        std::cout << "alpha" << alpha << std::endl;
//                        std::cout << "ai=" << ai << std::endl;
//                        std::cout << "bi=" << bi << std::endl;
//                        std::cout << "ci=" << ci << std::endl;
//                        std::cout << "aj=" << aj << std::endl;
//                        std::cout << "bj=" << bj << std::endl;
//                        std::cout << "cj=" << cj << std::endl;
                        throw boost::enable_current_exception(integral_error()) << err_description(__FILE__) << errno_code(ERR_IS_ZERO) << code_line(__LINE__);
                    }

                    if ((ai+alpha*bi).is_zero()) {
                        throw boost::enable_current_exception(integral_error()) << err_description(__FILE__) << errno_code(ERR_IS_ZERO) << code_line(__LINE__);
                    }

                    mpfr_class tmp1 = abs(ai+alpha*bi);
                    pr = logarithm(tmp1)*(4*power(ai,3)*bj*(bj*(power(A1,3)*(3*B1*bi-ci)-power(A2,3)*(3*B2*bi-ci))-2*(power(A1,3)-power(A2,3))*bi*cj)+6*ai*ai*bi*(A1*A1*(B1*bj-cj)*(2*bj*ci-bi*(3*B1*bj-cj))-A2*A2*(B2*bj-cj)*(2*bj*ci-bi*(3*B2*bj-cj)))+12*ai*bi*bi*(A1*(B1*bi-ci)*power(B1*bj-cj,2)-A2*(B2*bi-ci)*power(B2*bj-cj,2))-3*(power(A1,4)-power(A2,4))*power(ai,4)*bj*bj+power(bi,3)*(4*ci*(-3*(B1*B1-B2*B2)*bj*cj+(power(B1,3)-power(B2,3))*bj*bj+3*(B1-B2)*cj*cj)-bi*(-8*(power(B1,3)-power(B2,3))*bj*cj+3*(power(B1,4)-power(B2,4))*bj*bj+6*(B1*B1-B2*B2)*cj*cj)))/(12*bi*bi*power(aj*bi-ai*bj,2));
                    mpfr_class tmp2 = abs(aj+alpha*bj);
                    pr += logarithm(tmp2)*(4*power(aj,3)*bi*(bj*(3*(power(A1,4)-power(A2,4))*ai+2*(power(A1,3)*(3*B1*bi-ci)-power(A2,3)*(3*B2*bi-ci)))-4*(power(A1,3)-power(A2,3))*bi*cj)+6*aj*aj*(2*ai*bj*(2*(power(A1,3)-power(A2,3))*bi*cj-bj*(power(A1,3)*(3*B1*bi-ci)-power(A2,3)*(3*B2*bi-ci)))+bi*(A1*A1*(B1*bj-cj)*(2*bj*ci-bi*(3*B1*bj-cj))-A2*A2*(B2*bj-cj)*(2*bj*ci-bi*(3*B2*bj-cj))))+12*ai*aj*bj*(A2*A2*(B2*bj-cj)*(2*bj*ci-bi*(3*B2*bj-cj))-A1*A1*(B1*bj-cj)*(2*bj*ci-bi*(3*B1*bj-cj)))+bj*bj*(bi*(bi*(-8*(power(B1,3)-power(B2,3))*bj*cj+3*(power(B1,4)-power(B2,4))*bj*bj+6*(B1*B1-B2*B2)*cj*cj)-4*ci*(-3*(B1*B1-B2*B2)*bj*cj+(power(B1,3)-power(B2,3))*bj*bj+3*(B1-B2)*cj*cj))-12*ai*(A1*(B1*bi-ci)*power(B1*bj-cj,2)-A2*(B2*bi-ci)*power(B2*bj-cj,2)))-9*(power(A1,4)-power(A2,4))*power(aj,4)*bi*bi)/(12*bj*bj*power(aj*bi-ai*bj,2));
                    pr += (4*power(aj,3)*(bj*(power(A1,3)*(ci-3*B1*bi)+power(A2,3)*(3*B2*bi-ci))+2*(power(A1,3)-power(A2,3))*bi*cj)+6*aj*aj*(A1*A1*(B1*bj-cj)*(bi*(3*B1*bj-cj)-2*bj*ci)+A2*A2*(B2*bj-cj)*(bi*(cj-3*B2*bj)+2*bj*ci))+12*aj*bj*(A2*(B2*bi-ci)*power(cj-B2*bj,2)-A1*(B1*bi-ci)*power(cj-B1*bj,2))+3*(power(A1,4)-power(A2,4))*power(aj,4)*bi+bj*bj*(4*ci*(3*(B1*B1-B2*B2)*bj*cj+(power(B2,3)-power(B1,3))*bj*bj+3*(B2-B1)*cj*cj)+bi*(8*(power(B2,3)-power(B1,3))*bj*cj+3*(power(B1,4)-power(B2,4))*bj*bj+6*(B1*B1-B2*B2)*cj*cj)))/(12*bj*bj*(ai*bj-aj*bi)*(aj+alpha*bj));
                    pr += alpha*(aj*(bj*(3*(power(A1,4)-power(A2,4))*ai+4*power(A1,3)*(ci-3*B1*bi)+4*power(A2,3)*(3*B2*bi-ci))+6*(power(A1,4)-power(A2,4))*aj*bi+8*(power(A1,3)-power(A2,3))*bi*cj))/(12*bi*bj*(aj+alpha*bj));
                    pr += alpha*alpha*(2*bj*(3*(power(A1,4)-power(A2,4))*ai+4*power(A1,3)*(ci-3*B1*bi)+4*power(A2,3)*(3*B2*bi-ci))+9*(power(A1,4)-power(A2,4))*aj*bi+16*(power(A1,3)-power(A2,3))*bi*cj)/(24*bi*(aj+alpha*bj));
                    pr += power(alpha,3)*((power(A2,4)-power(A1,4))*bj)/(8*(aj+alpha*bj));
                }
            }
        }
        
        return pr;
    }
    
    inline bool lines_are_equal(const mpfr_class& a1, const mpfr_class& b1, const mpfr_class& c1,
                                const mpfr_class& a2, const mpfr_class& b2, const mpfr_class& c2) {
        if ((a1 == 0 && a2 != 0) || (a1 != 0 && a2 == 0))
            return false;
        
        if ((b1 == 0 && b2 != 0) || (b1 != 0 && b2 == 0))
            return false;
        
        if ((c1 == 0 && c2 != 0) || (c1 != 0 && c2 == 0))
            return false;
        
        return a2*b1 == a1*b2 && b2*c1 == b1*c2 && a2*c1 == a1*c2;
    }
    
    void volume4d(const mpfr_class (&alpha)[2],
                  const mpfr_class (&A)[2], const mpfr_class (&B)[2],
                  const mpfr_class (&a)[4], const mpfr_class (&b)[4], const mpfr_class (&c)[4],
                  mpfr_class* result) {
        mpfr_class X1X3X3(0);
        mpfr_class X3X1X1(0);
        if (!lines_are_equal(a[0], b[0], c[0], a[2], b[2], c[2])) {
//            std::cout << std::setprecision(20)
//                      << a[0] << " " << b[0] << " " << c[0] << std::endl
//                      << a[2] << " " << b[2] << " " << c[2] << std::endl;
            mpfr_class tmp1 = XiXjXj(alpha[1], A, B, a[0], b[0], c[0], a[2], b[2], c[2]);
            mpfr_class tmp2 = XiXjXj(alpha[0], A, B, a[0], b[0], c[0], a[2], b[2], c[2]);
//            std::cout << "X1X3X3=" << tmp1 << " - " << tmp2 << " = " << tmp1-tmp2 << std::endl;
            X1X3X3 = tmp1 - tmp2;
            mpfr_class tmp3 = XiXjXj(alpha[1], A, B, a[2], b[2], c[2], a[0], b[0], c[0]);
            mpfr_class tmp4 = XiXjXj(alpha[0], A, B, a[2], b[2], c[2], a[0], b[0], c[0]);
//            std::cout << "X3X1X1=" << tmp3 << " - " << tmp4 << " = " << tmp3-tmp4 << std::endl;
            X3X1X1 = tmp3 - tmp4;
        }
        mpfr_class X2X3X3(0);
        mpfr_class X3X2X2(0);
        if (!lines_are_equal(a[1], b[1], c[1], a[2], b[2], c[2])) {
//            std::cout << a[1] << " " << b[1] << " " << c[1] << std::endl
//                      << a[2] << " " << b[2] << " " << c[2] << std::endl;
            mpfr_class tmp1 = XiXjXj(alpha[1], A, B, a[1], b[1], c[1], a[2], b[2], c[2]);
            mpfr_class tmp2 = XiXjXj(alpha[0], A, B, a[1], b[1], c[1], a[2], b[2], c[2]);
//            std::cout << "X2X3X3=" << tmp1 << " - " << tmp2 << " = " << tmp1-tmp2 << std::endl;
            X2X3X3 = tmp1 - tmp2;
            mpfr_class tmp3 = XiXjXj(alpha[1], A, B, a[2], b[2], c[2], a[1], b[1], c[1]);
            mpfr_class tmp4 = XiXjXj(alpha[0], A, B, a[2], b[2], c[2], a[1], b[1], c[1]);
//            std::cout << "X3X2X2=" << tmp3 << " - " << tmp4 << " = " << tmp3-tmp4 << std::endl;
            X3X2X2 = tmp3 - tmp4;
        }
        mpfr_class X1X4X4(0);
        mpfr_class X4X1X1(0);
        if (!lines_are_equal(a[0], b[0], c[0], a[3], b[3], c[3])) {
//            std::cout << a[0] << " " << b[0] << " " << c[0] << std::endl
//                      << a[3] << " " << b[3] << " " << c[3] << std::endl;
            mpfr_class tmp1 = XiXjXj(alpha[1], A, B, a[0], b[0], c[0], a[3], b[3], c[3]);
            mpfr_class tmp2 = XiXjXj(alpha[0], A, B, a[0], b[0], c[0], a[3], b[3], c[3]);
//            std::cout << "X1X4X4=" << tmp1 << " - " << tmp2 << " = " << tmp1-tmp2 << std::endl;
            X1X4X4 = tmp1 - tmp2;
            mpfr_class tmp3 = XiXjXj(alpha[1], A, B, a[3], b[3], c[3], a[0], b[0], c[0]);
            mpfr_class tmp4 = XiXjXj(alpha[0], A, B, a[3], b[3], c[3], a[0], b[0], c[0]);
//            std::cout << "X4X1X1=" << tmp3 << " - " << tmp4 << " = " << tmp3-tmp4 << std::endl;
            X4X1X1 = tmp3 - tmp4;
        }
        mpfr_class X2X4X4(0);
        mpfr_class X4X2X2(0);
        if (!lines_are_equal(a[1], b[1], c[1], a[3], b[3], c[3])) {
//            std::cout << a[1] << " " << b[1] << " " << c[1] << std::endl
//                      << a[3] << " " << b[3] << " " << c[3] << std::endl;
            mpfr_class tmp1 = XiXjXj(alpha[1], A, B, a[1], b[1], c[1], a[3], b[3], c[3]);
            mpfr_class tmp2 = XiXjXj(alpha[0], A, B, a[1], b[1], c[1], a[3], b[3], c[3]);
//            std::cout << "X2X4X4=" << tmp1 << " - " << tmp2 << " = " << tmp1-tmp2 << std::endl;
            X2X4X4 = tmp1 - tmp2;
            mpfr_class tmp3 = XiXjXj(alpha[1], A, B, a[3], b[3], c[3], a[1], b[1], c[1]);
            mpfr_class tmp4 = XiXjXj(alpha[0], A, B, a[3], b[3], c[3], a[1], b[1], c[1]);
//            std::cout << "X4X2X2=" << tmp3 << " - " << tmp4 << " = " << tmp3-tmp4 << std::endl;
            X4X2X2 = tmp3 - tmp4;
        }

        mpfr_class tmp_alpha;
        if (alpha[0] < alpha[1])
            tmp_alpha = (alpha[0]+alpha[1])/2;
        else
            tmp_alpha = alpha[0] + (alpha[0]-alpha[1])/2;
        
        mpfr_class beta = (A[0]*tmp_alpha+B[0]+A[1]*tmp_alpha+B[1])/2;
        mpfr_class x1 = (b[0]*beta-c[0])/(a[0]+tmp_alpha*b[0]);
        mpfr_class x4 = (b[3]*beta-c[3])/(a[3]+tmp_alpha*b[3]);
        
//        std::cout << std::endl << result << std::endl << std::endl << std::endl;

        if (x1 < x4)
            *result += (X1X3X3-X3X1X1+X3X2X2-X2X3X3+X4X1X1-X1X4X4+X2X4X4-X4X2X2)/2;
        else
            *result -= (X1X3X3-X3X1X1+X3X2X2-X2X3X3+X4X1X1-X1X4X4+X2X4X4-X4X2X2)/2;
    }

    void volume4d(const Cell& cell,
                  const mpfr_class (&a)[4],
                  const mpfr_class (&b)[4],
                  const mpfr_class (&c)[4],
                  mpfr_class* result) {
//          std::cout << cell << std::endl;

        list<IntegralTrapezoid> trapezoids;
        cell.decompose(&trapezoids);

        IntegralTrapezoid* t_left = NULL;
        IntegralTrapezoid* t_right = NULL;
        for (list<IntegralTrapezoid>::iterator it = trapezoids.begin(); it != trapezoids.end(); ++it) {
            if (it->minus_infinity()) {
                t_left = &(*it);
                continue;
            }

            if (it->plus_infinity()) {
                t_right = &(*it);
                continue;
            }

//            std::cout << "left_alpha " << it->alpha_min().exact() << " right_alpha " << it->alpha_max().exact()
//                      << " bottom y=" << it->bottom().slope().exact() << " x + " << it->bottom().y_intersept().exact()
//                      << " top y=" << it->top().slope().exact() << " x + " << it->top().y_intersept().exact() << std::endl;

            mpfr_class alpha[2] = { mpfr_class(it->alpha_min()), mpfr_class(it->alpha_max()) };
            mpfr_class A[2] = { mpfr_class(it->bottom().slope()), mpfr_class(it->top().slope()) };
            mpfr_class B[2] = { mpfr_class(it->bottom().y_intersept()), mpfr_class(it->top().y_intersept()) };

//            std::cout << "left_alpha " << alpha[0] << " right_alpha " << alpha[1]
//                      << " bottom y=" << A[0] << " x + " << B[0]
//                      << " top y=" << A[1] << " x + " << B[1] << std::endl;

            mpfr_class pr(0);
            volume4d(alpha, A, B, a, b, c, &pr);
            if (pr < -PRECISION)
                throw boost::enable_current_exception(integral_error()) << err_description("prob =" + boost::to_string(pr.get_d()) + "\n in" + __FILE__) << code_line(__LINE__) << errno_code(ERR_LESS_THAN_ZERO);
            *result += pr;

//                std::cout << "trapezoid value " << v << std::endl;
        }
        
        if (t_left != NULL && t_right != NULL) {
            
//                std::cout << "left_alpha " << t_right->alpha_min() << " right_alpha " << t_left->alpha_max() <<
//                    " bottom y=" << t_left->bottom().slope() << " x + " << t_left->bottom().y_intersept() <<
//                    " top y=" << t_left->top().slope() << " x + " << t_left->top().y_intersept() << std::endl;
            
            mpfr_class alpha[2] = { mpfr_class(t_right->alpha_min()), mpfr_class(t_left->alpha_max()) };
            mpfr_class A[2] = { mpfr_class(t_right->bottom().slope()), mpfr_class(t_right->top().slope()) };
            mpfr_class B[2] = { mpfr_class(t_right->bottom().y_intersept()), mpfr_class(t_right->top().y_intersept()) };

            mpfr_class pr(0);
            volume4d(alpha, A, B, a, b, c, &pr);
            if (pr < -PRECISION)
                throw boost::enable_current_exception(integral_error()) << err_description(__FILE__) << code_line(__LINE__) << errno_code(ERR_LESS_THAN_ZERO);
            *result += pr;
            
//                std::cout << "trapezoid value " << v << std::endl;
        } else if (t_left != NULL || t_right != NULL) {
            throw boost::enable_current_exception(integral_error()) << err_description(__FILE__) << code_line(__LINE__);
        }
    }
    
    void volume4d(const Segment (&segs)[4], const Obstacles& obstacles, mpfr_class* result) {
        set<Line, bool (*)(const Line&, const Line&)> lines(operator<);
        set<Point> points;
        InexactSegment segs_i[4];
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 2; ++j)
                points.insert(segs[i].vertex(j));

            segs_i[i] = InexactSegment(InexactPoint(CGAL::to_double(segs[i].vertex(0).x()), CGAL::to_double(segs[i].vertex(0).y())),
                                       InexactPoint(CGAL::to_double(segs[i].vertex(1).x()), CGAL::to_double(segs[i].vertex(1).y())));
        }

        for (set<Point>::iterator it = points.begin(); it != points.end(); ++it) {
            for (set<Point>::iterator jt = it; jt != points.end(); ++jt) {
                if (it == jt)
                    continue;

                InexactPoint p1_i(CGAL::to_double(it->x()), CGAL::to_double(it->y()));
                InexactPoint p2_i(CGAL::to_double(jt->x()), CGAL::to_double(jt->y()));
                InexactLine l_i(p1_i, p2_i);
//                std::cout << p1 << " " << p2 << std::endl;
                if ( intersects_in_order(l_i, segs_i) ) {
//                    std::cout << p1_i << " " << p2_i << std::endl;
                    lines.insert(Line(*it, *jt));
                }
            }
        }

        if (lines.size() < 3)
            return;

//        CGAL::IO::Mode old_mode = CGAL::get_mode(std::cout);
//        CGAL::set_pretty_mode(std::cout);
//        std::cout << std::setprecision(20) << segs[0] << "\t" << segs[1] << std::endl;
//        std::cout << std::setprecision(20) << segs[2] << "\t" << segs[3] << std::endl << std::endl;
//        CGAL::set_mode(std::cout, old_mode);
//
//            for (Line l : lines) {
//                if (l.is_vertical())
//                    std::cout << l << std::endl << std::endl;
//                else
//                    std::cout << "y=" << -l.a()/l.b() << " x + " << -l.c()/l.b() << std::endl << std::endl;
//            }
//            std::cout << std::endl << std::endl;
        
        Point separator;
        {
            bool flag_sep = false;
            for (int i = 0; i < 4; ++i) {
                Line l(segs[i].source(), segs[i].target());
                if (lines.find(l) == lines.end()) {
                    separator = segs[i].source()+2*(segs[i].target()-segs[i].source());
                    flag_sep = true;
                    break;
                }
            }

            if (!flag_sep)
                throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);

//            Point A, B, C;
//            if (segs[0].vertex(0) == segs[1].vertex(0)) {
//                A = segs[0].vertex(0); B = segs[0].vertex(1); C = segs[1].vertex(1);
//            } else if (segs[0].vertex(0) == segs[1].vertex(1)) {
//                A = segs[0].vertex(0); B = segs[0].vertex(1); C = segs[1].vertex(0);
//            } else if (segs[0].vertex(1) == segs[1].vertex(0)) {
//                A = segs[0].vertex(1); B = segs[0].vertex(0); C = segs[1].vertex(1);
//            } else {
//                A = segs[0].vertex(1); B = segs[0].vertex(0); C = segs[1].vertex(0);
//            }
//
//            separator = Point(3*A.x()-B.x()-C.x(), 3*A.y()-B.y()-C.y());
        }
        
//        std::cout << "separator " << separator << std::endl;
        list<Cell> cells;

        Cell original_cell;
        original_cell.Initialize(lines.begin(), lines.end(), separator);

        if (original_cell.type() == Cell::EMPTY)
            return;

        cells.push_back(original_cell);
        SimplePolygon ch;
        CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(ch));
        assert(ch.is_convex());
        for (list<Obstacle>::const_iterator ob_it = obstacles.obstacles().begin();
             ob_it != obstacles.obstacles().end(); ++ob_it) {
            if (cells.empty())
                return;

            if (!do_intersect_interiors(ob_it->obstacle, ch)) {
//                std::cout << obstacles.obstacles()[i] << std::endl << ch << std::endl;
                continue;
            }

//            std::cout << ob_it->obstacle << std::endl;

            list<Cell> tmp_cells;
            for (list<Cell>::iterator cell_it = cells.begin(); cell_it != cells.end(); ++cell_it) {
                list<Cell> intersections;
//                std::cout << "Intersect " << *cell_it << std::endl << "with " << ob_it->dual_cell << std::endl << std::endl;
                cell_it->intersect(ob_it->dual_cell, &intersections);
//                for (Cell c : intersections)
//                    std::cout << "Result " << c << std::endl;

                tmp_cells.insert(tmp_cells.end(), intersections.begin(), intersections.end());
            }

            swap(cells, tmp_cells);
        }

        mpfr_class a[4];
        mpfr_class b[4];
        mpfr_class c[4];
        for (int i = 0; i < 4; ++i) {
            a[i] = mpfr_class(segs[i].supporting_line().a());
            b[i] = mpfr_class(segs[i].supporting_line().b());
            c[i] = mpfr_class(segs[i].supporting_line().c());
        }

        for (list<Cell>::iterator cell_it = cells.begin(); cell_it != cells.end(); ++cell_it) {
//            std::cout << *cell_it << std::endl << std::endl;
            mpfr_class pr(0);
            volume4d(*cell_it, a, b, c, &pr);
            if (pr < -PRECISION)
                throw boost::enable_current_exception(integral_error()) << err_description(__FILE__) << code_line(__LINE__) << errno_code(ERR_LESS_THAN_ZERO);
            *result += pr;
//            std::cout << "=" << pr << std::endl << std::endl << std::endl;
        }
    }

    void volume4d(const SimplePolygon& p1, const SimplePolygon& p2,
                  const Obstacles& obstacles, mpfr_class* result) {
        if (p1 == p2) {
            *result += mpfr_class(CGAL::abs(p1.area()*p2.area()));
            return;
        }

        for (int i1 = 0; i1 < p1.size(); ++i1)
            for (int i2 = 0; i2 < p1.size(); ++i2) {
                if (i1 == i2)
                    continue;

                const Segment& s1 = p1.edge(i1);
                const Segment& s2 = p1.edge(i2);

                for (int i3 = 0; i3 < p2.size(); ++i3)
                    for (int i4 = 0; i4 < p2.size(); ++i4) {
                        if (i3 == i4)
                            continue;

                        const Segment& s3 = p2.edge(i3);
                        const Segment& s4 = p2.edge(i4);
                        if (s1 == s3 || s1 == s4 || s2 == s4)
                            continue;

                        const Segment segs[4] = {s1, s2, s3, s4};
                        mpfr_class pr(0);
                        volume4d(segs, obstacles, &pr);
                        if (pr < -PRECISION)
                            throw boost::enable_current_exception(integral_error()) << err_description(__FILE__) << code_line(__LINE__) << errno_code(ERR_LESS_THAN_ZERO);
                        *result += pr;
                    }
            }

//          std::cout << "polygon 1 " << p1 << std::endl << "polygon 2 " << p2 << std::endl;
//          std::cout << pr << " / " << p1.area()*p2.area() << std::endl;
    }

    void volume4d(const Triangle& t1, const Triangle& t2,
                  const Obstacles& obstacles, mpfr_class* result) {
        if (t1 == t2) {
            *result += mpfr_class(t1.area()*t2.area());
            return;
        }
        
        Segment segments1[3] = { Segment(t1.vertex(0), t1.vertex(1)),
                                 Segment(t1.vertex(1), t1.vertex(2)),
                                 Segment(t1.vertex(2), t1.vertex(0)) };
        Segment segments2[3] = { Segment(t2.vertex(0), t2.vertex(1)),
                                 Segment(t2.vertex(1), t2.vertex(2)),
                                 Segment(t2.vertex(2), t2.vertex(0)) };
        
        for (int i1 = 0; i1 < 3; ++i1)
            for (int i2 = 0; i2 < 3; ++i2) {
                Segment& s1 = segments1[i1];
                Segment& s2 = segments1[i2];
                if (s1 == s2)
                    continue;
                
                for (int i3 = 0; i3 < 3; ++i3)
                    for (int i4 = 0; i4 < 3; ++i4) {
                        Segment& s3 = segments2[i3];
                        Segment& s4 = segments2[i4];
                        if (s3 == s4 || s1 == s3 || s1 == s4 || s2 == s4)
                            continue;
                        
//                        std::cout << std::setprecision(10) << s1 << "\t" << s2 << "\t" << s3 << "\t" << s4 << std::endl;
                        const Segment segs[4] = {s1, s2, s3, s4};
                        mpfr_class pr(0);
                        volume4d(segs, obstacles, &pr);
                        if (pr < -PRECISION)
                            throw boost::enable_current_exception(integral_error()) << err_description(__FILE__) << code_line(__LINE__) << errno_code(ERR_LESS_THAN_ZERO);
                        *result += pr;
                    }
            }
        
//          std::cout << "triangle 1 " << t1 << std::endl << "triangle 2 " << t2 << std::endl;
//          std::cout << pr << std::endl << t1.area()*t2.area() << std::endl;
    }
    
}

using namespace visibility_intl;

mpfr_class volume(const SimplePolygon& p1, const SimplePolygon& p2, const Obstacles& obstacles) {
    list<Polygon> dif1, dif2, intersec;
    CGAL::difference(p1, p2, std::back_inserter(dif1));
    CGAL::difference(p2, p1, std::back_inserter(dif2));
    CGAL::intersection(p1, p2, std::back_inserter(intersec));

    list<SimplePolygon> convex_polys1;
    list<SimplePolygon> convex_polys2;

    for (list<Polygon>::iterator it = dif1.begin(); it != dif1.end(); ++it)
        partition(*it, &convex_polys1);
    for (list<Polygon>::iterator it = dif2.begin(); it != dif2.end(); ++it)
        partition(*it, &convex_polys2);
    for (list<Polygon>::iterator it = intersec.begin(); it != intersec.end(); ++it) {
        list<SimplePolygon> convex_polys;
        partition(*it, &convex_polys);
        convex_polys1.insert(convex_polys1.end(), convex_polys.begin(), convex_polys.end());
        convex_polys2.insert(convex_polys2.end(), convex_polys.begin(), convex_polys.end());
    }

    mpfr_class result(0);
//      int i = 0;
//      int total = triangles1.size()*triangles2.size();
    for (list<SimplePolygon>::iterator it1 = convex_polys1.begin(); it1 != convex_polys1.end(); ++it1)
        for (list<SimplePolygon>::iterator it2 = convex_polys2.begin(); it2 != convex_polys2.end(); ++it2) {
            mpfr_class pr(0);
//            std::cout << *it1 << std::endl << *it2 << std::endl;
            volume4d(*it1, *it2, obstacles, &pr);

            if (pr < -PRECISION)
                throw boost::enable_current_exception(integral_error()) << err_description(__FILE__) << code_line(__LINE__) << errno_code(ERR_LESS_THAN_ZERO);
//            std::cout << tr << "/" << it1->area()*it2->area() << std::endl;
            result += pr;
        }
//      std::cout << std::endl;
//      std::cout << result << std::endl;
    
    return result;
}

mpfr_class volume(const Segment (&segs)[4], const Obstacles& obstacles) {
    mpfr_class result(0);
    volume4d(segs, obstacles, &result);
    return result;
}
