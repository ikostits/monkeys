/*
 * File:   dual_representation.cc
 * Author: irina
 *
 * Created on May 7, 2014, 3:44 PM
 */

#include "dual_representation.h"

#include <list>
#include <set>
#include <vector>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/convex_hull_2.h>

#include "errors.h"

using std::list;
using std::set;
using std::swap;
using std::vector;

namespace {

    bool compare_left_bottom(const Point& p1, const Point& p2) {
        return p1.x() < p2.x() || (p1.x() == p2.x() && p1.y() < p2.y());
    }

    bool compare_left_top(const Point& p1, const Point& p2) {
        return p1.x() < p2.x() || (p1.x() == p2.x() && p1.y() > p2.y());
    }

    void build_hg_part(const vector<Point>& dual_pts,
                       const NonVerticalLine<Kernel>& begin_line,
                       const NonVerticalLine<Kernel>& end_line,
                       bool positive_direction,
                       HourGlassPart* out) {
        SimplePolygon ch;
        CGAL::convex_hull_2(dual_pts.begin(), dual_pts.end(), std::back_inserter(ch));
        assert(ch.is_convex());

        //std::cout << ch << std::endl;
        //std::cout << "begin line " << begin_line << std::endl;
        //std::cout << "end line " << end_line << std::endl;

        double x = positive_direction ? 1 : -1;
        out->begin_dir() = Direction(x, x*begin_line.slope());
        out->end_dir() = Direction(x, x*end_line.slope());

        SimplePolygon::Vertex_const_circulator circ = ch.vertices_circulator();
        while (!begin_line.has_on(*circ)) {
            ++circ;

            if (circ == ch.vertices_circulator())
                throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);
        }

        SimplePolygon::Vertex_const_circulator part_begin_circ = circ;
        out->part().push_back(*part_begin_circ);

        while (!end_line.has_on(*circ)) {
            ++circ;

            if (circ == part_begin_circ)
                throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);

            out->part().push_back(*circ);
        }

        out->correct();
    }

    void decompose_simple_polygon(const SimplePolygon& simple_cell, list<IntegralTrapezoid>* trapezoids) {
        assert(simple_cell.is_convex());

        vector<Point> lower_envelope, upper_envelope;
        get_lower_envelope(simple_cell, &lower_envelope);
        get_upper_envelope(simple_cell, &upper_envelope);

        int low = 1;
        int upp = 1;
        NT left_alpha = lower_envelope[0].x();
        NonVerticalLine<Kernel> bot(lower_envelope[0], lower_envelope[1]);
        NonVerticalLine<Kernel> top(upper_envelope[0], upper_envelope[1]);
        while (low < lower_envelope.size() || upp < upper_envelope.size()) {
            if (lower_envelope[low].x() < upper_envelope[upp].x()) {
                IntegralTrapezoid t(left_alpha, bot, top, lower_envelope[low].x());
                trapezoids->push_back(t);
                bot = NonVerticalLine<Kernel>(lower_envelope[low], lower_envelope[low+1]);
                left_alpha = lower_envelope[low].x();
                low++;
            } else if (lower_envelope[low].x() > upper_envelope[upp].x()) {
                IntegralTrapezoid t(left_alpha, bot, top, upper_envelope[upp].x());
                trapezoids->push_back(t);
                top = NonVerticalLine<Kernel>(upper_envelope[upp], upper_envelope[upp+1]);
                left_alpha = upper_envelope[upp].x();
                upp++;
            } else {
                IntegralTrapezoid t(left_alpha, bot, top, lower_envelope[low].x());
                trapezoids->push_back(t);

                if (low < lower_envelope.size()-1) {
                    bot = NonVerticalLine<Kernel>(lower_envelope[low], lower_envelope[low+1]);
                    top = NonVerticalLine<Kernel>(upper_envelope[upp], upper_envelope[upp+1]);
                }

                left_alpha = lower_envelope[low].x();
                low++;
                upp++;
            }
        }
    }

    enum PartType { WEST , EAST };

    void decompose_hour_glass_part(
                                   const PartType type, const HourGlassPart& hg,
                                   list<IntegralTrapezoid>* trapezoids) {
        //std::cout << hg <<std::endl;

        NonVerticalLine<Kernel> begin_line(hg.part().front(), hg.begin_dir());
        NonVerticalLine<Kernel> end_line(hg.part().back(), hg.end_dir());
        vector<Point>::const_iterator upper_envelope_begin;
        vector<Point>::const_iterator upper_envelope_end;
        vector<Point>::const_iterator lower_envelope_begin;
        vector<Point>::const_iterator lower_envelope_end;
        if (type == EAST) {
            upper_envelope_begin = hg.part().begin();
            upper_envelope_end = std::min_element(hg.part().begin(), hg.part().end(), compare_left_top);
            upper_envelope_end++;
            lower_envelope_begin = std::min_element(hg.part().begin(), hg.part().end(), compare_left_bottom);
            lower_envelope_end = hg.part().end();
        } else {
            upper_envelope_begin = std::max_element(hg.part().begin(), hg.part().end(), compare_left_bottom);
            upper_envelope_end = hg.part().end();
            lower_envelope_begin = hg.part().begin();
            lower_envelope_end = std::max_element(hg.part().begin(), hg.part().end(), compare_left_top);
            lower_envelope_end++;
        }

        vector<Point> lower_envelope, upper_envelope;
        lower_envelope.assign(lower_envelope_begin, lower_envelope_end);
        upper_envelope.assign(upper_envelope_begin, upper_envelope_end);
        std::reverse(upper_envelope.begin(), upper_envelope.end());

        NonVerticalLine<Kernel> bot;
        NonVerticalLine<Kernel> top;
        if (type == WEST) {
            bot = begin_line;
            top = end_line;
        } else {
            if (lower_envelope.size() > 1)
                bot = NonVerticalLine<Kernel>(lower_envelope[0], lower_envelope[1]);
            else
                bot = end_line;

            if (upper_envelope.size() > 1)
                top = NonVerticalLine<Kernel>(upper_envelope[0], upper_envelope[1]);
            else
                top = begin_line;
        }
        int low = 0;
        int upp = 0;
        NT left_alpha;
        if (lower_envelope[0].x() < upper_envelope[0].x()) {
            left_alpha = lower_envelope[0].x();
            low++;
        } else if (lower_envelope[0].x() > upper_envelope[0].x()) {
            left_alpha = upper_envelope[0].x();
            upp++;
        } else {
            left_alpha = lower_envelope[0].x();
            low++;
            upp++;
        }
        if (type == WEST) {
            IntegralTrapezoid t(bot, top, left_alpha);
            trapezoids->push_back(t);
            if (low > 0 && low < lower_envelope.size())
                bot = NonVerticalLine<Kernel>(lower_envelope[low-1], lower_envelope[low]);
            if (upp > 0 && upp < upper_envelope.size())
                top = NonVerticalLine<Kernel>(upper_envelope[upp-1], upper_envelope[upp]);
        }
        while (low < lower_envelope.size() || upp < upper_envelope.size()) {
            if (low < lower_envelope.size() &&
                (upp >= upper_envelope.size() ||
                 lower_envelope[low].x() < upper_envelope[upp].x())) {
                    IntegralTrapezoid t(left_alpha, bot, top, lower_envelope[low].x());
                    trapezoids->push_back(t);
                    if (low < lower_envelope.size() - 1) {
                        bot = NonVerticalLine<Kernel>(lower_envelope[low], lower_envelope[low+1]);
                    } else if (type == EAST) {
                        bot = end_line;
                    }
                    left_alpha = lower_envelope[low].x();
                    low++;
                } else if (upp < upper_envelope.size()&&
                           (low >= lower_envelope.size() ||
                            lower_envelope[low].x() > upper_envelope[upp].x())) {
                               IntegralTrapezoid t(left_alpha, bot, top, upper_envelope[upp].x());
                               trapezoids->push_back(t);
                               if (upp < upper_envelope.size() - 1) {
                                   top = NonVerticalLine<Kernel>(upper_envelope[upp], upper_envelope[upp+1]);
                               } else if (type == EAST) {
                                   top = begin_line;
                               }
                               left_alpha = upper_envelope[upp].x();
                               upp++;
                           } else {
                               IntegralTrapezoid t(left_alpha, bot, top, lower_envelope[low].x());
                               trapezoids->push_back(t);
                               if (low < lower_envelope.size()-1) {
                                   bot = NonVerticalLine<Kernel>(lower_envelope[low], lower_envelope[low+1]);
                               } else if (type == EAST) {
                                   bot = end_line;
                               }

                               if (upp < upper_envelope.size()-1) {
                                   top = NonVerticalLine<Kernel>(upper_envelope[upp], upper_envelope[upp+1]);
                               } else if (type == EAST) {
                                   top = begin_line;
                               }

                               left_alpha = lower_envelope[low].x();
                               low++;
                               upp++;
                           }
        }
        if (type == EAST) {
            IntegralTrapezoid t(left_alpha, bot, top);
            trapezoids->push_back(t);
        }
    }

    void decompose_hour_glass_shape(const HourGlassCell& hour_glass, list<IntegralTrapezoid>* trapezoids) {
        for (int i = 0; i < 2; ++i) {
            if (hour_glass.part(i).empty())
                continue;

            if (hour_glass.part(i).begin_dir().dx() > 0 && hour_glass.part(i).end_dir().dx() > 0) {
                decompose_hour_glass_part(EAST, hour_glass.part(i), trapezoids);
            } else if (hour_glass.part(i).begin_dir().dx() < 0 && hour_glass.part(i).end_dir().dx() < 0) {
                decompose_hour_glass_part(WEST, hour_glass.part(i), trapezoids);
            } else {
                throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);
            }
        }
    }

    bool within_hourglass_part(const HourGlassPart& hg, const Point& q) {
        Point begin = hg.begin_ray().second_point();
        Point end = hg.end_ray().second_point();

        if ( !CGAL::left_turn(begin, hg.part().front(), q) )
            return false;

        for (int i = 0; i < hg.part().size()-1; ++i)
            if ( !CGAL::left_turn(hg.part()[i], hg.part()[i+1], q) )
                return false;

        if ( !CGAL::left_turn(hg.part().back(), end, q) )
            return false;

        return true;
    }
    /*
     bool on_hourglass_part_boundary(const HourGlassPart& hg, const Point& q) {
     if (hg.begin_ray().has_on(q) || hg.end_ray().has_on(q))
     return true;

     for (int i = 0; i < hg.part().size()-1; ++i) {
     Segment s(hg.part()[i], hg.part()[i+1]);
     if ( s.has_on(q) )
     return true;
     }

     return false;
     }
     */
    void overlap_hourglass_part_with_box(const HourGlassPart& hg,
                                         const IsoRectangle& iso,
                                         SimplePolygon* out) {
        list<Point> pts(hg.part().begin(), hg.part().end());

        NonVerticalLine<Kernel> begin_line(hg.part().front(), hg.begin_dir());
        NonVerticalLine<Kernel> end_line(hg.part().back(), hg.end_dir());

        if (hg.part().front().x() < iso.xmax() && hg.begin_dir().dx() > 0) {
            NT y = begin_line.y_at_x(iso.xmax());
            if (y > iso.ymax())
                pts.push_back(Point(begin_line.x_at_y(iso.ymax()), iso.ymax()));
            else if (y < iso.ymin())
                pts.push_back(Point(begin_line.x_at_y(iso.ymin()), iso.ymin()));
            else
                pts.push_back(Point(iso.xmax(), begin_line.y_at_x(iso.xmax())));
        }

        if (hg.part().front().x() > iso.xmin() && hg.begin_dir().dx() < 0) {
            NT y = begin_line.y_at_x(iso.xmin());
            if (y > iso.ymax())
                pts.push_back(Point(begin_line.x_at_y(iso.ymax()), iso.ymax()));
            else if (y < iso.ymin())
                pts.push_back(Point(begin_line.x_at_y(iso.ymin()), iso.ymin()));
            else
                pts.push_back(Point(iso.xmin(), begin_line.y_at_x(iso.xmin())));
        }

        if (hg.part().back().x() < iso.xmax() && hg.end_dir().dx() > 0) {
            NT y = end_line.y_at_x(iso.xmax());
            if (y > iso.ymax())
                pts.push_back(Point(end_line.x_at_y(iso.ymax()), iso.ymax()));
            else if (y < iso.ymin())
                pts.push_back(Point(end_line.x_at_y(iso.ymin()), iso.ymin()));
            else
                pts.push_back(Point(iso.xmax(), end_line.y_at_x(iso.xmax())));
        }

        if (hg.part().back().x() > iso.xmin() && hg.end_dir().dx() < 0) {
            NT y = end_line.y_at_x(iso.xmin());
            if (y > iso.ymax())
                pts.push_back(Point(end_line.x_at_y(iso.ymax()), iso.ymax()));
            else if (y < iso.ymin())
                pts.push_back(Point(end_line.x_at_y(iso.ymin()), iso.ymin()));
            else
                pts.push_back(Point(iso.xmin(), end_line.y_at_x(iso.xmin())));
        }

        for (int i = 0; i < 4; ++i)
            if (within_hourglass_part(hg, iso[i]))
                pts.push_back(iso[i]);

        if (pts.size() < 3)
            return;

        //std::cout << bg::dsv(pts) << std::endl;
        CGAL::convex_hull_2(pts.begin(), pts.end(), std::back_inserter(*out));
        if (out->size() < 3)
            throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);

        //std::cout << bg::dsv(*out) << std::endl;
    }

    void intersect_polygon_with_hourglass_part(const SimplePolygon& cell,
                                               const HourGlassPart& hg,
                                               SimplePolygon* out) {
        assert(cell.is_convex());

        Bbox bbox = CGAL::bbox_2(hg.part().begin(), hg.part().end()) + cell.bbox();
        IsoRectangle iso(NT(bbox.xmin()-1), NT(bbox.ymin()-1), NT(bbox.xmax()+1), NT(bbox.ymax()+1));
        SimplePolygon hg_cell;
        overlap_hourglass_part_with_box(hg, iso, &hg_cell);
        assert(hg_cell.is_convex());

        intersect_convex_polygons(cell, hg_cell, out);
    }

    Cell::CellType intersect_hourglass_parts(const HourGlassPart& hg1, const HourGlassPart& hg2,
                                             SimplePolygon* out_cell, HourGlassPart* out_hg) {
        if (hg1.part().empty() || hg2.part().empty())
            return Cell::EMPTY;
//        std::cout << __LINE__ << ": "  << hg1 << std::endl << std::endl;
//        std::cout << __LINE__ << ": "  << hg2 << std::endl << std::endl;

        SimplePolygon poly1, poly2;

        Bbox bbox;
        {
            bbox = CGAL::bbox_2(hg1.part().begin(), hg1.part().end()) + CGAL::bbox_2(hg2.part().begin(), hg2.part().end());

            Ray rays[4] = { hg1.begin_ray(), hg1.end_ray(), hg2.begin_ray(), hg2.end_ray() };

            for (int i = 0; i < 2; ++i)
                for (int j = 2; j < 4; ++j) {
                    CGAL::Object obj = CGAL::intersection(rays[i], rays[j]);
                    Point p;
                    if (CGAL::assign(p, obj))
                        bbox = bbox + p.bbox();
                }
        }

        IsoRectangle iso(NT(bbox.xmin()-1), NT(bbox.ymin()-1), NT(bbox.xmax()+1), NT(bbox.ymax()+1));
        overlap_hourglass_part_with_box(hg1, iso, &poly1);
        overlap_hourglass_part_with_box(hg2, iso, &poly2);
        assert(poly1.is_convex() && poly2.is_convex());

//        std::cout << __LINE__ << ": "  << poly1 << std::endl;
//        std::cout << __LINE__ << ": "  << poly2 << std::endl;

        if (!CGAL::do_intersect(poly1, poly2))
            return Cell::EMPTY;

        SimplePolygon intersection;
        intersect_convex_polygons(poly1, poly2, &intersection);
        assert(intersection.is_convex());

        if (intersection.is_empty())
            return Cell::EMPTY;

//        std::cout << __LINE__ << ": "  << intersection << std::endl;

        vector<int> corners;
        for (int i = 0; i < intersection.size(); ++i) {
            if (intersection[i].x() == iso.xmin() ||
                intersection[i].x() == iso.xmax() ||
                intersection[i].y() == iso.ymin() ||
                intersection[i].y() == iso.ymax()) {
                corners.push_back(i);
            }
        }

        if (corners.empty()) {
            *out_cell = intersection;
            return Cell::SIMPLE;
        }

        if (corners.front() == 0 && corners.back() == intersection.size()-1) {
            int i = 0;
            int dif;
            do {
                dif = corners[i+1]-corners[i];
                corners[i] += intersection.size();
                i++;
            } while (i < corners.size()-1 && dif == 1);
        }

        //        Direction out_begin_dir, out_end_dir;
        //        if (abs(hg1.begin_dir().dy()) > abs(hg2.begin_dir().dy()))
        //            out_begin_dir = hg1.begin_dir();
        //        else
        //            out_begin_dir = hg2.begin_dir();
        //
        //        if (abs(hg1.end_dir().dy()) > abs(hg2.end_dir().dy()))
        //            out_end_dir = hg1.end_dir();
        //        else
        //            out_end_dir = hg2.end_dir();

        int i_begin = *std::max_element(corners.begin(), corners.end()) + 1;
        int i_end = *std::min_element(corners.begin(), corners.end());

        Vector begin_dir = intersection[(intersection.size()+i_begin-1)%intersection.size()] - intersection[i_begin%intersection.size()];
        Vector end_dir = intersection[i_end%intersection.size()] - intersection[(intersection.size()+i_end-1)%intersection.size()];
        Direction out_begin_dir(begin_dir.x()/abs(begin_dir.x()), begin_dir.y()/abs(begin_dir.x()));
        Direction out_end_dir(end_dir.x()/abs(end_dir.x()), end_dir.y()/abs(end_dir.x()));

        vector<Point> out_part;
        for (int i = i_begin%intersection.size(); i != i_end%intersection.size();
             i = (i+1)%intersection.size())
            out_part.push_back(intersection[i]);

        *out_hg = HourGlassPart(out_part, out_begin_dir, out_end_dir);

//        std::cout << __LINE__ << ": "  << *out_hg <<std::endl;

        return Cell::HOURGLASS;
    }

    void intersect_hourglass_shapes(const HourGlassCell& h1, const HourGlassCell& h2, list<Cell>* out) {
//        std::cout << __LINE__ << ": " << h1 << std::endl << std::endl << h2 << std::endl << std::endl << std::endl << std::endl;

        Cell::CellType t[2][2];
        SimplePolygon poly[2][2];
        HourGlassPart part[2][2];

        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) {
                t[i][j] = intersect_hourglass_parts(h1.part(i), h2.part(j), &poly[i][j], &part[i][j]);
                assert(poly[i][j].is_empty() || poly[i][j].is_convex());
            }

        if ((t[0][0] == Cell::SIMPLE && t[1][1] == Cell::HOURGLASS) ||
            (t[0][0] == Cell::HOURGLASS && t[1][1] == Cell::SIMPLE) ||
            (t[0][1] == Cell::SIMPLE && t[1][0] == Cell::HOURGLASS) ||
            (t[0][1] == Cell::HOURGLASS && t[1][0] == Cell::SIMPLE)) {
//            std::cout << h1 << std::endl << std::endl << h2 << std::endl << std::endl << std::endl << std::endl;
            throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);
        }

        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                if (t[i][j] == Cell::SIMPLE)
                    out->push_back(Cell(poly[i][j]));

        if (t[0][0] == Cell::HOURGLASS || t[1][1] == Cell::HOURGLASS) {
            Cell c(HourGlassCell(part[0][0], part[1][1]));
            out->push_back(c);
        }

        if (t[0][1] == Cell::HOURGLASS || t[1][0] == Cell::HOURGLASS) {
            Cell c(HourGlassCell(part[0][1], part[1][0]));
            out->push_back(c);
        }
    }
}

HourGlassPart::HourGlassPart() {
}

HourGlassPart::HourGlassPart(const vector<Point>& part,
                             const Direction& begin_dir,
                             const Direction& end_dir)
: part_(part), begin_dir_(begin_dir), end_dir_(end_dir) {
    correct();
}

HourGlassPart::~HourGlassPart() {
}

const vector<Point>& HourGlassPart::part() const {
    return part_;
}

vector<Point>& HourGlassPart::part() {
    return part_;
}

const Direction& HourGlassPart::begin_dir() const {
    return begin_dir_;
}

Direction& HourGlassPart::begin_dir() {
    return begin_dir_;
}

const Direction& HourGlassPart::end_dir() const {
    return end_dir_;
}

Direction& HourGlassPart::end_dir() {
    return end_dir_;
}

const Ray& HourGlassPart::begin_ray() const {
    return begin_ray_;
}

const Ray& HourGlassPart::end_ray() const {
    return end_ray_;
}

bool HourGlassPart::empty() const {
    return part_.empty();
}

void HourGlassPart::correct() {
    if (part_.empty())
        throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);

    begin_ray_ = Ray(part_.front(), begin_dir_);
    end_ray_ = Ray(part_.back(), end_dir_);

    SimplePolygon poly(part_.begin(), part_.end());
    poly.insert(poly.vertices_begin(), part_.front() + begin_dir_.to_vector());
    poly.push_back(part_.back() + end_dir_.to_vector());

    if (!poly.is_simple())
        throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);

    if (!poly.is_convex()) {
//        std::cout << __LINE__ << ": " << poly << std::endl;
        throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);
    }

    if (poly.is_empty())
        throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);

    if (begin_ray_ == end_ray_ || (begin_ray_.source() != end_ray_.source() && do_intersect(begin_ray_, end_ray_))) {
//        std::cout << __LINE__ << ": "  << *this << std::endl;
        throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);
    }

    if (poly.is_clockwise_oriented()) {
        std::reverse(part_.begin(), part_.end());
        Direction old_begin = begin_dir_;
        Direction old_end = end_dir_;
        swap(begin_dir_, end_dir_);
        swap(begin_ray_, end_ray_);
        // TODO : test that swap works and delete this
        if (begin_dir_ != old_end || end_dir_ != old_begin)
            throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);
    }
}

HourGlassCell::HourGlassCell() {
}

HourGlassCell::HourGlassCell(const HourGlassPart& first, const HourGlassPart& second) {
    part_[0] = first;
    part_[1] = second;
    correct();
}

HourGlassCell::~HourGlassCell() {
}

void HourGlassCell::correct() {
    if (!part_[0].empty())
        part_[0].correct();

    if (!part_[1].empty())
        part_[1].correct();

    if (!part_[0].empty() && !part_[1].empty()) {
        Ray rays[] = { part_[0].begin_ray(),
            part_[0].end_ray(),
            part_[1].begin_ray(),
            part_[1].end_ray() };

        for (int i = 0; i < 4; ++i)
            for (int j = i+1; j < 4; ++j)
                if (rays[i] == rays[j] ||
                    (rays[i].source() != rays[j].source() && do_intersect(rays[i], rays[j])))
                    throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);
    }
}

void HourGlassCell::intersect(const SimplePolygon& cell, list<SimplePolygon>* out) const {
    assert(cell.is_convex());

    for (int i = 0; i < 2; ++i) {
        if (part_[i].empty())
            continue;

        SimplePolygon convex_poly;
        intersect_polygon_with_hourglass_part(cell, part_[i], &convex_poly);
        assert(convex_poly.is_convex());

        if (convex_poly.size() >= 3)
            out->push_back(convex_poly);
    }
}

/*
 bool HourGlassCell::within(const Point& p) const {
 SimplePolygon r;
 Point plus(1.0, 1.0);
 Point minus(1.0, -1.0);
 r.push_back(p+plus);
 r.push_back(p+minus);
 r.push_back(p-plus);
 r.push_back(p-minus);
 bg::correct(r);
 vector<SimplePolygon> out;
 intersect(r, &out);
 BOOST_FOREACH(SimplePolygon poly, out) {
 if (bg::within(p, poly))
 return true;
 }

 return false;
 }

 void HourGlassCell::get_linestrings_for_print(vector<Linestring>* out) const {
 if (!bottom_.empty()) {
 Linestring s;
 NT x1 = bottom_begin_dir_.x();  NT y1 = bottom_begin_dir_.y();
 NT x2 = bottom_end_dir_.x();    NT y2 = bottom_end_dir_.y();
 s.push_back(bottom_.front()+bottom_begin_dir_*4/sqrt(x1*x1+y1*y1));
 s.insert(s.end(), bottom_.begin(), bottom_.end());
 s.push_back(bottom_.back()+bottom_end_dir_*4/sqrt(x2*x2+y2*y2));

 out->push_back(s);
 }

 if (!top_.empty()) {
 Linestring s;
 NT x1 = top_begin_dir_.x();  NT y1 = top_begin_dir_.y();
 NT x2 = top_end_dir_.x();    NT y2 = top_end_dir_.y();
 s.push_back(top_.front()+top_begin_dir_*4/sqrt(x1*x1+y1*y1));
 s.insert(s.end(), top_.begin(), top_.end());
 s.push_back(top_.back()+top_end_dir_*4/sqrt(x2*x2+y2*y2));

 out->push_back(s);
 }
 }
 */
void Cell::Initialize(const SimplePolygon& obstacle) {
    if (obstacle.is_empty()) {
        type_ = EMPTY;
        return;
    }

    assert(obstacle.is_convex());

    type_ = HOURGLASS;
    vector<Point> lower_envelope, upper_envelope;
    get_lower_envelope(obstacle, &lower_envelope);
    get_upper_envelope(obstacle, &upper_envelope);

    NonVerticalLine<Kernel> bottom_begin_line(lower_envelope.front());
    NonVerticalLine<Kernel> bottom_end_line(lower_envelope.back());
    Direction bottom_begin_dir(-1, -bottom_begin_line.slope());
    Direction bottom_end_dir(1, bottom_end_line.slope());
    vector<Point> bottom_part;
    for (int i = 0; i < lower_envelope.size()-1; ++i) {
        NonVerticalLine<Kernel> l(lower_envelope[i], lower_envelope[i+1]);
        Point p;
        l.dual(&p);
        bottom_part.push_back(p);
    }
    HourGlassPart bottom(bottom_part, bottom_begin_dir, bottom_end_dir);

    NonVerticalLine<Kernel> top_begin_line(upper_envelope.front());
    NonVerticalLine<Kernel> top_end_line(upper_envelope.back());
    Direction top_begin_dir(1, top_begin_line.slope());
    Direction top_end_dir(-1, -top_end_line.slope());
    vector<Point> top_part;
    for (int i = 0; i < upper_envelope.size()-1; ++i) {
        NonVerticalLine<Kernel> l(upper_envelope[i], upper_envelope[i+1]);
        Point p;
        l.dual(&p);
        top_part.push_back(p);
    }
    HourGlassPart top(top_part, top_begin_dir, top_end_dir);

    hour_glass_cell_ = HourGlassCell(bottom, top);
}

void Cell::Initialize(const list<Line>& vertical_lines,
                      const list<NonVerticalLine<Kernel> >& below_separator_lines,
                      const list<NonVerticalLine<Kernel> >& above_separator_lines,
                      const Point& separator) {
    typedef list<NonVerticalLine<Kernel> >::const_iterator NonVerticalLine_iterator;
    typedef list<Line>::const_iterator Line_iterator;

    vector<Point> intersections;
    for (NonVerticalLine_iterator it = above_separator_lines.begin(); it != above_separator_lines.end(); ++it) {
        for (NonVerticalLine_iterator jt = below_separator_lines.begin(); jt != below_separator_lines.end(); ++jt) {
            CGAL::Object obj = CGAL::intersection(*it, *jt);
            Point p;
            if (!CGAL::assign(p, obj))
                throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);

            intersections.push_back(p);
        }
    }

    for (Line_iterator it = vertical_lines.begin(); it != vertical_lines.end(); ++it) {
        for (NonVerticalLine_iterator jt = below_separator_lines.begin(); jt != below_separator_lines.end(); ++jt) {
            CGAL::Object obj = CGAL::intersection(*it, *jt);
            Point p;
            if (!CGAL::assign(p, obj))
                throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);

            intersections.push_back(p);
        }
        for (NonVerticalLine_iterator jt = above_separator_lines.begin(); jt != above_separator_lines.end(); ++jt) {
            CGAL::Object obj = CGAL::intersection(*it, *jt);
            Point p;
            if (!CGAL::assign(p, obj))
                throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);

            intersections.push_back(p);
        }
    }

    if (intersections.empty()) {
        set<Point> dual_pts;
        for (NonVerticalLine_iterator jt = below_separator_lines.begin(); jt != below_separator_lines.end(); ++jt) {
            Point p;
            jt->dual(&p);
            dual_pts.insert(p);
        }
        for (NonVerticalLine_iterator jt = above_separator_lines.begin(); jt != above_separator_lines.end(); ++jt) {
            Point p;
            jt->dual(&p);
            dual_pts.insert(p);
        }

        SimplePolygon simple_cell;
        CGAL::convex_hull_2(dual_pts.begin(), dual_pts.end(), std::back_inserter(simple_cell));

        if (is_degenerate(simple_cell)) {
            type_ = EMPTY;
        } else {
            simple_cell_ = simple_cell;
            type_ = SIMPLE;
        }

        return;
    } else if (intersections.size() == 1) {
        type_ = EMPTY;
        return;
    }

    type_ = HOURGLASS;
    vector<Point> left_dual_pts, right_dual_pts;
        for (NonVerticalLine_iterator jt = below_separator_lines.begin(); jt != below_separator_lines.end(); ++jt) {
        Point p;
        jt->dual(&p);
        if (separator < *(intersections.begin())) {
            right_dual_pts.push_back(p);
        } else {
            left_dual_pts.push_back(p);
        }
    }
        for (NonVerticalLine_iterator jt = above_separator_lines.begin(); jt != above_separator_lines.end(); ++jt) {
        Point p;
        jt->dual(&p);
        if (separator < *(intersections.begin())) {
            left_dual_pts.push_back(p);
        } else {
            right_dual_pts.push_back(p);
        }
    }

    SimplePolygon ch;
    CGAL::convex_hull_2(intersections.begin(), intersections.end(), std::back_inserter(ch));

    NonVerticalLine<Kernel> right_begin(*(ch.right_vertex()));// *std::min_element(intersections.begin(), intersections.end(), compare_left_bottom);
    NonVerticalLine<Kernel> right_end(*(ch.left_vertex()));// *std::max_element(intersections.begin(), intersections.end(), compare_left_bottom);
    NonVerticalLine<Kernel> left_begin(*(right_bottom_vertex(ch)));// *std::max_element(intersections.begin(), intersections.end(), compare_left_top);
    NonVerticalLine<Kernel> left_end(*(left_top_vertex(ch)));// *std::min_element(intersections.begin(), intersections.end(), compare_left_top);

//    std::cout << right_begin << std::endl;
//    std::cout << right_end << std::endl;
//    std::cout << left_begin << std::endl;
//    std::cout << left_end << std::endl;

    HourGlassPart *part_of_left_pts = &(hour_glass_cell_.part(0));
    HourGlassPart *part_of_right_pts = &(hour_glass_cell_.part(1));
//    if (right_begin.slope() < 0) {
//        part_of_left_pts = &(hour_glass_cell_.bottom());
//        part_of_right_pts = &(hour_glass_cell_.top());
//    } else {
//        part_of_left_pts = &(hour_glass_cell_.top());
//        part_of_right_pts = &(hour_glass_cell_.bottom());
//    }

    if (!left_dual_pts.empty())
        build_hg_part(left_dual_pts, left_begin, left_end, false, part_of_left_pts);

    if (!right_dual_pts.empty())
        build_hg_part(right_dual_pts, right_begin, right_end, true, part_of_right_pts);
}

void Cell::decompose(list<IntegralTrapezoid>* trapezoids) const {
    if (type_ == SIMPLE) {
        decompose_simple_polygon(simple_cell_, trapezoids);
    } else if (type_ == HOURGLASS) {
        decompose_hour_glass_shape(hour_glass_cell_, trapezoids);
    }
}

void Cell::intersect(const Cell& c, list<Cell>* out) const {
    if (type_ == SIMPLE && c.type_ == SIMPLE) {
        SimplePolygon intersection;
        intersect_convex_polygons(simple_cell_, c.simple_cell_, &intersection);

        if (!intersection.is_empty()) {
            out->push_back(Cell(intersection));
        }
    } else if (type_ == SIMPLE && c.type_ == HOURGLASS) {
        list<SimplePolygon> convex_polys;
        c.hour_glass_cell_.intersect(simple_cell_, &convex_polys);
//        for (int i = 0; i < intersections.size(); ++i) {
//            out->push_back(Cell(intersections[i]));
//        }
        out->insert(out->end(), convex_polys.begin(), convex_polys.end());
    } else if (type_ == HOURGLASS && c.type_ == SIMPLE) {
        list<SimplePolygon> convex_polys;
        hour_glass_cell_.intersect(c.simple_cell_, &convex_polys);
//        for (list<ConvexPolygon>::iterator it = intersections.begin(); it != intersections.end(); ++it) {
//            out->push_back(Cell(*it));
//        }
        out->insert(out->end(), convex_polys.begin(), convex_polys.end());
    } else {
        intersect_hourglass_shapes(hour_glass_cell_, c.hour_glass_cell_, out);
    }
}

std::ostream& operator<<(std::ostream& os, const HourGlassPart& p) {
    os << "begin: " << p.begin_dir() << std::endl << "part: ";
    for (int i = 0; i < p.part().size(); ++i)
        os << "(" << p.part()[i] << ")";
    os << std::endl << "end: " << p.end_dir() << std::endl
    << "( rays: " << p.begin_ray() << "; " << p.end_ray() << " )" << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const HourGlassCell& c) {
    os << "Part 1: " << std::endl << c.part(0) << std::endl;
    os << "Part 2: " << std::endl << c.part(1) << std::endl;
    return os;
}


std::ostream& operator<<(std::ostream& os, const Cell& c) {
    switch (c.type()) {
        case Cell::HOURGLASS :
            os << "HOURGLASS" << std::endl;
            os << c.hour_glass_cell() << std::endl;
            break;
        case Cell::SIMPLE :
            os << "SIMPLE" << std::endl;
            os << c.simple_cell() << std::endl;
            break;
        default :
            os << "EMPTY" << std::endl;
    }

    return os;
}

/*
 void swap(ConvexPolygon& c1, ConvexPolygon& c2) {
 swap(c1, c2);
 c1.correct();
 c2.correct();
 }

 void swap(HourGlassPart& c1, HourGlassPart& c2) {
 swap(c1.begin_dir_, c2.begin_dir_);
 swap(c1.end_dir_, c2.end_dir_);
 swap(c1.begin_ray_, c2.begin_ray_);
 swap(c1.end_ray_, c2.end_ray_);
 swap(c1.part_, c2.part_);
 }

 void swap(HourGlassCell& c1, HourGlassCell& c2) {
 swap(c1.bottom_, c2.bottom_);
 swap(c1.top_, c2.top_);
 }
 
 void swap(Cell& c1, Cell& c2) {
 swap(c1.type_, c2.type_);
 swap(c1.simple_cell_, c2.simple_cell_);
 swap(c1.hour_glass_cell_, c2.hour_glass_cell_);
 }
 */