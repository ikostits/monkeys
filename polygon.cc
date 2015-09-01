/*
 * File:   polygon.cc
 * Author: irina
 *
 * Created on May 2, 2014, 1:59 PM
 */

#include "polygon.h"

#include <cmath>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include "simple_svg_1.0.0.hpp"

using std::list;
using std::vector;
using std::string;

namespace {
    
    struct FaceInfo2
    {
        FaceInfo2() {}
        int nesting_level;
        bool in_domain(){
            return nesting_level%2 == 1;
        }
    };
    
    typedef CGAL::Triangulation_vertex_base_2<Kernel>                      Vb;
    typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,Kernel>    Fbb;
    typedef CGAL::Constrained_triangulation_face_base_2<Kernel,Fbb>        Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                    TDS;
    typedef CGAL::Exact_predicates_tag                                     Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag>  CDT;
    
    void mark_domains(CDT& ct,
                      CDT::Face_handle start,
                      int index,
                      std::list<CDT::Edge>& border ) {
        if (start->info().nesting_level != -1) {
            return;
        }
        
        std::list<CDT::Face_handle> queue;
        queue.push_back(start);
        while(!queue.empty()) {
            CDT::Face_handle fh = queue.front();
            queue.pop_front();
            if(fh->info().nesting_level == -1) {
                fh->info().nesting_level = index;
                for(int i = 0; i < 3; i++) {
                    CDT::Edge e(fh,i);
                    CDT::Face_handle n = fh->neighbor(i);
                    if(n->info().nesting_level == -1) {
                        if(ct.is_constrained(e))
                            border.push_back(e);
                        else
                            queue.push_back(n);
                    }
                }
            }
        }
    }
    
    //explore set of facets connected with non constrained edges,
    //and attribute to each such set a nesting level.
    //We start from facets incident to the infinite vertex, with a nesting
    //level of 0. Then we recursively consider the non-explored facets incident
    //to constrained edges bounding the former set and increase the nesting level by 1.
    //Facets in the domain are those with an odd nesting level.
    void mark_domains(CDT& cdt) {
        for (CDT::All_faces_iterator it = cdt.all_faces_begin();
             it != cdt.all_faces_end(); ++it) {
            it->info().nesting_level = -1;
        }
        
        std::list<CDT::Edge> border;
        mark_domains(cdt, cdt.infinite_face(), 0, border);
        while(!border.empty()) {
            CDT::Edge e = border.front();
            border.pop_front();
            CDT::Face_handle n = e.first->neighbor(e.second);
            if(n->info().nesting_level == -1) {
                mark_domains(cdt, n, e.first->info().nesting_level+1, border);
            }
        }
    }
    
    void insert_polygon(CDT& cdt, const SimplePolygon& polygon) {
        if (polygon.is_empty())
            return;
        
//        CDT::Vertex_handle v_prev = cdt.insert(*CGAL::cpp11::prev(polygon.vertices_end()));
        SimplePolygon::Vertex_iterator tmp_it = polygon.vertices_end();
        tmp_it--;
        CDT::Vertex_handle v_prev = cdt.insert(*tmp_it);
        for (SimplePolygon::Vertex_iterator vit = polygon.vertices_begin();
             vit != polygon.vertices_end(); ++vit) {
            CDT::Vertex_handle vh = cdt.insert(*vit);
            cdt.insert_constraint(vh,v_prev);
            v_prev=vh;
        }
    }

    void get_triangulation(const Polygon& poly, CDT* cdt) {
        insert_polygon(*cdt, poly.outer_boundary());
        for (Polygon::Hole_const_iterator hole = poly.holes_begin();
             hole != poly.holes_end(); ++hole)
            insert_polygon(*cdt, *hole);

        //Mark facets that are inside the domain bounded by the polygon
        mark_domains(*cdt);
    }

    bool merge(const SimplePolygon& p1, const SimplePolygon& p2, SimplePolygon* result) {
        for (SimplePolygon::Edge_const_iterator it = p1.edges_begin(); it != p1.edges_end(); ++it)
            for (SimplePolygon::Edge_const_iterator jt = p2.edges_begin(); jt != p2.edges_end(); ++jt) {
                if ((*it) == jt->opposite()) {
                    Polygon tmp;
                    if (CGAL::join(p1, p2, tmp) && tmp.outer_boundary().is_convex()) {
                        *result = tmp.outer_boundary();
                        return true;
                    }
                }
            }

        return false;
    }

    void merge(list<SimplePolygon> *convex_polys) {
        for (list<SimplePolygon>::iterator it = convex_polys->begin(); it != convex_polys->end(); ++it)
            for (list<SimplePolygon>::iterator jt = it; jt != convex_polys->end(); ++jt) {
                if (it == jt)
                    continue;

                SimplePolygon p;
                if (merge(*it, *jt, &p)) {
                    convex_polys->erase(it);
                    convex_polys->erase(jt);
                    convex_polys->push_back(p);
                    it = convex_polys->begin();
                    jt = convex_polys->begin();
                }
            }
    }
}
/*
 SimplePolygon::SimplePolygon() : CGAL::Polygon_2<Kernel>() {
 }
 
 SimplePolygon::~SimplePolygon() {
 }
 */
/*
 Polygon::Polygon() : CGAL::Polygon_with_holes_2<Kernel>() {
 }
 
 Polygon::Polygon(const SimplePolygon& p) : CGAL::Polygon_with_holes_2<Kernel>(p) {
 }
 
 Polygon::~Polygon() {
 }
 */

//ConvexPolygon::ConvexPolygon() : SimplePolygon()  {
//}
//
//ConvexPolygon::ConvexPolygon(const SimplePolygon& poly) : SimplePolygon(poly) {
//    this->correct();
//}
//
//ConvexPolygon::~ConvexPolygon() {
//}
//
//void ConvexPolygon::correct() {
//    if (!this->is_simple() || !this->is_convex() || this->is_empty())
//        throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);
//
//    if (this->is_clockwise_oriented())
//        this->reverse_orientation();
//    
//    //this->init_envelopes_ends();
//}

/*
 const ConvexPolygon::Vertex_const_circulator& ConvexPolygon::lower_envelope_begin() const {
 return lower_envelope_begin_;
 }
 
 const ConvexPolygon::Vertex_const_circulator& ConvexPolygon::upper_envelope_end() const {
 return upper_envelope_end_;
 }
 
 const ConvexPolygon::Vertex_const_circulator& ConvexPolygon::lower_envelope_end() const {
 return lower_envelope_end_;
 }
 
 const ConvexPolygon::Vertex_const_circulator& ConvexPolygon::upper_envelope_begin() const {
 return upper_envelope_begin_;
 }
 */

//ConvexPolygon::Vertex_const_iterator ConvexPolygon::left_top_vertex() const {
//    ConvexPolygon::Vertex_const_iterator left = left_vertex();
//    if (left == vertices_end())
//        return left;
//    
//    while(true) {
//        ConvexPolygon::Vertex_const_iterator it = left;
//        if (it == vertices_begin())
//            it = vertices_end();
//        
//        it--;
//        if (it->x() != left->x() ||
//            it->y() > left->y() ||
//            it == left_vertex())
//            return left;
//        
//        left = it;
//    }
//}
//
//ConvexPolygon::Vertex_const_iterator ConvexPolygon::right_bottom_vertex() const {
//    ConvexPolygon::Vertex_const_iterator right = right_vertex();
//    if (right == vertices_end())
//        return right;
//    
//    while(true) {
//        ConvexPolygon::Vertex_const_iterator it = right + 1;
//        if (it == vertices_end())
//            it = vertices_begin();
//        
//        if (it->x() != right->x() ||
//            it->y() < right->y() ||
//            it == right_vertex())
//            return right;
//        
//        right = it;
//    }
//}
//
//void ConvexPolygon::intersect(const ConvexPolygon& cell, ConvexPolygon* out) const {
//    list<Polygon> intersections;
//    CGAL::intersection(*this, cell, std::back_inserter(intersections));
//    if (intersections.empty())
//        return;
//    if (intersections.size() > 1 || intersections.front().has_holes())
//        throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);
//    
//    *out = ConvexPolygon(intersections.front().outer_boundary().vertices_begin(),
//                         intersections.front().outer_boundary().vertices_end());
//}
//
//void ConvexPolygon::get_lower_envelope(vector<Point>* out) const {
//    if (this->is_empty())
//        return;
//    
//    Vertex_const_iterator it = this->left_vertex();
//    
//    Point cur;
//    
//    do {
//        cur = *it;
//        out->push_back(cur);
//        
//        it++;
//        if (it == this->vertices_end())
//            it = this->vertices_begin();
//    } while (cur.x() < it->x());
//}
//
//void ConvexPolygon::get_upper_envelope(vector<Point>* out) const {
//    if (this->is_empty())
//        return;
//    
//    Vertex_const_iterator it = this->right_vertex();
//    
//    Point cur;
//    
//    do {
//        cur = *it;
//        out->push_back(cur);
//        
//        it++;
//        if (it == this->vertices_end())
//            it = this->vertices_begin();
//    } while (cur.x() > it->x());
//    
//    std::reverse(out->begin(), out->end());
//}
//
//bool ConvexPolygon::is_degenerate() const {
//    for (int i = 0; i < this->size()-2; ++i)
//        if (!CGAL::collinear(this->vertex(i), this->vertex(i+1), this->vertex(i+2)))
//            return false;
//    
//    return true;
//}
//
//namespace CGAL {
//    bool do_intersect(const ConvexPolygon& p1, const ConvexPolygon& p2) {
//        return do_intersect(static_cast<const SimplePolygon&>(p1),
//                            static_cast<const SimplePolygon&>(p2));
//    }
//}

SimplePolygon::Vertex_const_iterator left_top_vertex(const SimplePolygon& convex_poly) {
    assert(convex_poly.is_empty() || convex_poly.is_convex());

    SimplePolygon::Vertex_const_iterator left = convex_poly.left_vertex();
    if (left == convex_poly.vertices_end())
        return left;

    while(true) {
        SimplePolygon::Vertex_const_iterator it = left;
        if (it == convex_poly.vertices_begin())
            it = convex_poly.vertices_end();

        it--;
        if (it->x() != left->x() ||
            it->y() > left->y() ||
            it == convex_poly.left_vertex())
            return left;

        left = it;
    }
}

SimplePolygon::Vertex_const_iterator right_bottom_vertex(const SimplePolygon& convex_poly) {
    assert(convex_poly.is_empty() || convex_poly.is_convex());

    SimplePolygon::Vertex_const_iterator right = convex_poly.right_vertex();
    if (right == convex_poly.vertices_end())
        return right;

    while(true) {
        SimplePolygon::Vertex_const_iterator it = right + 1;
        if (it == convex_poly.vertices_end())
            it = convex_poly.vertices_begin();

        if (it->x() != right->x() ||
            it->y() < right->y() ||
            it == convex_poly.right_vertex())
            return right;

        right = it;
    }
}

void intersect_convex_polygons(const SimplePolygon& poly1, const SimplePolygon& poly2, SimplePolygon* out) {
    if (poly1.is_empty() || poly2.is_empty())
        return;

    assert(poly1.is_convex() && poly2.is_convex());

    list<Polygon> convex_polys;
    CGAL::intersection(poly1, poly2, std::back_inserter(convex_polys));
    if (convex_polys.empty())
        return;
    if (convex_polys.size() > 1 || convex_polys.front().has_holes())
        throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);

    *out = SimplePolygon(convex_polys.front().outer_boundary().vertices_begin(),
                         convex_polys.front().outer_boundary().vertices_end());
}

void get_lower_envelope(const SimplePolygon& convex_poly, std::vector<Point>* out) {
    if (convex_poly.is_empty())
        return;

    assert(convex_poly.is_convex());

    SimplePolygon::Vertex_const_iterator it = convex_poly.left_vertex();

    Point cur;

    do {
        cur = *it;
        out->push_back(cur);

        it++;
        if (it == convex_poly.vertices_end())
            it = convex_poly.vertices_begin();
    } while (cur.x() < it->x());
}

void get_upper_envelope(const SimplePolygon& convex_poly, std::vector<Point>* out) {
    if (convex_poly.is_empty())
        return;

    assert(convex_poly.is_convex());

    SimplePolygon::Vertex_const_iterator it = convex_poly.right_vertex();

    Point cur;

    do {
        cur = *it;
        out->push_back(cur);

        it++;
        if (it == convex_poly.vertices_end())
            it = convex_poly.vertices_begin();
    } while (cur.x() > it->x());

    std::reverse(out->begin(), out->end());
}

bool is_degenerate(const SimplePolygon& convex_poly) {
    for (int i = 0; i < convex_poly.size()-2; ++i)
        if (!CGAL::collinear(convex_poly.vertex(i), convex_poly.vertex(i+1), convex_poly.vertex(i+2)))
            return false;

    return true;
}

bool do_intersect_interiors(const SimplePolygon& poly1, const SimplePolygon& poly2) {
    if (!CGAL::do_intersect(poly1, poly2))
        return false;

    SimplePolygon inters;
    intersect_convex_polygons(poly1, poly2, &inters);
    if (is_degenerate(inters))
        return false;

    return true;
}

/*
 void ConvexPolygon::init_envelopes_ends() {
 Vertex_const_circulator start = this->vertices_circulator();
 
 lower_envelope_begin_ = start;
 upper_envelope_end_ = start;
 lower_envelope_end_ = start;
 upper_envelope_begin_ = start;
 
 Vertex_const_circulator cur = start; cur++;
 
 while (cur != start) {
 Point p = *cur;
 
 if (p < *lower_envelope_begin_)
 lower_envelope_begin_ = cur;
 
 if (p.x() < upper_envelope_end_->x() ||
 (p.x() == upper_envelope_end_->x() && p.y() > upper_envelope_end_->y()))
 upper_envelope_end_ = cur;
 
 if (p.x() > lower_envelope_end_->x() ||
 (p.x() == lower_envelope_end_->x() && p.y() < lower_envelope_end_->y()))
 lower_envelope_end_ = cur;
 
 if (p > *upper_envelope_begin_)
 upper_envelope_begin_ = cur;
 
 cur++;
 }
 }
*/

void triangulate(const Polygon& poly, list<Triangle>* triangles) {
    CDT cdt;
    
    get_triangulation(poly, &cdt);
    
    for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
         fit != cdt.finite_faces_end(); ++fit) {
        if (fit->info().in_domain()) {
            CDT::Face_handle face = fit;
            triangles->push_back(cdt.triangle(face));
        }
    }
}

void partition(const Polygon& poly, list<SimplePolygon>* convex_polys) {
//    CGAL::Partition_traits_2<Kernel> partition_traits;
//    CGAL::Partition_is_valid_traits_2<CGAL::Partition_traits_2<Kernel>, CGAL::Is_convex_2<CGAL::Partition_traits_2<Kernel> > > validity_traits;
    if (poly.number_of_holes() == 0 && poly.outer_boundary().is_convex()) {
        convex_polys->push_back(poly.outer_boundary());
        return;
    }

    list<CGAL::Partition_traits_2<Kernel>::Polygon_2> partition_polys;
    if (poly.number_of_holes() == 0) {
        CGAL::optimal_convex_partition_2(poly.outer_boundary().vertices_begin(),
                                         poly.outer_boundary().vertices_end(),
                                         std::back_inserter(partition_polys));
        if (!CGAL::partition_is_valid_2(poly.outer_boundary().vertices_begin(),
                                        poly.outer_boundary().vertices_end(),
                                        partition_polys.begin(),
                                        partition_polys.end()))
            throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);

        for (list<CGAL::Partition_traits_2<Kernel>::Polygon_2>::iterator it = partition_polys.begin();
             it != partition_polys.end(); ++it) {
            convex_polys->push_back(SimplePolygon(it->vertices_begin(), it->vertices_end()));
        }
    } else {
        CDT cdt;

        get_triangulation(poly, &cdt);

        list<SimplePolygon> partition_polys;
        for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
             fit != cdt.finite_faces_end(); ++fit) {
            if (fit->info().in_domain()) {
                CDT::Face_handle face = fit;
                Triangle t = cdt.triangle(face);
                SimplePolygon p;
                p.push_back(t.vertex(0));
                p.push_back(t.vertex(1));
                p.push_back(t.vertex(2));
                partition_polys.push_back(p);
            }
        }

        merge(&partition_polys);
        convex_polys->insert(convex_polys->end(), partition_polys.begin(), partition_polys.end());
    }
}

void generate_polygon(double center_x, double center_y, double r_min, double r_max, SimplePolygon* out) {
    out->clear();

    int n_poly = 1+ceil(M_PI/acos(r_min/r_max));

    if (n_poly < 3)
        throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);

    double rand_start_angle = rand___.get_double(0, M_PI_4);
    for (int j = 0; j < n_poly; ++j) {
        double x(r_max*cos(rand_start_angle + 2.0*j*M_PI/n_poly));
        double y(r_max*sin(rand_start_angle + 2.0*j*M_PI/n_poly));
//        Point pt(center_x+x, center_y+y);
        Point pt(center_x+round(100*x)/100, center_y+round(100*y)/100);
        out->push_back(pt);
    }
}

NT area(const Polygon& poly) {
    NT s = poly.outer_boundary().area();
    for (Polygon::Hole_const_iterator hole = poly.holes_begin(); hole != poly.holes_end(); ++hole)
        s += hole->area();
    return s;
}

void print_to_file(const SimplePolygon& poly, std::string file_name) {
    svg::Dimensions dimensions(100, 100);
    svg::Document doc(file_name, svg::Layout(dimensions, svg::Layout::BottomLeft));

    Bbox box = poly.bbox();
    double scale = 100.0/std::max(CGAL::to_double(box.xmax()-box.xmin()),
                                  CGAL::to_double(box.ymax()-box.ymin()));
    svg::Polygon svg_poly(svg::Color(200, 160, 220), svg::Stroke(.5, svg::Color(150, 160, 200)));
    for (SimplePolygon::Vertex_const_iterator it = poly.vertices_begin();
         it != poly.vertices_end(); ++it)
        svg_poly << svg::Point(scale*CGAL::to_double(it->x()-box.xmin()),
                               scale*CGAL::to_double(it->y()-box.ymin()));
    doc << svg_poly;
    doc.save();
}

void print_to_file(const Polygon& poly, std::string file_name) {
    svg::Dimensions dimensions(100, 100);
    svg::Document doc(file_name, svg::Layout(dimensions, svg::Layout::BottomLeft));

    Bbox box = CGAL::bbox_2(poly.outer_boundary().vertices_begin(), poly.outer_boundary().vertices_end());
    double scale = 100.0/std::max(CGAL::to_double(box.xmax()-box.xmin()),
                                  CGAL::to_double(box.ymax()-box.ymin()));
    svg::Polygon svg_poly(svg::Color(200, 160, 220), svg::Stroke(.5, svg::Color(150, 160, 200)));
    for (SimplePolygon::Vertex_const_iterator it = poly.outer_boundary().vertices_begin();
         it != poly.outer_boundary().vertices_end(); ++it)
        svg_poly << svg::Point(scale*CGAL::to_double(it->x()-box.xmin()),
                               scale*CGAL::to_double(it->y()-box.ymin()));
    doc << svg_poly;
    
    for (Polygon::Hole_const_iterator it = poly.holes_begin(); it != poly.holes_end(); ++it) {
        svg::Polygon svg_hole(svg::Color(svg::Color::White), svg::Stroke(.5, svg::Color(150, 160, 200)));
        for (SimplePolygon::Vertex_const_iterator jt = it->vertices_begin(); jt != it->vertices_end(); ++jt)
            svg_hole << svg::Point(scale*CGAL::to_double(jt->x()-box.xmin()),
                                   scale*CGAL::to_double(jt->y()-box.ymin()));
        doc << svg_hole;
    }
    
    
    doc.save();
}

void print_to_file(const SimplePolygon& poly1, const SimplePolygon& poly2, const list<SimplePolygon>& obs, std::string file_name) {
    svg::Dimensions dimensions(1000, 1000);
    svg::Document doc(file_name, svg::Layout(dimensions, svg::Layout::BottomLeft));

    Bbox box(-50, -50, 50, 50);
    double scale = 1000.0/std::max(CGAL::to_double(box.xmax()-box.xmin()),
                                  CGAL::to_double(box.ymax()-box.ymin()));
    svg::Polygon svg_poly1(svg::Color(200, 160, 220), svg::Stroke(.5, svg::Color(150, 160, 200)));
    for (SimplePolygon::Vertex_const_iterator it = poly1.vertices_begin(); it != poly1.vertices_end(); ++it)
        svg_poly1 << svg::Point(scale*CGAL::to_double(it->x()-box.xmin()),
                                scale*CGAL::to_double(it->y()-box.ymin()));
    doc << svg_poly1;

    svg::Polygon svg_poly2(svg::Color(160, 200, 220), svg::Stroke(.5, svg::Color(110, 200, 200)));
    for (SimplePolygon::Vertex_const_iterator it = poly2.vertices_begin(); it != poly2.vertices_end(); ++it)
        svg_poly2 << svg::Point(scale*CGAL::to_double(it->x()-box.xmin()),
                                scale*CGAL::to_double(it->y()-box.ymin()));
    doc << svg_poly2;

    for (list<SimplePolygon>::const_iterator it = obs.begin(); it != obs.end(); ++it) {
        svg::Polygon svg_hole(svg::Color(133,94,66), svg::Stroke(.5, svg::Color(133,94,66)));
        for (SimplePolygon::Vertex_const_iterator jt = it->vertices_begin(); jt != it->vertices_end(); ++jt)
            svg_hole << svg::Point(scale*CGAL::to_double(jt->x()-box.xmin()),
                                   scale*CGAL::to_double(jt->y()-box.ymin()));
        doc << svg_hole;
    }

    doc.save();
}

void print_to_file(const list<SimplePolygon>& poly1, const list<SimplePolygon>& poly2, const list<SimplePolygon>& obs, std::string file_name) {
    svg::Dimensions dimensions(1000, 1000);
    svg::Document doc(file_name, svg::Layout(dimensions, svg::Layout::BottomLeft));

    Bbox box(-50, -50, 50, 50);
    double scale = 1000.0/std::max(CGAL::to_double(box.xmax()-box.xmin()),
                                   CGAL::to_double(box.ymax()-box.ymin()));

    for (list<SimplePolygon>::const_iterator it = poly1.begin(); it != poly1.end(); ++it) {
        svg::Polygon svg_poly1(svg::Fill(svg::Color(153, 204, 0), 0.75), svg::Stroke(.5, svg::Color(153, 204, 0)));
        for (SimplePolygon::Vertex_const_iterator jt = it->vertices_begin(); jt != it->vertices_end(); ++jt)
            svg_poly1 << svg::Point(scale*CGAL::to_double(jt->x()-box.xmin()),
                                    scale*CGAL::to_double(jt->y()-box.ymin()));
        doc << svg_poly1;
    }

    for (list<SimplePolygon>::const_iterator it = poly2.begin(); it != poly2.end(); ++it) {
        svg::Polygon svg_poly2(svg::Fill(svg::Color(0, 204, 153), 0.75), svg::Stroke(.5, svg::Color(0, 204, 153)));
        for (SimplePolygon::Vertex_const_iterator jt = it->vertices_begin(); jt != it->vertices_end(); ++jt)
            svg_poly2 << svg::Point(scale*CGAL::to_double(jt->x()-box.xmin()),
                                    scale*CGAL::to_double(jt->y()-box.ymin()));
        doc << svg_poly2;
    }

    for (list<SimplePolygon>::const_iterator it = obs.begin(); it != obs.end(); ++it) {
        svg::Polygon svg_hole(svg::Color(133,94,66), svg::Stroke(.5, svg::Color(133,94,66)));
        for (SimplePolygon::Vertex_const_iterator jt = it->vertices_begin(); jt != it->vertices_end(); ++jt)
            svg_hole << svg::Point(scale*CGAL::to_double(jt->x()-box.xmin()),
                                   scale*CGAL::to_double(jt->y()-box.ymin()));
        doc << svg_hole;
    }
    
    doc.save();
}

void print_to_file(double x1, double y1,
                   double x2, double y2,
                   double sigma,
                   const list<ProbPolygon>& poly1,
                   const list<ProbPolygon>& poly2,
                   double opacity,
                   const list<SimplePolygon>& obs,
                   string file_name) {
    svg::Dimensions dimensions(1000, 1000);
    svg::Document doc(file_name, svg::Layout(dimensions, svg::Layout::BottomLeft));

    Bbox box(-50, -50, 50, 50);
    if (!obs.empty()) {
        box = obs.front().bbox();
        for (list<SimplePolygon>::const_iterator it = obs.begin(); it != obs.end(); ++it)
            box = box + it->bbox();
    }
    for (list<ProbPolygon>::const_iterator it = poly1.begin(); it != poly1.end(); ++it) {
        box = box + it->poly.bbox();
    }
    for (list<ProbPolygon>::const_iterator it = poly2.begin(); it != poly2.end(); ++it) {
        box = box + it->poly.bbox();
    }
    double scale = 1000.0/std::max(CGAL::to_double(box.xmax()-box.xmin()),
                                   CGAL::to_double(box.ymax()-box.ymin()));

    for (list<ProbPolygon>::const_iterator it = poly1.begin(); it != poly1.end(); ++it) {
        svg::Polygon svg_poly1(svg::Fill(svg::Color(204, 0, 153), opacity), svg::Stroke(.5, svg::Color(svg::Color::Transparent)));
        for (SimplePolygon::Vertex_const_iterator jt = it->poly.vertices_begin(); jt != it->poly.vertices_end(); ++jt)
            svg_poly1 << svg::Point(scale*CGAL::to_double(jt->x()-box.xmin()),
                                    scale*CGAL::to_double(jt->y()-box.ymin()));
        doc << svg_poly1;
    }

    for (list<ProbPolygon>::const_iterator it = poly2.begin(); it != poly2.end(); ++it) {
        svg::Polygon svg_poly2(svg::Fill(svg::Color(0, 204, 153), opacity), svg::Stroke(.5, svg::Color(svg::Color::Transparent)));
        for (SimplePolygon::Vertex_const_iterator jt = it->poly.vertices_begin(); jt != it->poly.vertices_end(); ++jt)
            svg_poly2 << svg::Point(scale*CGAL::to_double(jt->x()-box.xmin()),
                                    scale*CGAL::to_double(jt->y()-box.ymin()));
        doc << svg_poly2;
    }

    for (list<SimplePolygon>::const_iterator it = obs.begin(); it != obs.end(); ++it) {
        svg::Polygon svg_hole(svg::Color(133,94,66), svg::Stroke(.5, svg::Color(133,94,66)));
        for (SimplePolygon::Vertex_const_iterator jt = it->vertices_begin(); jt != it->vertices_end(); ++jt)
            svg_hole << svg::Point(scale*CGAL::to_double(jt->x()-box.xmin()),
                                   scale*CGAL::to_double(jt->y()-box.ymin()));
        doc << svg_hole;
    }

    svg::Circle svg_circle1(svg::Point(scale*(x1-box.xmin()), scale*(y1-box.ymin())), scale*sigma*2, svg::Fill(svg::Color(svg::Color::Transparent)), svg::Stroke(2, svg::Color(0, 204, 153)));
    doc << svg_circle1;
    svg::Circle svg_circle2(svg::Point(scale*(x2-box.xmin()), scale*(y2-box.ymin())), scale*sigma*2, svg::Fill(svg::Color(svg::Color::Transparent)), svg::Stroke(2, svg::Color(204, 0, 153)));
    doc << svg_circle2;

    doc.save();
}
