/* 
 * File:   dual_representation.h
 * Author: irina
 *
 * Created on May 7, 2014, 3:44 PM
 */

#ifndef DUAL_REPRESENTATION_H
#define	DUAL_REPRESENTATION_H

#include <vector>

#include "geometry.h"
#include "polygon.h"

class HourGlassPart {
  //friend void swap(HourGlassPart& c1, HourGlassPart& c2);
 public:
  HourGlassPart();
  HourGlassPart(const std::vector<Point>& part,
                const Direction& begin_dir,
                const Direction& end_dir);
  virtual ~HourGlassPart();

  const std::vector<Point>& part() const;
  std::vector<Point>& part();
  const Direction& begin_dir() const;
  Direction& begin_dir();
  const Direction& end_dir() const;
  Direction& end_dir();
  const Ray& begin_ray() const;
  const Ray& end_ray() const;
  
  bool empty() const;

  void correct();
 private:
  std::vector<Point> part_;
  Direction begin_dir_;
  Direction end_dir_;

  Ray begin_ray_;
  Ray end_ray_;
};

class HourGlassCell {
 //public:
  //friend void swap(HourGlassCell& c1, HourGlassCell& c2);
 public:
  HourGlassCell();
  HourGlassCell(const HourGlassPart& first, const HourGlassPart& second);
  virtual ~HourGlassCell();

  void correct();

  const HourGlassPart& part(int index) const { return part_[index]; }
  HourGlassPart& part(int index) { return part_[index]; }

  // gets the outside free space of the hourglass intersected with poly
  void intersect(const SimplePolygon& cell, std::list<SimplePolygon>* out) const;
 private:
  HourGlassPart part_[2];
};

class IntegralTrapezoid {
 public:
  IntegralTrapezoid(NT alpha_min,
                    const NonVerticalLine<Kernel>& bottom, const NonVerticalLine<Kernel>& top,
                    NT alpha_max)
        : minus_infinity_(false), plus_infinity_(false),
          alpha_min_(alpha_min), alpha_max_(alpha_max),
          bottom_(bottom), top_(top) {}
  IntegralTrapezoid(const NonVerticalLine<Kernel>& bottom, const NonVerticalLine<Kernel>& top, NT alpha_max)
        : minus_infinity_(true), plus_infinity_(false), alpha_max_(alpha_max),
          bottom_(bottom), top_(top) {}
  IntegralTrapezoid(NT alpha_min, const NonVerticalLine<Kernel>& bottom, const NonVerticalLine<Kernel>& top)
        : minus_infinity_(false), plus_infinity_(true), alpha_min_(alpha_min),
          bottom_(bottom), top_(top) {}
  virtual ~IntegralTrapezoid() {}
  NT alpha_min() const { return alpha_min_; }
  NT alpha_max() const { return alpha_max_; }
  bool minus_infinity() const { return minus_infinity_; }
  bool plus_infinity() const { return plus_infinity_; }
  const NonVerticalLine<Kernel>& bottom() const { return bottom_; }
  const NonVerticalLine<Kernel>& top() const { return top_; }
 private:
  bool minus_infinity_;
  bool plus_infinity_;
  NT alpha_min_;
  NT alpha_max_;
  NonVerticalLine<Kernel> bottom_;
  NonVerticalLine<Kernel> top_;
};

class Cell {
  //friend void swap(Cell& c1, Cell& c2);
 public:
  enum CellType { SIMPLE, HOURGLASS, EMPTY };

  Cell() : type_(EMPTY) {}
  Cell(const SimplePolygon& cell) : type_(SIMPLE), simple_cell_(cell) {}
  Cell(const HourGlassCell& hg) : type_(HOURGLASS), hour_glass_cell_(hg) {}
  virtual ~Cell() {}

  void Initialize(const SimplePolygon& obstacle);

    template<class LineIterator>
    void Initialize(LineIterator first, LineIterator last, const Point& separator) {
        std::list<Line> vertical_lines;
        std::list<NonVerticalLine<Kernel> > below_separator_lines;
        std::list<NonVerticalLine<Kernel> > above_separator_lines;

        for (LineIterator it = first; it != last; ++it) {
            if (it->is_vertical()) {
                vertical_lines.push_back(*it);
                continue;
            }

            NonVerticalLine<Kernel> tmp(it->point(0), it->point(1));
            if (tmp.y_at_x(separator.x()) < separator.y()) {
                below_separator_lines.push_back(tmp);
            } else {
                above_separator_lines.push_back(tmp);
            }
        }

        Initialize(vertical_lines, below_separator_lines, above_separator_lines, separator);
    }

  void decompose(std::list<IntegralTrapezoid>* trapezoids) const;
  void intersect(const Cell& c, std::list<Cell>* out) const;

  const CellType& type() const { return type_; }
  CellType& type() { return type_; }
  const SimplePolygon& simple_cell() const { return simple_cell_; }
  SimplePolygon& simple_cell() { return simple_cell_; }
  const HourGlassCell& hour_glass_cell() const { return hour_glass_cell_; }
  HourGlassCell& hour_glass_cell() { return hour_glass_cell_; }
private:
    void Initialize(const std::list<Line>& vertical_lines,
                    const std::list<NonVerticalLine<Kernel> >& below_separator_lines,
                    const std::list<NonVerticalLine<Kernel> >& above_separator_lines,
                    const Point& separator);
 private:
  CellType type_;
  SimplePolygon simple_cell_;
  HourGlassCell hour_glass_cell_;
};

std::ostream& operator<<(std::ostream& os, const HourGlassPart& p);
std::ostream& operator<<(std::ostream& os, const HourGlassCell& c);
std::ostream& operator<<(std::ostream& os, const Cell& c);

#endif	/* DUAL_REPRESENTATION_H */

