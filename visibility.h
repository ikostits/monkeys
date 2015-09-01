/* 
 * File:   visibility.h
 * Author: irina
 *
 * Created on May 7, 2014, 1:15 PM
 */

#ifndef VISIBILITY_H
#define	VISIBILITY_H

#include "geometry.h"
#include "polygon.h"
#include "precision.h"

class Obstacles;

mpfr_class volume(const SimplePolygon& p1, const SimplePolygon& p2, const Obstacles& obstacles);
mpfr_class volume(const Segment (&segs)[4], const Obstacles& obstacles);

namespace visibility_intl {
bool intersects_in_order(const Line& l, const Segment (&segs)[4]);
bool intersects_in_order(const InexactLine& l, const InexactSegment (&segs)[4]);
}
#endif	/* VISIBILITY_H */

