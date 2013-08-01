#ifndef __POINT_H__
#define __POINT_H__

#include <fstream>
using namespace std;

// coordinate
class Point{
public:
	Point();
	Point(long x, long y, long z);
	void set(long, long, long);
	friend bool operator == (Point & a, Point & b);
	friend bool operator != (Point & a, Point & b);
	friend ostream & operator << (ostream & , const Point & );
	long x,y,z;
};

inline void Point::set(long _x, long _y, long _z){
	x=_x;
	y=_y;
	z=_z;
}

#endif

