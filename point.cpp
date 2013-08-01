#include "point.h"

Point::Point():x(0),y(0),z(0){}

Point::Point(long _x,long _y,long _z):x(_x),y(_y),z(_z){}

bool operator == (Point & a, Point & b){
	if (a.x == b.x &&
	    a.y == b.y &&
	    a.z == b.z)
		return true;
	return false;
}

bool operator != (Point & a, Point & b){
	return !(a==b);
}

ostream & operator << (ostream & os, const Point & pt){
	os<<"("<<pt.x<<","<<pt.y<<","<<pt.z<<")";
	return os;
}
