#ifndef __NODE_H__
#define __NODE_H__

#include <string>
#include <algorithm>
#include "global.h"
#include "point.h"
#include "net.h"
#include <vector>
#include <iostream>
//#include "block.h"
using namespace std;

class Net;
class Circuit;
// a Node in the network
class Node{
public:
	// member functions
	Node();
	~Node();
	Node(string name, Point _pt, int flag=-1, double v=0.0);
	Node(const Node & nd);
	Node & operator = (const Node & nd);

	int isS() const;
	bool is_ground() const;

	double get_value() const;
	void set_value(double v);

	bool is_mergeable() const;

	friend ostream & operator << (ostream & os, const Node & node);
	friend class Circuit;
	friend class Block;
	friend class Parser;

	////////////////////////////////////////////////////////////
	// member variables
	string name;		// node name
	Point pt;		// coordinate
	vector<Net *> nbr_vec;

	long rid;		// id in rep_list
	// to record whether this node is
	// a boundary node or internal node
	// flag_bd = 1, bd node, else internal node
	int flag_bd;
	// if =1, internal bd node,
	// if =0, general internal node.
	int internal_bd;
	// stores whether 
	string ckt_name;

private:
	double value;		// voltage
	// flag = 1 --> X
	// flag = 2 --> Y
	// flag = 3 --> Z
	int flag;		// mark if the node is an X
	Node * rep;		// representative, if somewhere is short-circuit

	Node * end[4];		// south - north (west-east) ends
	double eqvr[4];		// equivalent resisotrs
};      	

inline int Node::isS() const{return flag;}
// use a tricky way to speed up
inline bool Node::is_ground() const{return name == "0";}
inline double Node::get_value() const{return value;}
inline void Node::set_value(double v){value = v;}

#endif
