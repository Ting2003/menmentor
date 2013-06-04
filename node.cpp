#include <iomanip>
#include "global.h"
#include "node.h"
using namespace std;

// empty constructor
Node::Node():name(""),rid(0),
	value(0.0), flag(-1), rep(NULL){
	nbr_vec.clear();
	for(int i=0;i<4;i++){
		eqvr[i]=0.0;
		end[i]=this;
	}
	rid = -1;
	flag_bd = 0;
	internal_bd = 0;
}

Node::~Node (){
	nbr_vec.clear();
	delete [] rep;
	delete [] end;
	delete [] eqvr;
}

Node::Node(string n, int x, double v): 
	name(n), rid(0), 
	value(v), flag(x), rep(NULL) {
	nbr_vec.clear();
	for(int i=0;i<4;i++){
		eqvr[i]=0.0;
		end[i]=this;
	}
	rid = -1;
	flag_bd = 0;
	internal_bd = 0;
}

Node::Node(const Node & nd){
	name = nd.name;
	pt = nd.pt;
	rid = nd.rid;
	value = nd.value;
	flag = nd.flag;
	//rep = nd.rep;
	rep = NULL;
	nbr_vec.clear();
	for(int i=0;i<4;i++){
		eqvr[i]=0.0;
		end[i]=this;
	}
	flag_bd = nd.flag_bd;
	internal_bd = nd.internal_bd;
}

Node & Node::operator = (const Node & nd){
	(*this) = Node(nd);
	return *this;
}
#if 0
// Ting
// get node's nbr from dir direction
// template<class T>
Node * Node::get_nbr_node(Node *node, DIRECTION dir) const{
	Node * node_nbr =node->nbr[dir]->ab[0];
	if(node_nbr->name == node->name){
		node_nbr = node->nbr[dir]->ab[1];
	}
	return node_nbr;	
}

// Ting
// get nodes' nbr from dir direction
// template<class T>
Node Node::get_nbr_node(Node node, DIRECTION dir) const{
	Node node_nbr = *node.nbr[dir]->ab[0];
	if(node_nbr.name == node.name){
		node_nbr = *node.nbr[dir]->ab[1];
	}
	return node_nbr;
}
#endif

ostream & operator << (ostream & os, const Node & node){
	//os<<node.name<<node.pt<<"="<<scientific<<node.value;
	//os<<" rep="<<node.rep->name;
	//Net * net = node.nbr[TOP];
	os<<setw(OUTPUT_WIDTH_STRING)<<node.name
	  //<<" -> "<<setw(OUTPUT_WIDTH_STRING)<<node.rep->name
	  //<<" top="<<(net==NULL?"NULL":net->name)
	  <<setw(OUTPUT_WIDTH_FLOAT)<<scientific<<node.value;
	return os;
}
