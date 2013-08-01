#ifndef __CIRCUIT_H__
#define __CIRCUIT_H__

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
//#include <tr1/unordered_map>
#include <map>
#include <list>
#include <cmath>
#include "cholmod.h"
#include "global.h"
#include "node.h"
#include "net.h"
#include "vec.h"
#include "triplet.h"
#include "block.h"
#include "transient.h"
#include "mpi_class.h"
#include <algorithm>
#include "sp_graph_table.h"
#include "sp_node.h"


using namespace std;
// using namespace std::tr1;

typedef vector<double> DoubleVector;
typedef vector<Net *> NetPtrVector;
typedef list<Net *> NetPtrList;
typedef vector<Node *> NodePtrVector;
typedef NetPtrVector NetList;

// functor of translating Node * to void *
/*namespace std{ 
	namespace tr1{
		template<> struct hash< Node * >{
			size_t operator()( const Node * x ) const {
				return hash< const char* >()( (char*) x );
			}
		};
	}
}*/

class Circuit{
public:
	Circuit(string name="");
	~Circuit();
	
	void pre_release_circuit();
	void clean_explore();
	void check_sys() const;
	friend class Block;
	// can be written as inline to speed up
	Node * get_node(string name);
	Net * get_net(string name);
	string get_name() const;

	// add a node into nodelist
	bool add_node(Node * nd);
	bool add_node_bd(int &count, Node * node);
	bool add_node_inter(int &count, Node * node);
	
	void add_node_inter_bd(Node *nd_0, Node *nd_1, double *geo_line, int &bid);

	// add a net into netset
	bool add_net(Net * net);

	bool has_node(string name) const;
	bool has_net(string name) const;

	// sort nodes according to predefined order
	void sort_nodes();
	void sort_bd_nodes(int &my_id);
	void sort_internal_nodes(int &my_id);

	// solve for node voltage
	void solve(int &my_id, int &num_procs, MPI_CLASS &mpi_class, Tran &tran);
	void block_solve(int my_id);
	void block_solve_Jacobi(int my_id);
	double find_diff(int my_id);
	
	static void set_parameters(double, double, size_t);
	static void get_parameters(double&, double&, size_t&);

	friend ostream & operator << (ostream & os, const Circuit & ckt);
	friend class Parser;

	// C style output
	void print();
	void print_matlab(Matrix A);
	void print_matrix(Matrix A);
	void print_rhs();
	void print_solution();

	// mpi related variables
	// member
	// ******* processor 0  variable ********
	// b_new_info is the global bd solution array
	double *bd_x_g;
	double *internal_x_g;	
	// ****** other processor *******

	// solution array for each processor
	double *bd_x;	
	double *internal_x;

	// stores boundary nodes value
	int *bd_base;
	int *internal_base;

	// stores 4 boundary base of each processor
	// into processor 0
	int *bd_base_gd;
	int *internal_base_gd;

	// stores the base for receiving bd_x_g
	int *bd_base_g;
	int *internal_base_g;
	
	// bd_size_g store the total bd_size of each block
	int *bd_size_g;
	int *internal_size_g;

	// stores 4 boundary size
	int *bd_dd_size;
	int *internal_dd_size;
	// stores toal block's 4 boundary size
	int *bd_dd_size_g;
	int *internal_dd_size_g;
	
	// total boundary node size
	int bd_size;
	int internal_size;
	// stores the whole grid size
	// into processor 0
	int total_size;
	int total_internal_size;

	int total_blocks;
	
	// block_size is either 0 or 1
	int block_size;

	// 8 boundary nodelist
	NodePtrVector bd_nodelist_sw;
	NodePtrVector bd_nodelist_s;
	NodePtrVector bd_nodelist_se;
	NodePtrVector bd_nodelist_w;
	NodePtrVector bd_nodelist_e;
	NodePtrVector bd_nodelist_nw;
	NodePtrVector bd_nodelist_n;
	NodePtrVector bd_nodelist_ne;

	// 8 internal nodelist	
	NodePtrVector internal_nodelist_sw;
	NodePtrVector internal_nodelist_s;
	NodePtrVector internal_nodelist_se;
	NodePtrVector internal_nodelist_w;
	NodePtrVector internal_nodelist_e;
	NodePtrVector internal_nodelist_nw;
	NodePtrVector internal_nodelist_n;
	NodePtrVector internal_nodelist_ne;

private:
	// member functions

	bool solve_IT(int &my_id, int&num_procs, MPI_CLASS &mpi_class, Tran &tran);
	void decomp_matrix(int &my_id, Matrix &A);
	
	// initialize things before solve_iteration
	void solve_init(int &my_id);
	void solve_init_LU();

	// updates nodes value in each iteration
	double solve_iteration(int &my_id, int &iter, int&num_procs, MPI_CLASS &mpi_class);

	void block_init(int &my_id, MPI_CLASS &mpi_class);
	void build_bd_netlist();
	void update_geometry(int my_id, MPI_CLASS & mpi_class);
	void assign_block_nodes(int my_id);
	void assign_block_nets(int my_id);

	// methods of stamping the matrix
	
	void check_matrix(Matrix &A);
	void current_tr(Net *net, double &time);

	void set_eq_induc(Tran &tran);
	void set_eq_capac(Tran &tran);
	void release_tr_nodes(Tran &tran);

	double solve_iteration_tr(int &my_id, int &iter,
		int&num_procs, MPI_CLASS &mpi_class);

	void solve_DC(int &num_blocks, int &my_id, MPI_CLASS &mpi_class);
	bool solve_tr_step(int &num_procs, int &my_id, MPI_CLASS &mpi_class);
	void solve_tr(Tran &tran, int &my_id);
	
	// ******** end of transient ****
	void stamp_boundary_matrix();
	void stamp_boundary_net(Net * net);
	void stamp_block_resistor(int &my_id, Net *net, Matrix &A);
	void stamp_block_current(int &my_id, Net * net, MPI_CLASS &mpi_class);
	void stamp_inductance_dc(Matrix & A, Net * net, int &my_id);
	void stamp_block_VDD(int &my_id, Net * net, Matrix &A);

	void boundary_init(int &my_id, int &num_procs);
	
	void internal_init(int &my_id, int &num_procs);

	void assign_bd_array(int &my_id);
	void assign_bd_array_dir(int &base, NodePtrVector &list, int &my_id);

	void reset_bd_array(int &my_id);
	void reset_replist(int &my_id);
	
	void assign_bd_base(int &my_id);
	void assign_bd_dd_size(int &my_id);
	void assign_internal_base(int &my_id);

	void assign_bd_internal_array(int &my_id);

	void assign_bd_internal_array_dir(int &base, 
		NodePtrVector &list, double *internal_x, int &my_id);

	void reorder_bd_x_g(MPI_CLASS &mpi_class);

	//  ******* method for PCG method  ********
	// solve circuit with preconditioned pcg method

	// ****** function for transient *******
	void release_ckt_nodes(Tran &tran);
	void print_ckt_nodes(Tran &tran);
	void save_ckt_nodes_to_tr(Tran &tran);
	void link_tr_nodes(Tran &tran);
	void link_ckt_nodes(Tran &tran, int &my_id);
	void save_tr_nodes(Tran &tran, double *x);
	void save_ckt_nodes(Tran &tran);//, double *x);

	void print_tr_nodes(Tran &tran);

	void set_len_per_block();
	void find_block_size (MPI_CLASS &mpi_class);

	double modify_voltage(int &my_id, Block &block_info, double* x_old);

	void set_type(CIRCUIT_TYPE type){circuit_type = type;};
	// ************* functions and members for thread **********
	// set s_col_FFS and FBS
	// ********* sparse vectors ******
	Path_Graph pg;
	void update_node_set_bx();           	     void push_bd_nodes(Path_Graph &pg, int &my_id);

	void push_bd_net_nodes();
	void push_bd_nodes_one_set(Path_Graph &pg, int&my_id, NodePtrVector internal_set);

	void get_samples();
	bool check_diverge() const;
	
	vector<Node_TR_PRINT> ckt_nodes;
	// ************** member variables *******************
	NodePtrVector nodelist;		// a set of nodes
	NodePtrVector replist;		// a set of representative nodes
	NodePtrVector mergelist;	// nodes for merging
	vector<long> x_list;
	vector<long> y_list;
	map<long, int> x_list_bd_map;
	map<long, int> y_list_bd_map;
	map<long, int> x_list_nd_map;
	map<long, int> y_list_nd_map;

	// bd_netlist does not belong to net_set
	NetPtrVector bd_netlist;	
	NetList net_set[NUM_NET_TYPE];// should be the same as size of NET_TYPE
	// defines the net direction in layers
	// mapping from name to Node object pointer
	//unordered_
	map<string, Node*> map_node;

	// mapping from Net pointer to their index in netlist
	//unordered_
	map<Net*, size_t> net_id;

	// circuit name
	string name;

	// blocks
	vector<Block*> block_vec;
	float x_min, y_min, x_max, y_max;
	long lx, ly, ux, uy;

	// control variables
	static double EPSILON;
	// static double OMEGA;
	static double OVERLAP_RATIO;
	static size_t MAX_BLOCK_NODES;
	// static int MODE; // 0 = IT, 1 = LU
	static int NUM_BLOCKS_X;
	static int NUM_BLOCKS_Y;
	static int DEBUG;
	int num_blocks;

	CIRCUIT_TYPE circuit_type;

	NodePtrVector sample;

	double VDD;
};

// adds a node into nodelist
inline bool Circuit::add_node(Node * node){
	node->flag_bd = 0;
	nodelist.push_back(node);
	map_node[node->name] = node;
	return true;
}

// adds a node into nodelist
inline bool Circuit::add_node_bd(int &count, Node * node){
	if(count ==1){
		bd_nodelist_s.push_back(node);
	}
	else if(count==2){
		bd_nodelist_n.push_back(node);
	}
	else if(count==3){
		bd_nodelist_w.push_back(node);
	}
	else if(count==4){
		bd_nodelist_e.push_back(node);
	}
	// node is boundary node
	node->flag_bd = 1;
	map_node[node->name] = node;
	return true;
}

// adds a node into nodelist
inline bool Circuit::add_node_inter(int &count, Node * node){
	if(count ==1){
		internal_nodelist_s.push_back(node);
	}
	else if(count==2){
		internal_nodelist_n.push_back(node);
	}
	else if(count==3){
		internal_nodelist_w.push_back(node);
	}
	else if(count==4){
		internal_nodelist_e.push_back(node);
	}
	// node is internal_bd boundary node
	node->internal_bd = 1;
	return true;
}

// adds a net into netset
inline bool Circuit::add_net(Net * net){
	// has at least one bd node, then belongs to bd net
	if(net->ab[0]->flag_bd ==1 || net->ab[1]->flag_bd ==1)
		net->flag_bd = 1;
	else net->flag_bd = 0;

	if( net->type == RESISTOR )
		net_id[net] = net_set[net->type].size();
	net_set[net->type].push_back(net);
	return true;
}

// fina a node by name
inline bool Circuit::has_node(string name) const{
	if( map_node.find(name) != map_node.end() ) return true;
	return false;
}

// get a node by name
inline Node * Circuit::get_node(string name){
	map<string, Node*>::const_iterator it = map_node.find(name);
	if( it != map_node.end() ) return it->second;
	else return NULL;
}
// bool compare_node_ptr(const Node *a, const Node *b);
ostream & operator << (ostream & os, const NodePtrVector & nodelist);
ostream & operator << (ostream & os, const NetList & nets);
//ostream & operator << (ostream & os, const vector<Block > & block_info);

#endif
