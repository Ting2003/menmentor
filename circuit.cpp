#include <fstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <utility>
#include <cassert>
#include <vector>
#include "cholmod.h"
//#include "umfpack.h"
#include "circuit.h"
#include "util.h"
#include "algebra.h"
#include "node.h"

#include "mpi.h"
using namespace std;

double Circuit::EPSILON = 1e-5;
size_t Circuit::MAX_BLOCK_NODES =100000;//5500;
double Circuit::OMEGA = 1.2;
double Circuit::OVERLAP_RATIO = 0;
int    Circuit::MODE = 0;
const int MAX_ITERATION = 100000000;
const int SAMPLE_INTERVAL = 5;
const size_t SAMPLE_NUM_NODE = 10;
const double MERGE_RATIO = 0.3;
int Circuit::NUM_BLOCKS_X = 1;
int Circuit::NUM_BLOCKS_Y = 1;
int Circuit::DEBUG=1;

//////////////////////////////////////////////////////////////////////////
// Constructor and utility functions goes here

// constructor of Circuit class, name is optional
Circuit::Circuit(string _name):name(_name),
	x_min(INFTY),y_min(INFTY),x_max(0),y_max(0),
	circuit_type(UNKNOWN), VDD(0.0){
	// add ground node
	Node * gnd = new Node(string("0"), Point(-1,-1,-1));
	gnd->rep = gnd;
	this->add_node(gnd);

	// mpi relate	
	bd_x_g=NULL;
	internal_x_g = NULL;
	
	bd_x = NULL;
	internal_x = NULL;
	bd_base = NULL;
	internal_base = NULL;
	bd_base_gd = NULL;
	internal_base_gd = NULL;
	bd_base_g = NULL;
	internal_base_g = NULL;
	bd_size_g = NULL;
	internal_size_g = NULL;
	bd_dd_size = NULL;
	internal_dd_size = NULL;
	bd_dd_size_g = NULL;
	internal_dd_size_g = NULL;

	bd_size = 0;
	internal_size = 0;
	block_size=0;
}

// relase resources in pre partition 
void Circuit::pre_release_circuit(){
	// delete node
	for(size_t i=0;i<nodelist.size();i++){
		if(!nodelist[i]->is_ground()){
			// clog<<"node: "<<*nodelist[i]<<endl;
		 	free(nodelist[i]);
		// clog<<"finish deleting i: "<<i<<" th node. "<<endl;
		}
	}
	// clog<<"after release nodelist. "<<endl;
	// nodelist.clear();
	// delete nets
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		for(it=ns.begin();it!=ns.end();++it)
			free(*it);
	}
	// clog<<"after release netlist. "<<endl;
}

// Trick: do not release memory to increase runtime
Circuit::~Circuit(){
	// delete node
	for(size_t i=0;i<nodelist.size();i++){
		nodelist[i]->nbr_vec.clear();
		 delete nodelist[i];
	}
	// delete node
	for(size_t i=0;i<replist.size();i++) delete replist[i];
	// delete node
	for(size_t i=0;i<mergelist.size();i++) delete mergelist[i];
	nodelist.clear();
	replist.clear();
	mergelist.clear();
	map_node.clear();
	block_vec.clear();
	for(size_t i=0;i<bd_netlist.size();i++)
    		delete bd_netlist[i];
    	bd_netlist.clear();

	// delete bd nodes
	bd_nodelist_sw.clear();
	bd_nodelist_s.clear();
	bd_nodelist_se.clear();
	bd_nodelist_w.clear();
	bd_nodelist_e.clear();
	bd_nodelist_nw.clear();
	bd_nodelist_n.clear();
	bd_nodelist_ne.clear();

	internal_nodelist_sw.clear();
	internal_nodelist_s.clear();
	internal_nodelist_se.clear();
	internal_nodelist_w.clear();
	internal_nodelist_e.clear();
	internal_nodelist_nw.clear();
	internal_nodelist_n.clear();
	internal_nodelist_ne.clear();
	
	// delete nets
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		for(it=ns.begin();it!=ns.end();++it)
			delete *it;
	}

	// free Li, Lx and so on
	/*delete [] Lx;
	delete [] Lp;
	delete [] Lnz;
	delete [] Li;*/
	// mpi related variables
	delete [] bd_x_g;
	delete [] internal_x_g;
	delete [] bd_x;
	delete [] internal_x;
	delete [] bd_base;
	delete [] internal_base;
	delete [] bd_base_gd;
	delete [] internal_base_gd;
	delete [] bd_base_g;
	delete [] internal_base_g;
	delete [] bd_size_g;
	delete [] internal_size_g;
	delete [] bd_dd_size;
	delete [] internal_dd_size;
	delete [] bd_dd_size_g;
	delete [] internal_dd_size_g;
}

void Circuit::check_sys() const{
	clog<<"**** CHECKING SYSTEM ENVIRONMENT ****"<<endl;
	clog<<"* int size     = "<< sizeof(int)<<endl;
	clog<<"* long size    = "<< sizeof(long)<<endl;
	clog<<"* size_t size  = "<< sizeof(size_t)<<endl;
	//clog<<"* UF_long size = "<< sizeof(UF_long)<<endl;
	clog<<"* Max nodelist = "<<(size_t)nodelist.max_size()<<endl;
	clog<<"****            END              ****"<<endl<<endl;
}
// sort the nodes according to their coordinate 
// sort nodelist
void Circuit::sort_nodes(){
	sort(nodelist.begin(), nodelist.end(), compare_node_ptr);
	// update node id mapping, 
	// NOTE: ground node will be the last
}

// sort 4 boudnary nodelist
void Circuit::sort_bd_nodes(int &my_id){
	sort(bd_nodelist_sw.begin(), bd_nodelist_sw.end(), 
			compare_node_ptr);

	sort(bd_nodelist_s.begin(), bd_nodelist_s.end(), 
			compare_node_ptr);
	
	sort(bd_nodelist_se.begin(), bd_nodelist_se.end(), 
			compare_node_ptr);

	sort(bd_nodelist_w.begin(), bd_nodelist_w.end(), 
			compare_node_ptr);

	sort(bd_nodelist_e.begin(), bd_nodelist_e.end(), 
			compare_node_ptr);

	sort(bd_nodelist_nw.begin(), bd_nodelist_nw.end(), 
			compare_node_ptr);
	
	sort(bd_nodelist_n.begin(), bd_nodelist_n.end(), 
			compare_node_ptr);

	sort(bd_nodelist_ne.begin(), bd_nodelist_ne.end(), 
			compare_node_ptr);
}

// sort 4 boudnary nodelist
void Circuit::sort_internal_nodes(int &my_id){
	sort(internal_nodelist_sw.begin(), 
		internal_nodelist_sw.end(), compare_node_ptr);

	sort(internal_nodelist_s.begin(), 
		internal_nodelist_s.end(), compare_node_ptr);
	
	sort(internal_nodelist_se.begin(), 
		internal_nodelist_se.end(), compare_node_ptr);

	sort(internal_nodelist_w.begin(), 
		internal_nodelist_w.end(), compare_node_ptr);
	
	sort(internal_nodelist_e.begin(), 
		internal_nodelist_e.end(), compare_node_ptr);

	sort(internal_nodelist_nw.begin(), 
		internal_nodelist_nw.end(), compare_node_ptr);
	
	sort(internal_nodelist_n.begin(), 
		internal_nodelist_n.end(), compare_node_ptr);

	sort(internal_nodelist_ne.begin(), 
		internal_nodelist_ne.end(), compare_node_ptr);
}

string Circuit::get_name() const{return this->name;}

ostream & operator << (ostream & os, const NodePtrVector & nodelist){
	for(size_t i=0;i<nodelist.size();i++)
		os<<*nodelist[i]<<endl;
	return os;
}

ostream & operator << (ostream & os, const NetList & nets){
	NetList::const_iterator it;
	for(it=nets.begin();it!=nets.end();++it)
		if( (*it) != NULL ) os<<**it<<endl;
	return os;
}

ostream & operator << (ostream & os, const Circuit & ckt){
	os<<"Circuit ["<<ckt.name<<"] info:"<<endl;

	os<<"==== Nodes ===="<<endl;
	os<<ckt.nodelist;

	os<<"==== Reps  ===="<<endl;
	os<<ckt.replist;

	os<<"==== Nets  ===="<<endl;
	os<<ckt.net_set[RESISTOR];

	return os;
}

void Circuit::print(){
	// uncomment this if want to output to a file
	//freopen("output.txt","w",stdout);

	// don't output ground node
	for(size_t i=0;i<nodelist.size()-1;i++){
		printf("%s  %.5e\n", nodelist[i]->name.c_str(), 
				nodelist[i]->value);
	}
}

void Circuit::print_matlab(Matrix A){
	// uncomment this if want to output to a file
	//freopen("output.txt","w",stdout);

	// don't output ground node
	for(size_t i=0;i<A.size();i++){
		printf("%ld %ld %.5e\n", A.Ti[i]+1, A.Tj[i]+1, A.Tx[i]);
	}
}

void Circuit::print_matrix(Matrix A){
	bool DEBUG_flag = true;
	// uncomment this if want to output to a file
	stringstream ss;
	ss<<"OUTPUT_TRI/b2r_"<<name<<"_A.txt";
	clog<<"print file name: "<<ss.str()<<endl;
	// # build matrix merge map
	pair<long, long> index_pair;
	map<pair<long, long>, double> matrix_map;
	A.merge_matrix(matrix_map);
	
	size_t count = 0;
	map<pair<long, long>, double>::iterator it;
	for(it = matrix_map.begin(); it != matrix_map.end(); it++){
		count++;
	}

	//A.merge();
	FILE *f;
	f = fopen(ss.str().c_str(),"w");

# if 1
	// print title line
	fprintf(f, "\%\%MatrixMarket matrix coordinate real general\n");
	size_t num_nnz = (count-A.get_row())*2 + A.get_row();
	// print 1st line
	fprintf(f, "%d %d %d\n", A.get_row(), A.get_row(), num_nnz);
# endif
#if 0
	fprintf(f, "\%%MatrixMarket matrix coordinate real symmetric\n");
	fprintf(f, "%d %d %d\n", A.get_row(), A.get_row(), count);
# endif
	for(it = matrix_map.begin(); it != matrix_map.end(); it++){
		fprintf(f, "%d %d %.10e\n", it->first.first+1, it->first.second+1, 
			it->second);
# if 1
		/*if(it->first.first == it->first.second)
			fprintf(f, "%d %d %.10e\n", it->first.second+1, 
				it->first.first+1, it->second+1e-5);
		else
			fprintf(f, "%d %d %.10e\n", it->first.first+1, it->first.second+1, 
				it->second);
		*/

		if(it->first.first != it->first.second)
			fprintf(f, "%d %d %.10e\n", it->first.second+1, 
				it->first.first+1, it->second);
# endif

	}
#if 0	// don't output ground node
	for(size_t i=0;i<A.size();i++){
		fprintf(f, "%d %d %.10e\n", A.Ti[i]+1, A.Tj[i]+1, A.Tx[i]);
# if DEBUG_flag
		if(A.Ti[i] != A.Tj[i])
			fprintf(f, "%d %d %.10e\n", A.Tj[i]+1, A.Ti[i]+1, A.Tx[i]);
# endif
	}
#endif
	fclose(f);
	matrix_map.clear();
}

void Circuit::print_rhs(){
	// uncomment this if want to output to a file
	stringstream ss;
	ss<<"OUTPUT_TRI/b2r_"<<name<<"_b.txt";
	clog<<"print file name: "<<ss.str()<<endl;
	FILE *f;
	f = fopen(ss.str().c_str(),"w");
	// print title line
	fprintf(f, "\%\%MatrixMarket matrix array real general\n");
	/*size_t count = 0;
	for(size_t i=0;i<block_info.count;i++)
		if(block_info.bnewp[i] !=0)
			count++;
	clog<<"count: "<<count<<endl;*/
	// print 1st line
	fprintf(f, "%d %d\n", block_vec[0]->count, 1);

	// don't output ground node
	for(size_t i=0;i<block_vec[0]->count;i++){
		fprintf(f, "%.10e\n", block_vec[0]->bnewp[i]);
	}
	fclose(f);
}

void Circuit::print_solution(){
	// uncomment this if want to output to a file
	stringstream ss;
	ss<<"OUTPUT_TRI/b2r_"<<name<<"_x.txt";
	clog<<"print file name: "<<ss.str()<<endl;
	FILE *f;
	f = fopen(ss.str().c_str(),"w");
	// print title line
	fprintf(f, "\%\%MatrixMarket matrix array real general\n");
	// print 1st line
	fprintf(f, "%d %d %d\n", block_vec[0]->count, 1, block_vec[0]->count);

	// don't output ground node
	for(size_t i=0;i<block_vec[0]->count;i++){
		fprintf(f, "%d %.10e\n", i+1, block_vec[0]->xp[i]);
	}
	fclose(f);
}

///////////////////////////////////////////////////////////////////////////////
// Computation Functions

// initialization before solving the circuit
// 1. sort the nodes
// 2. set node representatives
// 3. find node in which block, update count
// 4. get representative lists
void Circuit::solve_init(int &my_id){
	sort_nodes();
	// if(my_id==0)
		// clog<<"total num_nodes: "<<nodelist.size()<<endl;
	sort_bd_nodes(my_id);
	sort_internal_nodes(my_id);

	// assign VDD value of circuit
	int type = VOLTAGE;
	NetList & ns = net_set[type];
	for(size_t i=0;i<ns.size();i++){
		VDD = ns[i]->value;
		break;		
	}
	/*if(my_id==3){
		clog<<"sw: "<<internal_nodelist_sw<<endl;
		clog<<"s: "<<internal_nodelist_s<<endl;
		clog<<"se: "<<internal_nodelist_se<<endl;
		clog<<"w: "<<internal_nodelist_w<<endl;
		clog<<"e: "<<internal_nodelist_e<<endl;
		clog<<"nw: "<<internal_nodelist_nw<<endl;
		clog<<"n: "<<internal_nodelist_n<<endl;
		clog<<"ne: "<<internal_nodelist_ne<<endl;
		clog<<endl;
	}

	if(my_id==3){
		clog<<"sw: "<<bd_nodelist_sw<<endl;
		clog<<"s: "<<bd_nodelist_s<<endl;
		clog<<"se: "<<bd_nodelist_se<<endl;
		clog<<"w: "<<bd_nodelist_w<<endl;
		clog<<"e: "<<bd_nodelist_e<<endl;
		clog<<"nw: "<<bd_nodelist_nw<<endl;
		clog<<"n: "<<bd_nodelist_n<<endl;
		clog<<"ne: "<<bd_nodelist_ne<<endl;
	}*/

	size_t size = nodelist.size() - 1;
	Node * p = NULL;
	size_t nr = 0;
	size_t i=0;
	for(i=0, nr=0;i<size;i++){
		p=nodelist[i];
		/*Net * net = p->nbr[TOP];	
		// test short circuit
		if( p->isS() !=Y && // Y must be representative 
		    net != NULL &&
		    fzero(net->value) ){
			// TODO: ensure ab[1] is not p itself
			assert( net->ab[1] != p );
			p->rep = net->ab[1]->rep;
		} // else the representative is itself
		*/
		// push the representatives into list
		if( p->rep == p ) {
			replist.push_back(p);
			//rep_id[p] = nr; // save the id
			p->rid = nr;
			++nr;
		}	
	}// end of for i

	/*if(my_id==0){
		cout<<endl;
		for(int i=0;i<nodelist.size()-1;i++)
			cout <<*nodelist[i]<<" "<<nodelist[i]->rid<<endl;
		//clog<<endl;
	}*/
	

	// if(nr >=0)
		// block_info.count = nr;

	size_t n_nodes = nodelist.size();
	size_t n_reps  = replist.size();
	
	/*clog<<"replist    "<<n_reps <<endl;
	clog<<"nodelist   "<<n_nodes<<endl;
	clog<<"ratio =    "<<ratio  <<endl;*/

	net_id.clear();
	//clog<<my_id<<" "<<block_info.count<<endl;
}

// build up block info
// 1. Find block divide point
// 2. Set each node in replist into block
// 3. Compute block size
// 4. Insert boundary netlist into map
void Circuit::block_init(int &my_id, MPI_CLASS &mpi_class){
	// assign nodes into blocks and sort
	assign_block_nodes(my_id);
	assign_block_nets(my_id);

	for(size_t i=0;i<block_vec.size();i++){
		block_vec[i]->allocate_resource();
		// copy_node_voltages_block();
		block_vec[i]->stamp_matrix(my_id, mpi_class);
	}
	
}
#if 0
// stamp the nets by sets, block version
// *NOTE* at the same time insert the net into boundary netlist
void Circuit::stamp_block_matrix(int &my_id, Matrix &A, MPI_CLASS &mpi_class){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		switch(type){
		case RESISTOR:
			if(my_id==0)
				clog<<"resis net. "<<ns.size()<<endl;
			for(it=ns.begin();it!=ns.end();++it){
				Net * net = *it;
				if( net == NULL ) continue;
				assert( fzero(net->value) == false );
				stamp_block_resistor(my_id, *it, A);
			}
			break;
		case CURRENT:
			for(it=ns.begin();it!=ns.end();++it){
				stamp_block_current(my_id, (*it), mpi_class);
			}
			break;
		case VOLTAGE:
			if(my_id==0)
				clog<<"VDD net: "<<ns.size()<<endl;
			for(it=ns.begin();it!=ns.end();++it){
				if( fzero((*it)->value)  && 
				    !(*it)->ab[0]->is_ground() &&
				    !(*it)->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_block_VDD(my_id,(*it), A);
			}
			break;
		case CAPACITANCE:
			break;
		case INDUCTANCE:
			if(my_id==0)
				clog<<"induc net: "<<ns.size()<<endl;
			for(size_t i=0;i<ns.size();i++)
				stamp_inductance_dc(A, ns[i], my_id);
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
	//if(my_id==0)
		//clog<<A<<endl;
	/*if(my_id==0){
		for(int i=0;i<10;i++)
			clog<<"b origin: "<<i<<" "<<block_info.bp[i]<<endl;
	}*/
	 make_A_symmetric(block_info.bp, my_id);
	
	A.set_row(block_info.count);
	//if(my_id==0)
		//check_matrix(A);
		//cout<<"before CK_decomp. "<<endl;
	if(block_info.count >0){
		block_info.CK_decomp(A, cm);
		//if(cm->status ==1)
			//clog<<" non SPD: "<<my_id<<" "<<cm->status<<endl;
	}
	//if(my_id==0)
		//clog<<"after CK_decomp. "<<endl;
}

// stamp the nets by sets, block version
// *NOTE* at the same time insert the net into boundary netlist
void Circuit::stamp_block_matrix_tr(int &my_id, Matrix &A, MPI_CLASS &mpi_class, Tran &tran){	
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		switch(type){
		case RESISTOR:
			for(it=ns.begin();it!=ns.end();++it){
				Net * net = *it;
				if( net == NULL ) continue;
				assert( fzero(net->value) == false );
				stamp_block_resistor_tr(my_id, *it, A);
			}
			break;
		case CURRENT:
			break;
		case VOLTAGE:
			for(it=ns.begin();it!=ns.end();++it){
				if( fzero((*it)->value)  && 
				    !(*it)->ab[0]->is_ground() &&
				    !(*it)->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_block_VDD_tr(my_id,(*it), A);
			}
			break;
		case CAPACITANCE:
			for(size_t i=0;i<ns.size();i++)
				stamp_capacitance_tr(A, ns[i], tran, my_id);
			break;
		case INDUCTANCE:
			for(size_t i=0;i<ns.size();i++)
				stamp_inductance_tr(A, ns[i], tran, my_id);
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
	//if(my_id==0)
		//clog<<A<<endl;
	/*if(my_id==0){
		for(int i=0;i<10;i++)
			clog<<"b origin: "<<i<<" "<<block_info.bp[i]<<endl;
	}*/	
}
#endif

void Circuit::solve(int &my_id, int&num_procs, MPI_CLASS &mpi_class, Tran &tran){
	// each block is solved by IT
	solve_IT(my_id, num_procs, mpi_class, tran);
	//clog<<my_id<<" finish solve: "<<endl;
}

// solve Circuit
bool Circuit::solve_IT(int &my_id, int&num_procs, MPI_CLASS &mpi_class, Tran &tran){
	double time=0;
	double t1, t2;
	
	/*cm = &c;
	cholmod_start(cm);
	cm->print = 5;*/

	total_blocks = mpi_class.X_BLOCKS *mpi_class.Y_BLOCKS;

	// did not find any `X' node
	if( circuit_type == UNKNOWN )
		circuit_type = C4;
	
	if(mpi_class.block_size>0){
		solve_init(my_id);
	}
	// update bd info of circuit and block vec
	update_geometry(my_id, mpi_class);

	// build boundary netlist of circuit
	build_bd_netlist();

	block_init(my_id, mpi_class);
	boundary_init(my_id, num_procs);
	internal_init(my_id, num_procs);
	
	bool successful = false;
	/*if(my_id==0){
		for(size_t i=0;i<block_vec.size();i++){
		clog<<"block: "<<i<<endl;
		clog<<"DC matrix: "<<block_vec[i]->A<<endl;
		}
	}*/
	//get_voltages_from_block_LU_sol();
	solve_DC(num_procs, my_id, mpi_class);
 #if 0	
	print_matrix(block_vec[0]->A);
	print_rhs();
	print_solution();
 #endif
	// cout<<nodelist;
	// if(my_id==0)
		// cout<<nodelist<<endl;
		/*for(size_t i=0;i<block_vec.size();i++){
			for(size_t j=0;j<block_vec[i]->count;j++)
				cout<<"b: "<<*block_vec[i]->replist[j]<<endl;
		}*/
	// return true;
	// then sync
	MPI_Barrier(MPI_COMM_WORLD);
	// return 0;
//#if 0
	for(size_t i=0;i<block_vec.size();i++){
		block_vec[i]->reset_array(block_vec[i]->bp);
		block_vec[i]->reset_array(block_vec[i]->bnewp);
	}
	/***** solve tran *********/
	// link transient nodes
	link_ckt_nodes(tran, my_id);
	for(size_t i=0;i<block_vec.size();i++){
		// block_vec[i]->flag_ck = 0;
		block_vec[i]->stamp_matrix_tr(my_id, mpi_class, tran);	
		block_vec[i]->make_A_symmetric_tr(my_id, tran);	   
		block_vec[i]->stamp_current_tr(my_id, time);
		// if(my_id==0)
			// clog<<block_vec[i]->A<<endl;

   		block_vec[i]->CK_decomp();
//#if 0   
		block_vec[i]->clear_A();
		block_vec[i]->build_id_map();
//#endif
		// bnewp = bp
		block_vec[i]->copy_vec(block_vec[i]->bnewp,
				block_vec[i]->bp);
	} 
  
   // set Geq for induc and capac
   set_eq_induc(tran);
   set_eq_capac(tran);
   
   // if(my_id==0)
	//   clog<<"before modify_rhs_tr_0. "<<endl;
   // already push back cap and induc into set_x and b
   for(size_t i=0;i<block_vec.size();i++){
   	block_vec[i]->modify_rhs_tr_0(block_vec[i]->bnewp, block_vec[i]->xp, my_id);

	// push node_set_b and part of set_x
  	block_vec[i]->push_nd_set_bx(tran);
	// push the boundary nodes of blocks
	block_vec[i]->push_bd_nets();
   }
   
   // push boundary nodes circuit's blocks
   push_bd_nodes(pg, my_id);
   // push_bd_net_nodes();

  int sum_n_ffs = 0;
  int sum_n_fbs = 0;
  int sum_n = 0;
  int block_flag = 0; 
   // get path_b, path_x, len_path_b, len_path_x
  for(size_t i=0;i<block_vec.size();i++){   
  	block_vec[i]->build_path_graph();
	sum_n_ffs += block_vec[i]->len_path_b;
	sum_n_fbs += block_vec[i]->len_path_x;
	sum_n += block_vec[i]->count;
	block_flag = 1;
	block_vec[i]->find_super();
	// block_vec[i]->test_path_super();
	// block_vec[i]->solve_eq_sp(block_vec[i]->xp, block_vec[i]->bnewp);
  }
  int root_sum_n_ffs = 0;
  int root_sum_n_fbs = 0;
  int root_sum_n = 0;  
  int root_block_flag = 0;

  MPI_Reduce(&sum_n_ffs, &root_sum_n_ffs, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum_n_fbs, &root_sum_n_fbs, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum_n, &root_sum_n, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&block_flag, &root_block_flag, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(my_id==0){
	double avg_n_ffs = 1.0 * root_sum_n_ffs / root_block_flag;
	double avg_n_fbs = 1.0 * root_sum_n_fbs / root_block_flag;
	double avg_n = 1.0 * root_sum_n / root_block_flag;
	clog<<"root_block_flag, avg_n_ffs, fbs, n: "<<root_block_flag<<" "<<avg_n_ffs<<" "<<avg_n_fbs<<" "<<avg_n<<endl;
  }
      /*if(my_id==0)
	   cout<<nodelist<<endl;
   for(size_t i=0;i<block_vec.size();i++){
	   for(size_t j=0;j<block_vec[i]->count;j++)
		   cout<<"b: "<<*block_vec[i]->replist[j]<<endl;
   }*/
  
   solve_tr_step(num_procs, my_id, mpi_class);
#if 0
   if(my_id==0){
	cout<<endl<<" first time step sol: "<<endl;
	for(size_t i=0;i<replist.size();i++){
		cout<<"i, nd: "<<i<<" "<<replist[i]->name<<" "<<replist[i]->pt<<" "<<replist[i]->value<<endl;
	}
	// cout<<nodelist<<endl;
   }
#endif
   //save_tr_nodes(tran, xp);
   // for(size_t i=0;i<block_vec.size();i++)
	save_ckt_nodes(tran);//, block_vec[i]->xp);

   time += tran.step_t;
   MPI_Barrier(MPI_COMM_WORLD);

   // return 0;
   int iter = 0;
   clock_t tr_ts, tr_te;
   //if(my_id==0)
	  // clog<<"after first time step. "<<endl;
   //for(; time <= tran.tot_t; time += tran.step_t){
   while(time <= tran.tot_t){// && iter < 1){
	if(iter ==0)
		tr_ts = clock();
	// bnewp[i] = bp[i];
	for(size_t i=0;i<block_vec.size();i++){
		block_vec[i]->copy_vec(block_vec[i]->bnewp, block_vec[i]->bp);

      		block_vec[i]->stamp_current_tr_1(time);
      		// get the new bnewp
      		block_vec[i]->modify_rhs_tr(block_vec[i]->bnewp, block_vec[i]->xp);
	}

      // if(my_id==0)
	  //    clog<<" ===== step: ===== "<<my_id<<" "<<time<<endl;
      // need to add bcast function for processors
      // need to be modified into block version
      //solve_eq_sp(block_info.xp, block_info.bnewp);
      // solution stored in block_info.xp
      //if(replist.size()>0)
      solve_tr_step(num_procs, my_id, mpi_class);

      // if(my_id==0)
	  // cout<<endl<<" "<<nodelist<<endl;
      //save_tr_nodes(tran, xp);

      // for(size_t i=0;i<block_vec.size();i++)
      	save_ckt_nodes(tran);//, block_vec[i]->xp);
      time += tran.step_t;
      // sync in the end of each time step
      MPI_Barrier(MPI_COMM_WORLD);
      /*if(iter ==0){
		tr_te = clock();
		clog<<"time to solve one step is: "<<1.0*(tr_te - tr_ts)/CLOCKS_PER_SEC<<endl;
	}*/
      //clog<<"after time-t step barrier. "<<my_id<<" "<<time<<endl;
      iter ++;
   }

   save_ckt_nodes_to_tr(tran);
   release_ckt_nodes(tran);
   for(size_t i=0;i<block_vec.size();i++){
	block_vec[i]->delete_paths();
   }
// #endif
#if 0
	/////////// release resources
	if(block_info.count > 0)
		block_info.free_block_cholmod(cm);
	//if(my_id==0) clog<<"free block info. "<<endl;
	cholmod_finish(cm);
	//clog<<"cholmod finish. "<<my_id<<endl;
#endif
        // MPI_Barrier(MPI_COMM_WORLD);
	return successful;
}
// solve blocks with mpi: multi-core
// One iteration during solving the circuit, for any block B:
// 1. update the righthand-side of the matrix of B
// 2. solve the matrix
// 3. update node voltages
// 4. track the maximum error of solution
double Circuit::solve_iteration_tr(int &my_id, int &iter,
		int&num_procs, MPI_CLASS &mpi_class){	
	double diff = .0;
	double diff_root=0;

	// 0 rank cpu will scatter all bd valuesfrom bd_x_g to bd_x
	if(iter >0){
		MPI_Scatterv(bd_x_g, bd_size_g, 
			bd_base_g, MPI_DOUBLE, bd_x, bd_size, 
			MPI_DOUBLE, 0, MPI_COMM_WORLD);

		assign_bd_array(my_id);
	}
	if(iter == 0){

		// reset_bd_array(my_id);
		// reset_replist(my_id);
		for(size_t i=0;i<block_vec.size();i++){
			block_vec[i]->copy_array(block_vec[i]->bnewp_temp, block_vec[i]->bnewp);
		}
	}
	
	for(size_t i=0;i<block_vec.size();i++){

		block_vec[i]->update_rhs(
			block_vec[i]->bnewp_temp, 
			block_vec[i]->bnewp, my_id);
	}

	if(iter==0){
		for(size_t i=0;i<block_vec.size();i++){
			block_vec[i]->reset_array(block_vec[i]->x_old);
		}
	}
	else{
		// x_old stores old solution
		for(size_t i=0;i<block_vec.size();i++){
			block_vec[i]->copy_array(block_vec[i]->x_old, block_vec[i]->xp);
		}
	}

	// new rhs store in bnewp
	for(size_t i=0;i<block_vec.size();i++){
		/*if(my_id==0)
			clog<<endl<<"before update rhs. "<<i<<endl;*/
		//block_info.solve_CK(cm);
		if(block_vec[i]->count >0){
			block_vec[i]->solve_CK_tr();
			/*if(my_id==0)
				for(size_t j=0;j<block_vec[i]->count;j++)
					clog<<"j, xp: "<<j<<" "<<block_vec[i]->xp[j]<<endl;*/
		}
		//solve_eq_sp(block_info.xp, block_info.bnewp);
			//block_info.xp = static_cast<double *>(block_info.x_ck->x);
		
		double local_diff = 
			block_vec[i]->modify_voltage(my_id);
		if(local_diff > diff)
			diff = local_diff;

	}
	/*if(my_id==1)
		for(int i=0;i<replist.size();i++)
			clog<<"xp: "<<i<<" "<<*replist[i]<<" "<<endl;*/

	// clog<<"diff: "<<my_id<<" "<<diff<<endl;
	assign_bd_internal_array(my_id);

	// 0 rank cpu will gather all the solution 
	// from bd_x to bd_x_g
	MPI_Gatherv(internal_x, internal_size, MPI_DOUBLE, 
		internal_x_g, internal_size_g, 
		internal_base_g, MPI_DOUBLE, 0, 
		MPI_COMM_WORLD);
	
	// reorder boundary array according to nbrs
	if(my_id==0){
		reorder_bd_x_g(mpi_class);
	}

	MPI_Reduce(&diff, &diff_root, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&diff_root, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return diff_root;
}
// solve blocks with mpi: multi-core
// One iteration during solving the circuit, for any block B:
// 1. update the righthand-side of the matrix of B
// 2. solve the matrix
// 3. update node voltages
// 4. track the maximum error of solution
double Circuit::solve_iteration(int &my_id, int &iter,
		int&num_procs, MPI_CLASS &mpi_class){	
	double diff = .0;
	double diff_root=0;

	// 0 rank cpu will scatter all bd valuesfrom bd_x_g to bd_x
	MPI_Scatterv(bd_x_g, bd_size_g, 
			bd_base_g, MPI_DOUBLE, bd_x, bd_size, 
			MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
	assign_bd_array(my_id);	

	// save old values first
	for(size_t i=0;i<block_vec.size();i++){
		block_vec[i]->copy_array(block_vec[i]->x_old, block_vec[i]->xp);
	}
	// solve_blocks until converge
	block_solve(my_id);
	//block_solve_Jacobi(my_id);

	diff = find_diff(my_id);
	// if(my_id==0)
		//cout<<"find diff: "<<diff<<endl;
	MPI_Barrier(MPI_COMM_WORLD);
	// new rhs store in bnewp and solve
	/*for(size_t i=0; i < block_vec.size();i++){
		

		block_vec[i]->solve_CK_DC(my_id);
		double local_diff = 
			block_vec[i]->modify_voltage(my_id);
		//if(my_id==0)
			//clog<<"local_diff: "<<local_diff<<endl;
		if(local_diff > diff)
			diff = local_diff;
	}*/
	// if(my_id==0)
		//clog<<"diff: "<<diff<<endl;
	assign_bd_internal_array(my_id);
	// 0 rank cpu will gather all the solution from bd_x
	// to bd_x_g
	MPI_Gatherv(internal_x, internal_size, MPI_DOUBLE, 
		internal_x_g, internal_size_g, 
		internal_base_g, MPI_DOUBLE, 0, 
		MPI_COMM_WORLD);
	
	// reorder boundary array according to nbrs
	if(my_id==0){
		reorder_bd_x_g(mpi_class);
	}
	
	MPI_Reduce(&diff, &diff_root, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Bcast(&diff_root, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);	
	
	// if(my_id==0) clog<<"iter, diff: "<<iter<<" "<<diff_root<<endl;
	return diff_root;
}


void Circuit:: release_ckt_nodes(Tran &tran){
   for(size_t j=0;j<ckt_nodes.size();j++){
         ckt_nodes[j].node = NULL;
   }
}

void Circuit::get_parameters(
		double & epsilon,
		double & omega,
		double & overlap_ratio,
		size_t & max_block_nodes,
		int & mode){
	epsilon		= EPSILON;
	omega		= OMEGA;
	overlap_ratio	= OVERLAP_RATIO; 
	max_block_nodes	= MAX_BLOCK_NODES;
	mode		= MODE;
}

// default values of these parameters are at the begining of this file
void Circuit::set_parameters(
		double epsilon, 
		double omega, 
		double overlap_ratio,
		size_t max_block_nodes,
		int mode){
	EPSILON		= epsilon;
	OMEGA		= omega;
	OVERLAP_RATIO 	= overlap_ratio;
	MAX_BLOCK_NODES	= max_block_nodes;
	MODE		= mode;
}

// choose an appropriate omega for the circuit s.t.
// - node size (use replist)
// - type (c4 or wb)
// - number of layers
void Circuit::select_omega(){
	double omega=OMEGA;
	size_t num_nodes = replist.size();
	//size_t num_layers = layers.size();
	if( num_nodes < 0.05e6 )
		omega=1.0;
	else if (num_nodes < 0.2e6 )
		omega = 1.1;
	else if (num_nodes < 0.3e6 )
		omega = 1.2;
	else if (num_nodes < 0.5e6 )
		omega = 1.3;
	else if (num_nodes < 1.2e6 )
		omega = 1.4;
	else
		omega = 1.5;

	//if( circuit_type == WB && num_layers >= 8 ) omega += 0.2;

	if( circuit_type == C4 ) omega += 0.1;

	if( name == "GND" && num_nodes < 1.2e6) omega -= 0.1;

	if( omega >= 1.6 ) omega = 1.6;
	if( omega <= 1.0 ) omega = 1.0;

	OMEGA = omega;
}

// Randomly choose a number of sample nodes to monitor
void Circuit::get_samples(){
	size_t num_nodes = replist.size();
	srand(time(NULL));
	while(sample.size()<SAMPLE_NUM_NODE){
		int id = rand() % num_nodes;
		sample.push_back(replist[id]);
	}
}

bool Circuit::check_diverge() const{
	for(size_t i=0;i<SAMPLE_NUM_NODE;i++){
		double x = sample[i]->value;
		if(VDD > 0){
			if( x < 0.0 || x > VDD ) return true;
		}
		else
			if(x<0.0) return true;
	}
	return false;
}

void Circuit::boundary_init(int &my_id, int &num_procs){
	assign_bd_base(my_id);

	assign_bd_dd_size(my_id);

	// allocate boundary_nodelist size
	bd_dd_size_g = new int [8*num_procs];	
	
	//if(my_id==0)
		//clog<<"before gathering bd_dd_size. "<<my_id<<endl;
	MPI_Gather(bd_dd_size, 8, MPI_INT, bd_dd_size_g, 
			8, MPI_INT, 0, MPI_COMM_WORLD);

	
	bd_x = new double[bd_size];

	// assign bd_size_g value
	bd_size_g = new int[num_procs];
	MPI_Gather(&bd_size, 1, MPI_INT, bd_size_g, 1, MPI_INT,
			0, MPI_COMM_WORLD);

	total_size = 0;
	if(my_id ==0){
		for(int i=0;i<total_blocks;i++)
			total_size += bd_size_g[i];
	}

	bd_x_g = new double[total_size];
	bd_base_g = new int [num_procs];
	//if(my_id != 0) return;
	int base = 0;
	for(int i=0;i<num_procs;i++){
		bd_base_g[i] = base;
		base += bd_size_g[i];
	}
	bd_base_gd = new int [8*num_procs];
	for(int i=0;i<8*num_procs;i++){
		bd_base_gd[i] = 0;
	}
}

void Circuit::assign_bd_base(int &my_id){
	// boundary base along e, w, n, s directions
	bd_base = new int [8];

	// assign 4 base values for array bd_x
	int  base = 0;
	bd_base[0] = base;
	base += bd_nodelist_sw.size();
	bd_base[1] = base;

	base += bd_nodelist_s.size();
	bd_base[2] = base;

	base += bd_nodelist_se.size();
	bd_base[3] = base;

	base += bd_nodelist_w.size();
	bd_base[4] = base;

	base += bd_nodelist_e.size();
	bd_base[5] = base;

	base += bd_nodelist_nw.size();
	bd_base[6] = base;

	base += bd_nodelist_n.size();
	bd_base[7] = base;
}

void Circuit::assign_bd_dd_size(int &my_id){
	bd_dd_size = new int [8];

	bd_dd_size[0] = bd_nodelist_sw.size();
	bd_dd_size[1] = bd_nodelist_s.size();
	bd_dd_size[2] = bd_nodelist_se.size();
	bd_dd_size[3] = bd_nodelist_w.size();	
	bd_dd_size[4] = bd_nodelist_e.size();
	bd_dd_size[5] = bd_nodelist_nw.size();
	bd_dd_size[6] = bd_nodelist_n.size();
	bd_dd_size[7] = bd_nodelist_ne.size();

	bd_size = 0;
	for(int i=0;i<8;i++)
		bd_size += bd_dd_size[i];
}

void Circuit::assign_internal_base(int &my_id){
	// boundary base along e, w, n, s directions
	internal_base = new int [8];

	// assign 4 base values for array bd_x
	int  base = 0;
	internal_base[0] = base;
	base += internal_nodelist_sw.size();
	internal_base[1] = base;
	base += internal_nodelist_s.size();
	internal_base[2] = base;
	base += internal_nodelist_se.size();
	internal_base[3] = base;
	base += internal_nodelist_w.size();
	internal_base[4] = base;
	base += internal_nodelist_e.size();
	internal_base[5] = base;
	base += internal_nodelist_nw.size();
	internal_base[6] = base;
	base += internal_nodelist_n.size();
	internal_base[7] = base;
}

// assign bd internal arrays
void Circuit::assign_bd_internal_array(int &my_id){
	//clog<<"sw: "<<endl;
	assign_bd_internal_array_dir(internal_base[0], 
		internal_nodelist_sw, internal_x, my_id);
	//clog<<"s: "<<endl;
	assign_bd_internal_array_dir(internal_base[1], 
		internal_nodelist_s, internal_x, my_id);
	//clog<<"se: "<<endl;
	assign_bd_internal_array_dir(internal_base[2], 
		internal_nodelist_se, internal_x, my_id);
	//clog<<"w: "<<endl;
	assign_bd_internal_array_dir(internal_base[3], 
		internal_nodelist_w, internal_x, my_id);
	//clog<<"e: "<<endl;
	assign_bd_internal_array_dir(internal_base[4], 
		internal_nodelist_e, internal_x, my_id);
	//clog<<"nw: "<<endl;
	assign_bd_internal_array_dir(internal_base[5], 
		internal_nodelist_nw, internal_x, my_id);
	//clog<<"n: "<<endl;
	assign_bd_internal_array_dir(internal_base[6], 
		internal_nodelist_n, internal_x, my_id);
	//clog<<"ne: "<<endl;
	assign_bd_internal_array_dir(internal_base[7], 
		internal_nodelist_ne, internal_x, my_id);
}

// assign 4 boundary internal nodes value, store
// them in array bd_x
void Circuit::assign_bd_internal_array_dir(int &base, NodePtrVector & list, double *internal_x, int &my_id){
	Node *nd;
	for(size_t i=0;i<list.size();i++){
		nd = list[i]->rep;
		internal_x[base+i] = nd->value;
		//if(my_id==0)
		//clog<<base+i<<" "<<*nd<<endl;
	}
}

// assign 4 boundary nodes value from receiving
// array bd_x
void Circuit::assign_bd_array(int &my_id){
	int base =0;
	base = bd_base[0];
	assign_bd_array_dir(base, bd_nodelist_sw, my_id);
	
	base = bd_base[1];
	assign_bd_array_dir(base, bd_nodelist_s, my_id);	
	
	base = bd_base[2];
	assign_bd_array_dir(base, bd_nodelist_se, my_id);	

	base = bd_base[3];
	assign_bd_array_dir(base, bd_nodelist_w, my_id);

	base = bd_base[4];
	assign_bd_array_dir(base, bd_nodelist_e, my_id);

	base = bd_base[5];
	assign_bd_array_dir(base, bd_nodelist_nw, my_id);

	base = bd_base[6];
	assign_bd_array_dir(base, bd_nodelist_n, my_id);

	base = bd_base[7];
	assign_bd_array_dir(base, bd_nodelist_ne, my_id);
}

// reset 4 boundary nodes value to 0
void Circuit::reset_bd_array(int &my_id){
	size_t n = bd_nodelist_sw.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_sw[i]->rep->value = 0;
	
	n = bd_nodelist_s.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_s[i]->rep->value = 0;

	n = bd_nodelist_se.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_se[i]->rep->value = 0;

	n = bd_nodelist_w.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_w[i]->rep->value = 0;

	n = bd_nodelist_e.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_e[i]->rep->value = 0;

	n = bd_nodelist_nw.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_nw[i]->rep->value = 0;

	n = bd_nodelist_n.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_n[i]->rep->value = 0;

	n = bd_nodelist_ne.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_ne[i]->rep->value = 0;
}

void Circuit::reset_replist(int &my_id){
	for(size_t i=0;i<replist.size();i++){
		if(replist[i]->isS()!=Y)
			replist[i]->value = 0;
	}
}

// input: internal_x
// output: bd_x
// func: assign internal_x to corresponding bd_x
void Circuit::reorder_bd_x_g(MPI_CLASS &mpi_class){
	// total_blocks to reorder
	// each block boundary is with e, w, n, s bd
	int X_BLOCKS = mpi_class.X_BLOCKS;
	int Y_BLOCKS = mpi_class.Y_BLOCKS;

	int base_p = 0; // stores current block base
	int base_q = 0; // stores nbr block base
	int base_glo_p = 0;
	int base_glo_q = 0;
	int by = 0;
	int bx = 0;
	int bid_nbr = 0;
	int size = 0;

	for(int i=0;i<total_blocks;i++){
		by = i / X_BLOCKS;
		bx = i % X_BLOCKS;
		//clog<<"reorder for "<<i<<" th block."<<endl<<endl;
		// compute block base for i
		base_glo_p = bd_base_g[i];

		// south west nbr dir
		base_p = base_glo_p+bd_base_gd[8*i];
		if(by>=1 && bx>=1){
			bid_nbr = i - X_BLOCKS - 1;
			//clog<<"east bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr+7];
			size = bd_dd_size_g[8*i];
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}
		
		// south nbr dir
		base_p = base_glo_p+bd_base_gd[8*i+1];
		if(by>=1){
			bid_nbr = i - X_BLOCKS;
			//clog<<"south bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr+6];
			size = bd_dd_size_g[8*i+1];
			//clog<<size<<endl;
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}

		// south east nbr dir
		base_p = base_glo_p+bd_base_gd[8*i+2];
		if(by>=1 && bx<X_BLOCKS-1){
			bid_nbr = i - X_BLOCKS + 1;
			//clog<<"east bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr+5];
			size = bd_dd_size_g[8*i+2];
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}

		// west nbr dir
		base_p = base_glo_p+bd_base_gd[8*i+3];
		if(bx>=1){
			bid_nbr = i - 1;
			//clog<<"west bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr+4];
			size = bd_dd_size_g[8*i+3];
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}

		// east nbr dir
		base_p = base_glo_p+bd_base_gd[8*i+4];
		//clog<<"base_glo_p: "<<base_glo_p<<" base_p: "<<base_p<<endl;
		if(bx<X_BLOCKS-1){
			bid_nbr = i+1;
			//clog<<"east bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			//clog<<"base_glo_q: "<<base_glo_q<<" ";
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr+3];
			//clog<<"base_q: "<<base_q<<endl;
			size = bd_dd_size_g[8*i+4];
			//clog<<"size: "<<size<<endl;
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}
		
		// north west nbr dir
		base_p = base_glo_p+bd_base_gd[8*i+5];
		if(by<Y_BLOCKS-1 && bx>=1){
			bid_nbr = i + X_BLOCKS - 1;
			//clog<<"notrhwest bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr+2];
			size = bd_dd_size_g[8*i+5];
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}

		// north nbr dir
		base_p = base_glo_p+bd_base_gd[8*i+6];
		if(by<Y_BLOCKS-1){
			bid_nbr = i + X_BLOCKS;
			//clog<<"north bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr+1];
			size = bd_dd_size_g[8*i+6];
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}

		// north east nbr dir
		base_p = base_glo_p+bd_base_gd[8*i+7];
		if(by<Y_BLOCKS-1 && bx<X_BLOCKS-1){
			bid_nbr = i + X_BLOCKS + 1;
			//clog<<"north east bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr];
			size = bd_dd_size_g[8*i+7];
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}
	}
}

void Circuit::internal_init(int &my_id, int &num_procs){
	assign_internal_base(my_id);

	// allocate boundary_nodelist size
	internal_dd_size = new int [8];
	internal_dd_size[0] = internal_nodelist_sw.size();
	internal_dd_size[1] = internal_nodelist_s.size();
	internal_dd_size[2] = internal_nodelist_se.size();
	internal_dd_size[3] = internal_nodelist_w.size();	
	internal_dd_size[4] = internal_nodelist_e.size();
	internal_dd_size[5] = internal_nodelist_nw.size();
	internal_dd_size[6] = internal_nodelist_n.size();
	internal_dd_size[7] = internal_nodelist_ne.size();

	internal_dd_size_g = new int [8*num_procs];
	
	MPI_Gather(internal_dd_size, 8, MPI_INT, 
		internal_dd_size_g, 8, MPI_INT, 0, 
		MPI_COMM_WORLD);

	internal_size = internal_nodelist_sw.size()+
		  internal_nodelist_s.size()+
		  internal_nodelist_se.size()+
		  internal_nodelist_w.size()+
		  internal_nodelist_e.size()+
		  internal_nodelist_nw.size()+
		  internal_nodelist_n.size()+
		  internal_nodelist_ne.size();

	internal_x = new double[internal_size];

	// assign bd_size_g value
	internal_size_g = new int[num_procs];
	MPI_Gather(&internal_size, 1, MPI_INT, 
			internal_size_g, 1, MPI_INT,
			0, MPI_COMM_WORLD);
	
	total_internal_size = 0;
	if(my_id ==0){
		for(int i=0;i<total_blocks;i++)
			total_internal_size += internal_size_g[i];
	}

	internal_x_g = new double[total_internal_size];
	internal_base_g = new int [num_procs];
	//if(my_id != 0) return;
	int base = 0;
	for(int i=0;i<num_procs;i++){
		internal_base_g[i] = base;
		base += internal_size_g[i];
	}
	internal_base_gd = new int [8*num_procs];
	for(int i=0;i<8*num_procs;i++){
		internal_base_gd[i] = 0;
	}
}

// assign dir boundary node value
void Circuit::assign_bd_array_dir(int &base, NodePtrVector &list, int &my_id){
	Node *p;
	for(size_t i=0;i<list.size();i++){
		p = list[i]->rep;
		p->value = bd_x[base+i];
		// if(my_id==0)
			// clog<<"bd_array: "<<*p<<endl;
	}
}

// assign value back to transient nodes
void Circuit:: save_ckt_nodes(Tran &tran){//, double *x){
   size_t id=0;
   for(size_t j=0;j<ckt_nodes.size();j++){
	 //cout<<"nodes: "<<ckt_nodes[j].node->name<<endl;
         // id = ckt_nodes[j].node->rep->rid;
	 //cout<<"value: "<<x[id]<<endl;
         // ckt_nodes[j].value.push_back(x[id]);
	 ckt_nodes[j].value.push_back(ckt_nodes[j].node->rep->value);

	   // if(my_id==0)
		//   cout<<"ckt_nodes, value: "<<*ckt_nodes[j]<<" "<<ckt_nodes[j].node->rep->value;
      }
}

void Circuit::save_ckt_nodes_to_tr(Tran &tran){
	int index = 0;
	double time = 0;
	int j=0;
	double value = 0;
	for(size_t i=0;i<ckt_nodes.size();i++){
		index = ckt_nodes[i].flag;
		//cout<<"ckt_nodes, index: "<<ckt_nodes[i].node->name
			//<<" "<<ckt_nodes[i].flag<<endl;
		//tran.nodes[index].node = ckt_nodes[i];
		j=0;
		for(time = 0; time < tran.tot_t; time+=tran.step_t){
			value = ckt_nodes[i].value[j];
			//cout<<"time, value: "<<time<<" "<<value<<endl;
			tran.nodes[index].value.push_back(value);
			j++;
		}
	}	
}

// link transient nodes with nodelist
void Circuit:: link_ckt_nodes(Tran &tran, int &my_id){
      for(size_t j=0;j<tran.nodes.size();j++){
	tran.nodes[j].node = get_node(tran.nodes[j].name);
	// if(tran.nodes[j].node)
	// clog<<"tran nodes: "<<*tran.nodes[j].node<<endl;
	}
   // clog<<"tran.node size: "<<tran.nodes.size()<<endl;
   Node_TR_PRINT nodes_temp;
   for(size_t i=0;i<nodelist.size();i++){
      for(size_t j=0;j<tran.nodes.size();j++){
         if(nodelist[i]->name == tran.nodes[j].name){
	    // record the index in tran.nodes
	    nodes_temp.flag = j;
	    // if(my_id==0)
	    	// clog<<"tran.nodes, index: "<<
		//nodelist[i]->name<<" "<<nodes_temp.flag<<endl; 
	    nodes_temp.node = nodelist[i];
	    ckt_nodes.push_back(nodes_temp);
            // record the id for ckt_nodes
            break;
         }
      }
   }
}

// decide transient step current values
void Circuit::current_tr(Net *net, double &time){
	double slope = 0;
	double Tr = net->tr[3]; // Tr
	double PW = Tr + net->tr[5]; // PW
	double Tf = PW + net->tr[4]; // Tf
	double t_temp = time - net->tr[2]; // TD
	double t = fmod(t_temp, net->tr[6]); // Period
	if(time <= net->tr[2])// TD
		net->value = net->tr[0];// V1
	else if(t > 0 && t<= Tr){
		slope = (net->tr[1] - net->tr[0]) / 
			(net->tr[3]);
		net->value = net->tr[0] + t*slope;
	}
	else if(t > Tr && t<= PW)
		net->value = net->tr[1];
	else if(t>PW && t<=Tf){
		slope = (net->tr[0]-net->tr[1])/(net->tr[4]);
		net->value = net->tr[1] + slope*(t-PW);
	}
	else
		net->value = net->tr[0];
	//return current;
}

#if 0
// add Ieq into rhs
// Ieq = i(t) + 2*C / delta_t *v(t)
void Circuit::modify_rhs_c_tr_0(Net *net, double * rhs, double *x, int &my_id){
	double i_t = 0;
	double temp = 0;
	double Ieq = 0;
	//if(my_id==0)
	//clog<<"c net: "<<*net<<" net->flag_bd: "<<net->flag_bd<< endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
        // nk point to Z node
        if(nk->isS() != Z){
		swap<Node *>(nk, nl);
		swap<Node*>(net->ab[0], net->ab[1]);
	}
	//if(my_id==0)
	//clog<<"nk, nl: "<<*nk<<" "<<*nl<<endl;
	size_t k = nk->rid;
	size_t l = nl->rid;
	//if(my_id==0)
	//clog<<"k, l: "<<k<<" "<<l<<" "<<nk->flag_bd<<" "<<nl->flag_bd<<endl;

	Net *r = nk->nbr[TOP];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to Z node
	if(a->isS()!=Z) {
		swap<Node *>(a, b);
		swap<Node*>(r->ab[0], r->ab[1]);
	}
	//clog<<"a, b: "<<*a<<" "<<*b<<endl;

	size_t id_a = a->rid;
	size_t id_b = b->rid;
	//i_t = (b->value - a->value) / r->value;
	i_t = (x[id_b] - x[id_a]) / r->value;
	//if(b->value != x[id_b] || a->value != x[id_a])
	   //cout<<"a, b, x_a, x_b: "<<a->value<<" "<<b->value<<" "<<
	     //x[id_a]<<" "<<x[id_b]<<endl;
	//clog<<"i_t: "<<i_t<<endl;
	//temp = 2*net->value / tran.step_t * 
		//(nk->value - nl->value);
       
        // push 2 nodes into node_set_x
        //clog<<*nk<<" "<<k<<endl;
//#if 0
        pg.node_set_x.push_back(k);
        if(!nl->is_ground()) {
              //clog<<*nl<<" "<<l<<endl;
           pg.node_set_x.push_back(l);
        }
        else if(!b->is_ground()){
              //clog<<*b<<" "<<id_b<<endl;
           pg.node_set_x.push_back(id_b);
        }
//#endif
	//if(my_id==0)
	//	clog<<"before calc temp. "<<endl;
	if(nk->is_ground())
	 //temp = 2*net->value/tran.step_t*(0-x[l]);
	 temp = net->value *(-x[l]);
        else if(nl->is_ground()){
         //temp = 2*net->value/tran.step_t *(x[k]);
	 temp = net->value *x[k];
        }
        else
         //temp = 2*net->value/tran.step_t *(x[k] - x[l]);
	 temp = net->value *(x[k]-x[l]);
	//if(nk->value != x[k] || nl->value != x[l])
	   //cout<<"k, l, x_k, x_l: "<<nk->value<<" "<<nl->value<<" "<<
	     //x[k]<<" "<<x[l]<<endl;
	//clog<<"nk-nl "<<(nk->value - nl->value)<<" "<<2*net->value/tran.step_t<<" "<<temp<<endl;
	
	Ieq  = (i_t + temp);
	//if(my_id==0)
	//clog<< "Ieq is: "<<Ieq<<endl;
	//clog<<"Geq is: "<<2*net->value / tran.step_t<<endl;
	if(!nk->is_ground()&& nk->isS()!=Y){
		 rhs[k] += Ieq;	// for VDD circuit
		 //if(my_id==1)
		    // clog<<k<<" "<<*nk<<" rhs +: "<<rhs[k]<<endl;
	}
	if(!nl->is_ground()&& nl->isS()!=Y){
		 rhs[l] += -Ieq; 
		 //if(my_id==1)
		    // clog<<l<<" "<<*nl<<" rhs +: "<<rhs[l]<<endl;
	}
	//if(my_id==0)
	//	clog<<"finish 1 net. "<<endl;
}

// add Ieq into rhs
// Ieq = i(t) + 2*C / delta_t *v(t)
void Circuit::modify_rhs_c_tr(Net *net, double * rhs, double *x){
	double temp = 0;
	//clog<<"c net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
        
	// nk point to Z node
	size_t k = nk->rid;
	size_t l = nl->rid;

	Net *r = nk->nbr[TOP];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to Z node

	size_t id_a = a->rid;
	size_t id_b = b->rid;
	double i_t = (x[id_b] - x[id_a]) / r->value;
	
	if(nk->is_ground())
	 temp = net->value *(-x[l]);
        else if(nl->is_ground())
	 temp = net->value *x[k];
        else
	 temp = net->value *(x[k]-x[l]);
	
	double Ieq  = i_t + temp;
	if(!nk->is_ground()&& nk->isS()!=Y){
		 rhs[k] += Ieq;	// for VDD circuit
	}
	if(!nl->is_ground()&& nl->isS()!=Y){
		 rhs[l] -= Ieq; 
	}
}
#endif

void Circuit::set_eq_induc(Tran &tran){
	NetPtrVector &ns = net_set[INDUCTANCE];
	for(size_t i=0;i<ns.size();i++)
		ns[i]->value = tran.step_t /(2*ns[i]->value);
}

void Circuit::set_eq_capac(Tran &tran){
	NetPtrVector &ns = net_set[CAPACITANCE];
	for(size_t i=0;i<ns.size();i++)
		ns[i]->value = 2*ns[i]->value/tran.step_t;
}
#if 0
// update rhs by transient nets
void Circuit::modify_rhs_tr_0(double * b, double *x, int &my_id){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type ==CAPACITANCE){	
			for(size_t i=0;i<ns.size();i++){
				modify_rhs_c_tr_0(ns[i], b, x, my_id);
			}
		}
		else if(type == INDUCTANCE){
			for(size_t i=0;i<ns.size();i++){
				modify_rhs_l_tr_0(ns[i], b, x, my_id);	
			}
		}
	}
}

// update rhs by transient nets
void Circuit::modify_rhs_tr(double * b, double *x){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type ==CAPACITANCE){
			for(size_t i=0;i<ns.size();i++)
				modify_rhs_c_tr(ns[i], b, x);
		}
		else if(type == INDUCTANCE){
			for(size_t i=0;i<ns.size();i++){
				modify_rhs_l_tr(ns[i], b, x);
			}
		}
	}
}
#endif
void Circuit::parse_path_table(){
   // build up nodelist info
      Node_G *node;
      for(size_t i=0;i<replist.size();i++){
         node = new Node_G();
         node->value = i;
         pg.nodelist.push_back(node);
      }
}

// push the 8 set of internal bd nodes into node_set_x
void Circuit::push_bd_nodes(Path_Graph &pg, int&my_id){
	int temp = 0;
	// name of the internal_nodelist set
	NodePtrVector internal_set;

	// sw direction
	internal_set = internal_nodelist_sw;
	push_bd_nodes_one_set(pg, my_id, internal_set);
	// s direction
	internal_set = internal_nodelist_s;
	push_bd_nodes_one_set(pg, my_id, internal_set);
	// se direction
	internal_set = internal_nodelist_se;
	push_bd_nodes_one_set(pg, my_id, internal_set);
	// w direction
	internal_set = internal_nodelist_w;
	push_bd_nodes_one_set(pg, my_id, internal_set);
	// e direction
	internal_set = internal_nodelist_e;
	push_bd_nodes_one_set(pg, my_id, internal_set);
	// nw direction
	internal_set = internal_nodelist_nw;
	push_bd_nodes_one_set(pg, my_id, internal_set);
	// n direction
	internal_set = internal_nodelist_n;
	push_bd_nodes_one_set(pg, my_id, internal_set);
	// ne direction
	internal_set = internal_nodelist_ne;
	push_bd_nodes_one_set(pg, my_id, internal_set);
}

// solve transient version
bool Circuit::solve_tr_step(int &num_procs, int &my_id, MPI_CLASS &mpi_class){
	int iter = 0;	
	double diff=0;
	bool successful = false;
	double time=0;
	double t1, t2;		

	t1= MPI_Wtime();
	while( iter < MAX_ITERATION ){
		diff = solve_iteration_tr(my_id, iter, num_procs, mpi_class);
		iter++;
		// if(my_id ==0)
			// clog<<"iter, diff: "<<iter<<" "<<diff<<endl;
		if( diff < EPSILON ){
			successful = true;
			break;
		}
	}
	// if(my_id==0)
		// clog<<"after solve tr. "<<endl;
	// copy the values from replist to nodelist
	for(size_t i=0;i<nodelist.size()-1;i++){
		Node *nd = nodelist[i];
		if(nd->name != nd->rep->name)
			nd->value = nd->rep->value;
	}
	//if(my_id ==0)
		//clog<<"iter: "<<iter<<endl;

	t2 = MPI_Wtime();
	time = t2-t1;
	return successful;
}

void Circuit::solve_DC(int &num_procs, int &my_id, MPI_CLASS &mpi_class){
	// stores 4 boundary base into bd_base_gd
	MPI_Gather(bd_base, 8, MPI_INT, bd_base_gd, 8, MPI_INT,
		0, MPI_COMM_WORLD);
	
	MPI_Gather(internal_base, 8, MPI_INT, internal_base_gd,
		8, MPI_INT, 0, MPI_COMM_WORLD);
		
	int iter = 0;	
	double diff=0;
	bool successful = false;

	// before iteration, copy boundary nodes value to corresponding blocks
	assign_bd_internal_array(my_id);
	MPI_Gatherv(internal_x, internal_size, 
		MPI_DOUBLE, 
		internal_x_g, internal_size_g, 
		internal_base_g, MPI_DOUBLE, 0, 
		MPI_COMM_WORLD);
	
	// reorder boundary array according to nbrs
	if(my_id==0)	reorder_bd_x_g(mpi_class);

	double time=0;
	clock_t start_e, end_e;
	start_e = clock();
	double t1= MPI_Wtime();
	while( iter < MAX_ITERATION ){
		diff = solve_iteration(my_id, iter, num_procs, mpi_class);
		iter++;
		// if(my_id ==0)
			// clog<<"iter, diff: "<<iter<<" "<<diff<<endl;
		if( diff < EPSILON ){
			successful = true;
			break;
		}
	}
	end_e = clock();
	double t2 = MPI_Wtime();
	time = t2-t1;
	double c_time = 1.0*(end_e-start_e)/CLOCKS_PER_SEC;
	// copy the values from replist to nodelist
	for(size_t i=0;i<nodelist.size()-1;i++){
		Node *nd = nodelist[i];
		if(nd->name != nd->rep->name)
			nd->value = nd->rep->value;
	}
	if(my_id==0){
		clog<<"DC iter: "<<iter<<endl;
		clog<<"solve DC cost (mpi): "<<time<<" secs."<<endl;
		clog<<"solve DC cost (c_time): "<<1.0*(end_e-start_e)/CLOCKS_PER_SEC<<" secs. "<<endl;
	}
	// get_voltages_from_block_LU_sol();
	// if(my_id==0)
		// clog<<nodelist<<endl;
}

// check if matrix is SPD or not
void Circuit::check_matrix(Matrix &A){
	size_t count = 0;
	A.merge();
	clog<<A.get_row()<<endl;
	for(size_t i=0;i<A.size();i++){
		cout<<A.Ti[i]+1<<" "<<A.Tj[i]+1<<" "<<A.Tx[i]<<endl;
		//if(A.Ti[i]!=A.Tj[i])
			//cout<<A.Tj[i]+1<<" "<<A.Ti[i]+1<<" "<<A.Tx[i]<<endl;
		//if(A.Ti[i]==A.Tj[i] && A.Tx[i]<10)
			//cout<<"diagonal:  "<<A.Ti[i]<<" "<<A.Tj[i]<<" "<<A.Tx[i]<<endl;
	}
	return;
	
	for(size_t i=0;i<A.size();i++){
		size_t row = A.Ti[i];
		size_t col = A.Tj[i];
		if(A.Ti[i] != A.Tj[i]){// && count<1){
			count++;
			if(A.Tx[i]>=0)
			cout<<"pos: "<<A.Tx[i]<<endl;
			double row_sum=0;
			double col_sum=0;
			for(size_t j=0;j<A.size();j++){
				if(A.Ti[j]==row && A.Tj[j]==row)
					row_sum += A.Tx[j];
				if(A.Ti[j]==col && A.Tj[j]==col)
					col_sum += A.Tx[j];

			}
			clog<<"row / col sum: "<<row_sum<<" "<<col_sum<<endl;
			if(A.Tx[i]+row_sum<0){
				cout<<"error data: "<<A.Tx[i]<<" "<<row_sum<<" "<<A.Ti[i]<<" "<<A.Tj[i]<<endl;
			}
			if(A.Tx[i]+col_sum<0){
				cout<<"error data: "<<A.Tx[i]<<" "<<col_sum<<" "<<A.Ti[i]<<" "<<A.Tj[i]<<endl;
			}
		}
	}
}

// update block 4 corners
void Circuit::update_geometry(int my_id, MPI_CLASS &mpi_class){
	// compute the geometrical information for the blocks
	x_min = mpi_class.block_geo[0];
	y_min = mpi_class.block_geo[1];
	x_max = mpi_class.block_geo[2];
	y_max = mpi_class.block_geo[3];
	// clog<<"ckt bd: "<<x_min<<" "<<y_min<<" "<<x_max<<" "<<y_max<<endl;

	num_blocks = NUM_BLOCKS_X * NUM_BLOCKS_Y;
	if(my_id==0)
		clog<<"total num_blocks for one core: "<<num_blocks<<endl;
	block_vec.clear();
	for(int i=0; i<num_blocks;i++){
		Block *temp_block = new Block();
		block_vec.push_back(temp_block);
	}

	double delta_x = (x_max - x_min) / NUM_BLOCKS_X;
	double delta_y = (y_max - y_min) / NUM_BLOCKS_Y;
	double base_x = x_min;
	double base_y = y_min;
	double new_base_x;
	double new_base_y;
	// range for a block is lx <= x < ux and
	// ly <= y < uy
	for(size_t j=0;j<NUM_BLOCKS_Y;j++){
		new_base_y = base_y + j*delta_y;
		for(size_t i=0;i<NUM_BLOCKS_X;i++){
			size_t id_block = j*NUM_BLOCKS_X + i;
			new_base_x = base_x + i*delta_x;
			block_vec[id_block]->lx = new_base_x;
			block_vec[id_block]->ly = new_base_y;  
			block_vec[id_block]->ux = new_base_x + delta_x;
			block_vec[id_block]->uy = new_base_y + delta_y;

			// if(my_id==0)
				// clog<<block_vec[id_block]->lx<<" "<<block_vec[id_block]->ux<<" "<<block_vec[id_block]->ly<<" "<<block_vec[id_block]->uy<<endl;
		}
	}
}

// assign circuit nodes into blocks
void Circuit::assign_block_nodes(int my_id){
	Node *nd = NULL;
	for(size_t i=0;i<nodelist.size();i++)
		if(nodelist[i]->is_ground()){
			nd = nodelist[i];
			break;
		}
	for(size_t j=0;j<block_vec.size();j++){
		block_vec[j]->nd_GND = nd;
	}

	if(my_id==0) clog<<"replist.size(): "<<replist.size()<<endl;
	for(size_t i=0;i<replist.size();i++){
		nd = replist[i];
		for(size_t j=0;j<block_vec.size();j++){
			// if a node belongs to 
			// some block
			if(block_vec[j]->node_in_block(nd)){
				block_vec[j]->replist.push_back(nd);
			}
		}
	}

	// sort internal nodes of blocks
	for(size_t i=0;i<block_vec.size();i++){
		block_vec[i]->count = block_vec[i]->replist.size();
		block_vec[i]->sort_nodes();
		/*if(my_id==0)
			for(size_t j=0;j<block_vec[i]->replist.size();j++)
				cout<<"j, node: "<<*block_vec[i]->replist[j]<<" "<<j<<endl;*/
		block_vec[i]->build_nd_IdMap();
	}
}

// assign circuit nets into blocks
void Circuit::assign_block_nets(int my_id){
	Net *net;
	Node *na, *nb;
	int net_flag = 0;
	
	for(int j=0;j<block_vec.size();j++){
		for(int type = 0; type < NUM_NET_TYPE;type++)
			block_vec[j]->net_set[type].clear();
		block_vec[j]->bd_netlist.clear();
	}

	int N_blocks = block_vec.size();
	// first handle internal nets of ckt
	for(int type= 0;type < NUM_NET_TYPE; type++){	
		NetList & ns = net_set[type];
		for(size_t i=0;i<ns.size();i++){
			net = ns[i];	
			net_flag = 0;
			
			for(size_t j=0;j<N_blocks;j++){
				net_flag = block_vec[j]->net_in_block(net);
				
				if(net_flag == 0)
					continue;
				if(net_flag == 2)
					block_vec[j]->net_set[type].push_back(net);

				else if(net_flag ==1){
					block_vec[j]->bd_netlist.push_back(net);
				}
			}
		}
	}
	int type = RESISTOR;
	if(my_id==0) clog<<"bd_netlist: "<<bd_netlist.size()<<endl;
	// then handle boundary nets of ckt
	for(size_t i=0;i<bd_netlist.size();i++){
		net = bd_netlist[i];
		net_flag = 0;
		for(int j=0;j<block_vec.size();j++){
			net_flag = block_vec[j]->net_in_block(net);

			// if(my_id==0)
				// cout<<"block, net, flag: "<<i<<" "<<j<<" "<<*net<<" "<<net_flag<<endl;
			if(net_flag == 2)
				block_vec[j]->net_set[type].push_back(net);
			else if(net_flag ==1)
				block_vec[j]->bd_netlist.push_back(net);
		}
	}
}

// build boundary netlist for circuit
void Circuit::build_bd_netlist(){
	int type = RESISTOR;
	NetPtrVector & ns = net_set[type];
	Net *net = NULL;
	for(size_t i=0;i<ns.size();i++){
		net = ns[i];
		if(net->flag_bd ==1)
			bd_netlist.push_back(net);
	}
}

void Circuit::block_solve_Jacobi(int my_id){
	double iter = 0;
	double diff = 0;
	while(iter ==0 || (diff > EPSILON) && iter <1){
		diff = 0;
		// first update bnewp for all the blocks
		for(size_t i=0;i<block_vec.size();i++){
			block_vec[i]->update_rhs(block_vec[i]->bnewp, block_vec[i]->bp, my_id);	
		}
		// then solve all the blocks with Jacobi
		// assign all the values into nodes
		// new rhs store in bnewp and solve
		for(size_t i=0; i < block_vec.size();i++){
			/*if(my_id==0){
				cout<<"solving block: "<<i<<endl;
				cout<<"bd lx, ly, ux, ly: "<<block_vec[i]->lx<<" "<<block_vec[i]->ly<<" "<<block_vec[i]->ux<<" "<<block_vec[i]->uy<<endl;
			}*/
			block_vec[i]->solve_CK_DC_Jacobi(my_id);
		}
		for(size_t i=0; i < block_vec.size();i++){
			// if(my_id==0)
				// clog<<"modify block: "<<i<<" "<<block_vec[i]->lx<<" "<<block_vec[i]->ly<<" "<<block_vec[i]->ux<<" "<<block_vec[i]->uy<<endl;
			double local_diff = 
				block_vec[i]->modify_voltage(my_id);
			// if(my_id==0)
				//clog<<"my_id, i, local_diff: "<<my_id<<" "<<i<<" "<<local_diff<<endl;
			if(local_diff > diff)
				diff = local_diff;
		}
		// if(my_id==0)
			// clog<<"my_id, iter, diff: "<<my_id<<" "<<iter<<" "<<diff<<endl;
		iter ++;
		//if(my_id==0)
			//cout<<"solution: "<<nodelist<<endl;
	}
}

void Circuit::block_solve(int my_id){
	double iter = 0;
	double diff = 0;
	while(iter ==0 || (diff > EPSILON) && iter <1){
		diff = 0;
		// new rhs store in bnewp and solve
		for(size_t i=0; i < block_vec.size();i++){
			/*if(my_id==0){
				cout<<"solving block: "<<i<<endl;
				cout<<"bd lx, ly, ux, ly: "<<block_vec[i]->lx<<" "<<block_vec[i]->ly<<" "<<block_vec[i]->ux<<" "<<block_vec[i]->uy<<endl;
			}*/
			block_vec[i]->solve_CK_DC(my_id);
			// if(my_id==0)
				// clog<<"modify block: "<<i<<" "<<block_vec[i]->lx<<" "<<block_vec[i]->ly<<" "<<block_vec[i]->ux<<" "<<block_vec[i]->uy<<endl;
			double local_diff = 
				block_vec[i]->modify_voltage(my_id);
			// if(my_id==0)
				// clog<<"my_id, i, local_diff: "<<my_id<<" "<<i<<" "<<local_diff<<endl;
			if(local_diff > diff)
				diff = local_diff;
		}
		// if(my_id==0)
			// clog<<"my_id, iter, diff: "<<my_id<<" "<<iter<<" "<<diff<<endl;
		iter ++;
		//if(my_id==0)
			//cout<<"solution: "<<nodelist<<endl;
	}
}

double Circuit::find_diff(int my_id){
	double diff = 0;
	for(size_t i=0;i<block_vec.size();i++){
		for(size_t j=0;j<block_vec[i]->count;j++){
			double local_diff = fabs(block_vec[i]->xp[j] - block_vec[i]->x_old[j]);
			if(local_diff > diff)
				diff = local_diff;
		}
	}
	// return the max vol diff of 2 iters
	return diff;
}

// clean the variables for exploring
void Circuit::clean_explore(){
	x_list.clear();
	y_list.clear();
	x_list_bd_map.clear();
	y_list_bd_map.clear();
	x_list_nd_map.clear();
	y_list_nd_map.clear();
	lx = -1;
	ly = -1;
	ux = -1;
	uy = -1;
	// clog<<"before pre release circuit. "<<endl;
	// clean nodelist and netlist
	pre_release_circuit();	
	// clog<<"after pre release circuit. "<<endl;
}

void Circuit::push_bd_nodes_one_set(Path_Graph &pg, int&my_id, NodePtrVector internal_set){
	Node *nd;
	// s direction
	size_t n_s = internal_set.size();
	// if(my_id==0)
	// cout<<"my_id, ns size: "<<my_id<<" "<<n_s<<endl;
	for(size_t i=0;i<n_s;i++){
		nd = internal_set[i]->rep;
		for(int j=0;j<block_vec.size();j++){
			if(!block_vec[j]->node_in_block(nd))
				continue;
			// if(my_id==0)
			// cout<<"block j push nd: "<<j<<" "<<*nd<<" "<<nd->pt<<endl;
			block_vec[j]->push_nd_pg_x(nd);
		}
	}
}

void Circuit::push_bd_net_nodes(){
	size_t n_s = bd_netlist.size();
	Net *net;
	Node *a = NULL; 
	Node *b = NULL;
	for(size_t i=0;i<n_s;i++){
		net = bd_netlist[i];
		if(net == NULL) 
			continue;
		a = net->ab[0]->rep;
		b = net->ab[1]->rep;
		for(int j=0;j<block_vec.size();j++){
			if(block_vec[j]->node_in_block(a))
				block_vec[j]->push_nd_pg_x(a);

			if(block_vec[j]->node_in_block(b))
				block_vec[j]->push_nd_pg_x(b);
		}
	}
}
