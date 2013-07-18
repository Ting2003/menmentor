// ----------------------------------------------------------------//
// Filename : block.cpp
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// implementation of block class
// ----------------------------------------------------------------//
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added opeartor !=
//   * added this log

#include <cassert>
#include "cholmod.h"
#include "block.h"
#include "node.h"
#include "util.h"
//#include "umfpack.h"

Block::Block(size_t _count):
        b_ck(NULL),
	b_new_ck(NULL),
	x_old(NULL),
	//bp(NULL),
	//bnewp(NULL),
	//xp(NULL),
	x_ck(NULL),	
	count(_count),
	lx(-1.0), ly(-1.0),
	ux(-1.0), uy(-1.0){
	
	Lx = NULL;
	Li = NULL;
	Lp = NULL;
	Lnz = NULL;
}

Block::~Block(){
    /* for(size_t i=0;i<boundary_netlist.size();i++)
    	delete boundary_netlist[i];
    boundary_netlist.clear();
    */
    
    A.clear();
    nd_IdMap.clear();
    // free Li, Lx and so on
    delete [] Lx;
    delete [] Lp;
    delete [] Lnz;
    delete [] Li;

    replist.clear();
    bd_netlist.clear();
    delete [] net_set;
    delete [] x_old;
    delete [] xp;
    delete [] bnewp;
    delete [] bp;
    delete [] bnewp_temp;
}

void Block::free_block_cholmod(cholmod_common *cm){
    cholmod_free_factor(&L, cm);
    cholmod_free_dense(&b_ck, cm);
    cholmod_free_dense(&b_new_ck, cm);
    cholmod_free_dense(&x_ck, cm);
    cholmod_free_dense(&bnew_temp, cm);
}

void Block::CK_decomp(Matrix & A, cholmod_common *cm){
	Algebra::CK_decomp(A, L, cm);
}

void Block::solve_CK_tr(){
	solve_eq_sp(xp, bnewp_temp);
#if 0
	x_ck = cholmod_solve(CHOLMOD_A, L, bnew_temp, cm);
	xp = static_cast<double *>(x_ck->x);
#endif
	// for(size_t i=0;i<count;i++)
		// clog<<"i, xp: "<<i<<" "<<xp[i]<<endl;
	//cholmod_solve_new(CHOLMOD_A, L, b_new_ck, x_ck, cm);
}

// start cm and allocate resources
void Block::allocate_resource(){
	cm = &c;
	cholmod_start(cm);
	cm->print = 5;

	if( count == 0 ) return;

	b_ck = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
	x_ck = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
	bp = static_cast<double*>(b_ck->x);
	xp = static_cast<double*>(x_ck->x);
	b_new_ck = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
	bnewp = static_cast<double*>(b_new_ck->x);
	bnew_temp = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
	bnewp_temp = static_cast<double*>(bnew_temp->x);

	x_old = new double [count];
}

// update rhs of each block with its boundary netlist
void Block::update_rhs(double *bnewp, double *bp, int &my_id){
	size_t k=0, l=0;

	int temp = 0;
	// bnewp  = bp
	copy_array(bnewp, bp);
	//b_new = b;
	/*for(size_t i=0;i<count;i++){
		// bnewp[i] = bp[i];
		if(my_id==temp && bnewp[i]!=bp[i])
			clog<<"i, bp: "<<i<<" "<<bp[i]<<" "<<bnewp[i]<<endl;
	}*/

	//if;(my_id==0)
		//clog<<"bd_netlist.size(): "<<bd_netlist.size()<<endl;
	// for each net in this block
	for(size_t i=0;i<bd_netlist.size();i++){
		Net * net = bd_netlist[i];
		// if(my_id==temp)
			// cout<<endl<<"net: "<<*net<<endl;
		double G = 1.0/net->value;

		Node * a = net->ab[0]->rep;
		Node * b = net->ab[1]->rep;
		// if(my_id==0)
			//cout<<"a, b: "<<*a<<" "<<*b<<endl;
	
		// if a is inside block
		//if(a->flag_bd == 0){
		if(node_in_block(a)){
			k = nd_IdMap[a];//a->rid;
			if(a->isS()!=Y){
				/*if(my_id==temp){
					cout<<"a inside block: "<<k<<" "<<*a<<" "<<G<<" "<<*b<<endl;
				}*/
				bnewp[k] += G * b->value;
				/*if(my_id==temp)
					cout<<"k, bnewp: "<<k<<" "<<bnewp[k]<<endl;*/
			}
		}
		else if(node_in_block(b)){
		//if(b->flag_bd ==0){
			l = nd_IdMap[b];//b->rid;
			if(b->isS()!=Y){
				/*if(my_id==temp){
					cout<<"b inside block: "<<l<<" "<<*b<<" "<<G<<" "<<*a<<endl;
				}*/
				bnewp[l] += G * a->value;
				/*if(my_id==temp)
					cout<<"l, bnewp: "<<l<<" "<<bnewp[l]<<endl;*/
			}
		}
	} // end of for i
}

/////////////////////////////////////////////////////////////////
// methods for BlockInfo
bool compare_node_ptr(const Node * a, const Node * b){
	if( a->is_ground() ) return false;
	if (b->is_ground() ) return true;

	if( a->pt.y == b->pt.y ){
		if( a->pt.x == b->pt.x ){
			if( a->pt.z == b->pt.z ){
				return (a->isS() > b->isS());
			}
			else{
				return (a->pt.z > b->pt.z);// top down
			}
		}
		else
			return ( a->pt.x < b->pt.x );
	}
	else
		return (a->pt.y < b->pt.y);
}

void Block::sort_nodes(){
	sort(replist.begin(), replist.end(), compare_node_ptr);
}

// judge whether a node is within a block
bool Block::node_in_block(Node *nd){
	long x = nd->pt.x;
	long y = nd->pt.y;
	// clog<<"block bd: "<<lx<<" "<<ly<<" "<<ux<<" "<<uy<<endl;
	// if a node belongs to some block
	if(x>=lx && x <=ux && 
		y>=ly && y<=uy){
		return true;
	}
	return false;
}

// judge whether a net is within a block
// 2: internal net of a block
// 1: bundary net of a block
// 0: outside net of a block
int Block::net_in_block(Net *net){
	Node *na, *nb;
	na = net->ab[0];
	nb = net->ab[1];
	bool flag_a = false;
	bool flag_b = false;

	if(!na->is_ground()){
		flag_a = node_in_block(na);
	}
	if(!nb->is_ground()){
		flag_b = node_in_block(nb);
	}	

	if(na->is_ground() && flag_b == true)
		return 2;
	if(nb->is_ground() && flag_a == true)
		return 2;
	if(flag_a == true && flag_b == true)
		return 2;
	if(flag_a == true || flag_b == true)
		return 1;
	return 0;
}

// 1. copy node voltages from the circuit to a Vec
//    from = true then copy circuit to x
//    else copy from x to circuit
// 2. map block voltage into global_index
void Block::copy_node_voltages_block(){
	size_t id;
	// copy node voltages from nodelist
	for(size_t i=0;i<replist.size();i++){
		Node *node = replist[i];
		id = nd_IdMap[node];//node->rid;
		xp[id] = replist[i]->value;
	}
}

void Block::stamp_matrix(int &my_id, MPI_CLASS &mpi_class){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		switch(type){
		case RESISTOR:
			// if(my_id==0)
				// clog<<"resis net. "<<ns.size()<<endl;
			for(it=ns.begin();it!=ns.end();++it){
				Net * net = *it;
				if( net == NULL ) continue;
				assert( fzero(net->value) == false );
				stamp_resistor(my_id, *it);
			}
			break;
		case CURRENT:
			// clog<<"current net. "<<ns.size()<<endl;
			for(it=ns.begin();it!=ns.end();++it){
				stamp_current(my_id, (*it), mpi_class);
			}
			break;
		case VOLTAGE:
			// if(my_id==0)
				// clog<<"VDD net: "<<ns.size()<<endl;
			for(it=ns.begin();it!=ns.end();++it){
				if( fzero((*it)->value)  && 
				    !(*it)->ab[0]->is_ground() &&
				    !(*it)->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_VDD(my_id,(*it));
			}
			break;
		case CAPACITANCE:
			break;
		case INDUCTANCE:
			// if(my_id==0)
				// clog<<"induc net: "<<ns.size()<<endl;
			for(size_t i=0;i<ns.size();i++)
				stamp_inductance_dc(ns[i], my_id);
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
	// if(my_id==0)
		// clog<<"bd_netlist.size: "<<bd_netlist.size()<<endl;
	// then stamp boundary nets
	for(size_t i=0;i<bd_netlist.size();i++){
		stamp_bd_net(my_id, bd_netlist[i]);
	}
	//if(my_id==0)
		//clog<<A<<endl;
	/*if(my_id==0){
		for(int i=0;i<10;i++)
			clog<<"b origin: "<<i<<" "<<block_info.bp[i]<<endl;
	}*/
	make_A_symmetric(bp, my_id);
	
	A.set_row(count);
	/*if(my_id==0){
		//check_matrix(A);
		cout<<"before CK_decomp. "<<endl;
		cout<<A<<endl;
	}*/
	if(count >0){
		CK_decomp(A, cm);
		// A.clear();
		//if(cm->status ==1)
			//clog<<" non SPD: "<<my_id<<" "<<cm->status<<endl;
	}
	//if(my_id==0)
		//clog<<"after CK_decomp. "<<endl;
}

// =========== stamp block version of matrix =======
// internal resistor net
void Block::stamp_resistor(int &my_id, Net * net){
	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
	
	double G;	
	G = 1./net->value;
	// if(nd[0]->isS() == Y || nd[1]->isS() == Y)
		// clog<<endl<<"resis net: "<<*net<<endl;	
	int count = 0;
	for(size_t j=0;j<2;j++){
		Node *nk = nd[j], *nl = nd[1-j];
		// if(my_id==0&& j==0)
			// clog<<"resis net: "<<*net<<endl;
		// if( nk->isS()!=Y && nl->isS()!=Y) 
		{
			size_t k1 = nd_IdMap[nk];//nk->rid;
			size_t l1 = nd_IdMap[nl];//nl->rid;

			// if(nk->isS() == Y || nl->isS() == Y)
				// clog<<"nk, nl: "<<*nk<<" "<<k1<<" "<<*nl<<" "<<l1<<" "<<nk->isS() <<" "<<nl->isS()<<endl;
			// search nk's nbr, skip insert for vol and induc nbr net
			if(nk->isS()!=Y && !nk->is_ground()){
				bool flag = false;
				for(size_t i=0;i<nk->nbr_vec.size();i++){
					Net *net = nk->nbr_vec[i];
					if(net->type == INDUCTANCE || 
					   net->type == VOLTAGE){
					   flag = true;
					   break;
					}
				}
				// if no inductance or voltage nbr net, stamp
				if(flag == false){
					// if(nk->isS() == Y || nl->isS() == Y)
						// clog<<"push: ("<<k1<<","<<k1<<","<<G<<")"<<endl;
					A.push_back(k1, k1, G);
				}

				if( !nl->is_ground()
					&& l1 < k1){
					bool flag_nl = false;				
					for(size_t i=0;i<nl->nbr_vec.size();i++){
						Net *net = nl->nbr_vec[i];
						if(net->type == INDUCTANCE || 
								net->type == VOLTAGE){
							flag_nl = true;
							break;
						}
					}
					if(flag_nl == false){
						A.push_back(k1,l1,-G);
						// if(nk->isS() == Y || nl->isS() == Y)
							// clog<<"push: ("<<k1<<","<<l1<<","<<-G<<endl;
					}
				}
			}
		}
	}// end of for j		
}

// only bottom layer has current net, no need to check Voltage net
void Block::stamp_current(int &my_id, Net * net, MPI_CLASS &mpi_class){
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;

	// if(my_id==0)
		// clog<<"cur net: "<<*net<<endl;
	// only stamp for internal node
	if( !nk->is_ground() && nk->isS()!=Y && node_in_block(nk)){//nk->flag_bd ==0) {
		size_t k = nd_IdMap[nk];//nk->rid;

		bp[k] += -net->value;

		// if(my_id==0) clog<<"bk: "<<*nk<<" "<<k<<" "<<bp[k]<<endl;
		//pk[k] += -net->value;
	}
	if( !nl->is_ground() && nl->isS()!=Y && node_in_block(nl)){//nl->flag_bd ==0) {
		size_t l = nd_IdMap[nl];//nl->rid;

		bp[l] += net->value;
		// if(my_id==0) clog<<"bk: "<<*nl<<" "<<l<<" "<<bp[l]<<endl;
		//pl[l] +=  net->value;
	}
}

void Block::stamp_VDD(int &my_id, Net * net){
	// find the non-ground node
	Node * X = net->ab[0]->rep;
	//if(my_id==0) clog<<"net: "<<*net<<endl;
	if( X->is_ground() ) X = net->ab[1]->rep;

	if(!node_in_block(X)) return;
	// if(X->rep->flag_bd ==1) return;
	// do stamping for internal node
	long id =nd_IdMap[X];//X->rep->rid;

	// if(my_id==0) clog<<"VDD net: "<<*net<<endl;
		// clog<<"id ("<<id<<","<<id<<",1)"<<endl;
	A.push_back(id, id, 1.0);
	
	bool flag = false;
	for(size_t i=0;i<X->nbr_vec.size();i++){
		Net *net_temp = X->nbr_vec[i];
		if(net_temp->type == CURRENT){
			flag = true;
			break;	
		}
	}
	// if vol node connects to current net
	if(flag == true){
		bp[id] = net->value;
		// clog<<"id=, bp: "<<id<<" "<<bp[id]<<endl;
	}
	else{
		bp[id] += net->value;	
		// clog<<"id+, bp: "<<id<<" "<<bp[id]<<endl;
	}
	// if(my_id==0) clog<<"id, bp: "<<id<<", "<<bp[id]<<endl;
}

// all cores stamp dc inductance
void Block::stamp_inductance_dc(Net * net, int &my_id){
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nd_IdMap[nk];//nk->rid;
	size_t l = nd_IdMap[nl];// nl->rid;
	G = 1./net->value;
	// if(my_id==0)
		// clog<<"induc net: "<<*net<<endl;
	if( nk->isS()!=Y && !nk->is_ground()&&node_in_block(nk)){//nk->flag_bd==0) {
		A.push_back(k,k, 1);
		// if(my_id==0)
			// clog<<"push ("<<*nk<<" "<<k<<" "<<k<<" 1)"<<endl;
		// general stamping
		if(!nl->is_ground())
		// A.push_back(k,l,-1);
		// make it symmetric
			bp[k] = bp[l];
		//clog<<"("<<k<<" "<<k<<" "<<1<<")"<<endl;
		//clog<<"("<<k<<" "<<l<<" "<<-1<<")"<<endl;
	}

	if( nl->isS() !=Y && !nl->is_ground() &&node_in_block(nl)){//nl->flag_bd==0) {
		A.push_back(l,l, 1);
		// if(my_id==0)
			// clog<<"push ("<<*nl<<" "<<l<<" "<<l<<" 1)"<<endl;
		if(!nk->is_ground())
		// general stamping
		// A.push_back(l,k,-1);
		bp[l] = bp[k];

		//clog<<"("<<l<<" "<<l<<" "<<1<<")"<<endl;
		//clog<<"("<<l<<" "<<k<<" "<<-1<<")"<<endl;
	}
}

// search for the nodes connects to inductance and voltage sources
void Block::make_A_symmetric(double *b, int &my_id){
	int type = VOLTAGE;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p=NULL, *q=NULL, *r =NULL;

	for(size_t i=0;i<ns.size();i++){
	   Net *net = ns[i];
           if( net == NULL ) continue;
	   Node *na = net->ab[0]->rep;
	   if(na->isS()!=Y)
		na = net->ab[1]->rep;

	   p = na; 
	   // node p points to Y node
	   for(size_t j=0;j<na->nbr_vec.size();j++){
		Net *net_nbr = na->nbr_vec[j];
		if(net_nbr->type != RESISTOR) 
			continue;
		// modify rhs of resistor net
		if(net_nbr->ab[0]->rep->name == na->name){
			q = net_nbr->ab[1]->rep;
		}
		else if(net_nbr->ab[1]->rep->name == na->name){
			q = net_nbr->ab[0]->rep;
		}

	   	/*if(my_id==0){
		   clog<<"modify rhs net: "<<*net_nbr<<endl;
		   clog<<"p, q: "<<*p<<" "<<*q<<" "<<endl;
	   	}*/
           	size_t id = nd_IdMap[q];// q->rid;
           	double G = 1.0 / net_nbr->value;

          	// clog<<"id, old_b: "<<id<<" "<<b[id]<<endl; 
           	b[id] += p->value * G;
	   	// if(my_id==0)
		   // clog<<"id, new b: "<<id<<" "<<b[id]<<endl;
	    }
        }
}

void Block::copy_array(double *x_old, double *xp){
	for(size_t j=0;j<count;j++)
		x_old[j] = xp[j];	
}

// update rhs and solve CK
void Block::solve_CK_DC_Jacobi(int my_id){
	if(count<=0)
		return;
	x_ck = cholmod_solve(CHOLMOD_A, L, b_new_ck, cm);
	// solve_CK(cm);
	xp = static_cast<double *>(x_ck->x);
	/*if(my_id==0)
	for(int i=0;i<replist.size();i++)
		if(xp[i]!= replist[i]->value)
		cout<<"xp: "<<i<<" "<<*replist[i]<<" "<<xp[i]<<endl;*/
		
}


// update rhs and solve CK
void Block::solve_CK_DC(int my_id){
	// new rhs store in bnewp and solve
	//double *bnewp_old;
	//bnewp_old = new double [count];
	//copy_array(bnewp_old, bnewp);
	update_rhs(bnewp, bp, my_id);
	/*for(size_t i=0;i<count;i++)
		if(bnewp[i]!=bp[i])
			cout<<"i, bold, bnew: "<<i<<" "<<bnewp[i]<<" "<<bp[i]<<endl;
	delete [] bnewp_old;*/
	// x_old = xp
	// copy_array(x_old, xp);

	if(count<=0)
		return;
	x_ck = cholmod_solve(CHOLMOD_A, L, b_new_ck, cm);
	// solve_CK(cm);
	xp = static_cast<double *>(x_ck->x);
	/*if(my_id==0)
	for(int i=0;i<replist.size();i++)
		// if(xp[i]!= replist[i]->value)
		cout<<"replist, xp: "<<*replist[i]<<" "<<xp[i]<<endl;//, xp:  "<<*replist[i]<<" "<<bnewp[i]<<" "<<xp[i]<<endl;*/
		
}

double Block::modify_voltage(int &my_id){
	double max_diff = 0.0;
	//if(get_name()=="VDDA") OMEGA = 1.0;
	//else OMEGA = 1.15;
	double OMEGA = 1.0;
	for(size_t i=0;i<count;i++){
		//xp[i] = (1-OMEGA)*x_old[i] + OMEGA*
			// xp[i];
		double vol_old = replist[i]->value;
		// update block nodes value
		replist[i]->value = xp[i];
		double diff = fabs(replist[i]->value - vol_old);
		// double diff = fabs(replist[i]->value - x_old[i]);
		if( diff > max_diff ) max_diff = diff;
		//if(my_id==0 && diff >= 0.54)
			// clog<<"0 id: rep, old, diff: "<<lx<<" "<<ly<<" "<<ux<<" "<<uy<<" "<<vol_old<<" "<<*replist[i]<<" "<<max_diff<<endl;
	}
	return max_diff;
}

// build block node id map
void Block::build_nd_IdMap(){
	pair<Node *, size_t> id_pair;
	for(size_t i=0;i<replist.size();i++){
		id_pair.first = replist[i];
		id_pair.second = i;
		nd_IdMap.insert(id_pair);
	}
}

void Block::stamp_matrix_tr(int &my_id, MPI_CLASS &mpi_class, Tran &tran){
	// cout<<A<<endl;	
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		switch(type){
		case RESISTOR:
			/*for(it=ns.begin();it!=ns.end();++it){
				Net * net = *it;
				if( net == NULL ) continue;
				assert( fzero(net->value) == false );
				stamp_resistor_tr(my_id, *it);
			}*/
			break;
		case CURRENT:
			break;
		case VOLTAGE:
			for(it=ns.begin();it!=ns.end();++it){
				if( fzero((*it)->value)  && 
				    !(*it)->ab[0]->is_ground() &&
				    !(*it)->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_VDD_tr(my_id,(*it));
			}
			break;
		case CAPACITANCE:
			for(size_t i=0;i<ns.size();i++)
				stamp_capacitance_tr(ns[i], tran, my_id);
			break;
		case INDUCTANCE:
			// for(size_t i=0;i<ns.size();i++)
				// stamp_inductance_tr(ns[i], tran, my_id);
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

// only stamp resistor node connected to inductance
void Block::stamp_resistor_tr(int &my_id, Net * net){
	// no boundary nets here
	// if(net->flag_bd ==1)
		// return;

	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
	double G;	
	G = 1./net->value;

	for(size_t j=0;j<2;j++){
		Node *nk = nd[j], *nl = nd[1-j];

		// skip resistor layer net
		if(nk->isS()!=Y && nl->isS()!=Y) continue;
		if(nk->isS()!=Y) continue;

		size_t k1 = nd_IdMap[nk];//nk->rid;
		size_t l1 = nd_IdMap[nl];//nl->rid;

		// if(my_id==0) clog<<"net: "<<*net<<endl;
		if( !nk->is_ground()){
			bool flag = false;
			for(size_t i=0;i<nk->nbr_vec.size();i++){
				Net *nbr_net = nk->nbr_vec[i];
				if(nbr_net->type == INDUCTANCE){
					flag = true;
					break;
				}
			}
			if(flag == true)	
				A.push_back(k1,k1, G);
		}
		// also need to push back connection
		if(!nl->is_ground()){
			bool flag = false;
			for(size_t i=0;i<nk->nbr_vec.size();i++){
				Net *nbr_net = nk->nbr_vec[i];
				if(nbr_net->type == INDUCTANCE){
					flag = true;
					break;
				}
			}
			if(flag == false){ 
				if(l1 < k1){
					//if(my_id==0) clog<<"push ("<<k1<<","<<l1<<","<<-G<<")"<<endl;
					A.push_back(k1,l1,-G);
				}
				else if(l1>k1){
					//if(my_id==0) clog<<"push ("<<l1<<","<<k1<<","<<-G<<")"<<endl;
					A.push_back(l1, k1, -G);
				}
			}
		}
	}// end of for j	
}

// stamp a voltage source
void Block::stamp_VDD_tr(int &my_id, Net * net){
	//if(my_id==0)
		// clog<<"VDD net: "<<*net<<endl;
	// find the non-ground node
	Node * X = net->ab[0]->rep;
	if( X->is_ground() ) X = net->ab[1]->rep;

	if(!node_in_block(X)) return;
	// if(X->rep->flag_bd ==1) return;
	// do stamping for internal node
	long id = nd_IdMap[X];//X->rep->rid;
	// A.push_back(id, id, 1.0);
	bool flag = false;
	for(size_t i=0;i<X->nbr_vec.size();i++){
		Net *nbr_net = X->nbr_vec[i];
		if(nbr_net->type == CURRENT){
			flag = true;
			break;
		}
	}
	if(flag == true){
		// this node connects to a VDD and a current
		// the current should be stamped
		//assert( feqn(1.0, q[id]) ); 
		assert( feqn(1.0, bp[id]) );
		bp[id] = net->value;
		/*if(my_id==0)
			clog<<"id, bp: "<<id<<" "<<bp[id]<<endl;*/
		//q[id] = net->value;	    // modify it
	}
	else{
		bp[id] += net->value;
		/*if(my_id==0)
			clog<<"id, bp: "<<id<<" "<<bp[id]<<endl;*/
		//q[id] += net->value;
	}
}

// stamp capacitance Geq = 2C/delta_t
void Block::stamp_capacitance_tr(Net *net, Tran &tran, int &my_id){
	/*if(my_id==0)
		clog<<"cap net: "<<*net<<endl;*/
	double Geq = 0;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nd_IdMap[nk];//nk->rid;
	size_t l = nd_IdMap[nl];//nl->rid;
	// Geq = 2*C / delta_t
	Geq = (2*net->value) / tran.step_t;
	//net->value = Geq;
	//clog<<"C delta_t Geq: "<<net->value<<" "<<tran.step_t<<" "<<Geq<<endl;
	// Ieq = i(t) + 2*C / delta_t * v(t)

	if( nk->isS()!=Y  && !nk->is_ground() && node_in_block(nk)){//nk->flag_bd==0) {
		A.push_back(k,k, Geq);
		/*if(my_id==0)
			clog<<"("<<k<<" "<<k<<" "<<Geq<<")"<<endl;*/
		if(!nl->is_ground()&& k > l){
			A.push_back(k,l,-Geq);
			/*if(my_id==0)
				clog<<"("<<k<<" "<<l<<" "<<-Geq<<")"<<endl;*/
		}
	}

	if( nl->isS() !=Y && !nl->is_ground() && node_in_block(nl)){//nl->flag_bd==0) {
		A.push_back(l,l, Geq);
		/*if(my_id==0)
			clog<<"("<<l<<" "<<l<<" "<<Geq<<")"<<endl;*/
		if(!nk->is_ground()&& l > k){
			A.push_back(l,k,-Geq);
			/*if(my_id==0)
				clog<<"("<<l<<" "<<k<<" "<<-Geq<<")"<<endl;*/
		}
	}
}

// stamp inductance Geq = delta_t/(2L)
void Block::stamp_inductance_tr(Net * net, Tran &tran, int &my_id){
	/*if(my_id==0)
		clog<<"induc net: "<<*net<<endl;*/
	double Geq = 0;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nd_IdMap[nk];//nk->rid;
	size_t l = nd_IdMap[nl];//nl->rid;
	// Geq = delta_t / (2*L)
	Geq = tran.step_t / (2*net->value);
	//net->value = Geq;

	if( nk->isS()!=Y  && !nk->is_ground() && node_in_block(nk)){//nk->flag_bd==0) {
		// -1 is to clear formal inserted 1 at (k,k)
		A.push_back(k,k, Geq-1);
		/*if(my_id==0)
			clog<<"("<<k<<" "<<k<<" "<<Geq-1<<")"<<endl;*/
		//clog<<nl->isS()<<endl;
		if(!nl->is_ground()&& nl->isS()!=Y && k>l){
			A.push_back(k,l,-Geq);
		        /*if(my_id==0)
			clog<<"("<<k<<" "<<l<<" "<<-Geq<<")"<<endl;*/
		}
	}

	if( nl->isS() !=Y && !nl->is_ground() && node_in_block(nl)){//nl->flag_bd==0) {
		// -1 is to clear formal inserted 1 at (l,l)
		A.push_back(l,l, Geq-1);
		/*if(my_id==0)
			clog<<"("<<l<<" "<<l<<" "<<Geq-1<<")"<<endl;*/
		if(!nk->is_ground() && nk->isS()!=Y && l>k){
			A.push_back(l,k,-Geq);
			/*if(my_id==0)
				clog<<"("<<l<<" "<<k<<" "<<-Geq<<")"<<endl;*/
		}
	}
}

// reset b and bnew
void Block::reset_array(double *bp){
	for(size_t i=0;i<count;i++){
		bp[i] = 0;
	}
}

// make A symmetric for tran: resis around voltage sources
void Block::make_A_symmetric_tr(int &my_id, Tran &tran){
	int type = VOLTAGE;
	NetList & ns = net_set[type];
	NetList::iterator it;
	// Node *p=NULL, *q=NULL, *r =NULL;
	Net *net = NULL;
	Node *nb = NULL;

	for(it=ns.begin();it!=ns.end();it++){
           if( (*it) == NULL ) continue;
	   // clog<<"modify net: "<<*(*it)<<endl;
	   Node *na = (*it)->ab[0]->rep;
	   if(na->isS()!=Y){
		na= (*it)->ab[1]->rep;
	   }
	   for(size_t i=0;i<na->nbr_vec.size();i++){
		net = na->nbr_vec[i];
		if(net->type != RESISTOR)
			continue;
		nb = net->ab[0]->rep;
		if(nb->name == na->name)
			nb = net->ab[1]->rep;
	   	// na to Y, nb to resistor node
	   	// clog<<"na, nb: "<<*na<<" "<<*nb<<endl;
           	size_t id = nd_IdMap[nb];
	   	double G = 1.0/net->value;
           
           	bp[id] += xp[nd_IdMap[na]]*G;

	   	// if(my_id==0)
		   // clog<<"id, bp: "<<id<<" "<<bp[id]<<endl;
	   }
	}
}

// stamp transient current values into rhs
void Block::stamp_current_tr(int &my_id, double &time){
	NetPtrVector & ns = net_set[CURRENT];
	for(size_t i=0;i<ns.size();i++)
		stamp_current_tr_net(ns[i], time, my_id);
}

void Block::stamp_current_tr_net(Net * net, double &time, int &my_id){
	current_tr(net, time);
	/*if(my_id==0){
		clog<<"net: "<<*net<<endl;
		// clog<<"current: "<<current<<endl;
	}*/
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground()&& nk->isS()!=Y) { 
		size_t k = nd_IdMap[nk];//nk->rid;
		//clog<<"node, rid: "<<*nk<<" "<<k<<endl;
		bp[k] += -net->value;//current;
		// if(my_id==0)
			// clog<<"time, k, b: "<<time<<" "<<k<<" "<<bp[k]<<endl;
	}
	if( !nl->is_ground() && nl->isS()!=Y) {
		size_t l = nd_IdMap[nl];//nl->rid;
		//clog<<"node, rid: "<<*nl<<" "<<l<<endl;
		bp[l] +=  net->value;// current;
		// if(my_id==0)
			// clog<<"time, l, b: "<<time<<" "<<l<<" "<<bp[l]<<endl;
	}
}

void Block::clear_A(){
	A.clear();
}

void Block::CK_decomp(){
	Algebra::CK_decomp(A, L, cm);
}

void Block::copy_vec(double *bnewp, double *bp){
	for(size_t i=0;i<count;i++){
		bnewp[i] = bp[i];
	}
}

// update rhs by transient nets
void Block::modify_rhs_tr_0(double * b, double *x, int &my_id){
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

// stamp transient current values into rhs
void Block::stamp_current_tr_1(double &time){
	NetPtrVector & ns = net_set[CURRENT];
	for(size_t i=0;i<ns.size();i++)
		stamp_current_tr_net_1(ns[i], time);
}

void Block::stamp_current_tr_net_1(Net * net, double &time){
	double diff = 0;
	double current = net->value;
	current_tr(net, time);
	//if(time / 1e-11 == 22)
	//cout<<"time, co, cn: "<<time<<" "<<current<<" "<<net->value<<endl;
	// only stamps when net got a different current
	if(current != net->value){
		diff = net->value - current;
		
		//cout<<"time, old, new, diff: "<<time <<" "<<current<<" "<<net->value<<" "<<diff<<endl;
		//cout<<"net: "<<*net;
		//clog<<"current: "<<current<<endl;
		Node * nk = net->ab[0]->rep;
		Node * nl = net->ab[1]->rep;
		if( !nk->is_ground()&& nk->isS()!=Y) { 
			size_t k = nd_IdMap[nk];//nk->rid;
			//clog<<"node, rid: "<<*nk<<" "<<k<<endl;
			//clog<<"time, k, b bef: "<<time<<" "<<k<<" "<<b[k]<<endl;
			bnewp[k] += -diff;//current;
			bp[k] = bnewp[k];
			//clog<<"time, k, b: "<<time <<" "<<k<<" "<<b[k]<<endl;
		}
		if( !nl->is_ground() && nl->isS()!=Y) {
			size_t l = nd_IdMap[nl];//nl->rid;
			//clog<<"time, l, b bef: "<<time<<" "<<l<<" "<<b[l]<<endl;
			//clog<<"node, rid: "<<*nl<<" "<<l<<endl;
			bnewp[l] +=  diff;// current;
			bp[l] = bnewp[l];
			//clog<<"time, l, b: "<<time<<" "<<l<<" "<<b[l]<<endl;
		}
	}
}

// update rhs by transient nets
void Block::modify_rhs_tr(double * b, double *x){
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

// add Ieq into rhs
// Ieq = i(t) + 2*C / delta_t *v(t)
void Block::modify_rhs_c_tr(Net *net, double * rhs, double *x){
	double temp = 0;
	//clog<<"c net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
        if(nk->is_ground()){
		nk = net->ab[1]->rep;
		nl = net->ab[0]->rep;
	} 
	// nk point to Z node
	size_t k = nd_IdMap[nk];//nk->rid;
	size_t l = nd_IdMap[nl];//nl->rid;

	Net * nbr_resis;
	for(size_t i=0;i<nk->nbr_vec.size();i++){
		if(nk->nbr_vec[i]->type == RESISTOR){
			nbr_resis = nk->nbr_vec[i];
			break;
		}
	}
	Net *r = nbr_resis; //nk->nbr[TOP];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to Z node

	size_t id_a = nd_IdMap[a];//a->rid;
	size_t id_b = nd_IdMap[b];//b->rid;
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

// add Ieq into rhs
// Ieq = i(t) + delta_t / (2*L) *v(t)
void Block::modify_rhs_l_tr(Net *net, double *rhs, double *x){
	//clog<<"l net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
	// nk point to X node
	//if(nk->isS() !=X){ 
		//swap<Node*>(nk, nl);
	//}
	size_t k = nd_IdMap[nk];//nk->rid;
	size_t l = nd_IdMap[nl];//nl->rid;
	double Ieq = 0;

	double i_t = 0;
	double temp = 0;
	//temp = tran.step_t / (2*net->value) * 
		//(nl->value - nk->value);
	//temp = tran.step_t / (2*net->value)*(x[l] - x[k]);
	temp = net->value *(x[l] - x[k]);	
	//if(nk->value != x[k] || nl->value != x[l])
	   //clog<<"k, l, x_k, x_l: "<<nk->value<<" "<<nl->value<<" "<<
	     //x[k]<<" "<<x[l]<<endl;

	//clog<<"delta_t/2L, nl-nk, temp: "<<tran.step_t / (2*net->value)<<" "<<(nl->value-nk->value)<<" "<<temp<<endl;
	Net *nbr_resis;
	for(size_t i=0;i<nk->nbr_vec.size();i++){
		Net *nbr_net = nk->nbr_vec[i];
		if(nbr_net->type == RESISTOR){
			nbr_resis = nbr_net;
			break;
		}
	}

	Net *r = nbr_resis;//nk->nbr[BOTTOM];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to X node
	//if(a->isS()!=X) {
		//swap<Node*>(a, b);
	//}
	size_t id_a = nd_IdMap[a];//a->rid;
	size_t id_b = nd_IdMap[b];//b->rid;
	i_t = (x[id_a] - x[id_b]) / r->value;
	//i_t = (a->value - b->value) / r->value;
        //if(b->value != x[id_b] || a->value != x[id_a])
	   //clog<<"a, b, x_a, x_b: "<<a->value<<" "<<b->value<<" "<<
	     //x[id_a]<<" "<<x[id_b]<<endl;

	//clog<<"resiste r: "<<*r<<endl;
	//clog<<*a<<" "<<*b<<endl;
	//clog<<"a, b, r, i_t: "<<a->value<<" "<<b->value<<" "<<
		//r->value<<" "<<i_t<<endl;
       
        // push inductance nodes into node_set_x
        //clog<<*nk<<" "<<k<<endl;
        //clog<<*b<<" "<<id_b<<endl;
 #if 0
        if(iter==0){
           pg.node_set_x.push_back(k);
           pg.node_set_x.push_back(id_b);
        }
#endif
	Ieq  = i_t + temp;
	//clog<<"Ieq: "<<Ieq<<endl;
	if(nk->isS() !=Y && !nk->is_ground()){
		 rhs[k] += Ieq; // VDD circuit
		//clog<<*nk<<" "<<rhs[k]<<endl;
	}
	if(nl->isS()!=Y && !nl->is_ground()){
		 rhs[l] += -Ieq; // VDD circuit
		//clog<<*nl<<" "<<rhs[l]<<endl;
	}
}

// decide transient step current values
void Block::current_tr(Net *net, double &time){
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

// add Ieq into rhs
// Ieq = i(t) + 2*C / delta_t *v(t)
void Block::modify_rhs_c_tr_0(Net *net, double * rhs, double *x, int &my_id){
	double i_t = 0;
	double temp = 0;
	double Ieq = 0;
	// if(my_id==0)
	//clog<<"c net: "<<*net<<" net->flag_bd: "<<net->flag_bd<< endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
	// nk points to non-zero node
	if(nk->is_ground()){
		nk = net->ab[1]->rep;
		nl = net->ab[0]->rep;
	}
#if 0
        // nk point to Z node
        if(nk->isS() != Z){
		swap<Node *>(nk, nl);
		swap<Node*>(net->ab[0], net->ab[1]);
	}
#endif
	// if(my_id==0)
	//clog<<"nk, nl: "<<*nk<<" "<<*nl<<endl;
	size_t k = nd_IdMap[nk];//nk->rid;
	size_t l = nd_IdMap[nl];//nl->rid;
	//if(my_id==0)
	//clog<<"k, l: "<<k<<" "<<l<<" "<<nk->flag_bd<<" "<<nl->flag_bd<<endl;
	Net *nbr_resis = NULL;
	for(size_t i=0;i<nk->nbr_vec.size();i++){
		Net *nbr_net = nk->nbr_vec[i];
		if(nbr_net->type == RESISTOR){
			nbr_resis = nbr_net;
			// clog<<"nbr_resis: "<<*nbr_resis<<endl;
			break;
		}
	}

	Net *r = nbr_resis; //nk->nbr[TOP];
	// if(r!=NULL)
	//clog<<"net_r: "<<*r<<endl;
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to Z node
	if(a->isS()!=Z) {
		swap<Node *>(a, b);
		swap<Node*>(r->ab[0], r->ab[1]);
	}
	// clog<<"a, b: "<<*a<<endl;//" "<<*b<<endl;

	size_t id_a = nd_IdMap[a];//a->rid;
	size_t id_b = nd_IdMap[b];//b->rid;
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
	// if(my_id==0)
		//clog<<"finish 1 net. "<<endl;
}

// add Ieq into rhs
// Ieq = i(t) + delta_t / (2*L) *v(t)
void Block::modify_rhs_l_tr_0(Net *net, double *rhs, double *x, int &my_id){
	//clog<<"l net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
	// nk point to X node
	if(nk->isS() !=X){ 
		swap<Node*>(nk, nl);
		swap<Node*>(net->ab[0], net->ab[1]);
	}
	size_t k = nd_IdMap[nk];//]nk->rid;
	size_t l = nd_IdMap[nl];//nl->rid;
	double Ieq = 0;

	double i_t = 0;
	double temp = 0;
	//temp = tran.step_t / (2*net->value) * 
		//(nl->value - nk->value);
	//temp = tran.step_t / (2*net->value)*(x[l] - x[k]);
	temp = net->value *(x[l] - x[k]);
	//if(nk->value != x[k] || nl->value != x[l])
	   //clog<<"k, l, x_k, x_l: "<<nk->value<<" "<<nl->value<<" "<<
	     //x[k]<<" "<<x[l]<<endl;

	//clog<<"delta_t/2L, nl-nk, temp: "<<tran.step_t / (2*net->value)<<" "<<(nl->value-nk->value)<<" "<<temp<<endl;
	Net *nbr_resis;
	for(size_t i=0;i<nk->nbr_vec.size();i++){
		Net *nbr_net = nk->nbr_vec[i];
		if(nbr_net->type == RESISTOR){
			nbr_resis = nbr_net;
			break;
		}
	}

	Net *r = nbr_resis;//nk->nbr[BOTTOM];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to X node
	if(a->isS()!=X) {
		swap<Node*>(a, b);
		swap<Node*>(r->ab[0], r->ab[1]);
	}
	size_t id_a = nd_IdMap[a];//a->rid;
	size_t id_b = nd_IdMap[b];//b->rid;
	i_t = (x[id_a] - x[id_b]) / r->value;
	//i_t = (a->value - b->value) / r->value;
        //if(b->value != x[id_b] || a->value != x[id_a])
	   //clog<<"a, b, x_a, x_b: "<<a->value<<" "<<b->value<<" "<<
	     //x[id_a]<<" "<<x[id_b]<<endl;

	//clog<<"resiste r: "<<*r<<endl;
	//clog<<*a<<" "<<*b<<endl;
	//clog<<"a, b, r, i_t: "<<a->value<<" "<<b->value<<" "<<
		//r->value<<" "<<i_t<<endl;
       
        // push inductance nodes into node_set_x
        //clog<<*nk<<" "<<k<<endl;
        //clog<<*b<<" "<<id_b<<endl;
// #if 0
        pg.node_set_x.push_back(k);
        pg.node_set_x.push_back(id_b);
//#endif
	Ieq  = i_t + temp;
	if(nk->isS() !=Y && !nk->is_ground()){
		 rhs[k] += Ieq; // VDD circuit
		// if(my_id==0)
			// cout<<"k, Ieq, rhs: "<<k<<" "<<Ieq<<" "<<rhs[k]<<endl;
	}
	if(nl->isS()!=Y && !nl->is_ground()){
		 rhs[l] += -Ieq; // VDD circuit

		// if(my_id==0)
			// cout<<"l, Ieq, rhs: "<<l<<" "<<Ieq<<" "<<rhs[l]<<endl;
	}
}

void Block::stamp_bd_net(int my_id, Net *net){
	// if(net->type != RESISTOR)
		// return;

	Node *na= NULL;
	Node *nb = NULL;

	na = net->ab[0]->rep;
	nb = net->ab[1]->rep;

	size_t id;
	if(node_in_block(na)){
		id = nd_IdMap[na];
		// if(my_id==0)
			// clog<<"bd net: "<<*net<<" "<<*na<<" "<<id<<endl;

		A.push_back(id, id, 1.0/net->value);
		// if(my_id==0 && id ==567)
			//clog<<"bd net: "<<id<<" "<<id<<" "<<1.0/net->value<<endl;
	}
	else if(node_in_block(nb)){
		id = nd_IdMap[nb];

		// if(my_id==0)
			// clog<<"bd net: "<<*net<<" "<<*nb<<" "<<id<<endl;
		A.push_back(id, id, 1.0/net->value);
		// if(my_id==0 && id ==567)
			//clog<<"bd net: "<<id<<" "<<id<<" "<<1.0/net->value<<endl;
	}
}

void Block::build_id_map(){
   double *temp;
   int *id_map;
   id_map = new int [count];
   cholmod_build_id_map(CHOLMOD_A, L, cm, id_map);

   temp = new double [count];
   // then substitute all the nodes rid
   for(size_t i=0;i<count;i++){
	int id = id_map[i];
	replist[id]->rid = i;
	temp[i] = bp[i];
   }

   for(size_t i=0;i<count;i++)
	bp[i] = temp[id_map[i]];
   for(size_t i=0;i<count;i++)
        temp[i] = xp[i];
   for(size_t i=0;i<count;i++)
        xp[i] = temp[id_map[i]];
   delete [] temp;
   delete [] id_map;
}

void Block::push_nd_set_bx(Tran &tran){
   //if(flag_ck ==1){
      // push rhs node into node_set b
      for(size_t i=0;i<count;i++){
         if(bnewp[i] !=0)
            pg.node_set_b.push_back(i);
     // } 

      // push back all nodes in output list
      vector<size_t>::iterator it;
      size_t id;
      for(size_t i=0;i<tran.nodes.size();i++){
         if(tran.nodes[i].node == NULL) continue;
	 // only handles inside block nd
	 if(!node_in_block(tran.nodes[i].node->rep))
		continue;
         if(!tran.nodes[i].node->rep->is_ground()){
            id = tran.nodes[i].node->rep->rid;
            it = find(pg.node_set_x.begin(), pg.node_set_x.end(), id);
            if(it == pg.node_set_x.end()){
               pg.node_set_x.push_back(id);
            }
         } 
      }

      // then push back boundary nodes into node_se_x
      // push_bd_nodes(pg, my_id); 

      // get path_b, path_x, len_path_b, len_path_x
      // build_path_graph();
      //clog<<"len_b, x to n: "<<len_path_b<<" "<<
        // len_path_x<<" "<<count<<endl;
   }
   // clog<<"build up path time: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;   
}
 
void Block::build_path_graph(){
   build_FFS_path();

   build_FBS_path();

   // only keep the 2 paths, switch from List into array
   len_path_b = pg.path_FFS.get_size();
   len_path_x = pg.path_FBS.get_size();

   path_b = new int[len_path_b];
   path_x = new int [len_path_x];
   
   Node_G *nd;
   nd = pg.path_FFS.first;
   for(int i=0;i<len_path_b;i++){
      path_b[i] = nd->value;
      if(nd->next != NULL)
         nd = nd->next;
   }
   pg.path_FFS.destroy_list();

   nd = pg.path_FBS.first;
   for(int i=0;i<len_path_x;i++){
      path_x[i] = nd->value;
      if(nd->next != NULL)
         nd = nd->next;
   }
   pg.path_FBS.destroy_list();

   pg.nodelist.clear();
   pg.node_set_b.clear();
   pg.node_set_x.clear();
}

void Block::build_FFS_path(){
   parse_path_table(); 

   set_up_path_table();

   find_path(pg.node_set_b, pg.path_FFS);

   pg.path_FFS.assign_size();

   for(size_t i=0;i<replist.size();i++)
      pg.nodelist[i]->flag = 0;
}

void Block::build_FBS_path(){
  pg.nodelist.clear();
  parse_path_table();
  set_up_path_table();
  find_path(pg.node_set_x, pg.path_FBS);
   pg.path_FBS.assign_size();
}

void Block::parse_path_table(){
   // build up nodelist info
      Node_G *node;
      for(size_t i=0;i<replist.size();i++){
         node = new Node_G();
         node->value = i;
         pg.nodelist.push_back(node);
      }
}

void Block::set_up_path_table(){
   size_t n = L->n;
   //int *Lp, *Li, *Lnz;
   int p, lnz, s, e;
   //Lp = static_cast<int *> (L->p);
   //Lnz = (int *)L->nz;
   //Li = static_cast<int *> (L->i);
   for(size_t i=0;i<n;i++){
      p = Lp[i];
      lnz = Lnz[i];

      s = Li[p];
      e = s;
      if(lnz >1) 
         e = Li[p+1];

      if(s<e)
         pg.nodelist[s]->next = pg.nodelist[e];
   }
}

bool compare_Node_G(const Node_G *nd_1, const Node_G *nd_2){
   return (nd_1->value < nd_2->value);
 }

void Block::find_path(vector<size_t> &node_set, List_G &path){
   Node_G* ne = pg.nodelist[pg.nodelist.size()-1];
   vector <Node_G *> insert_list;
   sort(node_set.begin(), node_set.end());
   if(node_set.size() == 0) return;
   
   // build up the first path start with id = min 
   int id = node_set[0];
   do{
      path.add_node(pg.nodelist[id]);
      pg.nodelist[id]->flag =1;
      if(pg.nodelist[id]->next == NULL) break;
      pg.nodelist[id] = pg.nodelist[id]->next;
   }while(pg.nodelist[id]->value != ne->value);
   path.add_node(ne);

   for(size_t i=0; i<node_set.size();i++){
      int id = node_set[i];
      if(pg.nodelist[id]->flag == 1) continue;
      // stops at first place where flag = 0
      // clog<<"stops at i: "<<i<<endl;
      do{
         if(pg.nodelist[id]->flag ==0){
            insert_list.push_back(pg.nodelist[id]);
            pg.nodelist[id]->flag =1;
         }
         if(pg.nodelist[id]->next == NULL || 
           pg.nodelist[id]->next->flag ==1)
            break;
         // else clog<<"next node: "<<*pg.nodelist[id]->next;
         pg.nodelist[id] = pg.nodelist[id]->next;
      }while(pg.nodelist[id]->value != ne->value); 
   }

   //clog<<"insert_list.size: "<<insert_list.size()<<endl;
   sort(insert_list.begin(), insert_list.end(), compare_Node_G);
   //for(int i=0;i<insert_list.size();i++)
      //clog<<"i, insert: "<<i<<" "<<*insert_list[i]<<endl;

   //clog<<"path: "<<&path<<endl;
   // p is the old pointer to the list
   // will be updated into new one
   Node_G *q=NULL;
   Node_G *p=NULL;
//#if 0
   for(size_t k=0;k<insert_list.size();k++){
      if(k ==0) p = path.first;
      else p = q;
      q = path.insert_node(insert_list[k], p);
   }
//#endif

   //clog<<"path: "<<&path<<endl;
   //clog<<endl;
   insert_list.clear();
}

// add a nodd intp pg.node_set_x
void Block::push_nd_pg_x(Node *nd){
	vector<size_t>::iterator it;
	size_t id;
	id = nd->rid;
	it = find(pg.node_set_x.begin(), pg.node_set_x.end(), id);
	if(it == pg.node_set_x.end())
		pg.node_set_x.push_back(id);
}

 // find super node columns for path_b and path_x
void Block::find_super(){
    s_col_FFS = new int [len_path_b];
    s_col_FBS = new int [len_path_x];
    int p, lnz;
    int j, k;
    // FFS loop
    for(k=0;k<len_path_b;){
       j = path_b[k];
       p = Lp[j];
       lnz = Lnz[j];
       if (lnz < 4 || path_b[k+1]!=j+1 || lnz != Lnz [j+1] + 1
         || Li [p+1] != j+1){
          s_col_FFS[k]= 1;
          k++;
       }
       else if (path_b[k+2]!= j+2 || lnz != Lnz [j+2] + 2
         || Li [p+2] != j+2){
         s_col_FFS[k]=2;
         k+=2;
       }
       else{
          s_col_FFS[k]=3;
          k+=3;
       }
    }
    //FBS loop
    for(k=len_path_x-1;k>=0;){
       j = path_x[k];
       p = Lp[j];
       lnz = Lnz[j];
       if (j < 4 || path_x[k-1]!=j-1||lnz != Lnz [j-1] - 1
         || Li [Lp [j-1]+1] != j){
         s_col_FBS[k]=1;
         k--;
       }
       else if (path_x[k-2] != j-2 ||lnz != Lnz [j-2]-2 ||
         Li [Lp [j-2]+2] != j){
         s_col_FBS[k]=2;
         k-=2;
       }
       else{
          s_col_FBS[k]=3;
          k-=3;
       }
    }
}
 
void Block::solve_eq_sp(double *X, double *bnewp){
    int p, q, r, lnz, pend;
    int j, k, n = L->n ;
    for(int i=0;i<n;i++){
       X[i] = bnewp[i];
    }
    // FFS solve
    for(k=0; k < len_path_b;){
       j = path_b[k];
       //for (j = 0 ; j < n ; ){
       /* get the start, end, and length of column j */
       p = Lp [j] ;
       lnz = Lnz [j] ;
       pend = p + lnz ;
 
       if (s_col_FFS[k]==1)//lnz < 4 || lnz != Lnz [j+1] + 1 || Li [p+1] != j+1)
       {
 
          /* -------------------------------------------------------------- */
          /* solve with a single column of L */
          /* -------------------------------------------------------------- */
          double y = X [j] ;
          if(L->is_ll == true){
             X[j] /= Lx [p] ;
          }
          for (p++ ; p < pend ; p++)
          {
             X [Li [p]] -= Lx [p] * y ;
          }
          k++ ;  /* advance to next column of L */
 
       }
       else if (s_col_FFS[k]==2)//lnz != Lnz [j+2] + 2 || Li [p+2] != j+2)
       {
          {
             double y [2] ;
             q = Lp [j+1] ;
             if(L->is_ll == true){
                y [0] = X [j] / Lx [p] ;
                y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
                X [j  ] = y [0] ;
                X [j+1] = y [1] ;
             }
 
             else{
                y [0] = X [j] ;
                y [1] = X [j+1] - Lx [p+1] * y [0] ;
                X [j+1] = y [1] ;
             }
             for (p += 2, q++ ; p < pend ; p++, q++)
             {
                X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] ;
             }
 
          }
          k += 2 ;           /* advance to next column of L */
 
       }
       else
       {
          {
             double y [3] ;
             q = Lp [j+1] ;
             r = Lp [j+2] ;
             //#ifdef LL
             if(L->is_ll == true){
                y [0] = X [j] / Lx [p] ;
                y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
                y [2] = (X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1]) / Lx [r] ;
                X [j  ] = y [0] ;
                X [j+1] = y [1] ;
                X [j+2] = y [2] ;
             }
 
             else{
                y [0] = X [j] ;
                y [1] = X [j+1] - Lx [p+1] * y [0] ;
                y [2] = X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1] ;
                X [j+1] = y [1] ;
                X [j+2] = y [2] ;
             }
             for (p += 3, q += 2, r++ ; p < pend ; p++, q++, r++)
             {
                X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] + Lx [r] * y [2] ;
             }
          }
             // marched to next 2 columns of L
          k += 3;
       }
    }
    // FBS solve
    for(k = len_path_x - 1; k >=0;){
       j = path_x[k];
       //for(j = n-1; j >= 0; ){
 
       /* get the start, end, and length of column j */
       p = Lp [j] ;
       lnz = Lnz [j] ;
       pend = p + lnz ;
 
       /* find a chain of supernodes (up to j, j-1, and j-2) */
 
        if (s_col_FBS[k]==1)//j < 4 || lnz != Lnz [j-1] - 1 || Li [Lp [j-1]+1] != j)
       {
 
          /* -------------------------------------------------------------- */
          /* solve with a single column of L */
          /* -------------------------------------------------------------- */
 
          double d = Lx [p] ;
          if(L->is_ll == false)
             X[j] /= d ;
          for (p++ ; p < pend ; p++)
          {
             X[j] -= Lx [p] * X [Li [p]] ;
          }
          if(L->is_ll == true)
             X [j] /=  d ;
          k--;
       }
       else if (s_col_FBS[k]==2)//lnz != Lnz [j-2]-2 || Li [Lp [j-2]+2] != j)
       {
          {
             double y [2], t ;
             q = Lp [j-1] ;
             double d [2] ;
             d [0] = Lx [p] ;
             d [1] = Lx [q] ;
             t = Lx [q+1] ;
             if(L->is_ll == false){
                y [0] = X [j  ] / d [0] ;
                y [1] = X [j-1] / d [1] ;
             }
             else{
                y [0] = X [j  ] ;
                y [1] = X [j-1] ;
             }
             for (p++, q += 2 ; p < pend ; p++, q++)
             {
                int i = Li [p] ;
                y [0] -= Lx [p] * X [i] ;
                y [1] -= Lx [q] * X [i] ;
             }
             if(L->is_ll == true){
                y [0] /= d [0] ;
                y [1] = (y [1] - t * y [0]) / d [1] ;
             }
             else
                y [1] -= t * y [0] ;
             X [j  ] = y [0] ;
             X [j-1] = y [1] ;
          }
          k -= 2;
       }
       else
       {
          {
             double y [3], t [3] ;
             q = Lp [j-1] ;
             r = Lp [j-2] ;
             double d [3] ;
             d [0] = Lx [p] ;
             d [1] = Lx [q] ;
             d [2] = Lx [r] ;
             t [0] = Lx [q+1] ;
             t [1] = Lx [r+1] ;
             t [2] = Lx [r+2] ;
             if(L->is_ll == false){
                y [0] = X [j]   / d [0] ;
                y [1] = X [j-1] / d [1] ;
                y [2] = X [j-2] / d [2] ;
             }
             else{
                y [0] = X [j] ;
                y [1] = X [j-1] ;
                y [2] = X [j-2] ;
             }
             for (p++, q += 2, r += 3 ; p < pend ; p++, q++, r++)
             {
                int i = Li [p] ;
                y [0] -= Lx [p] * X [i] ;
                y [1] -= Lx [q] * X [i] ;
                y [2] -= Lx [r] * X [i] ;
             }
             if(L->is_ll == true){
                y [0] /= d [0] ;
                y [1] = (y [1] - t [0] * y [0]) / d [1] ;
                y [2] = (y [2] - t [2] * y [0] - t [1] * y [1]) / d [2] ;
             }
             else{
                y [1] -= t [0] * y [0] ;
                y [2] -= t [2] * y [0] + t [1] * y [1] ;
             }
             X [j-2] = y [2] ;
             X [j-1] = y [1] ;
             X [j  ] = y [0] ;
          }
          k -= 3;
       }
    }
}

void Block::delete_paths(){	
   delete [] s_col_FFS;
   delete [] s_col_FBS;
}
