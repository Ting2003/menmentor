// ----------------------------------------------------------------//
// Filename : parser.cpp
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// implementation file of parser.h
// ----------------------------------------------------------------//
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include "util.h"
#include "parser.h"
#include "mpi.h"
using namespace std;

// store the pointer to circuits
Parser::Parser(vector<Circuit*> * ckts):p_ckts(ckts){
}

Parser::~Parser(){ 
	// x_list.clear();
	// y_list.clear();
}

// node234_2_3_4 
// name_z_x_y
void Parser::extract_node(char * str, Node & nd, char * coord){
	//static Node gnd(string("0"), Point(-1,-1,-1));
	if( str[0] == '0' ) {
		nd.name="0";
		//nd.pt.set(-1,-1,-1);
		return;
	}

	long z, y, x;
	int flag = -1;
	char * chs;
	char * saveptr;
	char l[MAX_BUF];
	strcpy(l, coord);
	const char * sep = "(_ )";

	// node name
	chs = strtok_r(l, sep, &saveptr);
	string node_name;
	// node_name.append(chs);
	node_name.append(str);
	// for transient, 'Y' is the VDD source node
	/*if( chs[0] == 'X' || chs[0]== 'Y' || chs[0]== 'Z' ){
		flag = chs[0]-'X';
		chs = strtok_r(NULL, sep, &saveptr);
	}*/

	z = atol(chs);
	chs = strtok_r(NULL, sep, &saveptr);
	x = atol(chs);
	chs = strtok_r(NULL, sep, &saveptr);
	y = atol(chs);
	nd.name.assign(node_name);
	nd.pt.set(x, y, z);
	nd.flag = flag;
	nd.rid = -1;
}

// given a line, extract net and node information
void Parser::insert_net_node(char * line, int &my_id, MPI_CLASS &mpi_class){
	char *chs, *saveptr;
	const char* sep = " (),\n";

	static char sname[MAX_BUF];
	static char sa[MAX_BUF];
	static char sb[MAX_BUF];
	
	static char star[MAX_BUF];
	static char coord1[MAX_BUF];
	static char coord2[MAX_BUF];
	static char line_s[MAX_BUF];
	static char star_check[MAX_BUF];
	static Node nd[2];
	Node * nd_ptr[2];	// this will be set to the two nodes found
	double value;
	int count=0;
	
	char *chs_1;
	char *saveptr_1;
	const char *sep_1 = "_";
	bool pulse_flag = false;
	// find grid boundary x and y
	sscanf(line, "%s %s %s", sname,sa,sb);
	// if(my_id==1)
		// clog<<"block 1 line, sa, sb: "<<line<<" "<<sa<<" "<<sb<<endl;
	if(sa[0] == '0' || sb[0] == '0'){
		// clog<<"line: "<<line<<endl;
// #if 0
	  if(sname[0] == 'I' || sname[0] == 'i'){
		strcpy(line_s, line);
		// clog<<"current line: "<<line<<endl;
		chs = strtok_r(line_s, sep, &saveptr);
		strcpy(sname, chs);
		// clog<<"chs, sname: "<<chs<<" "<<sname<<endl;
		chs = strtok_r(NULL, sep, &saveptr);
		strcpy(sa, chs);
		chs = strtok_r(NULL, sep, &saveptr);
		strcpy(sb, chs);
		chs = strtok_r(NULL, sep, &saveptr);
		value = atof(chs);	
		// clog<<"sname, sa, sb, value: "<<sname<<" "<<sa<<" "<<sb<<" "<<value<<endl;
		// now skip the pulse current if any exists
		while(chs !=NULL){
			chs = strtok_r(NULL, sep, &saveptr);
			if(chs == NULL) break;
			// clog<<"chs: "<<chs<<endl;
			strcpy(star_check, chs);
			if(chs[0] == 'P' || chs[0] == 'p')
				pulse_flag = true;
			// read the coordinate
			if(chs[0] == '*'){
				chs = strtok_r(NULL, sep, &saveptr);
				if(chs == NULL) break;
				// clog<<"coord: chs: "<<chs<<endl;
				strcpy(coord1, chs);
				strcpy(coord2, coord1);
				break;
			}
		}
		// clog<<"finish one line. "<<endl;
	   }else{
		sscanf(line, "%s %s %s %lf %s %s", sname,sa,sb, &value, star, coord1);
		// copy string
		strcpy(coord2, coord1);
		// clog<<"coord1, coord2: "<<coord1<<" "<<coord2<<endl;
		// coord2 = coord1;
	  }
// #endif
	}
	else{
		// clog<<"regular line: "<<line<<endl; 
		sscanf(line, "%s %s %s %lf %s %s %s", sname,sa,sb, &value, star, coord1, coord2);
	}
// #if 0
	// net type
	chs_1 = strtok_r(sname, sep_1, &saveptr_1);
	// ckt name
	chs_1 = strtok_r(NULL, sep_1, &saveptr_1);
	string ckt_name(chs_1);
	// ckt_name.append(chs_1);
	
	extract_node(sa, nd[0], coord1);
	extract_node(sb, nd[1], coord2);

	int ckt_id = 0;
 	for(size_t i=0;i<(*p_ckts).size();i++){
		if((*p_ckts)[i]->name == ckt_name){
			ckt_id = i;
			break;
		}
	}
	Circuit * ckt = (*p_ckts)[ckt_id];

	// clog<<"ckt, line: "<<ckt->name<<" "<<line<<" "<<endl;	
	for(int i=0;i<2;i++){
		if ( nd[i].is_ground() ){
			nd_ptr[i] = ckt->nodelist[0]; // ground node
		}
		else if ( (nd_ptr[i] = ckt->get_node(nd[i].name) ) == NULL ){
			// create new node and insert
			nd_ptr[i] = new Node(nd[i]); // copy constructor
			
			nd_ptr[i]->ckt_name = ckt_name;
			nd_ptr[i]->rep = nd_ptr[i];  // set rep to be itself
			count = cpr_nd_block(nd_ptr[i], mpi_class.block_geo, my_id);
			if (count==1) // internal node
				ckt->add_node(nd_ptr[i]);
				// if(my_id==0)
					// clog<<"circuit add node: "<<*nd_ptr[i]<<endl;
			else{
				nd_ptr[i]->flag_bd = 1;
			}

			if( nd_ptr[i]->isS()==Y)	     // determine circuit type
				ckt->set_type(WB);
		}
	}

	// if there is a cap with it, need to modify
	if((line[0]=='r' || line[0] =='R') && 
		!(nd[0].pt.x == nd[1].pt.x && 
		nd[0].pt.y == nd[1].pt.y)){
		// add node into bd and internal vector
		add_node_inter(nd_ptr[0], nd_ptr[1], 
			mpi_class, ckt, my_id);
	}
	
	NET_TYPE net_type = RESISTOR;
	// find net type
	switch(sname[0]){
	case 'r': // resistor
	case 'R':
		net_type = RESISTOR;
		break;
	case 'v': // VDD
	case 'V':
		net_type = VOLTAGE;
		break;
	case 'i': // current
	case 'I':
		net_type = CURRENT;
		break;
	case 'c': // capacitance
	case 'C':
		net_type = CAPACITANCE;
		break;
	case 'l':
	case 'L':
		net_type = INDUCTANCE;
		break;
	default:
		report_exit("Invalid net type!\n");
		break;
	}

	
	// create a Net
	Net * net = new Net(net_type, value, nd_ptr[0], nd_ptr[1]);

	if(net_type == CURRENT && pulse_flag == true){
		// clog<<"to current net. "<<line<<endl;
		// clog<<"pulse_flag: "<<pulse_flag<<endl;
		net->tr = new double [7];
		// assign pulse paramter for pulse input
		chs = strtok_r(line, sep, &saveptr);
		for(int i=0;i<3;i++)
			chs = strtok_r(NULL, sep, &saveptr);
		if(chs != NULL){
			// clog<<"chs: "<<chs<<endl;
			chs = strtok_r(NULL, sep, &saveptr);
			}
		if(chs != NULL){
			// clog<<"chs again: "<<chs<<endl;
			chs = strtok_r(NULL, sep, &saveptr);
			// V1
			net->tr[0] = atof(chs);
			chs = strtok_r(NULL, sep, &saveptr);
			// V2
			net->tr[1] = atof(chs);
			chs = strtok_r(NULL, sep, &saveptr);
			// TD
			net->tr[2] = atof(chs);
			chs = strtok_r(NULL, sep, &saveptr);
			// Tr
			net->tr[3] = atof(chs);
			chs = strtok_r(NULL, sep, &saveptr);
			// Tf
			net->tr[4] = atof(chs);
			chs = strtok_r(NULL, sep, &saveptr);
			// PW
			net->tr[5] = atof(chs);
			chs = strtok_r(NULL, sep, &saveptr);
			// Period
			net->tr[6] = atof(chs);
		}
	}
	// trick: when the value of a resistor via is below a threshold,
	// treat it as a 0-voltage via
	//if( Circuit::MODE == (int)IT ) {
		// try_change_via(net);
	//}

	// insert this net into circuit
	ckt->add_net(net);
	// IMPORTANT: set the relationship between node and net
	// update node voltage if it is an X node
	// set node to be X node if it connects to a voltage source
	update_node(net);
}

// Given a net with its two nodes, update the connection information for thet two nodes
void Parser::update_node(Net * net){
	// first identify their connection type:
	// 1. horizontal/vertical   2. via/VDD 3. current
	//
	// swap *a and *b so that a is:
	// WEST   for horizontal
	// SOUTH  for vertical
	// BOTTOM for via / XVDD
	// ground node for CURRENT
	Node *a=net->ab[0], *b=net->ab[1];
	// assign isS() == Y
	if(net->type == VOLTAGE){
		if(a->is_ground()){
			b->flag== 1;
			b->value = net->value;
		}
		else{
			a->flag = 1;
			a->value = net->value;
		}
	}

	// clog<<"parsing net: "<<*net<<endl;
	if(a->is_ground()){
		b->nbr_vec.push_back(net);
		// clog<<"b add net: "<<*b<<endl;
	}
	else if(b->is_ground()){
		a->nbr_vec.push_back(net);
		// clog<<"a add net. "<<*a<<endl;
	}
	else {
		a->nbr_vec.push_back(net);
		b->nbr_vec.push_back(net);
		// clog<<"both a and b add net: "<<*a<<" "<<*b<<endl;
	}
	return;
	//cout<<"setting "<<net->name<<" nd1="<<nd1->name<<" nd2="<<nd2->name<<endl;	
}

// parse the file and create circuits
int Parser::create_circuits(vector<CKT_NAME> &ckt_name_vec){
	int n_circuit=0;

	string prev_ckt_name("");
	Circuit * p_last_circuit=NULL;
	// now read filename.info to create circuits (they are SORTED)
	for(size_t i=0;i<ckt_name_vec.size();i++){
		string name_string;
		for(int j=0;ckt_name_vec[i].name[j]!='\n';j++)
			name_string += ckt_name_vec[i].name[j];
		// name_string.append(ckt_name_vec[i].name);
		// compare with previous circuit name 
		//name_string.assign(name);
		if( prev_ckt_name == "" ||
		    name_string != prev_ckt_name ){
			Circuit * circuit = new Circuit(name_string);
			(*p_ckts).push_back(circuit);
			++n_circuit;
			prev_ckt_name = name_string;
			p_last_circuit = circuit;
		}
	}	
	return n_circuit;
}

// parse the file
// Note: the file will be parsed twice
// the first time is to find the network information
// and the second time is to create nodes
void Parser::parse(int &my_id, char * filename, MPI_CLASS &mpi_class, Tran &tran, int num_procs, bool partition_flag){	
	MPI_Datatype MPI_Vector;
	int count =1;
	int lengths = 10;
	MPI_Aint offsets = 0;
	MPI_Datatype types[2]={MPI_CHAR};
	MPI_Type_struct(count, &lengths, &offsets, types, 
			&MPI_Vector);
	int error = MPI_Type_commit(&MPI_Vector);

	// if(my_id==0) clog<<"create new type: "<<error<<endl;

	this->filename = filename;

	// processor 0 will extract layer info
	// and bcast it into other processor
	vector<CKT_NAME >ckt_name_vec;
	
	// clog<<"before extract ckt name: "<<endl;
	if(my_id==0){
		extract_ckt_name(my_id, ckt_name_vec, mpi_class, tran);
	}
	// broadcast info for transient 
	
	int ckt_name_size = ckt_name_vec.size();
	
	MPI_Bcast(&ckt_name_size, 1, MPI_INT, 0,
			MPI_COMM_WORLD);
	if(my_id!=0) ckt_name_vec.resize(ckt_name_size);

	// if(my_id==0) clog<<"before bcast mpi_vector"<<endl;
	MPI_Bcast(&ckt_name_vec[0], ckt_name_size, 
			MPI_Vector, 0, MPI_COMM_WORLD);

	// if(my_id==0) clog<<"after bcast mpi_vector"<<endl;
	/*clog<<my_id<<"==== "<<endl;
	for(size_t i=0;i<ckt_name_info.size();i++){
		clog<<ckt_name_info[i].name<<" "<<ckt_name_info[i].layer<<endl;
	}*/
	// first time parse:
	create_circuits(ckt_name_vec);
	if(my_id ==0){
		// clog<<"before pre partitioning. "<<endl;
		if(partition_flag == true){
			pre_partition(my_id, mpi_class, tran, num_procs);
		}
	}
	
	if(partition_flag)
		return;	
	if(my_id==0) clog<<"after pre partitioning."<<endl;

	build_block_geo(my_id, mpi_class, tran, num_procs);

	MPI_Barrier(MPI_COMM_WORLD);

	// temporary comment second parse	
	// if(my_id==0) clog<<"before second parse. "<<endl;
	second_parse(my_id, mpi_class, tran, num_procs);

	// if(my_id==0) clog<<"after second parse."<<endl;
}

void Parser::build_block_geo(int &my_id, MPI_CLASS &mpi_class, Tran &tran, int num_procs){
	// block info in cpu 0
	mpi_class.geo = new float[mpi_class.X_BLOCKS *mpi_class.Y_BLOCKS *4];
	// block info in other cpus
	mpi_class.block_geo = new float [4];

	// update block geometry
	if(my_id==0) set_block_geometry(mpi_class.geo, mpi_class);
	int total_blocks = 4 * mpi_class.X_BLOCKS * mpi_class.Y_BLOCKS;
	MPI_Bcast(mpi_class.geo, total_blocks, MPI_FLOAT, 0, MPI_COMM_WORLD);

	MPI_Bcast(&mpi_class.len_ovr_x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&mpi_class.len_ovr_y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// clog<<my_id<<" "<<mpi_class.len_ovr_x<<" "<<mpi_class.len_ovr_y<<endl;

	MPI_Scatter(mpi_class.geo, 4, MPI_FLOAT, 
		mpi_class.block_geo, 4, MPI_FLOAT, 
		0, MPI_COMM_WORLD);
	/*clog<<my_id<<" "<<mpi_class.block_geo[0]<<" "<<mpi_class.block_geo[1]<<" "<<
		mpi_class.block_geo[2]<<" "<<mpi_class.block_geo[3]<<endl;*/
	if(my_id==0){
		net_to_block(mpi_class.geo, mpi_class, tran, num_procs, my_id);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//if(my_id==3)
		//clog<<"lx, ly, ux, uy: "<<mpi_class.block_geo[0]<<" "<<mpi_class.block_geo[1]<<" "<<mpi_class.block_geo[2]<<" "<<mpi_class.block_geo[3]<<endl;
	// stores the geo boundary line for internal node:
	// e, w, n, s
	mpi_class.block_geo_origin = new float [4];
	
	mpi_class.set_geo_origin(mpi_class);
}

void Parser::second_parse(int &my_id, MPI_CLASS &mpi_class, Tran &tran, int num_procs){
	char  buff[100];
	FILE *f = NULL;

	int color = 0;
	int block_size = mpi_class.block_size;
	//if(block_size==0) return;
	//for(int i=0;i<block_size;i++){
		sprintf(buff, "./INPUT_FILE/netlist_%d%d.txt", color, my_id);
		f = fopen(buff, "r");
		if(f==NULL) report_exit("Input file not exist!\n");
	//}
	
	char line[MAX_BUF];
	char type;

	tran.nodes.clear();
	while( fgets(line, MAX_BUF,f)!=NULL){
		type = line[0];
		switch(type){
			case 'r': // resistor
			case 'R':
			case 'v': // VDD
			case 'V':
			case 'i': // current
			case 'I':
			case 'c':
			case 'C':
			case 'l':
			case 'L':
				insert_net_node(line, my_id, mpi_class);
				break;
			case '.': // parse tran nodes: need to write
				 block_parse_dots(line, tran, my_id); 
				 break;
			case '*': // comment
			case ' ':
			case '\n':
				break;
			default:
				printf("Unknown input line: ");
				report_exit(line);
				break;
		}
	}	

	fclose(f);
	// release map_node resource
	/*for(size_t i=0;i<(*p_ckts).size();i++){
		Circuit * ckt = (*p_ckts)[i];
		if(ckt->map_node.size()>0)
			ckt->map_node.clear();
	}*/
	
}// end of parse

// done by core 0
void Parser::parse_dot(char *line, Tran &tran){
	char *chs;
	char *saveptr;
	char sname[MAX_BUF];
	Node_TR_PRINT item;
	const char *sep = "= v() \n";
	switch(line[1]){
		case 't': // transient steps
			/*sscanf(line, "%s %lf %lf", sname, 
				&tran.step_t, &tran.tot_t);
			tran.isTran = 1; // do transient ana;*/
			//clog<<"step: "<<tran.step_t<<" tot: "<<tran.tot_t<<endl;
			break;
		case 'w': // output length
			/*chs = strtok_r(line, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			tran.length = atoi(chs);*/
			//clog<<"out len: "<<tran.length<<endl;
			break;
		case 'p': // print
			chs = strtok_r(line, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			while(chs != NULL){
				chs = strtok_r(NULL, sep, &saveptr);
				if(chs == NULL) break;
				item.name = chs;
				tran.nodes.push_back(item);
				// disribute nodes into cores
			};
			break;
		default: 
			break;
	}
}

int Parser::extract_ckt_name(int &my_id, 
		vector<CKT_NAME > &ckt_name_vec,
		MPI_CLASS &mpi_class,
		Tran &tran){
	char line[MAX_BUF];
	char word[MAX_BUF];
	string word_s;
	char name[10];
	static char sname[MAX_BUF];
	static char sa[MAX_BUF];
	static char sb[MAX_BUF];
	static char line_s[MAX_BUF];
	
	static char star[MAX_BUF];
	static char coord1[MAX_BUF];
	static char coord2[MAX_BUF];
	static char star_check[MAX_BUF];
	static Node nd[2];
	double value;
	int i=0;
	long x_max=-1;
	long x_min=-1;
	long y_max=-1;
	long y_min=-1;

	char *chs;
	char *saveptr;
	const char *sep = " (,)";
	// only processor 0 will extract layer info
	if(my_id!=0) return 0;

	FILE *f;
	f = fopen(filename, "r");
	if(f==NULL) report_exit("Input file not exist!\n");

	CKT_NAME ckt_name;

	// skip the first comment line
	fgets(line, MAX_BUF, f);

	while(fgets(line, MAX_BUF, f)!=NULL){
		//clog<<"line: "<<line<<endl;
		if(line[0]=='*'){
			// copy the entire line into stringstream
			stringstream ss;
			ss<< line;
			// if(my_id==0) clog<<"line: "<<line<<endl;
			int word_count = 0;
			while(ss.getline(word, 10, ' ')){
				if(word_count ==1 && (word[0] == 'v' || word[0] == 'V' ||
					word[0] == 'G' || word[0] == 'g')){
					//ckt_name(word);
					strcpy(ckt_name.name, word);
					ckt_name_vec.push_back(ckt_name);
					// clog<<"ckt_name: "<<ckt_name.name<<endl;	
				}
				if(word_count >=2) break;	
				word_count++;
			}
		}
		else if(line[0]!='.'){
			// find grid boundary x and y
			sscanf(line, "%s %s %s", sname,sa,sb);
			if(sa[0] == '0' || sb[0] == '0'){
				// clog<<"line: "<<line<<endl;
				// current nets
			  if(sname[0] == 'I' || sname[0] == 'i'){
				strcpy(line_s, line);
				// clog<<"current line: "<<line<<endl;
				chs = strtok_r(line_s, sep, &saveptr);
				strcpy(sname, chs);
				// clog<<"chs, sname: "<<chs<<" "<<sname<<endl;
				chs = strtok_r(NULL, sep, &saveptr);
				strcpy(sa, chs);
				chs = strtok_r(NULL, sep, &saveptr);
				strcpy(sb, chs);
				chs = strtok_r(NULL, sep, &saveptr);
				value = atof(chs);	
				// clog<<"sname, sa, sb, value: "<<sname<<" "<<sa<<" "<<sb<<" "<<value<<endl;
				// now skip the pulse current if any exists
				while(chs !=NULL){
					chs = strtok_r(NULL, sep, &saveptr);
					if(chs == NULL) break;
					// clog<<"chs: "<<chs<<endl;
					strcpy(star_check, chs);
					// read the coordinate
					if(chs[0] == '*'){
						chs = strtok_r(NULL, sep, &saveptr);
						if(chs == NULL) break;
						// clog<<"coord: chs: "<<chs<<endl;
						strcpy(coord1, chs);
						strcpy(coord2, coord1);
						break;
					}
				}
			  }else{
				sscanf(line, "%s %s %s %lf %s %s", sname,sa,sb, &value, 
					star, coord1);
				// copy string
				strcpy(coord2, coord1);
			  }
			}
			else 
				sscanf(line, "%s %s %s %lf %s %s %s", sname,sa,sb, &value, 
					star, coord1, coord2);
			// clog<<"line: "<<line<<endl;	
			//if( sa[0] != '0' )
				extract_node(sa, nd[0], coord1);

			//if( sb[0] != '0' ) 
				extract_node(sb, nd[1], coord2);
			for(i=0;i<2;i++){
				if(nd[i].pt.x > x_max) 
					x_max = nd[i].pt.x;
				if(x_min == -1) x_min = x_max;
				else if(nd[i].pt.x>0 && nd[i].pt.x <x_min)
					x_min = nd[i].pt.x;
				if(nd[i].pt.y > y_max)
					y_max = nd[i].pt.y;
				if(y_min == -1) y_min = y_max;
				else if(nd[i].pt.y>0 && nd[i].pt.y <y_min)
					y_min = nd[i].pt.y;
			}
		}else if(line[0] == '.'){
			// parse_dot by core 0
			parse_dot(line, tran);
		}
	}
	clog<<"finish extract coords from lines. "<<endl;
	mpi_class.x_max = x_max;
	mpi_class.x_min = x_min;
	mpi_class.y_max = y_max;
	mpi_class.y_min = y_min;
	clog<<"power grid bd: "<<x_min<<" "<<y_min<<" "<<x_max<<" "<<y_max<<endl;
	fclose(f);
	// sort resulting vector by the ckt name
	sort(ckt_name_vec);
	return 0;
}
// sort ckt according to its name and layer
// ckt name decrease, layer rising
bool Parser::sort(vector<CKT_NAME > &a){
	CKT_NAME tmp;
	// sort according to circuit name
	for(size_t i=0;i< a.size()-1;i++){
		int minIndex = i;
		for(size_t j = i+1;j< a.size();j++)
			if(strcmp(a[j].name, a[minIndex].name)>0)
				minIndex = j;
		if(minIndex !=i){
			tmp = a[i];
			a[i] = a[minIndex];
			a[minIndex] = tmp;
		}	
	}
	return true;
}

void Parser::set_block_geometry(float *geo, MPI_CLASS &mpi_class){
	double x, y;
	x = (double)(mpi_class.x_max-mpi_class.x_min+0.5)*1.0 
		/ mpi_class.X_BLOCKS;
	y = (double)(mpi_class.y_max-mpi_class.y_min+0.5)*1.0
		/ mpi_class.Y_BLOCKS;
	//if( fzero(x) ) x = 1.0;
	//if( fzero(y) ) y = 1.0;
	double len_per_block_x = x;
	double len_per_block_y = y;
	double len_ovr_x = x * mpi_class.overlap_ratio;
	double len_ovr_y = y * mpi_class.overlap_ratio;
	mpi_class.len_per_block_x = x;
	mpi_class.len_per_block_y = y;
	mpi_class.len_ovr_x = len_ovr_x;
	mpi_class.len_ovr_y = len_ovr_y;
	clog<<"len_per_x, len_per_y: "<<x<<" "<<y<<endl;
	clog<<"len_ovr_x, len_ovr_y: "<<len_ovr_x
		<<" / "<<len_ovr_y<<endl;

	size_t bid = 0;
	// update block 4 corners
	for(size_t y=0;y<mpi_class.Y_BLOCKS;y++){
		for(size_t x=0;x<mpi_class.X_BLOCKS;x++){
			bid = y * mpi_class.X_BLOCKS + x;
			// lx
			geo[4*bid] = mpi_class.x_min + x * len_per_block_x - len_ovr_x;
			// ly
			geo[4*bid+1] = mpi_class.y_min + y * len_per_block_y - len_ovr_y;
			// ux
			geo[4*bid+2] = mpi_class.x_min + (x+1) * len_per_block_x + len_ovr_x;
			// uy
			geo[4*bid+3] = mpi_class.y_min + (y+1) * len_per_block_y + len_ovr_y;
		}
	}
}

// done by processor 0
void Parser::net_to_block(float *geo, MPI_CLASS &mpi_class, Tran &tran, int num_procs, int &my_id){
	static char sname[MAX_BUF];
	static char sa[MAX_BUF];
	static char sb[MAX_BUF];
	static Node nd[2];
	double value;	
	
	int temp = 0;

	FILE *f;
	f = fopen(this->filename, "r");
	if(f==NULL) report_exit("Input file not exist!\n");

	// handles the long print sentence
	// int temp_buf = 10000000;
	char line[MAX_BUF];
	
	static char star[MAX_BUF];
	static char coord1[MAX_BUF];
	static char coord2[MAX_BUF];
	static char star_check[MAX_BUF];
	static char line_s[MAX_BUF];
	
	char *chs;
	char *saveptr;
	const char *sep = " (,)";

	int color = 0;
	vector<FILE *> of;
	int num_blocks  = mpi_class.X_BLOCKS * mpi_class.Y_BLOCKS;
	clog<<"num_blocks. "<<num_blocks<<endl;
	InitialOF(of, num_procs, color);//num_blocks, color);
	// clog<<"after initial ofs. "<<endl;

	int count_1 = 0, count_2 = 0;
	while( fgets(line, MAX_BUF, f)!=NULL){
		if(line[0]=='r' || line[0] =='R' ||
		   line[0]=='v' || line[0] =='V' ||
		   line[0]=='i' || line[0]=='I' ||
		   line[0]=='c' || line[0] == 'C' ||
		   line[0]=='l' || line[0] == 'L'){
			sscanf(line, "%s %s %s", sname,sa,sb);
			if(sa[0] == '0' || sb[0] == '0'){
				if(sname[0] == 'I' || sname[0] == 'i'){
				  // clog<<endl<<"current line: "<<line;
				  strcpy(line_s, line);
				  chs = strtok_r(line_s, sep, &saveptr); 
				  // clog<<"line again: "<<line<<endl;	
				  strcpy(sname, chs);
				  // clog<<"chs, sname: "<<chs<<" "<<sname<<endl;
				  chs = strtok_r(NULL, sep, &saveptr);
				  strcpy(sa, chs);
				  chs = strtok_r(NULL, sep, &saveptr);
				  strcpy(sb, chs);
				  chs = strtok_r(NULL, sep, &saveptr);
				  value = atof(chs);	
				  // clog<<"sname, sa, sb, value: "<<sname<<" "<<sa<<" "<<sb<<" "<<value<<endl;
				
				  // now skip the pulse current if any exists
				  while(chs !=NULL){
					chs = strtok_r(NULL, sep, &saveptr);
					if(chs == NULL) break;
					// clog<<"chs: "<<chs<<endl;
					// strcpy(star_check, chs);
					// read the coordinate
					if(chs[0] == '*'){
						chs = strtok_r(NULL, sep, &saveptr);
						if(chs == NULL) break;
						// clog<<"coord: chs: "<<chs<<endl;
						strcpy(coord1, chs);
						strcpy(coord2, coord1);
						break;
					}
				  }
				}else{	  
				  sscanf(line, "%s %s %s %lf %s %s", sname,sa,sb, &value, 
					star, coord1);
				  // copy string
				  strcpy(coord2, coord1);
				}
			}
			else 
				sscanf(line, "%s %s %s %lf %s %s %s", sname,sa,sb, &value, 
					star, coord1, coord2);

			extract_node(sa, nd[0], coord1);
			extract_node(sb, nd[1], coord2);
			for(int i=0;i<num_blocks;i++){
				// at least one node is inside block
				count_1 = cpr_nd_block(nd[0], geo, i);
				count_2 = cpr_nd_block(nd[1], geo, i);

				/*if(my_id==0){//&&nd[0].name == "n1_2024_174" && nd[1].name == "n1_2072_174"){
					clog<<"line: "<<line<<endl;
					clog<<count_1<<" "<<count_2<<endl;
				}*/
				if(nd[0].is_ground() && count_2 ==1){
					temp = fprintf(of[i], "%s", line);
					// clog<<"print line: "<<line<<endl;
				}
				else if(nd[1].is_ground() && count_1 ==1){
					temp = fprintf(of[i], "%s", line);
					// clog<<"print line: "<<line<<endl;
				}
				// write all voltage sources
				//else if(count_1 + count_2 ==2){
				else if(!nd[0].is_ground() && !nd[1].is_ground() 
						&& (count_1 + count_2 >=1)){
					temp = fprintf(of[i], "%s", line);
					//if(my_id==0)
					// clog<<"write: "<<temp<<" "<<i<<" "<<line<<endl;
				}
			}
		}
		else{
			// clog<<"special lines. "<<line<<endl;
			for(int i=0;i<num_procs;i++){
				fprintf(of[i], "%s", line);
			}
		}
	}
	// free(line);
	// finally print end file symbol (not need to)
	fclose(f);
	for(int i=0;i<num_procs;i++){
		fclose(of[i]);
	}
	of.clear();
}

int Parser::cpr_nd_block(Node &nd, float *geo, int &bid){
	if(nd.pt.x >= geo[4*bid] &&
	   nd.pt.x <= geo[4*bid+2] &&
	   nd.pt.y >= geo[4*bid+1] &&
	   nd.pt.y <= geo[4*bid+3]){
		return 1;
	}
	else return 0;
}

int Parser::cpr_nd_block(Node *nd, float *geo, int &bid){
	if(nd->pt.x >= geo[0] &&
	   nd->pt.x <= geo[2] &&
	   nd->pt.y >= geo[1] &&
	   nd->pt.y <= geo[3]){
		return 1;
	}
	else return 0;
}

int Parser::cpr_nd_block(Node *nd, float &lx, float &ly, float &ux, float &uy){
	if(nd->pt.x >= lx&&
	   nd->pt.x <= ux &&
	   nd->pt.y >= ly &&
	   nd->pt.y <= uy){
		return 1;
	}
	else {
		return 0;
	}
}

void Parser::add_node_inter(Node *nd_0, Node *nd_1, 
	MPI_CLASS &mpi_class, Circuit *ckt, int &my_id){
	int bx, by;
	int bid_nbr;
	float lx, ly, ux, uy;
	float lx_0, ly_0, ux_0, uy_0;
	int count_1, count_2;
	int count_10, count_20;

	by = my_id / mpi_class.X_BLOCKS;
	bx = my_id % mpi_class.X_BLOCKS;

	find_bound_line(my_id, mpi_class, lx_0, ly_0,
			ux_0, uy_0);
	
	count_10 = cpr_nd_block(nd_0, lx_0, ly_0, ux_0, uy_0);
	count_20 = cpr_nd_block(nd_1, lx_0, ly_0, ux_0, uy_0);
	
	// sw block
	if((by>=1 && bx>=1)){
		bid_nbr = my_id - mpi_class.X_BLOCKS - 1;
		
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_sw, ckt->internal_nodelist_sw, 
		nd_0, nd_1, count_10, count_20);
	}

	// s block
	if(by>=1){
		bid_nbr = my_id - mpi_class.X_BLOCKS;
		
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_s, ckt->internal_nodelist_s, 
		nd_0, nd_1, count_10, count_20);
	}
	
	// se block
	if((by>=1 && bx<mpi_class.X_BLOCKS-1)){
		bid_nbr = my_id - mpi_class.X_BLOCKS + 1;
	
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_se, ckt->internal_nodelist_se, 
		nd_0, nd_1, count_10, count_20);
	}

	// w block
	if(bx>=1){
		bid_nbr = my_id - 1;
		
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_w, ckt->internal_nodelist_w, 
		nd_0, nd_1, count_10, count_20);	
	}

	// e block
	if((bx<mpi_class.X_BLOCKS-1)){
		bid_nbr = my_id + 1;
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_e, ckt->internal_nodelist_e, 
		nd_0, nd_1, count_10, count_20);
	}

	// nw block
	if((by<mpi_class.Y_BLOCKS-1 && bx>=1)){
		bid_nbr = my_id + mpi_class.X_BLOCKS -1;
		
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_nw, ckt->internal_nodelist_nw, 
		nd_0, nd_1, count_10, count_20);
	}

	// n block
	if((by<mpi_class.Y_BLOCKS-1)){
		bid_nbr = my_id + mpi_class.X_BLOCKS;
		
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_n, ckt->internal_nodelist_n, 
		nd_0, nd_1, count_10, count_20);
	}

	// ne block
	if((by<mpi_class.Y_BLOCKS-1 && bx<mpi_class.X_BLOCKS-1)){
		bid_nbr = my_id + mpi_class.X_BLOCKS + 1;
		
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_ne, ckt->internal_nodelist_ne, 
		nd_0, nd_1, count_10, count_20);
	}
}

void Parser::InitialOF(vector<FILE *> & of, int &num_blocks, int &color){
	char  buff[100];
	FILE *f;
	for(int i=0;i<num_blocks;i++){
		sprintf(buff, "./INPUT_FILE/netlist_%d%d.txt", color, i);
		f = fopen(buff, "w");
		of.push_back(f);
	}
}

void Parser::insert_node_list(Node *nd_0, Node *nd_1, int &count_10, 
		int &count_20, int &count_1, int &count_2, NodePtrVector &list, bool &flag){
	// add into bd boundary list
	if(count_10 ==1 && count_20 ==0 && count_2 ==1){
		list.push_back(nd_1);
		if(flag == true)
			nd_1->internal_bd = 1;
	}
	else if(count_10 ==0 && count_20 ==1 && count_1 ==1){
		list.push_back(nd_0);
		if(flag == true)
			nd_0->internal_bd = 1;
	}
}

void Parser::find_bound_line(int &bid_nbr, MPI_CLASS &mpi_class, float &lx, float &ly,  float &ux, float &uy){
	lx = mpi_class.geo[4*bid_nbr];
	ux = mpi_class.geo[4*bid_nbr+2];
	ly = mpi_class.geo[4*bid_nbr+1];
	uy = mpi_class.geo[4*bid_nbr+3];
}

void Parser::insert_node_dir(int &bid_nbr, MPI_CLASS &mpi_class, NodePtrVector &bd_list, NodePtrVector &inter_list, Node *nd_0, Node *nd_1, int &count_10, int &count_20){
	float lx, ly, ux, uy;
	int count_1, count_2;
	bool flag = true;
	
	find_bound_line(bid_nbr, mpi_class, lx, ly,
			ux, uy);

	count_1 = cpr_nd_block(nd_0, lx, ly, ux, uy);
	count_2 = cpr_nd_block(nd_1, lx, ly, ux, uy);

	insert_node_list(nd_0, nd_1, count_1, 
			count_2, count_10, count_20, 
			inter_list, flag);
	flag = false;

	insert_node_list(nd_0, nd_1, count_10, 
			count_20, count_1, count_2, 
			bd_list, flag);
}

//every core parses in data simultaneously
// need to extract its corresponding transient nodes
// use the find 
void Parser::block_parse_dots(char *line, Tran &tran, int &my_id){
	char *chs;
	char *saveptr;
	char sname[MAX_BUF];
	Node_TR_PRINT item;
	map<string, Node*>::iterator it;
	// clear tran.nodes, especially for core 0
	const char *sep = "= v() \n";
	switch(line[1]){
		case 't': // transient steps
			sscanf(line, "%s %lf %lf", sname, 
				&tran.step_t, &tran.tot_t);
			tran.isTran = 1; // do transient ana;
			//clog<<"step: "<<tran.step_t<<" tot: "<<tran.tot_t<<endl;
			break;
		case 'w': // output length
			chs = strtok_r(line, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			tran.length = atoi(chs);
			//clog<<"out len: "<<tran.length<<endl;
			break;
		case 'p': // print
			// clog<<"line: "<<line<<endl;
			Node *nd_ptr;	
			chs = strtok_r(line, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			while(chs != NULL){
				chs = strtok_r(NULL, sep, &saveptr);
				if(chs == NULL) break;
				item.name = chs;
				// clog<<"name: "<<item.name<<endl;
				for(size_t i=0;i<(*p_ckts).size();i++){
					Circuit *ckt = (*p_ckts)[i];
					it = ckt->map_node.find(item.name);
					if(it!=ckt->map_node.end()){
						// nd_ptr = it->second;
						tran.nodes.push_back(item);
						// if(my_id==0)
						//clog<<"node: "<<*nd_ptr<<" "<<nd_ptr->pt<<endl;
						break;
					}
				}
				// disribute nodes into cores
			};
			break;
		default: 
			break;
	}
}

void Parser::write_print(Tran &tran, vector<FILE *> &of, MPI_CLASS &mpi_class, char *line){
// then write tran nodes into files
	char *chs;
	char *saveptr;
	char *saveptr_1;
	const char *sep = " _";
	const char *sep_1 = " v()V";
	string ndname;
	vector<string> name_vec;

	chs = strtok_r(line, sep_1, &saveptr_1);

	int num_blocks  = mpi_class.X_BLOCKS * mpi_class.Y_BLOCKS;
	// first print .tran to each file
	for(int j=0;j<num_blocks;j++)
		fprintf(of[j], "%s ", chs);
	chs = strtok_r(NULL, sep_1, &saveptr_1);
	for(int j=0;j<num_blocks;j++)
		fprintf(of[j], "%s ", chs);

	chs = strtok_r(NULL, sep_1, &saveptr_1);
	while(chs != NULL){
		// clog<<"chs: "<<chs<<endl;
		if(string(chs) != "\n")
			name_vec.push_back(string(chs));
		chs = strtok_r(NULL, sep_1, &saveptr_1);
	}
	// clog<<"name_vec, size(): "<<name_vec.size()<<endl;
	// then print the nodes
	for(size_t i=0;i<name_vec.size();i++){
		string tr_nd_name = name_vec[i];
		//clog<<"name: "<<tr_nd_name<<endl;
		char *p = &tr_nd_name[0];
		chs = strtok_r(p, sep, &saveptr);
		ndname = string(chs);
		// skip the X or Y part
		if(ndname[0]=='X' || ndname[0]=='Y') 
			chs = strtok_r(NULL, sep, &saveptr);
		// skip the nz
		chs = strtok_r(NULL, sep, &saveptr);
		// extract x coordinate
		double x = atoi(chs);
		chs = strtok_r(NULL, sep, &saveptr);
		// extract y coordinate
		double y = atoi(chs);

		// judge which blocks this (x,y) belongs to
		for(int j=0;j<num_blocks;j++){
			if(x >= mpi_class.geo[4*j] && x<= mpi_class.geo[4*j+2]){
				if(y >= mpi_class.geo[4*j+1]&& y <= mpi_class.geo[4*j+3]){
					//clog<<"name again: "<<name_vec[i]<<endl;
					//clog<<"belongs to block: "<<j<<" "<<x<<" "<<mpi_class.geo[4*j]<<" "<<mpi_class.geo[4*j+2]<<" "<<y<<" "<<mpi_class.geo[4*j+1]<<" "<<mpi_class.geo[4*j+3]<<endl;
					fprintf(of[j], "v(%s) ", name_vec[i].c_str());	
				}
			}
		}
	}
	for(int j=0;j<num_blocks;j++){
		fprintf(of[j], "\n");
	}
	name_vec.clear();
}
#if 0
// only fit for the same layer net
// map net into horizontal or vertical dir
bool Parser::map_res_net(Net*net){
	Node *a = net->ab[0];
	Node *b = net->ab[1];
	if(a->is_ground() || b->is_ground())
		return true;
	long xa = a->pt.x;
	long ya = a->pt.y;
		
	long xb = b->pt.x;
	long yb = b->pt.y;

	long delta_x = abs(xa-xb);
	long delta_y = abs(ya-yb);

	double tan = 1.0* delta_y / delta_x;
	// 45 degree
	double ref = 1.0* sqrt(2) / 2;
	// > 45 degree, goes to y then x dir
	if(tan - ref > 1e-5){
		// clog<<" > 45 degree. "<<endl;
		bool flag_y = map_net_y(net);
		if(flag_y == true) return true;
		bool flag_x = map_net_x(net);
		if(flag_x == true) return true;	
		// clog<<"not empty for both south and north. "<<endl;
	}
	// <= 45 degree, goes to x dir
	else{
		// clog<<" < 45 degree. "<<endl;
		bool flag_x = map_net_x(net);
		if(flag_x == true) 
			return true;
		bool flag_y = map_net_y(net);
		if(flag_y == true) 
			return true;	
		// clog<<"not empty for both south and north. "<<endl;	
	}
	clog<<"no spot to map the net. Fault. "<<endl;
	return false;
}

// project along y direction
bool Parser::map_net_y(Net *net){
	Node *a = net->ab[0];
	Node *b = net->ab[1];
	long xa = a->pt.x;
	long xb = b->pt.x;
	long ya = a->pt.y;
	long yb = b->pt.y;
	
	if(ya < yb && b->nbr[SOUTH] == NULL && a->nbr[NORTH] == NULL){
		b->set_nbr(SOUTH, net);
		a->set_nbr(NORTH, net);
		clog<<"SOUTH of: "<<*b<<" NORTH of "<<*a<<endl;
		return true;
	}
	else if(yb < ya && b->nbr[NORTH] == NULL && a->nbr[SOUTH] == NULL){
		b->set_nbr(NORTH, net);
		a->set_nbr(SOUTH, net);

		clog<<"NORTH of: "<<*b<<" SOUTH of "<<*a<<endl;
		return true;
	}
	return false;
}

// project along x direction
bool Parser::map_net_x(Net *net){
	Node *a = net->ab[0];
	Node *b = net->ab[1];
	long xa = a->pt.x;
	long xb = b->pt.x;
	long ya = a->pt.y;
	long yb = b->pt.y;
	
	if(xa < xb && a->nbr[EAST] == NULL && b->nbr[WEST] == NULL){
		a->set_nbr(EAST, net);
		b->set_nbr(WEST, net);

		clog<<"WEST of: "<<*b<<" EAST of "<<*a<<endl;
		return true;

	}
	else if(xb < xa && b->nbr[EAST] == NULL && a->nbr[WEST] == NULL){
		b->set_nbr(EAST, net);
		a->set_nbr(WEST, net);

		clog<<"WEST of: "<<*a<<" EAST of "<<*b<<endl;
		return true;

	}
}
#endif
 
void Parser::pre_partition(int my_id, MPI_CLASS &mpi_class, Tran &tran, int num_procs){
	char line[MAX_BUF];
	// processing original input file
	FILE *f = NULL;
	f = fopen(filename, "r");
	if(f==NULL) report_exit("Input file not exist!\n");	
	for(int i=0;i<(*p_ckts).size();i++){
		Circuit *ckt = (*p_ckts)[i];
		ckt->lx = -1;
		ckt->ux = -1;
		ckt->ly = -1;
		ckt->uy = -1; 
	}

	char type;
	while( fgets(line, MAX_BUF,f)!=NULL){
		type = line[0];
		//clog<<line<<endl;
		switch(type){
			case 'r': // resistor
			case 'R':
			case 'v': // VDD
			case 'V':
			case 'i': // current
			case 'I':
			case 'c':
			case 'C':
			case 'l':
			case 'L':
				pre_insert_net_node(line, my_id, mpi_class);
				break;
			case '.': // parse tran nodes: need to write
				 // block_parse_dots(line, tran, my_id); 
				 // break;
			case '*': // comment
			case ' ':
			case '\n':
				break;
			default:
				printf("Unknown input line: ");
				report_exit(line);
				break;
		}
	}	

	fclose(f);

	// release map_node resource
	for(size_t i=0;i<(*p_ckts).size();i++){
		Circuit * ckt = (*p_ckts)[i];
		if(ckt->map_node.size()>0)
			ckt->map_node.clear();	
	}
	build_x_y_list_map();	
	// then explore the partition with no overlap
	explore_partition(num_procs, mpi_class);
	if(my_id==0) clog<<"after explore_partition. "<<endl; 
	// now clean all the resources for exploring partitioning
	//clean_explore_partition(num_procs, mpi_class);
	// if(my_id==0) clog<<"after clean explore partition. "<<endl;
	// clog<<"after release ckt: "<<(*p_ckts).size()<<endl;	
}// end of parse

// given a line, extract net and node information
void Parser::pre_insert_net_node(char * line, int &my_id, MPI_CLASS &mpi_class){
	char *chs, *saveptr;
	const char* sep = " (),\n";

	static char sname[MAX_BUF];
	static char sa[MAX_BUF];
	static char sb[MAX_BUF];
	
	static char star[MAX_BUF];
	static char coord1[MAX_BUF];
	static char coord2[MAX_BUF];
	static Node nd[2];
	Node * nd_ptr[2];	// this will be set to the two nodes found
	double value;
	int count=0;
	
	char *chs_1;
	char *saveptr_1;
	const char *sep_1 = "_";
	// find grid boundary x and y
	sscanf(line, "%s %s %s", sname,sa,sb);
	// if(my_id==1)
		// clog<<"block 1 line, sa, sb: "<<line<<" "<<sa<<" "<<sb<<endl;
	if(sa[0] == '0' || sb[0] == '0'){
	// if(sa == "0" || sb == "0"){
		sscanf(line, "%s %s %s %lf %s %s", sname,sa,sb, &value, star, coord1);
		// copy string
		strcpy(coord2, coord1);
		// clog<<"coord1, coord2: "<<coord1<<" "<<coord2<<endl;
		// coord2 = coord1;
	}
	else 
		sscanf(line, "%s %s %s %lf %s %s %s", sname,sa,sb, &value, star, coord1, coord2);

	// net type
	chs_1 = strtok_r(sname, sep_1, &saveptr_1);
	// ckt name
	chs_1 = strtok_r(NULL, sep_1, &saveptr_1);
	string ckt_name(chs_1);
	// ckt_name.append(chs_1);
	
	extract_node(sa, nd[0], coord1);
	extract_node(sb, nd[1], coord2);
	//clog<<"after extract nodes. "<<endl;

	int ckt_id = 0;
 	for(size_t i=0;i<(*p_ckts).size();i++){
		if((*p_ckts)[i]->name == ckt_name){
			ckt_id = i;
			break;
		}
	}
	Circuit * ckt = (*p_ckts)[ckt_id];
	// cout<<endl<<"ckt, line: "<<ckt->name<<" "<<line;	
	for(int i=0;i<2;i++){
		if ( nd[i].is_ground() ){
			nd_ptr[i] = ckt->nodelist[0]; // ground node
		}
		else if ( (nd_ptr[i] = ckt->get_node(nd[i].name) ) == NULL ){
			// create new node and insert
			nd_ptr[i] = new Node(nd[i]); // copy constructor
			// extract the boundary for circuit
			if(nd_ptr[i]->pt.x > ckt->ux)
				ckt->ux = nd_ptr[i]->pt.x;
			if(ckt->lx == -1) ckt->lx = ckt->ux;
			else if(nd_ptr[i]->pt.x>0 && nd_ptr[i]->pt.x <ckt->lx)
				ckt->lx = nd_ptr[i]->pt.x;
			if(nd_ptr[i]->pt.y > ckt->uy)
				ckt->uy = nd_ptr[i]->pt.y;
			if(ckt->ly == -1) ckt->ly = ckt->uy;
			else if(nd_ptr[i]->pt.y>0 && nd_ptr[i]->pt.y <ckt->ly)
				ckt->ly = nd_ptr[i]->pt.y;

			
			nd_ptr[i]->ckt_name = ckt_name;
			nd_ptr[i]->rep = nd_ptr[i];  // set rep to be itself
			ckt->add_node(nd_ptr[i]);
			// cout<<"ckt, add_node: "<<ckt->get_name()<<" "<<*nd_ptr[i]<<endl;
		}
	}	
	
	NET_TYPE net_type = RESISTOR;
	// find net type
	switch(sname[0]){
	case 'r': // resistor
	case 'R':
		net_type = RESISTOR;
		break;
	case 'v': // VDD
	case 'V':
		net_type = VOLTAGE;
		break;
	case 'i': // current
	case 'I':
		net_type = CURRENT;
		break;
	case 'c': // capacitance
	case 'C':
		net_type = CAPACITANCE;
		break;
	case 'l':
	case 'L':
		net_type = INDUCTANCE;
		break;
	default:
		report_exit("Invalid net type!\n");
		break;
	}
	
	// create a Net
	Net * net = new Net(net_type, value, nd_ptr[0], nd_ptr[1]);

	// insert this net into circuit
	ckt->add_net(net);
}

// build the statistic info for x and y dir
void Parser::build_x_y_list_map(){
	Node *nd;
	// all the ckts together
	// now have the nodes and nets info
	for(size_t i=0;i<(*p_ckts).size();i++){
		Circuit *ckt = (*p_ckts)[i];
		for(size_t j=0;j<ckt->nodelist.size();j++){
			nd = ckt->nodelist[j];
			if(nd->is_ground())
				continue;
			ckt->x_list.push_back(nd->pt.x);
			ckt->y_list.push_back(nd->pt.y);
			// cout<<"i, node, x, y: "<<i<<" "<<*nd<<" "<<nd->pt<<" "<<nd->pt.x<<" "<<nd->pt.y<<endl;
		}
		// build ckt->x_list_nd_map
		for(size_t j=0;j<ckt->x_list.size();j++){
			ckt->x_list_nd_map[ckt->x_list[j]] ++;	
		}
		// build ckt->x_list_nd_map
		for(size_t j=0;j<ckt->y_list.size();j++){
			ckt->y_list_nd_map[ckt->y_list[j]] ++;	
		}

		std::sort(ckt->x_list.begin(), ckt->x_list.end());
		std::sort(ckt->y_list.begin(), ckt->y_list.end());
		vector<long>::iterator it_vec;
		it_vec = std::unique(ckt->x_list.begin(), ckt->x_list.end());
		ckt->x_list.resize(std::distance(ckt->x_list.begin(), it_vec));
	
		it_vec = std::unique(ckt->y_list.begin(), ckt->y_list.end());
		ckt->y_list.resize(std::distance(ckt->y_list.begin(), it_vec));

		// clog<<"x_list size: "<<x_list.size()<<endl;
		// clog<<"y_list size: "<<y_list.size()<<endl;	
		// now build x/y_list idmap
		map<long, int> x_list_id_map;
		map<long, int> y_list_id_map;
		for(size_t j=0;j<ckt->x_list.size();j++)
			x_list_id_map[ckt->x_list[j]] = j;
	
		for(size_t j=0;j<ckt->y_list.size();j++)
			y_list_id_map[ckt->y_list[j]] = j;
		
#if 0
		for(size_t i=0;i<x_list.size();i++){
			cout<<"i, x_list: "<<i<<" "<<x_list[i]<<endl;
		}
		for(size_t i=0;i<y_list.size();i++){
			cout<<"i, y_list: "<<i<<" "<<y_list[i]<<endl;
		}
#endif

// #if 0	
		// then accumulate the count
		// map<long, int> x_list_map;
		// map<long, int> y_list_map;
		Net *net;
		Node *na, *nb;
		long x1, x2, y1, y2;
		// only explore resistor net
		int type_r = RESISTOR;
		// for(size_t i=0;i<(*p_ckts).size();i++){
			// Circuit *ckt = (*p_ckts)[i];
		NetList &ns = ckt->net_set[type_r];
		// clog<<"ckt net.size: "<<ckt->get_name()<<" "<<ns.size()<<endl;
		for(size_t k=0;k<ns.size();k++){
			na = ns[k]->ab[0]->rep;
			nb = ns[k]->ab[1]->rep;
			if(na->is_ground() || nb->is_ground())
				continue;
			x1 = na->pt.x;
			x2 = nb->pt.x;
			// x1 point to smaller one
			if(x1 > x2){
				x1 = nb->pt.x;
				x2 = na->pt.x;
			}
			
			y1 = na->pt.y;
			y2 = nb->pt.y;
			if(y1 > y2){
				y1 = nb->pt.y;
				y2 = na->pt.y;
			}
			
			// clog<<"net: "<<*ns[k]<<endl;
			int id_1 = x_list_id_map[x1];
			int id_2 = x_list_id_map[x2];
			// clog<<"id_1, id_2: "<<id_1<<" "<<id_2<<endl;
			// map +1
			for(size_t l = id_1;l<=id_2;l++){
				// clog<<"i, x_list, x1, x2: "<<l<<" "<<x_list[l]<<" "<<x1<<" "<<x2<<endl;
				ckt->x_list_bd_map[ckt->x_list[l]]++;	
			}
			
			id_1 = y_list_id_map[y1];
			id_2 = y_list_id_map[y2];
			// map +1
			for(size_t l = id_1;l <= id_2;l++){
				ckt->y_list_bd_map[ckt->y_list[l]]++;	
			}
		}
#if 0		
		clog<<"start to output x and y list. "<<endl;
		// scan the nets to count bd nets 
		map<long, int>::iterator it;
		for(it = ckt->x_list_bd_map.begin();it != ckt->x_list_bd_map.end();it++){
			cout<<"ckt, x, count: "<<ckt->get_name()<<" "<<it->first<<" "<<it->second<<endl;	
		}
		cout<<endl;
		for(it = ckt->y_list_bd_map.begin();it != ckt->y_list_bd_map.end();it++){
			cout<<"ckt, y, count: "<<ckt->get_name()<<" "<<it->first<<" "<<it->second<<endl;	
		}
#endif
		x_list_id_map.clear();
		y_list_id_map.clear();
		// x_list.clear();
		// y_list.clear();
	}
// #endif
}

// explore the partitions with x_y_list_map
// no overlap here
void Parser::explore_partition(int num_procs, MPI_CLASS &mpi_class){
	int core_limit = 8;//num_procs;
	// assume the maximum cores are 8
	Core_x = core_limit;
	Core_y = core_limit;
	int num_cores=0;
	x_coord_vec.resize((*p_ckts).size());
	y_coord_vec.resize((*p_ckts).size());
	x_bd_vec.resize((*p_ckts).size());
	y_bd_vec.resize((*p_ckts).size());
	// first need to fix the number of partitions
	// based on number of boundary nets
	for(int i=0;i<Core_x;i++){
		//int j=1;{
		for(int j=0;j<Core_y;j++){
			num_cores = i*j;
			// generate the combination
			if(num_cores <=1 || num_cores > core_limit)
				continue;
			// for each of the partition, statistically calc the bd nets and so on	
			// x_list_bd_map, x_list_nd_map
			// explore_one_partition(i, j);	
			explore_one_partition_balance(i, j);	
		}
	}

	// find out whether worth doing x or y partition
	for(size_t i=0;i<(*p_ckts).size();i++){
		Circuit *ckt = (*p_ckts)[i];
		// compare the avg value of x_bd_vec and y_bd_vec
		double avg_x_bd = 0;
		for(size_t j=0;j<x_bd_vec[i].size();j++)
			avg_x_bd += x_bd_vec[i][j];
		if(x_bd_vec[i].size()!=0)
			avg_x_bd /= x_bd_vec[i].size();
		double avg_y_bd = 0;
		for(size_t j=0;j<y_bd_vec[i].size();j++)
			avg_y_bd += y_bd_vec[i][j];
		if(y_bd_vec[i].size()!=0)
			avg_y_bd /= y_bd_vec[i].size();
		cout<<"ckt, avg_x_bd and y_bd: "<<ckt->get_name()<<" "<<avg_x_bd<<" "<<avg_y_bd<<endl;
		if(avg_x_bd > 3*avg_y_bd){
			clog<<"only doing y partition. "<<endl;
			mpi_class.X_BLOCKS=1;
			mpi_class.Y_BLOCKS=num_procs;	
		}
		else if(avg_y_bd > 3*avg_x_bd){
			clog<<"only doing x partition. "<<endl;
			mpi_class.X_BLOCKS=num_procs;
			mpi_class.Y_BLOCKS=1;
		}
		else{
			clog<<"doing both x and y partition. "<<endl;
			clog<<"need to manually assign the partition. "<<endl;
		}
		// clog<<"mpi_class.X_BLOCKS and Y_BLOCKS: "<<mpi_class.X_BLOCKS<<" "<<mpi_class.Y_BLOCKS<<endl;
	}
}

// explore one partition: X x Y
// try to balance the number of nodes: +(-)10% of nodes of each block
void Parser::explore_one_partition_balance(int x_blocks, int y_blocks){
	// cout<<endl<<"par: X x Y: "<<x_blocks<<" "<<y_blocks<<endl;
	double len_per_x, len_per_y;
	size_t num_nodes = 0;
	size_t num_nodes_x = 0;
	size_t num_nodes_y = 0;
	// only compare x or y direction
	if(x_blocks > 1 && y_blocks > 1)
		return;
	for(int i=0;i<(*p_ckts).size();i++){
		Circuit *ckt = (*p_ckts)[i];
		num_nodes = ckt->nodelist.size()-1;
		num_nodes_x = num_nodes / x_blocks;
		num_nodes_y = num_nodes / y_blocks;
		// cout<<"ckt, num_nodes_x, num_nodes_y: "<<ckt->get_name()<<" "<<num_nodes_x<<" "<<num_nodes_y<<endl;
		size_t sum_nodes_x = 0;
		size_t sum_nodes_y = 0;
		if(x_blocks >1){
			double thresh_l = num_nodes_x * 0.5;
			double thresh_u = num_nodes_x * 1.5;
			// cout<<"thresh_x l / u: "<<thresh_l<<" "<<thresh_u<<endl;
			long min_bd = 0;
			long min_coord = 0;
			size_t min_sum_nodes_x = 0;
			size_t min_j = 0;
			int count = 0;
			for(size_t j=0;j<ckt->x_list.size();j++){
				long coord = ckt->x_list[j];
				sum_nodes_x += ckt->x_list_nd_map[coord];
				// cout<<"coord, num_nodes_x: "<<coord<<" "<<num_nodes_x<<endl;
				// start to sample the bd nodes
				if(sum_nodes_x >= thresh_l && sum_nodes_x <= thresh_u){
					// clog<<"current, min: "<<ckt->x_list_bd_map[coord]<<" "<<min_bd<<endl;
					if(min_bd == 0 || ckt->x_list_bd_map[coord] < min_bd){
						min_bd = ckt->x_list_bd_map[coord];
						min_coord = coord;
						min_sum_nodes_x = sum_nodes_x;
						min_j = j;
						// clog<<"min_bd, coord: "<<min_bd<<" "<<min_coord<<" "<<min_sum_nodes_x<<endl;
					}
				}
				if(sum_nodes_x > thresh_u){
					x_coord_vec[i].push_back(min_coord);
					x_bd_vec[i].push_back(min_bd);
					// cout<<"x / count, num_nodes, min_bd, min_coord: "<<count++<<" "<<min_sum_nodes_x<<" "<<min_bd<<" "<<min_coord<<endl;
					min_bd = 0;
					if(j == ckt->x_list.size()-1)
						break;
					// start to search next one
					j = min_j;
					sum_nodes_x = 0;
				}
			}
			// cout<<"x min_bd, min_coord: "<<min_bd<<" "<<min_coord<<endl;
		}
		if(y_blocks >1){
			double thresh_l = num_nodes_y * 0.9;
			double thresh_u = num_nodes_y * 1.1;
			// cout<<"thresh_y l / u: "<<thresh_l<<" "<<thresh_u<<endl;
			long min_bd = 0;
			long min_coord = 0;
			size_t min_sum_nodes_y = 0;
			size_t min_j=0;
			int count = 0;
			// cout<<"start to calculate y bd info. "<<y_blocks<<endl;
			for(size_t j=0;j<ckt->y_list.size();j++){
				long coord = ckt->y_list[j];
				sum_nodes_y += ckt->y_list_nd_map[coord];
				// cout<<"coord, num_nodes_x: "<<coord<<" "<<num_nodes_x<<endl;
				// start to sample the bd nodes
				if(sum_nodes_y >= thresh_l && sum_nodes_y <= thresh_u){
					// clog<<"current, min: "<<ckt->x_list_bd_map[coord]<<" "<<min_bd<<endl;
					if(min_bd == 0 || ckt->y_list_bd_map[coord] < min_bd){
						min_bd = ckt->y_list_bd_map[coord];
						min_coord = coord;
						min_sum_nodes_y = sum_nodes_y;
						min_j = j;
						// clog<<"min_bd, coord: "<<min_bd<<" "<<min_coord<<" "<<min_sum_nodes_x<<endl;
					}
				}
				if(sum_nodes_y > thresh_u){
					y_coord_vec[i].push_back(min_coord);
					y_bd_vec[i].push_back(min_bd);
					// cout<<"y / count, num_nodes, min_bd, min_coord: "<<count++<<" "<<min_sum_nodes_y<<" "<<min_bd<<" "<<min_coord<<endl;
					min_bd = 0;
					if(j == ckt->y_list.size()-1)
						break;
					// start to search next one
					j = min_j;
					sum_nodes_y = 0;	
				}
			}
			
			// clog<<"min_bd, min_coord: "<<min_bd<<" "<<min_coord<<endl;
		}
	}
}

// explore one partition: X x Y
void Parser::explore_one_partition(int x_blocks, int y_blocks){
	// clog<<"par: X x Y: "<<x_blocks<<" "<<y_blocks<<endl;
	double len_per_x, len_per_y;
	for(int i=0;i<1;i++){//(*p_ckts).size();i++){
		Circuit *ckt = (*p_ckts)[i];
		// clog<<"ckt->lx, ly, ux, uy: "<<ckt->lx<<" "<<ckt->ly<<" "<<ckt->ux<<" "<<ckt->uy<<endl;
		len_per_x = (ckt->ux - ckt->lx+0.5)*1.0 / x_blocks;
		len_per_y = (ckt->uy - ckt->ly+0.5)*1.0 / y_blocks;
		// clog<<"ckt, len_x / y: "<<ckt->get_name()<<" "<<len_per_x<<" "<<len_per_y<<endl;
		int sum_bd_nets = 0;
		
		// bd nets along x dir
		for(int j=0;j<x_blocks-1;j++){
			long bd_x = j*len_per_x + floor(len_per_x);
			// clog<<"x, bd: "<<bd_x<<endl;
			double min_dist = -1;
			int min_id = 0;
			for(int k=0;k<ckt->x_list.size();k++){
				double dist = fabs(bd_x - ckt->x_list[k]);
				if(min_dist == -1){
					min_dist = dist;
					min_id = ckt->x_list[k];
				}
				else if(dist < min_dist){
					min_dist = dist;
					min_id = ckt->x_list[k];
				}	
			}
			clog<<"min_dist, x_list: "<<min_id<<" "<<ckt->x_list_bd_map[min_id]<<endl;
			sum_bd_nets += ckt->x_list_bd_map[min_id];
			clog<<"ckt, par, x bd: "<<ckt->get_name()<<" "<<j<<" "<<ckt->x_list_bd_map[min_id]<<endl;
			clog<<"sum_bd nets: "<<sum_bd_nets<<endl;
		}
		// bd nets along y dir
		for(int j=0;j<y_blocks-1;j++){
			long bd_y = j*len_per_y + floor(len_per_y);
			// clog<<"x, bd: "<<bd_x<<endl;
			double min_dist = -1;
			int min_id = 0;
			for(int k=0;k<ckt->y_list.size();k++){
				double dist = fabs(bd_y - ckt->y_list[k]);
				if(min_dist == -1){
					min_dist = dist;
					min_id = ckt->y_list[k];
				}
				else if(dist < min_dist){
					min_dist = dist;
					min_id = ckt->y_list[k];
				}	
			}
			clog<<"min_dist, y_list: "<<min_id<<" "<<ckt->y_list_bd_map[min_id]<<endl;
			sum_bd_nets += ckt->y_list_bd_map[min_id];
			clog<<"ckt, par, y bd: "<<ckt->get_name()<<" "<<j<<" "<<ckt->y_list_bd_map[min_id]<<endl;
			// clog<<"sum_bd nets: "<<sum_bd_nets<<endl;
		}
	}
}

// clean all the resources for exploring the partitioning
void Parser::clean_explore_partition(int num_procs, MPI_CLASS &mpi_class){
	x_coord_vec.clear();
	y_coord_vec.clear();
	x_bd_vec.clear();
	y_bd_vec.clear();
	for(size_t i=0;i<(*p_ckts).size();i++){
		Circuit *ckt = (*p_ckts)[i];
		ckt->clean_explore();	
	}
	// clog<<"after release ckt resource"<<endl;
}
