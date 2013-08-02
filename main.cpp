#include "main.h"
#include "cholmod.h"
#include "mpi.h"
#include "circuit.h"
//#include <omp.h>

const char * usage="Usage: %s [-eorbILifl] benchmark\n\
    -e EPSILON\n\
    -r mpi overlap ratio\n\
    -b max block nodes\n\
    -i input file\n\
    -f output file\n\
    -l log file (default to screen)\n\
    -x mpi x block number\n\
    -y mpi y block number\n\
    -a explore partition(1) or solving(0)"
;

const char * usage2="Usage: %s -i input -a partition_flag -x X_BLOCK -y Y_BLOCK -f output\n";

int main(int argc, char * argv[]){
	int my_id;
	int num_procs;
	int c;
	int x_blocks = 0;
	int y_blocks = 0;
	float mpi_olap_ratio;
	double epsilon, omega, overlap_ratio;
	size_t max_block_nodes;
	//char * logfile="/dev/null";
	char * logfile=NULL;
	char * input=NULL, * output=NULL;
	bool partition_flag = false;
	bool input_flag = false, output_flag = false;
	Circuit::get_parameters(epsilon, overlap_ratio, 
			max_block_nodes);
	MPI_CLASS::get_parameters(x_blocks, y_blocks, mpi_olap_ratio);

	while( ( c = getopt(argc, argv, "i:a:f:e:r:b:l:x:y:L")) != -1 ){
		switch(c){
		case 'e':
			epsilon = atof(optarg);
			break;
		case 'r':
			mpi_olap_ratio = atof(optarg);
			break;
		case 'b':
			max_block_nodes = atof(optarg);
			break;
		case 'a':
			partition_flag = atoi(optarg);
			break;
		case 'l':
			logfile = optarg;
			break;
		case 'i':
			input = optarg;
			input_flag = true;
			break;
		case 'f':
			output = optarg;
			output_flag = true;
			break;
		case 'x':
			x_blocks = atoi(optarg);
			break;
		case 'y':
			y_blocks = atoi(optarg);
			break;
		case '?':
		default:
			fprintf(stderr, usage2, argv[0]);
			exit(EXIT_FAILURE);
		}
	}
	//if( argc == optind ) report_exit(usage2);
	if( !input_flag || ! output_flag ){
		fprintf(stderr, usage2, argv[0]);
		exit(EXIT_FAILURE);
	}
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	if(my_id==0) clog<<"num_procs: "<<num_procs<<endl;	

	open_logfile(logfile);
#if ENABLE_OUTPUT
	stringstream ss;
	ss<<output<<"_"<<my_id;
	if( freopen(ss.str().c_str(), "w", stdout) == NULL )
		report_exit("Ouptut file error\n");
#endif
	Circuit::set_parameters(epsilon,overlap_ratio, 
			max_block_nodes);
	MPI_CLASS::set_parameters(x_blocks, y_blocks, mpi_olap_ratio);
	if(my_id==0)
	clog<<"mpi.x_blocks, y_blocks, overlap_ratio: "<<MPI_CLASS::X_BLOCKS<<" "<<MPI_CLASS::Y_BLOCKS<<" "<<MPI_CLASS::overlap_ratio<<endl;
	// start to parfile
	vector<Circuit *> cktlist;
	MPI_CLASS mpi_class;
	// allocate class Tran for all cores
	Tran tran;
	if(my_id==0){
		mpi_class.start_task = new int [num_procs];
		mpi_class.end_task = new int [num_procs];
		mpi_class.tasks_n = new int [num_procs];
		mpi_class.MPI_Assign_Task(num_procs);
	}
	MPI_Scatter(mpi_class.tasks_n, 1, MPI_INT, 
		&mpi_class.block_size, 
		1, MPI_INT, 0, MPI_COMM_WORLD);
		
	Parser parser(&cktlist);
	clock_t t1,t2;
	t1=clock();
	if(my_id ==0)
		clog<<"readed par flag is: "<<partition_flag<<endl;
	// bool partition_flag = false;
	parser.parse(my_id, input, mpi_class, tran, num_procs, partition_flag);
	MPI_Barrier(MPI_COMM_WORLD);

	// after parsing, this mem can be released
	t2=clock();
	if(my_id==0) clog<<"Parse time="<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;

	if(partition_flag == true)
		return 0;	
	double mpi_t11, mpi_t12;
	mpi_t11 = MPI_Wtime();
	
	for(size_t i=0;i<cktlist.size();i++){
		Circuit * ckt = cktlist[i];
		if(my_id==0){
			clog<<"<======== solving: "<<ckt->get_name()<<" =========>"<<my_id<<endl;
		}
		ckt->solve(my_id, num_procs, mpi_class, 
				tran);	
		free(ckt);
	
		MPI_Barrier(MPI_COMM_WORLD);
	}
#if ENABLE_OUTPUT
	if(my_id==0)
	 	tran.print_tr_nodes();
#endif

	mpi_t12 = MPI_Wtime();
	
	// output a single ground node
	if(my_id==0){
		clog<<"solve using: "<<1.0*(mpi_t12-mpi_t11)<<endl;
	}
	// close_logfile();
	MPI_Barrier(MPI_COMM_WORLD);
	int err = MPI_Finalize();
	return 0;
}
