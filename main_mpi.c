#include "lab4_io.h"
#include "lab4_mpi.h"

#include <stdlib.h>
#include "mpi.h"

/*
	Arguments:
		arg: input filename (consist text, pattern_set)
*/

int main(int argc, char *argv[])
{
	if (argc < 2){
		printf("\nLess Arguments\n");
		return 0;
	}

	if (argc > 2){
		printf("\nTOO many Arguments\n");
		return 0;
	}

	//---------------------------------------------------------------------
	int 	n;				// length of text (input)
	char* 	text;			// text (input)
	int 	num_patterns;	// #patterns to be searched in the text (input)
	int* 	m_set;			// lengths of patterns in pattern_set (input)
	int* 	p_set;			// periods of patterns in pattern_set (input)
	char** 	pattern_set;	// set of patterns to be searched (input)
	int* 	match_counts;	// #match of pattern_i in text (to be computed)
	int* 	matches;		// set of all match of each pattern_i in text (to be computed)
	//---------------------------------------------------------------------

	double start_time, computation_time, total_time;
	int id;

	/*
		-- Pre-defined function --
		reads input text and patterns from input file and creats array text, m_set, p_set and pattern_set
	    see lab4_io.h for details
	*/
	read_data (argv[1], &n, &text, &num_patterns, &m_set, &p_set, &pattern_set);

	MPI_Init(&argc, &argv);

	start_time = MPI_Wtime();
	
	// /*
	// 	*****************************************************
	// 		TODO -- You must implement this two function
	// 	*****************************************************
	// */
	periodic_pattern_matching (
		n, 
		text, 
		num_patterns, 
		m_set, 
		p_set, 
		pattern_set, 
		&match_counts, 
		&matches);

	computation_time = MPI_Wtime() - start_time;
	
	MPI_Reduce(&computation_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	if (id == 0) {
		/*
		--Pre-defined function --
		checks for correctness of results computed by periodic_pattern_matching
		and outputs the results
		*/
		write_result(match_counts, matches, total_time);
	}

	MPI_Finalize();

	return 0;
}
