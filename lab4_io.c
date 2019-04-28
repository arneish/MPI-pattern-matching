#include "lab4_io.h"

void read_data (
        const char* input_filename, 
        int *n, 
        char **text, 
        int *num_patterns, 
        int **m_set, 
        int **p_set, 
        char ***pattern_set)
{
	FILE *fin = fopen(input_filename, "r");

	fscanf(fin, "%d%d", n, num_patterns);
	int np = *num_patterns;

	*text = (char*) malloc(sizeof(char)*(*n + 1));
	fscanf(fin, "%s", *text);

	*m_set = (int *) malloc(sizeof(int) * np);
	for (int i=0; i<np; i++) {
		fscanf(fin, "%d", (*m_set + i));
	}

	*p_set = (int *) malloc(sizeof(int) * np);
	for (int i=0; i<np; i++) {
		fscanf(fin, "%d", (*p_set + i));
	}

	*pattern_set = (char **) malloc(sizeof(char*) * np);
	for (int i=0; i<np; i++) {
		*(*pattern_set + i) = (char *) malloc(sizeof(char) * (*(*m_set+i) + 1));
		fscanf(fin, "%s", *(*pattern_set + i));
	}

	fclose(fin);
}

// Will contain output code
void write_result (
        int *match_counts, 
        int *matches, 
        double computation_time)
{
	printf("\nTime elapsed: %f\n", computation_time);
}

void format_checker (
        int num_patterns,
        int *match_counts, 
        int *matches)
{
	int index=0;
	for (int i=0; i<num_patterns; i++) {
		printf("pattern %d found at %d locations:\n", i, match_counts[i]);
		for (int j=0; j<match_counts[i]; j++) {
			printf("%d\t", matches[index]);
			index++;
		}
		printf("\n");
	}	
}
