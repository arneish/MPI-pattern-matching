#ifndef LAB4_MPI_H
#define LAB4_MPI_H


/*
To be implemented
Arguments:
    n     : integer containing number of characters (length) of text 
            in which patterns are to be searched (input)
    text  : character array containing text in which patterns are to be searched (input)
    num_patterns: integer containing #patterns to be matched in the text (input)
    m_set : integer array containing length of patterns (input)
            #elements in m_set = num_patterns
            --------------------------------------------------------------------------
            | len(pattern[0]) | len(pattern[1]) | ... | len(pattern[num_patterns-1]) |
            --------------------------------------------------------------------------
    p_set : integer array containing length of period in patterns (input)
            #elements in p_set = num_patterns
            -----------------------------------------------------------------------------------
            | period(pattern[0]) | period(pattern[1]) | ... | period(pattern[num_patterns-1]) |
            -----------------------------------------------------------------------------------
    pattern_set : array of character array containing patterns to be matchede (input)
            #elements in pattern_set = num_patterns
            -----------------------------------------------------------
            | pattern[0] | pattern[1] | ... | pattern[num_patterns-1] |
            -----------------------------------------------------------
    match_counts : 1D integer array containing number of matches of each pattern in text (output)
            #elements in ocuurance_count = num_patterns
            -----------------------------------------------------------------------------------------
            | #matches(pattern[0]) | #matches(pattern[1]) | ... | #matches(pattern[num_patterns-1]) |
            -----------------------------------------------------------------------------------------
    matches : 1D array of integers containing list of all matches (start index of) of pattern_i in text (output)
            consider index of text starting from 0 (not 1)
            #elements in matches = sum(match_counts)
            -------------------------------------------------------------------------------------------------
            | match(pattern[0])[0] | match(pattern[0])[1] | ... | match(pattern[0])[#matches(pattern[0])-1] |
            -------------------------------------------------------------------------------------------------
            | match(pattern[1])[0] | match(pattern[1])[1] | ... | match(pattern[1])[#matches(pattern[1])-1] |
            -------------------------------------------------------------------------------------------------
            | ... ... ... ... ... ... | match(pattern[num_patterns-1])[#matches(pattern[num_patterns-1])-1] |
            -------------------------------------------------------------------------------------------------
*/
void periodic_pattern_matching (
		int n, 
		char *text, 
		int num_patterns, 
		int *m_set, 
		int *p_set, 
		char **pattern_set, 
		int **match_counts, 
		int **matches);

#endif