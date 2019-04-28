#ifndef LAB4_IO_H
#define LAB4_IO_H

#include <stdio.h>
#include <malloc.h>

/*
    n     : integer containing number of characters (length) of text 
            in which patterns are to be searched
    text  : character array containing text in which patterns are to be searched
    num_patterns: integer containing #patterns to be searched in the text
    m_set : integer array containing length of patterns
            #elements in m_set = num_patterns
            --------------------------------------------------------------------------
            | len(pattern[0]) | len(pattern[1]) | ... | len(pattern[num_patterns-1]) |
            --------------------------------------------------------------------------
    p_set : integer array containing length of periods in patterns
            #elements in p_set = num_patterns
            -----------------------------------------------------------------------------------
            | period(pattern[0]) | period(pattern[1]) | ... | period(pattern[num_patterns-1]) |
            -----------------------------------------------------------------------------------
    pattern_set : array of character array containing patterns to be matched
            #elements in pattern_set = num_patterns
            -----------------------------------------------------------
            | pattern[0] | pattern[1] | ... | pattern[num_patterns-1] |
            -----------------------------------------------------------
*/
void read_data (
        const char* input_filename, 
        int *n, 
        char **text, 
        int *num_patterns, 
        int **m_set, 
        int **p_set, 
        char ***pattern_set);

/*
check correctess of periodic_pattern_matching
Arguments:
    match_counts : integer array containing number of matches of each pattern in text
           #elements in match_counts = num_patterns
            -----------------------------------------------------------------------------------------
            | #matches(pattern[0]) | #matches(pattern[1]) | ... | #matches(pattern[num_patterns-1]) |
            -----------------------------------------------------------------------------------------
    matches : 1D array of integers containing list of all matches (start index of) of pattern_i in text
            consider index of text starting from 0 (not 1)
            -------------------------------------------------------------------------------------------------
            | match(pattern[0])[0] | match(pattern[0])[1] | ... | match(pattern[0])[#matches(pattern[0])-1] |
            -------------------------------------------------------------------------------------------------
            | match(pattern[1])[0] | match(pattern[1])[1] | ... | match(pattern[1])[#matches(pattern[1])-1] |
            -------------------------------------------------------------------------------------------------
            | ... ... ... ... ... ... | match(pattern[num_patterns-1])[#matches(pattern[num_patterns-1])-1] |
            -------------------------------------------------------------------------------------------------
    computation_time : Time elapsed in computing periodic pattern matching
*/
void write_result (
        int *match_counts, 
        int *matches, 
        double computation_time);


/* 
Function to check the format of output code.
You can call it from main to check if the array returned from your code matches with our specifications
This is a dummy function. It is not the function that will be used for evaluation.
*/
void format_checker (
        int num_patterns,
        int *match_counts, 
        int *matches);

#endif