#include "lab4_mpi.h"
#include <iostream>
#include <malloc.h>
#include <cstring>
#include <string>
#include <cmath>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <unordered_set>
#include "mpi.h"

using namespace std;

#pragma GCC optimize("Ofast")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")

inline void duel(char *Z, int &n, char *Y, int &m, int *witY, int &len_witY, int &periodY, int &i, int &j)
{
	//assert((i >= 0) && (i < j) && (j < n) && ((j - i) < periodY));
	int k = witY[j - i];
	if (((j + k) >= n) || Z[j + k] != Y[k])
	{
		return;
	}
	else
	{
		i = j;
	}
}

void compute_witness(char *Y, int &m, int &periodY, int *witY)
{
	//length of witY is periodY
	bool found = false;
	witY[0] = 0;
	for (int i = 1; i < periodY; i++)
	{
		//compute witY[i]: find 'k'
		found = false;
		for (int k = 0; k <= m - 1 - i; k++)
		{
			if (Y[k] != Y[k + i])
			{
				witY[i] = k;
				found = true;
				break;
			}
		}
		//assert(found == true);
	}
}

inline bool brute_force_match(char *T, int &n, char *P, int &m, int pos_T)
{
	if (pos_T + m > n) //exceeding T's size
	{
		return false;
	}
	for (int i = pos_T; i < pos_T + m; i++)
	{
		if (T[i] != P[i - pos_T])
		{
			return false;
		}
	}
	return true;
}

void NPPM(char *T, int &n, char *P, int &m, int *witP, int &len_witP, int &periodP, vector<int> &all_match_pos) //Non-Periodic Pattern Matching
{
	//partition T[0:n-1] into n/(ceil(m/2)) blocks (each block handles ceil(m/2) consecutive chars of T)
	int block_size = ceil((double)m / 2), num_blocks = n / block_size;
	int *potential_positions = (int *)malloc(sizeof(int) * num_blocks);
	int start_index, end_index, winner;
	for (int bid = 0; bid < num_blocks; bid++)
	{
		//compute start and end for block "bid":
		start_index = bid * block_size;
		end_index = start_index + block_size - 1;
		winner = start_index;

		//O(m/2) duel:
		for (int j = start_index + 1; j <= end_index; j++)
		{
			duel(T, n, P, m, witP, len_witP, periodP, winner, j);
		}
		potential_positions[bid] = winner;
	}
	bool match;
	for (int bid = 0; bid < num_blocks; bid++)
	{
		//brute-force check for P in T at potential position:
		match = brute_force_match(T, n, P, m, potential_positions[bid]);
		if (match)
		{
			all_match_pos.push_back(potential_positions[bid]);
		}
	}
	free(potential_positions);
}

void PPM(char *text, int n, char *P, int m, int *witP, int len_witP, int periodP, vector<int> &matches_local) //Periodic Pattern Matching
{
	//assert(len_witP == periodP);
	//P_sub:
	int len_Psub = 2 * periodP - 1;
	
	char Psub[len_Psub + 1];
	memcpy(Psub, P, len_Psub);
	Psub[len_Psub] = '\0';

	//witP: witness array for P_sub
	compute_witness(Psub, len_Psub, periodP, witP);

	//computing pos_Psub;
	vector<int> matchpos_Psub;
	NPPM(text, n, Psub, len_Psub, witP, len_witP, periodP, matchpos_Psub);

	//P = u^k.v
	int k = m / periodP; 
	int len_v = m - k * periodP;
	int v_startindex = k*periodP;
	
	int *M = (int *)calloc(n, sizeof(int));
	bool IsPresent;
	int nextindex_text, i, j;
	for (const auto &Psub_index : matchpos_Psub)
	{
		//check if u^2.v is at text[Psub_index]
		IsPresent = true;
		nextindex_text = Psub_index + 2 * periodP - 1;
		if (nextindex_text == n || text[nextindex_text] != P[periodP - 1])
		{
			IsPresent = false;
		}
		else
		{
			//check if remainder is v if text exists at that index
			nextindex_text++;
			j = 0;
			while (j < len_v)
			{
				if (nextindex_text == n || text[nextindex_text] != P[v_startindex+j])//v[j])
				{
					IsPresent = false;
					break;
				}
				else
				{
					nextindex_text++;
					j++;
				}
			}
		}
		if (IsPresent)
		{
			M[Psub_index] = 1;
		}
	}
	//successfully populated M

	int num_blocks = ceil((double)n/periodP), last_index, current_index, sum_next=0, sum_current=0;
	for (i=0; i<periodP; i++)
	{
		last_index = i+(num_blocks-1)*periodP;
		if (last_index>=n)
		{
			last_index-=periodP;
		}
		current_index = last_index;
		while (current_index>=0)
		{
			if (M[current_index]==0)
			{
				sum_current=0;
			}
			else
			{
				sum_current=1+sum_next;
				if (sum_current>=k-1)
				{
					matches_local.push_back(current_index);
				}
			}
			sum_next = sum_current;
			current_index-=periodP;
		}
	}
	free(M);
}

void periodic_pattern_matching(
	int n,
	char *text,
	int num_patterns,
	int *m_set,
	int *p_set,
	char **pattern_set,
	int **match_counts,
	int **matches)
{
	/*Implementing parallel algorithm (on a single processor) for many patterns concurrently*/
	// high_resolution_clock::time_point t_begin, t_end, t1, t2, t3;
	// t_begin = high_resolution_clock::now();
	int num_procs, proc_id, rc;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

	//computed patternid_start and number of patterns for ALL procs:
	vector<pair<int, int>> allprocs_patternID(num_procs);
	int pattern_index = 0, patterns_per_proc = 0, num_patterns_remain = num_patterns, num_procs_remain = num_procs;
	for (int proc = 0; proc < num_procs; proc++)
	{
		patterns_per_proc = ceil((double)num_patterns_remain / num_procs_remain);
		allprocs_patternID[proc] = make_pair(pattern_index, patterns_per_proc);
		num_patterns_remain -= patterns_per_proc;
		num_procs_remain--;
		pattern_index += patterns_per_proc;
	}
	int start_index = allprocs_patternID[proc_id].first, assigned_size = allprocs_patternID[proc_id].second, end_index = start_index + assigned_size - 1;

	vector<int> match_counts_local(assigned_size, 0);
	vector<int> matches_local;
	int index = 0;
	for (int i = start_index; i <= end_index; i++)
	{
		// if (proc_id == 0)
		// {
		// 	printf("\r %d/%d", i, assigned_size);
		// 	fflush(stdout);
		// }
		char *pattern = *(pattern_set + i);
		int m = m_set[i], periodP = p_set[i], lenwitP = periodP;
		int *witP = (int *)malloc(sizeof(int) * periodP);
		int size_matches_local_prev = matches_local.size();
		PPM(text, n, pattern, m, witP, lenwitP, periodP, matches_local);
		int current_matches_count = matches_local.size() - size_matches_local_prev;
		match_counts_local[index] = current_matches_count;
		index++;
		free(witP);
	}

	//send match_counts_local array (size "assigned_size") to proc0:
	//send-buffer setup:
	int *sendbuf_first = nullptr;
	if (assigned_size)
	{
		sendbuf_first = &match_counts_local[0];
	}
	int sendcount_first = assigned_size;

	//receiver-buffer setup:
	int *recvcounts_first;
	int recvbuf_first_size;
	int *displs_first;
	if (proc_id == 0)
	{
		recvcounts_first = (int *)malloc(sizeof(int) * num_procs);
		for (int i = 0; i < num_procs; i++)
		{
			recvcounts_first[i] = allprocs_patternID[i].second; //assigned size
		}
		displs_first = (int *)malloc(sizeof(int) * num_procs);
		for (int i = 0; i < num_procs; i++)
		{
			displs_first[i] = (i > 0) ? (displs_first[i - 1] + recvcounts_first[i - 1]) : 0;
		}
		recvbuf_first_size = displs_first[num_procs - 1] + recvcounts_first[num_procs - 1];
		*match_counts = (int *)malloc(sizeof(int) * recvbuf_first_size);
	}

	//sending match-counts for the "assigned patterns"
	MPI_Gatherv(sendbuf_first, sendcount_first, MPI_INT, *match_counts, recvcounts_first, displs_first, MPI_INT, 0, MPI_COMM_WORLD);

	//Now, proc0 recv matches from other procs IN ORDER OF PATTERNS EXPLORED BY THE PROCS
	//other procs send matches data for patterns IN ORDER OF PATTERNS EXPLORED BY THE PROCS

	//send-buffer-setup:
	int *sendbuf_second = nullptr;
	if (matches_local.size() > 0)
	{
		sendbuf_second = &matches_local[0];
	}
	int sendcount_second = matches_local.size();

	//receiver-buffer setup:
	int *recvcounts_second;
	int recvbuf_second_size;
	int *displs_second;
	if (proc_id == 0)
	{
		recvcounts_second = (int *)calloc(num_procs, sizeof(int));
		for (int i = 0; i < num_procs; i++)
		{
			int pattern_index_start = allprocs_patternID[i].first, pattern_index_end = pattern_index_start + allprocs_patternID[i].second - 1;
			for (int j = pattern_index_start; j <= pattern_index_end; j++)
			{
				recvcounts_second[i] += (*match_counts)[j];
			}
		}
		displs_second = (int *)malloc(sizeof(int) * num_procs);
		for (int i = 0; i < num_procs; i++)
		{
			displs_second[i] = (i > 0) ? (displs_second[i - 1] + recvcounts_second[i - 1]) : 0;
		}
		recvbuf_second_size = displs_second[num_procs - 1] + recvcounts_second[num_procs - 1];
		*matches = (int *)malloc(sizeof(int) * recvbuf_second_size);
	}

	MPI_Gatherv(sendbuf_second, sendcount_second, MPI_INT, *matches, recvcounts_second, displs_second, MPI_INT, 0, MPI_COMM_WORLD);
	
	// if (proc_id==0)
	// {
	// 	sort((*matches), (*matches)+recvbuf_second_size);
	// 	for (int i=0; i<recvbuf_second_size; i++)
	// 	{
	// 		cout<<(*matches)[i]<<",";
	// 		assert((*matches)[i]>=0);
	// 	}
	// 	cout<<"\n";
	// }
	// t_end = high_resolution_clock::now();
	// duration<double> time_span = duration_cast<duration<double>>(t_end - t_begin);
	//printf("\nTOTAL TIME:%d\n", proc_id)
}
