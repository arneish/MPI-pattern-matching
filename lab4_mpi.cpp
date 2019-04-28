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
	//int i = min(i_, j_), j = max(i_, j_);
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
	//bool match = true;
	for (int i = pos_T; i < pos_T + m; i++)
	{
		if (T[i] != P[i - pos_T])
		{
			return false;
			// match = false;
			// break;
		}
	}
	return true;
}

void NPPM(char *T, int &n, char *P, int &m, int *witP, int &len_witP, int &periodP, vector<int> &all_match_pos) //Non-Periodic Pattern Matching
{
	//partition T[0:n-1] into n/(ceil(m/2)) blocks (each block handles ceil(m/2) consecutive chars of T)
	int block_size = ceil((double)m / 2), num_blocks = n / block_size;
	//int num_blocks = ceil((double)n / block_size);
	//int potential_positions[num_blocks];
	int *potential_positions = (int *)malloc(sizeof(int) * num_blocks);
	int start_index, end_index, winner;
	for (int bid = 0; bid < num_blocks; bid++)
	{
		//compute start and end for block "bid":
		start_index = bid * block_size;
		end_index = start_index + block_size - 1;
		//end_index = min(start_index + block_size - 1, n - 1);
		winner = start_index;

		//O(m/2) duel:
		for (int j = start_index + 1; j <= end_index; j++)
		{
			// int k=witP[j-winner];
			// if (((j+k)>=n)||T[j+k]!=P[k])
			// {

			// }
			// else
			// {
			// 	winner = j;
			// }
			duel(T, n, P, m, witP, len_witP, periodP, winner, j);
		}
		potential_positions[bid] = winner;
	}
	//all_match_pos.clear();
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
	//matchpos_Psub.reserve(ceil(2.0*n/m));
	NPPM(text, n, Psub, len_Psub, witP, len_witP, periodP, matchpos_Psub);

	//P = u^k.v
	// char u[periodP + 1];
	// memcpy(u, P, periodP);
	// u[periodP] = '\0';
	

	int k = m / periodP; //CHECKPOINT: what if m%periodP ==0?
	int len_v = m - k * periodP;

	int v_startindex = k*periodP;
	// char v[len_v + 1];
	// memcpy(v, P + k * periodP, len_v);
	// v[len_v] = '\0';
	
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

	// vector<int> *S;
	// int len_S, index;
	// int *sum_S;
	
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


	// S = new vector<int>();
	// for (i = 0; i < periodP; i++)
	// {
	// 	S->clear();
	// 	for (int j = i; j < n; j += periodP)
	// 	{
	// 		S->push_back(M[j]);
	// 	}
	// 	len_S = S->size();
	// 	//vector<int> *C = new vector<int> (len_S, 0);
	// 	sum_S = (int *)calloc(len_S, sizeof(int));
	// 	if (S->at(len_S - 1) == 1)
	// 	{
	// 		sum_S[len_S - 1] = 1;
	// 		if (1 >= k - 1)
	// 		{
	// 			index = i + (len_S - 1) * periodP; //i+l*p
	// 			matches_local.push_back(index);
	// 			//matching_indices->push_back(index);
	// 			//C->at(len_S-1)=1;
	// 		}
	// 	}
	// 	for (j = len_S - 2; j >= 0; j--)
	// 	{
	// 		if (S->at(j) == 0)
	// 			sum_S[j] = 0;
	// 		else
	// 		{
	// 			sum_S[j] = sum_S[j + 1] + 1;
	// 		}
	// 		if (sum_S[j] >= k - 1)
	// 		{
	// 			index = i + j * periodP;
	// 			matches_local.push_back(index);
	// 			//matching_indices->push_back(index);
	// 			//C->at(j)=1;
	// 		}
	// 	}
	// 	//C_all.push_back(C);
	// 	free(sum_S);
	// 	//free(S);
	// }
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
	//int total_alloted = 0;
	for (int proc = 0; proc < num_procs; proc++)
	{
		patterns_per_proc = ceil((double)num_patterns_remain / num_procs_remain);
		allprocs_patternID[proc] = make_pair(pattern_index, patterns_per_proc);
		//total_alloted += patterns_per_proc;
		num_patterns_remain -= patterns_per_proc;
		num_procs_remain--;
		pattern_index += patterns_per_proc;
	}
	//assert(total_alloted == num_patterns);
	int start_index = allprocs_patternID[proc_id].first, assigned_size = allprocs_patternID[proc_id].second, end_index = start_index + assigned_size - 1;

	vector<int> match_counts_local(assigned_size, 0);
	//int total_count = 0;
	vector<int> matches_local;
	//matches_local.reserve(1194419);
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
		//total_count += current_matches_count;
		index++;
		free(witP);
	}
	//assert(total_count == matches_local.size());

	//send match_counts_local array (size "assigned_size") to proc0:

	//send-buffer setup:
	int *sendbuf_first = nullptr;
	if (assigned_size)
	{
		sendbuf_first = &match_counts_local[0];
	}
	int sendcount_first = assigned_size;

	//receiver-buffer setup:
	//int *recvbuf_first;
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
		// for (int i=0; i<recvbuf_first_size; i++) //COMMENT OUT!
		// {
		// 	(*match_counts)[i]=-911;
		// }
	}

	//sending match-counts for the "assigned patterns"
	MPI_Gatherv(sendbuf_first, sendcount_first, MPI_INT, *match_counts, recvcounts_first, displs_first, MPI_INT, 0, MPI_COMM_WORLD);

	// if (proc_id==0)
	// {
	// 	for (int i=0; i<recvbuf_first_size; i++)
	// 	{
	// 		cout<<(*match_counts)[i]<<",";
	// 		assert((*match_counts)[i]>=0);
	// 	}
	// 	cout<<"\n";
	// }

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
	//int *recvbuf_second;
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
		//cout<<"recvbuf_second_size:"<<recvbuf_second_size<<"\n";
		*matches = (int *)malloc(sizeof(int) * recvbuf_second_size);
		// for (int i=0; i<recvbuf_second_size; i++) //COMMENT OUT!
		// {
		// 	(*matches)[i]=-786;
		// }
	}

	MPI_Gatherv(sendbuf_second, sendcount_second, MPI_INT, *matches, recvcounts_second, displs_second, MPI_INT, 0, MPI_COMM_WORLD);

	// if (proc_id==0)
	// {

	// 	cout<<"\n"<<recvbuf_second_size;
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
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
	//printf("\nTOTAL TIME:%d\n", proc_id);

	return;
}
	// vector<string> *pattern_vec = new vector<string>(num_patterns);
	// string T = "babababababaabab";
	// int lenT = T.size();
	// string P = "abababa";
	// int lenP = P.size();
	// int periodP = 2;
	// int *witP = (int *)malloc(sizeof(int) * periodP);
	// int lenwitP = 2;
	// PPM(T, lenT, P, lenP, witP, lenwitP, periodP);
/* rough
	// int numtasks, rank, next, prev, buf[2], tag1=1, tag2=2;
	// MPI_Request reqs[4];
	// MPI_Status stats[4];
	// MPI_Init(NULL, NULL);
	// MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	// prev = rank-1;
	// next = rank+1;
	// if (rank==0)
	// {
	// 	prev = numtasks-1;
	// }
	// if (rank==numtasks-1)
	// {
	// 	next = 0;
	// }
	// MPI_Irecv(&buf[0], 1, MPI_INT, prev, tag1, MPI_COMM_WORLD, &reqs[0]);
	// MPI_Irecv(&buf[1], 1, MPI_INT, next, tag2, MPI_COMM_WORLD, &reqs[1]);

	// MPI_Isend(&rank, 1, MPI_INT, prev, tag2, MPI_COMM_WORLD, &reqs[2]);
	// MPI_Isend(&rank, 1, MPI_INT, next, tag1, MPI_COMM_WORLD, &reqs[3]);

	// MPI_Waitall(4, reqs, stats);
	// printf("pID:%d :: %d, %d\n", rank, buf[0], n);
	// MPI_Finalize();

// 	int numtasks, rank, sendcount, recvcount, source;
//    float sendbuf[SIZE][SIZE] = {
//      {1.0, 2.0, 3.0, 4.0},
//      {5.0, 6.0, 7.0, 8.0},
//      {9.0, 10.0, 11.0, 12.0},
//      {13.0, 14.0, 15.0, 16.0}  };
//    float recvbuf[SIZE];

//    MPI_Init(NULL,NULL);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

//    if (numtasks == SIZE) {
//      // define source task and elements to send/receive, then perform collective scatter
//      source = 1;
//      sendcount = SIZE;
//      recvcount = SIZE;
//      MPI_Scatter(sendbuf,sendcount,MPI_FLOAT,recvbuf,recvcount,
//                  MPI_FLOAT,source,MPI_COMM_WORLD);

//      printf("rank= %d  Results: %f %f %f %f\n",rank,recvbuf[0],
//             recvbuf[1],recvbuf[2],recvbuf[3]);
//      }
//    else
//      printf("Must specify %d processors. Terminating.\n",SIZE);

//    MPI_Finalize();
*/