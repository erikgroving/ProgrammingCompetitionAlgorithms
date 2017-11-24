#include <stdio.h>
#include <stdlib.h>
#include "contest.h"
#include <time.h>
#include <math.h>
#define TIME
int main() {
	struct group* groups = NULL;
	struct group** in_range = (struct group**)malloc(sizeof(struct group*));
	struct ap* aps = NULL;
	struct wall* walls = NULL;
	int num_groups;
	int num_aps;
	int num_walls;
	int num_valid;
	int num_out;

	
	struct bucket b[BUCK * BUCK];
		
	#ifdef TIME
	struct timespec start,  end;
	clock_gettime(CLOCK_REALTIME, &start);
	#endif	


	// Parse and initialize
	parseInput(&groups, &aps, &walls, &num_groups, &num_aps, &num_walls);


	float incy;
	float incx;
	// find max and min for wall

	if (num_walls){
		int max_x = 0;
		int max_y = 0;
		for (int i = 0; i < num_walls; i++) {
			int wmax_x = max (walls[i].x1, walls[i].x2);
			max_x = max(max_x, wmax_x);
			int wmax_y = max (walls[i].y1, walls[i].y2);
			max_y = max(max_y, wmax_y);
		}
		
		incx = (float)max_x / BUCK;
		incy = (float)max_y / BUCK;
		
		makeBuckets(num_walls, walls, b, incx, incy);
	}



	// Calculate which access point are in range for a group of students
	groupsInRange(in_range, groups, b, aps, walls, num_groups, num_aps, num_walls, &num_valid, &num_out, incx, incy);

	
	

	
	// Create the graph with the groups that are in range!
	// and the access points
	struct edge** adj_list = NULL;
	int* degree = (int *)malloc(sizeof(int) * (num_valid + num_aps + 2));

	adj_list = createGraph(&adj_list, *in_range, aps, &degree, num_aps, num_valid);
	

	
	// Ford-Fulkerson on the graph!
	int vertices = num_valid + num_aps + 2;

	
	

	int max_flow = maxFlow(adj_list, vertices, degree, num_groups, num_walls);


	printf("%d %d\n", num_out, max_flow);

		
	/* it actually takes a few milliseconds to iterate over all the 
	linked lists and free the nodes..........
	WE ARE TERMINATING ANYWAY D:, THE MEMORY LINKED TO THE PROCESS
	IS FREED... DON'T YELL AT ME T___T */
	/*
	for (int i = 0; i < vertices; i++) {
		if (adj_list[i]) {
			free(adj_list[i]);
		}
	}
	
	for (int i = 0; i < num_groups; i++) {
		if (groups[i].head != NULL) {
			struct ap_ll* s1 = groups[i].head;
			while (s1 != NULL) {
				struct ap_ll* s2 = s1->next;
				free(s1);
				s1 = s2;
			}
		}
	}
	if(num_walls) {  
		for (int i = 0; i < BUCK * BUCK; i++) {

				if (b[i].size) {
					free(b[i].w);
				}

		}
	}
	free((*in_range));
	free(in_range);
	free(adj_list);
	free(degree);
	free(groups);
	free(aps);
	free(walls);*/
	#ifdef TIME
	clock_gettime(CLOCK_REALTIME, &end);
	
	double diff = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
	printf("TIME TO RUN: %lf\n", diff);
	#endif	

	return 0;
}


