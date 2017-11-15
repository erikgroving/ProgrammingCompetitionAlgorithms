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
	
	#ifdef TIME
	struct timespec start, end;
	clock_gettime(CLOCK_REALTIME, &start);
	#endif

	// Parse and initialize
	parseInput(&groups, &aps, &walls, &num_groups, &num_aps, &num_walls);

	
	// Calculate which access point are in range for a group of students
	groupsInRange(in_range, groups, aps, walls, num_groups, num_aps, num_walls, &num_valid, &num_out);


	
	// Create the graph with the groups that are in range!
	// and the access points
	struct edge** adj_list = NULL;
	int* degree = (int *)malloc(sizeof(int) * (num_valid + num_aps + 2));

	adj_list = createGraph(&adj_list, *in_range, aps, &degree, num_aps, num_valid);

	
	// Ford-Fulkerson on the graph!
	int vertices = num_valid + num_aps + 2;


	
	int max_flow = maxFlow(adj_list, vertices, degree);
	clock_gettime(CLOCK_REALTIME, &start);
	for (int i = 0; i < 3054322; i++) {
		malloc(sizeof(struct vll));
	}
	#ifdef TIME
	clock_gettime(CLOCK_REALTIME, &end);
	
	double diff = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
	printf("TIME TO RUN: %lf\n", diff);
	#endif	
	
	printf("%d %d\n", num_out, max_flow);
	
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

	free((*in_range));
	free(in_range);
	free(adj_list);
	free(degree);
	free(groups);
	free(aps);
	free(walls);
*/
	return 0;
}


