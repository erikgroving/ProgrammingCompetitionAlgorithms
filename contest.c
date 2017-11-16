#include "contest.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a : b
#define BY_DIST 0
#define BY_SIZE 1
#define DEBUG

// Parse the input
void parseInput(struct group** groups, struct ap** aps, struct wall** walls,
				int* num_groups, int* num_ap, int* num_walls) {
					
	if (scanf("%d %d %d", num_groups, num_ap, num_walls) != 3) {
		printf("Incorrect input.\n");
		exit(1);
	}

	*groups = (struct group*)malloc(*num_groups * sizeof(struct group));
	*aps = (struct ap*)malloc(*num_ap * sizeof(struct ap));
	*walls = (struct wall*)malloc(*num_walls * sizeof(struct wall));

	for (int i = 0; i < *num_groups; i++) {
		if (scanf ("%d %d %d", &(*groups)[i].x, &(*groups)[i].y, &(*groups)[i].size) != 3){
			printf("Incorrect input.\n");
			exit(1);
		}
		(*groups)[i].head = NULL;
		(*groups)[i].num_aps_in_range = 0;
		(*groups)[i].origin_dist_sq = (*groups)[i].x * (*groups)[i].x + 
										(*groups)[i].y * (*groups)[i].y;
	}

	for (int i = 0; i < *num_ap; i++) {
		if (scanf("%d %d %d %d", &(*aps)[i].x, &(*aps)[i].y, &(*aps)[i].r, &(*aps)[i].c) != 4) {
			printf("Incorrect input\n");
			exit(1);
		}
		(*aps)[i].num_g_in_range = 0;
		double ap_dist = (*aps)[i].x * (*aps)[i].x + 
										(*aps)[i].y * (*aps)[i].y;
		ap_dist = sqrt(ap_dist);
		double min, max;
		min = ap_dist - (*aps)[i].r;
		min *= abs(min);
		max = ap_dist + (*aps)[i].r;
		max *= max;
		(*aps)[i].min = min;		
		(*aps)[i].max = (max + 0.5);
	}
	
	apquickSort(aps, 0, *num_ap);
	gquickSort(groups, 0, *num_groups, BY_DIST);
	
	
	for (int i = 0; i < *num_walls; i++) {
		if (scanf ("%d %d %d %d", &(*walls)[i].x1,
									&(*walls)[i].y1, 
									&(*walls)[i].x2,
									&(*walls)[i].y2) != 4) {
			printf("Incorrect input.\n");
			exit(1);
		}
		(*walls)[i].md = max (&(*walls)[i].x1 * &(*walls)[i].x1 + 
								&(*walls)[i].y1 * &(*walls)[i].y1,
								&(*walls)[i].x2 * &(*walls)[i].x2 + 
								&(*walls)[i].y2 * &(*walls)[i].y2);
	}
}

#ifdef DEBUG
int num_checks = 0;
#endif
int inRange(struct group g, struct ap aps, struct wall* walls, int num_walls) {
	// first check if point g is with the box instead of circle to avoid expensive
	// multiplication operations
	
	int x_dist = abs(aps.x - g.x);
	int y_dist = abs(aps.y - g.y);
	
	if (y_dist > aps.r || x_dist > aps.r) {
		return FALSE;
	}
	
	

	int dist_sq = (x_dist) * (x_dist) + (y_dist) * (y_dist);
	if (dist_sq <= aps.r * aps.r) {
		// Check to see if there is a wall
		wall g_to_ap;
		g_to_ap.x1 = g.x;
		g_to_ap.y1 = g.y;
		g_to_ap.x2 = aps.x;
		g_to_ap.y2 = aps.y;
		for (int i = 0; i < num_walls; i++) {
			#ifdef DEBUG
			num_checks++;
			#endif
			if(lineIntersect(g_to_ap, walls[i])) {
				return FALSE;
			}
		}
		return TRUE;
	}
	return FALSE;
}

int lineIntersect(struct wall a, struct wall b) {
	return (boxIntersect(a, b) &&
		isLeftAndRight(a, b) &&
		isLeftAndRight(b, a));	
}

int isLeftAndRight(struct wall a, struct wall b) {
	double pos_p1 = (b.x2 - b.x1) * (a.y1 - b.y1) - (b.y2 - b.y1) * (a.x1 - b.x1);
	double pos_p2 = (b.x2 - b.x1) * (a.y2 - b.y1) - (b.y2 - b.y1) * (a.x2 - b.x1);
	return (pos_p1 * pos_p2 <= 0);
}

int boxIntersect(struct wall a, struct wall b) {
	
	return ((min(a.x1, a.x2)) <= (max(b.x2, b.x1)) &&
			(max(a.x1, a.x2)) >= (min(b.x1, b.x2)) &&
			(min(a.y1, a.y2)) <= (max(b.y2, b.y1)) &&
			(max(a.y1, a.y2)) >= (min(b.y1, b.y2)));
}

// Returns a group pointer that contains all the groups that are in range
// of an access point.
void groupsInRange(struct group** in_range, struct group* groups, struct ap* aps, struct wall* walls,
							int num_groups, int num_ap, int num_walls, int* num_valid, int* num_out) {
	*num_valid = 0;
	*num_out = 0;
	int grp = 0;
	
	for (int j = 0; j < num_ap; j++) {
		// find search starting point
		while (grp < num_groups && groups[grp].origin_dist_sq < aps[j].min) {
			grp++;
		}
		for (int i = grp; i < num_groups; i++) {
			if (groups[i].origin_dist_sq > aps[j].max) {				
				break;
			}
			if (inRange(groups[i], aps[j], walls, num_walls)) {
				*num_valid += (groups[i].num_aps_in_range == 0);
				groups[i].num_aps_in_range++;
				aps[j].num_g_in_range++;
			
				// Create a new node
				struct ap_ll* node = (struct ap_ll*)malloc(sizeof(struct ap_ll));
				node->next = NULL;
				node->ap_idx = j;
				// Insert into the list
				if (groups[i].head == NULL) {
					groups[i].head = node;
					groups[i].tail = node;
				}
				else {
					groups[i].tail->next = node;
					groups[i].tail = groups[i].tail->next;
				}
			}
		}
	}

	(*in_range) = (struct group*)malloc((*num_valid) * sizeof(struct group));
	int j = 0;
	for (int i = 0; i < num_groups; i++) {
		if (groups[i].num_aps_in_range) {
			(*in_range)[j] = groups[i];
			j++;
		}
		else {
			*num_out += groups[i].size;
		}
	}
	
	gquickSort(in_range, 0, *num_valid, BY_SIZE);

	return;
}


void apquickSort(struct ap** aps, int l, int r) {
	if (l >= r - 1)
		return;
	struct ap tmp;
	struct ap* a = (*aps);
	int rand_piv = (rand() % (r - l)) + l;
	struct ap pivot = a[rand_piv];
	int wall = l;
	a[rand_piv] = a[l];
	a[l] = pivot;
	for (int i = l + 1; i < r; i++) {
		if (pivot.min > a[i].min) {
			wall++;
			tmp = a[i];
			a[i] = a[wall];
			a[wall] = tmp;
		}
	}
	tmp = a[wall];
	a[wall] = pivot;
	a[l] = tmp;
	aps = &a;
	apquickSort(aps, l, wall);
	apquickSort(aps, wall + 1, r);
	
}

void gquickSort(struct group** in_range, int l, int r, int spec) {
	if (l >= r - 1)
		return;
	struct group tmp;
	struct group* a = (*in_range);
	int rand_piv = (rand() % (r - l)) + l;
	struct group pivot = a[rand_piv];
	int wall = l;
	a[rand_piv] = a[l];
	a[l] = pivot;
	
	if (spec == BY_SIZE) {
		for (int i = l + 1; i < r; i++) {
			if (pivot.size < a[i].size) {
				wall++;
				tmp = a[i];
				a[i] = a[wall];
				a[wall] = tmp;
			}
		}
	}
	else {
		for (int i = l + 1; i < r; i++) {
			if (pivot.origin_dist_sq > a[i].origin_dist_sq) {
				wall++;
				tmp = a[i];
				a[i] = a[wall];
				a[wall] = tmp;
			}			
		}
	}
	tmp = a[wall];
	a[wall] = pivot;
	a[l] = tmp;
	in_range = &a;
	gquickSort(in_range, l, wall, spec);
	gquickSort(in_range, wall + 1, r, spec);
	return;
}

struct edge** createGraph(struct edge*** a, struct group* in_range, struct ap* aps, int** degree, int num_aps, int num_valid) {
	struct edge*** g = a;
	struct edge** graph;
	int* edge_idx = (int*)malloc(sizeof(int) * num_aps);
	graph = (struct edge**)malloc((2 + num_valid + num_aps) * sizeof(struct edge*));
	memset(*degree, 0, sizeof(int) *(num_valid + num_aps + 2));
	
	
	for (int i = 0; i < num_aps; i++) {
		edge_idx[i] = 1;
	}
	
	// source node and terminal node
	graph[0] = (struct edge*)malloc(sizeof(struct edge) * num_valid);
	graph[1 + num_valid + num_aps] = (struct edge*)malloc(sizeof(struct edge) * num_aps);

	// create edges from source to all student groups
	for (int i = 0; i < num_valid; i++) {
		graph[0][i].dest = i + 1;
		graph[0][i].cap = in_range[i].size;
	}
	(*degree)[0] = num_valid;

	// create edges from all aps to terminal
	for (int i = 0; i < num_aps; i++) {
		graph[1 + num_valid + i] = 
			(struct edge*)malloc(sizeof(struct edge) * (aps[i].num_g_in_range + 1));
		
		graph[1 + num_valid + i][0].dest = 1 + num_valid + num_aps;
		graph[1 + num_valid + i][0].cap = aps[i].c;
		(*degree)[1 + num_valid + i] = 1 + aps[i].num_g_in_range;
	}

	// create edges from student groups to acces points
	for (int i = 0; i < num_valid; i++) {
		struct ap_ll* sentinel = in_range[i].head;
		graph[i + 1] = 
			(struct edge*)malloc(sizeof(struct edge) * 
			in_range[i].num_aps_in_range);

		for (int j = 0; j < in_range[i].num_aps_in_range; j++) {
			int ap = sentinel->ap_idx;
			int ap_v = ap + 1 + num_valid;
			graph[i + 1][j].dest = ap_v;
			
			graph[i + 1][j].cap = aps[ap].c;
			
			graph[ap_v][edge_idx[ap]].dest = i + 1;
			graph[ap_v][edge_idx[ap]].cap = 0;
			edge_idx[ap]++;
			sentinel = sentinel->next;

		}

		(*degree)[i + 1] = in_range[i].num_aps_in_range;

	}

	g = &graph;
	free(edge_idx);
	return (*g);
}

#ifdef DEBUG 
int nodes_popped = 0;
int nodes_pushed = 1;
int calls = 0;
#endif
int maxFlow(struct edge** adj_list, int vertices, int* degree, int num_groups, int num_walls) {
	int flow = 0;
	int prev_deg = 0;
	int end_deg = degree[0];
	int inc = ((num_walls << 1) >= num_groups) ? 20 : 1	;
	degree[0] = inc % (end_deg + 1);
	do {
		prev_deg = max(0, degree[0] - inc);

		degree[0] = min((degree[0] + inc), end_deg);

		#ifdef DEBUG
		while(findPath(&adj_list, vertices, degree, prev_deg, &flow, &end_deg)){
			calls++;
		} 
		#else
		while(findPath(&adj_list, vertices, degree, prev_deg, &flow, &end_deg));
		#endif
	} while(degree[0] != end_deg);
	#ifdef DEBUG
	printf("popped: %d\tpushed: %d\t ratio: %f\tcalls: %d\tpushes/call: %f\tpops/call: %f\npops+push/call: %f\n",
		nodes_popped, nodes_pushed, (float)nodes_pushed/nodes_popped, 
		calls, (float)nodes_pushed/calls, (float)nodes_popped/calls, (float)(nodes_popped+nodes_pushed)/calls);
	
	printf("num checks: %d\n", num_checks);
	#endif
	return flow;
}

int findPath(struct edge*** adj_list_tp, int vertices, int* degree, int prev_deg, int* flow, int* end_deg) {
	struct edge** adj_list = (*adj_list_tp);
	struct stack dfs;
	int* parent = (int *)malloc(sizeof(int) * vertices);
	memset(parent, -1, vertices * sizeof(int));
	parent[0] = -2;

	int term = vertices - 1;
	

	// Perform the depth first searchpath
	for (int m = degree[0] - 1; m >= prev_deg; m--) {
		int dest_pre = adj_list[0][m].dest;
		int cap_pre = adj_list[0][m].cap;
		dfs.head = NULL;
		if (parent[dest_pre] != -1 || cap_pre == 0) {
			continue;
		}
		#ifdef DEBUG
		nodes_pushed++;
		#endif
		struct vll* new_head = (struct vll*)malloc(sizeof(struct vll));
		new_head->v = adj_list[0][m].dest;
		new_head->next = NULL;
		parent[dest_pre] = 0;
		dfs.size = 1;
		dfs.head = new_head;
		while (dfs.size) {
			#ifdef DEBUG
			nodes_popped++;
			#endif
			// take the front of stack and pop
			int v = dfs.head->v;
			struct vll* tmp = dfs.head->next;

			free(dfs.head);
			dfs.head = tmp;
			dfs.size--;
		
			
			// prepend to the linked list all elements 
			// in range of current node
			for (int i = 0; i < degree[v]; i++) {
				int dest = adj_list[v][i].dest;
				int cap = adj_list[v][i].cap;
				
				if (parent[dest] != -1 || cap == 0) {
					continue;
				}
				if (dest == term) {
					int curr_node = dest;
					int bottleneck = INT_MAX;
					parent[curr_node] = v;
					// find the bottleneck
					while(curr_node != -2) {
						int from = parent[curr_node];
						if(from == -2) break;

						int to = curr_node;
						curr_node = parent[curr_node];
		
						// find path from parent to curr
						for (int k = 0; k < degree[from]; k++) {
							if(adj_list[from][k].dest == to) {
								bottleneck = min(bottleneck, adj_list[from][k].cap);
								break;
							}
						}
					}
					//subtract the bottleneck in forward flow
					// and add in backward flow
					curr_node = dest;
					while(parent[curr_node] != -2) {
						int from = parent[curr_node];
						int to = curr_node;
						// find path from parent to curr
						for (int k = 0; k < degree[from]; k++) {
							if(adj_list[from][k].dest == to) {
								if (parent[from] == -2 && adj_list[from][k].cap == bottleneck) {
									struct edge tmp = adj_list[from][k];
									adj_list[from][k] = adj_list[from][*end_deg - 1];
									adj_list[from][*end_deg - 1] = tmp;
									degree[0]--;
									(*end_deg)--;
								}
								else {
									adj_list[from][k].cap -= bottleneck;
									if (to != term) {
										for (int l = 0; l < degree[to]; l++) {
											if (adj_list[to][l].dest == from) {
												adj_list[to][l].cap += bottleneck;
												break;
											}
										}
									}
								}
								break;
							}
						}
						curr_node = parent[curr_node];
					}
					
					(*flow) += bottleneck;
					adj_list_tp = &adj_list;
					
					free(parent);
					if (dfs.head != NULL) {
						while (dfs.head->next != NULL) {
							struct vll* s1 = dfs.head;
							dfs.head = dfs.head->next;
							free(s1);
						}
						free(dfs.head);
					}
					return 1;
				}


				// push on the stack
				struct vll* node = (struct vll*)malloc(sizeof(struct vll));
				node->v = dest;
			
			
				// if empty, set head and tail
				/*if (dfs.size == 0) {
					node->next = NULL;
					dfs.head = node;
					dfs.tail = node;
				}*/
				// prepend if non-augmenting
				//else{ /*if (dest > v) */
					node->next = dfs.head;
					dfs.head = node;
					#ifdef DEBUG
					nodes_pushed++;
					#endif
				//}
				

				// append to tail if augmenting edge
				/*else {
					node->next = NULL;
					dfs.tail->next = node;
					dfs.tail = dfs.tail->next;
				}*/
				
				dfs.size++;
				parent[dest] = v;
				
				// check if the node has a connection to the terminal, if it does, we can stop looping			
				// or if the node's node has a connection to the terminal.
				int dest2 = adj_list[dest][0].dest;
				int cap2 = adj_list[dest][0].cap;
				if (dest2 == term && cap2 > 0) {
					i = degree[v];
				}
				
				
				
			}
		}
	}

	free(parent);
	return 0;
}