#include "contest.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#define BY_DIST 0
#define BY_SIZE 1
//#define DEBUG

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
		(*aps)[i].r_sq = (*aps)[i].r * (*aps)[i].r;
	}
	
	
	for (int i = 0; i < *num_walls; i++) {
		if (scanf ("%f %f %f %f", &(*walls)[i].x1,
									&(*walls)[i].y1, 
									&(*walls)[i].x2,
									&(*walls)[i].y2) != 4) {
			printf("Incorrect input.\n");
			exit(1);
		}
	}
}

void makeApBuckets(int num_aps, struct ap* aps, struct abucket a[(APBUCK << APSHIFT)], float incx, float incy) {
	for (int i = 0; i < APBUCK; i++) {
		for (int j = 0; j < APBUCK; j++) {
			a[(i << APSHIFT) + j].size = 0;
			a[(i << APSHIFT) + j].allocd = 16;
			a[(i << APSHIFT) + j].w = malloc(sizeof(int) * 16);
		}
	}	
	
	for (int k = 0; k < num_aps; k++) {
		float neg_x = aps[k].x - aps[k].r;
		float max_x = aps[k].x + aps[k].r;
		float neg_y = aps[k].y - aps[k].r;
		float max_y = aps[k].y + aps[k].r;
		float min_x = max(0, neg_x);
		float min_y = max(0, neg_y);
		
		int min_ax = min(min_x / incx, APBUCK - 1);
		int min_ay = min(min_y / incy, APBUCK - 1);	
		int max_ax = min(max_x / incx, APBUCK - 1);	
		int max_ay = min(max_y / incy, APBUCK - 1);	
		for (int i = min_ax; i <= max_ax; i++) {
			for (int j = min_ay; j <= max_ay; j++) {
				a[(i << APSHIFT) + j].w[a[(i << APSHIFT) + j].size] = k;
				a[(i << APSHIFT) + j].size++;
				if (a[(i << APSHIFT) + j].size == a[(i << APSHIFT) + j].allocd) {
					a[(i << APSHIFT) + j].allocd += a[(i << APSHIFT) + j].allocd;
					int* tmpw = malloc(sizeof(int) * a[(i << APSHIFT) + j].allocd);
					for (int k = 0; k < a[(i << APSHIFT) + j].size; k++) {
						tmpw[k] = a[(i << APSHIFT) + j].w[k];
					}
					free(a[(i << APSHIFT) + j].w);
					a[(i << APSHIFT) + j].w = tmpw;
				}
			}
		}
	}
	
	
}


void makeBuckets(int num_walls, struct wall* walls, struct bucket b[(BUCK << BSHIFT)], float incx, float incy) {
	// Making buckets
	for (int i = 0; i < BUCK; i++) {
		for (int j = 0; j < BUCK; j++) {
			b[(i << BSHIFT) + j].size = 0;
			struct wall tmp;
			tmp.x1 = (i * incx);
			tmp.y1 = (j * incy);
			tmp.x2 = tmp.x1 + incx;
			tmp.y2 = tmp.y1 + incy;
			b[(i << BSHIFT) + j].region = tmp;
			b[(i << BSHIFT) + j].allocd = 16;
			b[(i << BSHIFT) + j].w = malloc(sizeof(int) * 16);
		}
	}
	
	// Finding walls in each bucket
	for (int i = 0; i < BUCK; i++) {
		for (int j = 0; j < BUCK; j++) {
			struct bucket tmpb = b[(i << BSHIFT) + j];
			for (int k = 0; k < num_walls; k++) {
				if (boxIntersect(walls[k], tmpb.region)) {
					tmpb.w[tmpb.size] = k;
					tmpb.size++;
					if (tmpb.size == tmpb.allocd) {
						tmpb.allocd += tmpb.allocd;
						int* tmpw = malloc(sizeof(int) * tmpb.allocd);
						for (int k = 0; k < tmpb.size; k++) {
							tmpw[k] = tmpb.w[k];
						}
						free(tmpb.w);
						tmpb.w = tmpw;
					}
				}				
			}
			b[(i << BSHIFT) + j] = tmpb;
		}
	}
}





#ifdef DEBUG
int num_checks = 0;
int box_intersect = 0;
int num_block = 0;
int num_ok = 0;
#endif

int inRange(struct group g, struct ap aps, struct wall* walls, struct bucket b[(BUCK << BSHIFT)], int num_walls, float incx, float incy) {
	// first check if point g is with the box instead of circle to avoid expensive
	// multiplication operations
	
	int x_dist = abs(aps.x - g.x);
	int y_dist = abs(aps.y - g.y);
	
	if (y_dist > aps.r || x_dist > aps.r) {
		return FALSE;
	}
	
	

	int dist_sq = (x_dist) * (x_dist) + (y_dist) * (y_dist);

	if (dist_sq <= aps.r_sq && num_walls) {
		// Check to see if there is a wall
		struct wall g_to_ap;
		g_to_ap.x1 = g.x;
		g_to_ap.y1 = g.y;
		g_to_ap.x2 = aps.x;
		g_to_ap.y2 = aps.y;
		
		float max_y = (incy * BUCK);
		float max_x = (incx * BUCK);
					
		if (g_to_ap.x2 > max_x && g_to_ap.x1 > max_x) {
			return TRUE;
		}
		if (g_to_ap.y1 > max_y && g_to_ap.y2 > max_y) {
			return TRUE;
		}	
		
		// First, find the direction of the line!
		int x_dir = g_to_ap.x2 > g_to_ap.x1 ? 1 : -1;
		int y_dir = g_to_ap.y2 > g_to_ap.y1 ? 1 : -1;
		
	

		
		float slope = (g_to_ap.y2 - g_to_ap.y1) / (g_to_ap.x2 - g_to_ap.x1);
		float intercept = g_to_ap.y2 - (slope * g_to_ap.x2);

		
		// find the starting bucket
		int bx = min(g_to_ap.x1 / incx, BUCK - 1);
		int by = min(g_to_ap.y1 / incy, BUCK - 1);
		// then find the end bucket
		int final_bx = min(g_to_ap.x2 / incx, BUCK - 1);
		int final_by = min(g_to_ap.y2 / incy, BUCK - 1);

		int reached = 0;
		do {
			//printf("bx: %d\tby:%d\tfinalbx: %d\tfinal_by: %d\n", bx, by, final_bx, final_by);
			reached = (bx == final_bx && by == final_by);			
			
			// check current bucket
			struct bucket curr_b = b[(bx << BSHIFT) + by];
			
			for (int i = 0; i < curr_b.size; i++) {

				if (lineIntersect(g_to_ap, walls[curr_b.w[i]])) {		
					#ifdef DEBUG
					num_block++;
					#endif
					return FALSE;
				}
				#ifdef DEBUG
				num_ok++;
				#endif

			}
			if (bx == final_bx) {
				by += y_dir;
			}
			else if (by == final_by) {
				bx += x_dir;
			}
			else  {
				bx += x_dir;
				struct wall tmp;
				float inc = (x_dir == 1) ? incx : -incx;
				tmp.x1 = curr_b.region.x1 + inc;
				tmp.x2 = curr_b.region.x2 + inc;
				tmp.y1 = slope * tmp.x1 + intercept;
				tmp.y2 = slope * tmp.x2 + intercept;
				#ifdef DEBUG
				box_intersect++;
				#endif	
				if (!boxIntersect(tmp, curr_b.region)) {
					bx -= x_dir;
					by += y_dir;
				}
			}
		} while (!reached);
/*
		for (int i = 0; i < num_walls; i++) {
	
			#ifdef DEBUG
			num_checks++;
			#endif
			if(lineIntersect(g_to_ap, walls[i])) {
				//free(seen);
				return FALSE;
			}
		}*/
		return TRUE;
	}
	else if (dist_sq <= aps.r_sq ) {
		return TRUE;
	}
	return FALSE;
}

int lineIntersect(struct wall a, struct wall b) {
	return ( boxIntersect(a, b) && isLeftAndRight(a, b) &&
		isLeftAndRight(b, a));	
}

int isLeftAndRight(struct wall a, struct wall b) {
	double pos_p1 = (b.x2 - b.x1) * (a.y1 - b.y1) - (b.y2 - b.y1) * (a.x1 - b.x1);
	double pos_p2 = (b.x2 - b.x1) * (a.y2 - b.y1) - (b.y2 - b.y1) * (a.x2 - b.x1);
	return (pos_p1 * pos_p2 <= 0);
}

int boxIntersect(struct wall a, struct wall b) {
	
	return (min(a.x1, a.x2)) <= (max(b.x2, b.x1)) &&
			(max(a.x1, a.x2)) >= (min(b.x1, b.x2)) &&
			(min(a.y1, a.y2)) <= (max(b.y2, b.y1)) &&
			(max(a.y1, a.y2)) >= (min(b.y1, b.y2));
}

// Returns a group pointer that contains all the groups that are in range
// of an access point.
void groupsInRange(struct group** in_range, struct group* groups,  struct abucket a[(APBUCK << APSHIFT)], struct bucket b[(BUCK << BSHIFT)], struct ap* aps, struct wall* walls,
							int num_groups, int num_ap, int num_walls, int* num_valid, int* num_out, float incx, float incy, float ap_incx, float ap_incy) {
	*num_valid = 0;
	*num_out = 0;
	//int grp = 0;
	
	for (int i = 0; i < num_groups; i++) {
		// find search starting point
		/*while (grp < num_groups && groups[grp].origin_dist_sq < aps[j].min) {
			grp++;
		}
		for (int i = grp; i < num_groups; i++) {
			if (groups[i].origin_dist_sq > aps[j].max) {				
				break;
			}
			
		*/
		int ax = min((float)groups[i].x / ap_incx, APBUCK - 1);
		int ay = min((float)groups[i].y / ap_incy, APBUCK - 1);		
		for (int j = 0; j < a[(ax << APSHIFT) + ay].size; j++) {
			int ap_idx = a[(ax << APSHIFT) + ay].w[j];
			#ifdef DEBUG
			num_checks++;
			#endif
			if (inRange(groups[i], aps[ap_idx], walls, b, num_walls, incx, incy)) {
				*num_valid += (groups[i].num_aps_in_range == 0);
				groups[i].num_aps_in_range++;
				aps[ap_idx].num_g_in_range++;
			
				// Create a new node
				struct ap_ll* node = (struct ap_ll*)malloc(sizeof(struct ap_ll));
				node->ap_idx = ap_idx;
				node->next = NULL;
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
	
	//gquickSort(in_range, 0, *num_valid, BY_SIZE);

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
			graph[i + 1][j].mirror = &(graph[ap_v][edge_idx[ap]]);
			
			
			graph[ap_v][edge_idx[ap]].dest = i + 1;
			graph[ap_v][edge_idx[ap]].cap = 0;
			graph[ap_v][edge_idx[ap]].mirror = &(graph[i + 1][j]);

			edge_idx[ap]++;
			sentinel = sentinel->next;

		}

		

		(*degree)[i + 1] = in_range[i].num_aps_in_range;
	}
	/*
	for (int i = 0; i < num_valid; i++) {
		for (int j = 0; j < in_range[i].num_aps_in_range; j++) {
			int min_r = aps[graph[i + 1][j].dest - (1 + num_valid)].r;
			int min_idx = j;
			for (int k = j + 1; k < in_range[i].num_aps_in_range; k++) {
				if (aps[(graph[i + 1][k].dest - (1 + num_valid))].r < min_r) {
					min_idx = k;
					min_r = aps[graph[i + 1][k].dest - (1 + num_valid)].r;
				}
				steps++;
			}
			graph[i + 1][min_idx].mirror->mirror = &graph[i + 1][j];
			graph[i + 1][j].mirror->mirror = &graph[i + 1][min_idx];
			edge tmp = graph[i + 1][min_idx];
			graph[i + 1][min_idx]  = graph[i + 1][j];
			graph[i + 1][j] = tmp;
		}
	}*/

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
	int ap_start = degree[0] + 1;
	// Greedy ford fulkerson! No augmenting paths

	for(int i = 1; i <= degree[0]; i++) {
		int* group_cap = &(adj_list[0][i - 1].cap);
		//printf("node: %d\n", i);
		struct edge* a = adj_list[i];
		for (int j = 0; j < degree[i]; j++) {			
			int dest = a[j].dest;
			if (adj_list[dest][0].dest == vertices - 1 && adj_list[dest][0].cap) {
				int* ap_cap = &(adj_list[dest][0].cap);
				int bottleneck = min(*group_cap, *ap_cap);
				flow += bottleneck;
				//printf("g: %d ap: %d b: %d\n",*group_cap, *ap_cap, bottleneck);
				a[j].cap -= bottleneck;
				a[j].mirror->cap += bottleneck;
				
				(*ap_cap) -= bottleneck;
				(*group_cap) -= bottleneck;
				
				if (*group_cap == 0) {
					break;
				}				
			}
		}
	}
	
	for (int i = 0; i < degree[0]; i++) {
		if (!adj_list[0][i].cap) {
			degree[0]--;
			adj_list[0][i] = adj_list[0][degree[0]];
		}
	}


	int prev_deg = 0;
	int end_deg = degree[0];
	//int inc = 1;//((num_walls << 1) >= num_groups) ? 1 : 1	;
	degree[0] = 0;

	do {		
		prev_deg = degree[0] - 1;

		degree[0] = degree[0] + 1;
		#ifdef DEBUG
		while(findPath(&adj_list, vertices, ap_start, degree, prev_deg, &flow, &end_deg)){
			calls++;
		}
		calls++;
		#else
		while(findPath(&adj_list, vertices, ap_start, degree, prev_deg, &flow, &end_deg));
		#endif
	} while(degree[0] < end_deg);
	#ifdef DEBUG
	printf("popped: %d\npushed: %d\n ratio: %f\ncalls: %d\n\tpushes/call: %f\tpops/call: %f\npops+push/call: %f\n",
		nodes_popped, nodes_pushed, (float)nodes_pushed/nodes_popped, 
		calls, (float)nodes_pushed/calls, (float)nodes_popped/calls, (float)(nodes_popped+nodes_pushed)/calls);
	
	printf("box_intersects: %d\n", box_intersect);
	printf("CHECKS: %d\nnum_ok: %d\tnum_block: %d\n", num_checks, num_ok, num_block);
	#endif
	return flow;
}

int findPath(struct edge*** adj_list_tp, int vertices, int ap_start, int* degree, int prev_deg, int* flow, int* end_deg) {
	struct edge** adj_list = (*adj_list_tp);
	struct stack dfs;
	int* parent = (int *)malloc(sizeof(int) * vertices);
	memset(parent, -1, vertices * sizeof(int));
	parent[0] = -2;

	int term = vertices - 1;
	
	struct vll* nodes = (struct vll*)malloc(sizeof(struct vll) * vertices);
	int n_idx = 0;
	
	// Perform the depth first searchpath
	for (int m = degree[0] - 1; m > prev_deg; m--) {
		int dest_pre = adj_list[0][m].dest;
		int cap_pre = adj_list[0][m].cap;
		if (parent[dest_pre] != -1 || cap_pre == 0) {
			continue;
		}
		#ifdef DEBUG
		nodes_pushed++;
		#endif

		nodes[n_idx].v = adj_list[0][m].dest;
		nodes[n_idx].next = NULL;
		parent[dest_pre] = 0;
		dfs.size = 1;
		dfs.head = &nodes[n_idx];
		n_idx++;
		dfs.tail = dfs.head;
		while (dfs.size) {
			#ifdef DEBUG
			nodes_popped++;
			#endif
			// take the front of stack and pop
			int v = dfs.head->v;
			dfs.head = dfs.head->next;
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
									adj_list[from][k] = adj_list[from][*end_deg - 1];
									degree[0]--;
									prev_deg--;
									(*end_deg)--;
								}
								else {
									adj_list[from][k].cap -= bottleneck;
									if (from != 0 && to != term) {
										adj_list[from][k].mirror->cap += bottleneck;
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
					free(nodes);
					return 1;
				}




				
				
				// check if the node has a connection to the terminal, if it does, we can stop looping			
				// or if the node's node has a connection to the terminal.
				int dest2 = adj_list[dest][0].dest;
				int cap2 = adj_list[dest][0].cap;
				int found = (dest2 == term && cap2 > 0);
				if (found) {
					i = degree[v];
				}		

				// push on the stack
				nodes[n_idx].v = dest;

				if (!dfs.size) {
					nodes[n_idx].next = NULL;
					dfs.head = &nodes[n_idx];
					dfs.tail = dfs.head;
				}
				else if (found) {
					nodes[n_idx].next = dfs.head;
					dfs.head = &nodes[n_idx];
				}
				else {
					nodes[n_idx].next = NULL;
					dfs.tail->next = &nodes[n_idx];
					dfs.tail = dfs.tail->next;
				}
				n_idx++;
				#ifdef DEBUG
				nodes_pushed++;
				#endif
				
				parent[dest] = v;
				dfs.size++;				
			}
		}
	}

	free(parent);
	free(nodes);
	return 0;
}