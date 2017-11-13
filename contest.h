#ifndef CONTEST_H_
#define CONTEST_H_

#define TRUE 1
#define FALSE 0

typedef struct stack {	
	struct vll* head;
	struct vll* tail;
	int size;
} stack;

// vertex linked list
typedef struct vll {	
	struct vll* next;
	int v;
} vll;

typedef struct edge {
	int dest;
	int cap;
	int risk;
} edge;


typedef struct ap {	
	int x;
	int y;
	int r;
	int num_g_in_range;
	int min;
	int max;
	int c;
} ap;

typedef struct group {
	struct ap_ll* head;
	struct ap_ll* tail;
	int x;
	int y;
	int size;
	int num_aps_in_range;
	int origin_dist_sq;
} group;

typedef struct ap_ll {	
	struct ap_ll* next;
	int ap_idx;
} ap_ll;
	

typedef struct wall {
	int x1;
	int y1;
	int x2;
	int y2;
} wall;



void parseInput(struct group**, struct ap**, struct wall**, int* , int* , int*);
int inRange(struct group, struct ap, struct wall*, int);
int lineIntersect(struct wall, struct wall);
int isLeftAndRight(struct wall, struct wall);
int boxIntersect(struct wall, struct wall);
void groupsInRange(struct group**, struct group*, struct ap* , struct wall* ,
							int, int, int, int*, int*);
void apquickSort(struct ap**, int, int);
void gquickSort(struct group**, int, int, int);
struct edge** createGraph(struct edge***, struct group*, struct ap*, int**, int, int);
int maxFlow(struct edge**, int, int*);
int findPath(struct edge*** adj_list_tp, int vertices, int* degree, int* flow, int*) ;
#endif /* CONTEST_H_ */
