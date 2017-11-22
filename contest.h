#ifndef CONTEST_H_
#define CONTEST_H_

#define TRUE 1
#define FALSE 0
#define BUCK 8
#define BSHIFT 3
#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a : b

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
	struct edge* mirror;
} edge;


typedef struct ap {	
	int x;
	int y;
	int r;
	int r_sq;
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
	float x1;
	float y1;
	float x2;
	float y2;
} wall;
typedef struct bucket {
	int* w;
	int size;
	struct wall region;
}bucket;



void parseInput(struct group**, struct ap**, struct wall**, int* , int* , int*);
void makeBuckets(int num_walls, struct wall* walls, struct bucket buckets[BUCK * BUCK], float, float);
int inRange(struct group, struct ap, struct wall*, struct bucket b[BUCK * BUCK], int, float, float);
int lineIntersect(struct wall, struct wall);
int isLeftAndRight(struct wall, struct wall);
int boxIntersect(struct wall, struct wall);
void groupsInRange(struct group**, struct group*, struct bucket b[BUCK * BUCK], struct ap* , struct wall* ,
							int, int, int, int*, int*, float, float);
//void wquickSort(struct wall**, int, int);
void apquickSort(struct ap**, int, int);
void gquickSort(struct group**, int, int, int);
struct edge** createGraph(struct edge***, struct group*, struct ap*, int**, int, int);
int maxFlow(struct edge**, int, int*, int, int);
int findPath(struct edge*** adj_list_tp, int vertices, int, int* degree, int, int* flow, int*, int*) ;
#endif /* CONTEST_H_ */
