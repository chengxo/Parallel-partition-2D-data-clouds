#ifndef P4_TO_P8
#include <p4est.h>
#include <string.h>
#include <p4est_vtk.h>
#include <p4est_bits.h>
#include <p4est_base.h>
#include <p4est_algorithms.h>
#include <p4est_extended.h>
#include <p4est_iterate.h>
#include <p4est_io.h>
#include <p4est_communication.h>
#else
#include <p8est_vtk.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>
#include <p8est_io.h>
#include <p8est_communication.h>
#endif
#include <sc.h>
#include <sc_options.h>
#include <sc_mpi.h>

//*** make sure stdlib.h is included>

/** The resolution of the image data in powers of two. */
#define P4EST_STEP1_PATTERN_LEVEL 5
/** The dimension of the image data. */
#define P4EST_STEP1_PATTERN_LENGTH (1 << P4EST_STEP1_PATTERN_LEVEL)
#define BDIM 2
#define REFINE_THRES 1
static const int    plv = P4EST_STEP1_PATTERN_LEVEL;    /**< Shortcut */
static const int    ple = P4EST_STEP1_PATTERN_LENGTH;   /**< Shortcut */
#ifdef P4_TO_P8
static const p4est_qcoord_t eighth = P4EST_QUADRANT_LEN (3);
#endif

typedef struct bound_primitive
{
  int u[BDIM]; // global, u[0] <= u[1] (= if 1 point; u[1] = u[0] - 1 if empty, need to deal with u[1] = u[0] - 1 case)
}
bd_t;

typedef struct local_info
{
  double *pts; // array of points in the proc (1d, (x1, y1, x2, y2...))
  int *cumsum_num_pt; // cumsum of number of points in each proc; length is mpisize, first entry is 0.
  // this is true num, without *P4EST_DIM; 
}
local_info_t;

int log2_diy(int x){
  int ret = 0;
  while (x >>= 1){
    ret = ret + 1;
  }
  return ret;
}

int initial_refine_func(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
  return 1;
}

static int weight_func(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
  bd_t *data;
  data = (bd_t *) q->p.user_data;
  return (data->u[1] - data->u[0] + 1);
}

static int refine_func(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
  bd_t *data = (bd_t *)q->p.user_data;
  return (((data->u[1] - data->u[0] + 1)>(REFINE_THRES))? 1 : 0);
}

//Does the distribution of quadrants follow the mpi rank?


// is unit square [0, 1] * [0, 1] in p4est?
void get_local_range(int total_num, int my_rank, int my_size, double *range, char *filename){
  FILE *in_file = fopen(filename, "r");
  double curr_read[2];
  int offset = 0;//((total_num / my_size) * my_rank) * ((sizeof(double) + sizeof(int)) * P4EST_DIM);
  //printf("fseek\n");
  //fseek(in_file, offset, SEEK_SET);
  //printf("fscanf\n");
  fscanf(in_file, "%lf\n%lf\n", &curr_read[0], &curr_read[1]);
  //range 0, 1, 2, 3: left, right, bottom, top
  range[0] = curr_read[0];
  range[1] = curr_read[0];
  range[2] = curr_read[1];
  range[3] = curr_read[1];

  while(fscanf(in_file, "%lf\n%lf\n", &curr_read[0], &curr_read[1]) == 2){
    (range[0] > curr_read[0]) ? (range[0] = curr_read[0]) : (range[0] = range[0]);
    (range[1] < curr_read[0]) ? (range[1] = curr_read[0]) : (range[1] = range[1]);
    (range[2] > curr_read[1]) ? (range[2] = curr_read[1]) : (range[2] = range[2]);
    (range[3] < curr_read[1]) ? (range[3] = curr_read[1]) : (range[3] = range[3]);
  }
  fclose(in_file);
}


void transform_len_fromgridtodata(p4est_qcoord_t len, double *range, double *transformed_len){
  transformed_len[0] = (len * 1.0 / P4EST_ROOT_LEN) * (range[1] - range[0]);
  transformed_len[1] = (len * 1.0/ P4EST_ROOT_LEN) * (range[3] - range[2]);
}

void transform_loc_fromdatatogrid(double *ret, double *loc, double *range){
  ret[0] = ((loc[0] - range[0]) / (range[1] - range[0])) * (P4EST_ROOT_LEN); // confirm this choice... too large?
  ret[1] = ((loc[1] - range[2]) / (range[3] - range[2])) * (P4EST_ROOT_LEN);
}

void transform_loc_fromgridtodata(double *ret, p4est_qcoord_t *loc, double *range){
  ret[0] = ((loc[0] * 1.0 / P4EST_ROOT_LEN) * (range[1] - range[0])) + range[0];
  ret[1] = ((loc[1] * 1.0 / P4EST_ROOT_LEN) * (range[3] - range[2])) + range[2];
}

void transform_loc_fromgridtodata_double(double *ret, double *loc, double *range){
  ret[0] = ((loc[0] * 1.0 / P4EST_ROOT_LEN) * (range[1] - range[0])) + range[0];
  ret[1] = ((loc[1] * 1.0 / P4EST_ROOT_LEN) * (range[3] - range[2])) + range[2];
}


// check this... seems that printf is different from loc_xy...
int count_num(p4est_qcoord_t *loc, p4est_qcoord_t len, double *range, char *filename){
  FILE *in_file = fopen(filename, "r");
  double curr_read[2];
  double transformed_loc[2];
  double transformed_len[2];
  //printf("before transform, local_xy: %lf, %lf\n", loc[0], loc[1]);
  transform_loc_fromgridtodata(transformed_loc, loc, range);
  //printf("transformed_loc: %lf, %lf\n", transformed_loc[0], transformed_loc[1]);
  transform_len_fromgridtodata(len, range, transformed_len);
  //printf("transformed_len: %lf, %lf\n", transformed_len[0], transformed_len[1]);
  int count = 0;
  //printf("local_xy:%lf, %lf\n", loc[0], loc[1]);
  // printf("original loc: %lf, %lf\n", loc[0], loc[1]);
  // printf("transformed_loc: %lf, %lf\n", transformed_loc[0], transformed_loc[1]);
  // printf("transformed_len: %lf, %lf\n", transformed_len[0], transformed_len[1]);
  while(fscanf(in_file, "%lf\n%lf\n", &curr_read[0], &curr_read[1]) == 2){
    if ((curr_read[0] >= transformed_loc[0]) && (curr_read[0] < transformed_loc[0] + transformed_len[0]) && (curr_read[1] > transformed_loc[1]) && (curr_read[1] <= transformed_loc[1] + transformed_len[1])){
     // printf("loc%lf %lf cur_read in count_num is %lf, %lf\n", loc[0], loc[1], curr_read[0], curr_read[1]);
      count += 1;
    }
  }
  fclose(in_file);
  return count;
}

p4est_tree_t * p4est_tree_array_index_diy(sc_array_t * array, p4est_topidx_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p4est_tree_t));
  // if (it < 0){
  //   printf("it is %d\n", it);
  // }
  // if ( (size_t) it >= array->elem_count){
  //   printf("elem_count is %d\n", array->elem_count);
  // }
  P4EST_ASSERT (it >= 0 && (size_t) it < array->elem_count);
  return (p4est_tree_t *) (array -> array +  sizeof (p4est_tree_t) * (size_t) it);
}

p4est_quadrant_t * p4est_quadrant_array_index_diy(sc_array_t * array, size_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p4est_quadrant_t));
  P4EST_ASSERT (it < array->elem_count);
  return (p4est_quadrant_t *) (array -> array + sizeof (p4est_quadrant_t) * it);
}

p4est_quadrant_t * get_first_quadrant(p4est_t * p4est)
{
  p4est_topidx_t first_local_tree = p4est->first_local_tree;
  sc_array_t *trees = p4est->trees, *quadrants;
  p4est_tree_t *tree = p4est_tree_array_index_diy (trees, first_local_tree);
  // if (first_local_tree == -1){
  //   printf("get_first_quadrant_first_local_tree is -1\n");
  // }
  p4est_quadrant_t *q;
  quadrants = &(tree->quadrants);
  q = p4est_quadrant_array_index_diy (quadrants, 0);
  return q;
}


void get_pts(p4est_qcoord_t *loc, p4est_qcoord_t len, double *range, char *filename, double *pts, int num_pts){
  FILE *in_file = fopen(filename, "r");
  double curr_read[2];
  int count = 0;
  while((count < num_pts) &&  (fscanf(in_file, "%lf\n%lf\n", &curr_read[0], &curr_read[1]) == 2)){
    transform_loc_fromdatatogrid(curr_read, curr_read, range);
    if ((curr_read[0] >= loc[0]) && (curr_read[0] < loc[0] + len) && (curr_read[1] > loc[1]) && (curr_read[1] <= loc[1] + len)){
      pts[count * 2] = curr_read[0];
      pts[count * 2 + 1] = curr_read[1];
      count += 1;
    }
  }
  fclose(in_file);
}


// scan the points in a quadrant, order them in the original array (with a temp
// array in the intermediate steps), output the numbers of points in each children quadrant
void quadrant_scan_and_order(p4est_t * p4est, p4est_quadrant_t *q, int *children_num_pts)
{
  //children_num_pts: length 4 array, with allocated memory, initialized to 0.
  bd_t *data = (bd_t *)q->p.user_data;
  int local_num_pt = data->u[1] - data->u[0] + 1;
  double *temp_pts = (double *)malloc(local_num_pt * (sizeof(double)*P4EST_DIM));
  local_info_t *linfo = (local_info_t*) p4est->user_pointer;
  double *local_pts = linfo->pts;
  int offset;
  if (p4est->mpirank == 0){
    offset = 0;
  }
  else{
    offset = linfo->cumsum_num_pt[p4est->mpirank - 1];// ** take care, since the id start from 0.
  }
  p4est_qcoord_t half_length = P4EST_QUADRANT_LEN(q->level)/2;
  p4est_qcoord_t x_mid = q->x + half_length;
  p4est_qcoord_t y_mid = q->y + half_length;
  int start_id = (data->u[0] - offset);
  int end_id = (data->u[1] - offset);
  int cur_child_id;
  int cur_bulk_child_id;
  int *children_id = (int *)malloc(sizeof(int) * local_num_pt);
  int current_children_id[1<<P4EST_DIM];
  int cid;

  for (int i = start_id; i <= end_id; i++){
    cur_child_id = 0;
    (local_pts[i*P4EST_DIM] >= x_mid)? (cur_child_id += 1) : (cur_child_id = cur_child_id);
    (local_pts[i*P4EST_DIM+1] >= y_mid)? (cur_child_id += 2) : (cur_child_id = cur_child_id);
    children_num_pts[cur_child_id] += 1;
    children_id[i-start_id] = cur_child_id;
  }

  current_children_id[0] = 0; // get to know the starting id of each kind of children
  for (int i = 1; i < (1<<P4EST_DIM); i++){
    current_children_id[i] = children_num_pts[i-1] + current_children_id[i-1];
  }


  for (int i = start_id; i <= end_id; i++){
    cid = children_id[i-start_id];
    cur_bulk_child_id = current_children_id[cid];
    temp_pts[cur_bulk_child_id * P4EST_DIM] = local_pts[i*P4EST_DIM];
    temp_pts[cur_bulk_child_id * P4EST_DIM +1] = local_pts[i*P4EST_DIM+1];
    current_children_id[cid] += 1;
  }

  for (int i = 0; i < local_num_pt; i++){
    local_pts[(i+start_id)*P4EST_DIM] = temp_pts[i*P4EST_DIM];
    local_pts[(i+start_id)*P4EST_DIM+1] = temp_pts[i*P4EST_DIM+1];
  }
  
  free(children_id);
  free(temp_pts);
}



static void replace_quads (p4est_t * p4est, p4est_topidx_t which_tree,
                     int num_outgoing,
                     p4est_quadrant_t * outgoing[],
                     int num_incoming, p4est_quadrant_t * incoming[])
{
  int children_num_pts[1<<P4EST_DIM]={0};
  bd_t *parent_data, *child_data;
  int i;
  if (num_outgoing <= 1){
    parent_data = (bd_t*) outgoing[0]->p.user_data;
    quadrant_scan_and_order(p4est, outgoing[0], children_num_pts);

    for (i= 1; i < P4EST_CHILDREN; i++){
      children_num_pts[i] += children_num_pts[i-1]; // make children_num_pts a cumsum.
    }
    i = 0;
    child_data = (bd_t *)incoming[i]->p.user_data;
    child_data->u[0] = parent_data->u[0];
    for (i = 1; i < P4EST_CHILDREN; i++){
      child_data = (bd_t *)incoming[i]->p.user_data;
      child_data->u[0] = children_num_pts[i-1] + parent_data->u[0];
    }
    i = P4EST_CHILDREN - 1;
    child_data = (bd_t *)incoming[i]->p.user_data;
    child_data->u[1] = parent_data->u[1];
    for (i = 0; i < P4EST_CHILDREN - 1; i++){
      child_data = (bd_t *)incoming[i]->p.user_data;
      child_data->u[1] = children_num_pts[i] + parent_data->u[0] - 1;
    }
  }
}

// see function p4est_volume_iterate_simple in p4est_iterate.c
int get_num_each_proc(p4est_t * p4est)
{
  p4est_topidx_t first_local_tree = p4est->first_local_tree;
  sc_array_t *trees = p4est->trees, *quadrants;
  if (first_local_tree == -1){
    //printf("get_num_each_proc_first_local_tree is -1\n");
    return 0;
  }
  p4est_tree_t *tree= p4est_tree_array_index_diy(trees, first_local_tree);
  p4est_quadrant_t *q;
  bd_t *data;
  size_t n_quads;
  int bds[2];
  quadrants = &(tree->quadrants);
  n_quads = quadrants->elem_count;
  q = p4est_quadrant_array_index_diy(quadrants, 0);
  data = (bd_t *)q->p.user_data;
  bds[0] = data->u[0];
  q = p4est_quadrant_array_index_diy(quadrants, n_quads-1);
  data = (bd_t *)q->p.user_data;
  bds[1] = data->u[1];
  return (bds[1] - bds[0] + 1);
} 

int nothingtosend(int *old, int id){
  if (id == 0){
    if (old[id] == 0){
      return 1;
    }
  }
  else{
    if (old[id] == old[id-1]){
      return 1;
    }
  }
  return 0;
}

void find_send_indices(int *old, int *new, int id, int total_len, int *counts_send, int *displs_send)
{
  int lu, ru; // bound of indices
  int lv, rv; // bound of values
  if (nothingtosend(old, id)){
    return;
  }

  if (id == 0){
    lv = 0;
  }
  else {
    lv = old[id-1];
  }
  for (int i = 0; i < total_len; i++){
    if (new[i] > lv){
      lu = i;
      break;
    }
  }

  rv = old[id];
  for (int i = lu; i < total_len; i++){
    if (new[i] >= rv){
      ru = i;
      break;
    }
  }
  if (lu == ru){
    counts_send[lu] = (old[id] - lv) * P4EST_DIM;
    displs_send[lu] = 0;
  }
  else{
    counts_send[lu] = (new[lu] - lv) * P4EST_DIM;
    displs_send[lu] = 0;
    for (int i = lu + 1; i < ru; i++){
      counts_send[i] = (new[i] - new[i-1]) * P4EST_DIM;
      displs_send[i] = displs_send[i-1] + counts_send[i-1];
    }
    counts_send[ru] = (rv - new[ru-1]) * P4EST_DIM;
    displs_send[ru] =  displs_send[ru-1] + counts_send[ru-1];
  }

}

// get to know the mpi communication
void communications_among_procs(p4est_t * p4est)
{
  int loc_num;
  int my_size = p4est->mpisize;
  int *new_cumsum_num_pt = (int *) malloc(my_size * sizeof(int));
  loc_num = get_num_each_proc(p4est);
  double *newpts = (double *) malloc((loc_num * P4EST_DIM) * sizeof(double));
  // this assumes that the ordering between processors is z-order. It is related to the behavior of p4est_partition.
  sc_MPI_Allgather(&loc_num, 1, sc_MPI_INT, new_cumsum_num_pt, 1, sc_MPI_INT, p4est->mpicomm);
  for (int i = 1; i < my_size; i++){
    new_cumsum_num_pt[i] += new_cumsum_num_pt[i-1];
  }
  int my_rank = p4est->mpirank;
  int counts_send[my_size];
  int counts_recv[my_size];
  int displs_send[my_size];
  int displs_recv[my_size];
  memset(counts_send, 0, my_size * sizeof(int));
  memset(counts_recv, 0, my_size * sizeof(int));
  memset(displs_send, 0, my_size * sizeof(int));
  memset(displs_recv, 0, my_size * sizeof(int));
  local_info_t *linfo = (local_info_t*) p4est->user_pointer;
  int *old_cumsum_num_pt = linfo->cumsum_num_pt;
  find_send_indices(old_cumsum_num_pt, new_cumsum_num_pt, my_rank, my_size, counts_send, displs_send);
  // actually the following is recv
  find_send_indices(new_cumsum_num_pt, old_cumsum_num_pt, my_rank, my_size, counts_recv, displs_recv);
  MPI_Alltoallv(linfo->pts, counts_send, displs_send, MPI_DOUBLE, newpts, counts_recv, displs_recv, MPI_DOUBLE, p4est->mpicomm);
  free(linfo->pts);
  linfo->pts = newpts;
  free(linfo->cumsum_num_pt); // don't use array then free // ** notice memory leak management***
  linfo->cumsum_num_pt = new_cumsum_num_pt;
}

int find_max_level(p4est_t * p4est)
{
  p4est_topidx_t t;
  p4est_topidx_t first_local_tree = p4est->first_local_tree;
  // if (first_local_tree <0 ){
  //   printf("find_max_level_first_local_tree is %d\n", first_local_tree);
  // }
  p4est_topidx_t last_local_tree = p4est->last_local_tree;
  // if (last_local_tree <0){
  //   printf("find_max_level_last_local_tree is %d\n", last_local_tree);
  // }
  sc_array_t *trees = p4est->trees;
  p4est_tree_t *tree;
  size_t si, n_quads;
  sc_array_t *quadrants;
  p4est_quadrant_t *q;
  int local_max_level=0;
  int global_max_level=0;

  // if no quadrant, first_local_tree == -1 and last_local_tree == -2. so no loop.
  
  for (t = first_local_tree; t <= last_local_tree; t++){
    tree= p4est_tree_array_index_diy(trees, t);
    quadrants = &(tree->quadrants);
    n_quads = quadrants->elem_count;
    for (si= 0; si < n_quads; si++){
      q = p4est_quadrant_array_index_diy(quadrants, si);
      if (q->level > local_max_level){
        local_max_level = q->level;
      }
    }
  }
  sc_MPI_Allreduce(&local_max_level, &global_max_level, 1, sc_MPI_INT, sc_MPI_MAX, p4est->mpicomm);
  return global_max_level;
}

void init_initial_condition (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
  int my_rank = p4est->mpirank;
  local_info_t *linfo = (local_info_t*) p4est->user_pointer;
  bd_t *bd = (bd_t *)q->p.user_data;
  if (my_rank == 0){
    bd->u[0] = 0;
  }
  else{
    bd->u[0] = linfo->cumsum_num_pt[my_rank-1]; // not + 1 !
  }
  bd->u[1] = linfo->cumsum_num_pt[my_rank] - 1; // previously, not -1.
}

double max_dev_from_equal(int *num, int len)
{ // Note that num is a cumsum.
  double average = 0.0;
  double cur;
  average = (num[len-1]*1.0)/len;
  double dev = fabs(num[0]*1.0 - average); // need math.h
  for (int i = 1; i < len; i++){
    cur = fabs((num[i] - num[i-1])*1.0-average);
    if (cur > dev){
      dev = cur;
    }
  }
  return dev;
}

void printpoints(p4est_t * p4est, double * global_range)
{
  local_info_t *linfo = (local_info_t*) p4est->user_pointer;
  double *local_pts = linfo->pts;
  int my_rank = p4est->mpirank;
  int *local_cumsum = linfo->cumsum_num_pt;
  int num;
  double * output_pts = (double *) malloc(sizeof(double) * P4EST_DIM);
  printf("printpoints: rank = %d, worldsize = %d", my_rank, p4est->mpisize);
  for(int i = 0; i < p4est->mpisize; i++){
    printf("%d th cumsum = %d ", i, local_cumsum[i]);
  }
  if (my_rank==0){
    num = local_cumsum[my_rank];
  }
  else{
    num = local_cumsum[my_rank] - local_cumsum[my_rank-1];
  }
  for (int i = 0; i < num; i++){
    transform_loc_fromgridtodata_double(output_pts, local_pts + (i*P4EST_DIM), global_range);
    printf("%lf, %lf;", output_pts[0], output_pts[1]);
  }
  printf("\n");
}

void printfirstlocaltree(p4est_t * p4est)
{
  printf("rank = %d, firstlocaltree = %d\n", p4est->mpirank, p4est->first_local_tree);
}



/** The main function of the step1 example program.
 *
 * It creates a connectivity and forest, refines it, and writes a VTK file.
 */
int main (int argc, char **argv)
{
  int                 mpiret;
  int                 recursive, partforcoarsen, allowed_level;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  int world_size, world_rank;
  int min_level;
  local_info_t initial_info;
  local_info_t *linfo;

  int total_num_pt = atoi(argv[2]);
  int local_num_pt;
  double dev;
  double tol;
  double margin = 0.1;

  double *range = (double*) malloc(4*sizeof(double));
  double *global_range = (double*) malloc(4*sizeof(double));
  p4est_quadrant_t *quadrant;
  p4est_qcoord_t local_xy[2];
  double *local_pts;
  double *points;
  int *locstored_cumsum_num_pt;
  FILE *fp;
  char outputfile[30] = "output";
  char globaloutput[30]="global";
  int max_level;
  double output_pts[2];

  /* Initialize MPI; see sc_mpi.h.
   * If configure --enable-mpi is given these are true MPI calls.
   * Else these are dummy functions that simulate a single-processor run. */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  

  /* These functions are optional.  If called they store the MPI rank as a
   * static variable so subsequent global p4est log messages are only issued
   * from processor zero.  Here we turn off most of the logging; see sc.h. */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_PRODUCTION);
  P4EST_GLOBAL_PRODUCTIONF
    ("This is the p4est %dD %s partition\n",
     P4EST_DIM, P4EST_STRING);

  /* Create a forest that consists of just one quadtree/octree.
   * This file is compiled for both 2D and 3D: the macro P4_TO_P8 can be
   * checked to execute dimension-dependent code. */
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_unitsquare ();
#else
  conn = p8est_connectivity_new_unitcube ();
#endif

  // *** need allocate and initialize initial_info from file reading ***
  sc_MPI_Comm_size(mpicomm, &world_size);
  sc_MPI_Comm_rank(mpicomm, &world_rank);
  char mpirankstring[30]={'\0'};
  sprintf(mpirankstring,"%d", world_rank);
  sprintf(mpirankstring+strlen(mpirankstring),"%s", argv[1]);
  strcat(globaloutput, argv[1]);
  min_level = log2_diy(world_size)/P4EST_DIM; // make sure the mpi_size is set to be good
  tol = 2.0;//SC_MAX(2.0, ((total_num_pt * 1.0) / world_size) * 0.1); // in sc.h
  locstored_cumsum_num_pt = (int *) malloc(sizeof(int) * world_size);

  /* Create a forest that is not refined; it consists of the root octant. */
  // ***read p4est_new_ext to make sure we are doing the correct thing here.***
  // might need to first initialize some null to first know the x, y values of the grids, and then read the corresponding
  // data, and then initialize correspondingly (both quadrant and p4est).
  p4est = p4est_new(mpicomm, conn, sizeof(bd_t), NULL, NULL);
  for (int i = 0; i < min_level; i++){
    p4est_refine(p4est, 0, initial_refine_func, NULL);
  }
  p4est_partition(p4est, 0, NULL);
  recursive = 0;
  partforcoarsen = 0; 
  allowed_level = P4EST_QMAXLEVEL;

  // *** write a function to make sure every proc has one and only one quadrant
  // might be interesting to also print the locs, to see the relation with mpi_rank

  // read data
  
  get_local_range(total_num_pt, world_rank, world_size, range, argv[1]);
  sc_MPI_Allreduce(range, global_range, 1, sc_MPI_DOUBLE, sc_MPI_MIN, mpicomm);
  sc_MPI_Allreduce(range+1, global_range+1, 1, sc_MPI_DOUBLE, sc_MPI_MAX, mpicomm);
  sc_MPI_Allreduce(range+2, global_range+2, 1, sc_MPI_DOUBLE, sc_MPI_MIN, mpicomm);
  sc_MPI_Allreduce(range+3, global_range+3, 1, sc_MPI_DOUBLE, sc_MPI_MAX, mpicomm);
  for (int i = 0; i < P4EST_DIM; i++){
    global_range[2*i] = global_range[2*i] - margin;
    global_range[2*i+1] = global_range[2*i+1] + margin;
  }
  //printf("rank = %d, global_range= %lf, %lf, %lf, %lf\n", p4est->mpirank, global_range[0], global_range[1], global_range[2], global_range[3]);


  quadrant = get_first_quadrant(p4est);
  local_xy[0] = quadrant->x;
  local_xy[1] = quadrant->y;
  //printf("rank = %d, local_xy= %d, %d\n", p4est->mpirank, local_xy[0], local_xy[1]);
  local_num_pt = count_num(local_xy, P4EST_QUADRANT_LEN(quadrant->level), global_range, argv[1]);
  //printf("local_num_pt= %d\n", local_num_pt);
  local_pts = (double *)malloc((local_num_pt*P4EST_DIM)*sizeof(double));
  get_pts(local_xy, P4EST_QUADRANT_LEN(quadrant->level), global_range, argv[1], local_pts, local_num_pt);
  //printf("rank = %d, local_pts= %lf, %lf\n", p4est->mpirank, local_pts[0], local_pts[1]);
  initial_info.pts = local_pts;
  sc_MPI_Allgather(&local_num_pt, 1, sc_MPI_INT, locstored_cumsum_num_pt, 1, sc_MPI_INT, mpicomm);
  for(int i = 1; i < world_size; i++){
    locstored_cumsum_num_pt[i] += locstored_cumsum_num_pt[i-1];
  }
  initial_info.cumsum_num_pt = locstored_cumsum_num_pt;
  p4est->user_pointer = (void*) &initial_info;
  //printf("rank = %d, local_num_pt= %d, locstored_cumsum_num_pt= %d\n", p4est->mpirank, local_num_pt, locstored_cumsum_num_pt[p4est->mpirank]);

  // personal thinking of mempool: helps to manage memory. For what we do ourselves, remember to free the memory.
  // initialize the quadrant. see p4est_algorithms.c p4est_quadrant_init_data

  init_initial_condition(p4est, p4est->first_local_tree, quadrant);

  
  
  double tt = sc_MPI_Wtime();  
  // *** need to set the stopping criterion *** //
  // ** e.g. consider although still needs to be refined, it is balanced already
  linfo = (local_info_t *) p4est->user_pointer;
  dev = max_dev_from_equal(linfo->cumsum_num_pt, world_size);
  // if (p4est->mpirank ==0){
  //   printf("initial dev = %lf\n", dev);
  //   printpoints(p4est, global_range);
  // }
  max_level = find_max_level(p4est);
  //int flag = 0;
  int flag = 0;
  if (p4est->mpirank == 0){
    printf("tol = %lf\n initial dev = %lf\n", tol, dev);
    // for (int i =0; i<p4est->mpisize; i++){
    //   printf("locstored_cumsum_num_pt[%d] = %d\n", i, locstored_cumsum_num_pt[i]);
    // }
  }
  while((dev > tol) && (max_level < P4EST_QMAXLEVEL)){ // consider max_refine reached
    //printf("refine again\n");
    // if (p4est->mpirank==3){
    //   printf("I am 3 here!\n");
    // }
    //printpoints(p4est, global_range);
    //printfirstlocaltree(p4est);
    p4est_refine_ext (p4est, recursive, allowed_level, refine_func, NULL, replace_quads);
    //printf("after refine");
    //printfirstlocaltree(p4est);
    p4est_partition(p4est, partforcoarsen, weight_func);
    //printf("after partition");
    //printfirstlocaltree(p4est);
    communications_among_procs(p4est);
    max_level = find_max_level(p4est);
    dev = max_dev_from_equal(linfo->cumsum_num_pt, world_size);
    if (p4est->mpirank == 0){
      printf("flag = %d, dev = %lf\n", flag, dev);
    }
    flag = flag + 1;
    // if (flag > 5){
    //   break;
    // }
  }
  tt = sc_MPI_Wtime() - tt;
  if (p4est->mpirank == 0){
    printf("total number of refines: %d\n", flag);
    printf("total time: %lf s\n", tt);
  }

  linfo = (local_info_t *) p4est->user_pointer;

  // write the output
  // need stdio.h, string.h
  // when storing the points, we scale them; now we need to transform them back.
  fp = fopen(strcat(outputfile, mpirankstring), "w");
  if (fp == NULL){
    printf("could not open file");
  }
  local_num_pt = linfo->cumsum_num_pt[world_rank];
  if (world_rank > 0){
    local_num_pt -= linfo->cumsum_num_pt[world_rank-1];
  }
  points = linfo -> pts;
  for (int i = 0; i < local_num_pt ; i++){
    transform_loc_fromgridtodata_double(output_pts, points + (i*P4EST_DIM), global_range);
    //printf("rank %d , %lf, %lf\n", world_rank, output_pts[0], output_pts[1]);
    fprintf(fp, "%lf\n%lf\n", output_pts[0], output_pts[1]);
  }
  fclose(fp);
  if (world_rank == 0){
    fp = fopen(globaloutput, "w");
    if (fp == NULL){
      printf("could not open file");
    }
    for (int i = 0; i < world_size; i++){
      //printf("%d", linfo->cumsum_num_pt[i]);
      fprintf(fp, "%d ", linfo->cumsum_num_pt[i]);
    }
    fclose(fp);
  }
  // *** need to free the memory that allocated by ourselves and wasn't freed yet. ***
  /* Destroy the p4est and the connectivity structure. */
  free(range);
  free(global_range);
  free(linfo->cumsum_num_pt);
  free(linfo->pts);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  /* Verify that allocations internal to p4est and sc do not leak memory.
   * This should be called if sc_init () has been called earlier. */
  sc_finalize ();

  /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
