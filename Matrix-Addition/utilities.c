/*******************************************************************************
  Title          : utilities.c
  Author         : Stewart Weiss
  Created on     : February 10, 2014
  Description    : Various functions used for MPI and non-MPI programs
  Purpose        : 
  Build with     : mpicc -c utilities.c
 
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "utilities.h"

#define PROMPT_MSG         1
#define RESPONSE_MSG       2


/******************************************************************************/
/** get_type() returns a constant representing MPI_Datatype argument
 *  @param    MPI_Datatype  t  
 *  @return   integer constant
 */
int get_type (MPI_Datatype t) 
{
    if ( ( t == MPI_BYTE) ||
         ( t == MPI_CHAR) ||
         ( t == MPI_SIGNED_CHAR) ||
         ( t == MPI_UNSIGNED_CHAR) )
        return CHAR;

    if ( t == MPI_DOUBLE )
        return DOUBLE;
    if (t == MPI_FLOAT) 
        return FLOAT;

    if ( ( t == MPI_INT) ||
         ( t == MPI_UNSIGNED) )
        return INT;

    if ( ( t == MPI_LONG) ||
         ( t == MPI_UNSIGNED_LONG) )
        return LONG;

    if ( ( t == MPI_LONG_LONG) ||
         ( t == MPI_LONG_LONG_INT) ||
         ( t == MPI_UNSIGNED_LONG_LONG) )
        return LLONG;
    
   return -1;
}



/******************************************************************************/
/** get_size() returns the number of bytes in MPI_Datatype argument
 *  @param    MPI_Datatype  t  
 *  @return   number of bytes in t
 */
int get_size (MPI_Datatype t) 
{
    if ( ( t == MPI_BYTE) ||
         ( t == MPI_CHAR) ||
         ( t == MPI_SIGNED_CHAR) ||
         ( t == MPI_UNSIGNED_CHAR) )
        return sizeof(char);

    if ( t == MPI_SHORT )
        return sizeof(short int);

    if ( t == MPI_DOUBLE )
        return sizeof(double);
    if (t == MPI_FLOAT) 
        return sizeof(float);

    if ( ( t == MPI_INT) ||
         ( t == MPI_UNSIGNED) )
        return sizeof(int);

    if ( ( t == MPI_LONG) ||
         ( t == MPI_UNSIGNED_LONG) )
        return sizeof(long int);

    if ( ( t == MPI_LONG_LONG) ||
         ( t == MPI_LONG_LONG_INT) ||
         ( t == MPI_UNSIGNED_LONG_LONG) )
        return sizeof(long long int);
    
   return -1;
}

/******************************************************************************/
/** terminate() prints an error message and terminates calling process 
 *  @param   int id        [IN] rank of calling process
 *  @param   char* message [IN] error message to print
 *  @post    A message is printed on standard output only if id == 0
             but the caller is always terminated with MPI_Finalize.
 */
void terminate (
   int   id,            /* IN - Process rank */
   char *error_message) /* IN - Message to print */
{
   if ( 0 == id ) {
      printf ("%s", error_message);
      fflush (stdout);
   }
   MPI_Finalize();
   exit (0);
}

/******************************************************************************/

/* owner(r,p,n) is the rank of the process that owns element r */
int owner( int row, int num_procs, int total_elements )
{
    return ( num_procs * (row+1) -1 ) / total_elements;
}

/******************************************************************************/

int size_of_block( int id, int ntotal_elements, int p )
{
    return ( ( ( id + 1) * ntotal_elements ) / p ) - 
           ( ( id *      ntotal_elements ) / p );
}

/******************************************************************************/
void alloc_matrix( 
        int     nrows,          /* number of rows in matrix                   */
        int     ncols,          /* number of columns in matrix                */
        size_t  element_size,   /* number of bytes per matrix element         */
        void  **matrix_storage, /* address of linear storage array for matrix */
        void ***matrix,         /* address of start of matrix                 */
        int    *errvalue)       /* return code for error, if any              */
{
    int   i;
    void *ptr_to_row_in_storage; /* pointer to a place in linear storage array
                                    where a row begins                        */
    void **matrix_row_start;     /* address of a 2D matrix row start pointer
                                    e.g., address of (*matrix)[row]           */
    size_t total_bytes;          /* amount of memory to allocate              */

    //printf("alloc_matrix called with r=%d,c=%d,e=%d\n",nrows, ncols, element_size);

    total_bytes = nrows * ncols * element_size;

    /* Step 1: Allocate an array of nrows * ncols * element_size bytes  */  
    *matrix_storage = malloc(total_bytes);
    if ( NULL == *matrix_storage ) {
        /* malloc failed, so set error code and quit */
        *errvalue = MALLOC_ERROR;
        return;
    }

    memset(*matrix_storage, 0, total_bytes );

    /* Step 2: To create the 2D matrix, first allocate an array of nrows 
       void* pointers */   
    *matrix = malloc (nrows * sizeof(void*));
    if ( NULL == *matrix ) {
        /* malloc failed, so set error code and quit */
        *errvalue = MALLOC_ERROR;
        return;
    }


    /* Step 3: (The hard part) We need to put the addresses into the
       pointers of the 2D matrix that correspond to the starts of rows
       in the linear storage array. The offset of each row in linear storage 
       is a multiple of (ncols * element_size) bytes.  So we initialize
       ptr_to_row_in_storage to the start of the linear storage array and
       add (ncols * element_size) for each new row start.
       The pointers in the array of pointers to rows are of type void* 
       so an increment operation on one of them advances it to the next pointer.
       Therefore, we can initialize matrix_row_start to the start of the 
       array of pointers, and auto-increment it to advance it.
    */

    /* Get address of start of array of pointers to linear storage, 
       which is the address of first pointer, (*matrix)[0]   */
    matrix_row_start = (void*) &(*matrix[0]);

    /* Get address of start of linear storage array */
    ptr_to_row_in_storage = (void*) *matrix_storage;

    /* For each matrix pointer, *matrix[i], i = 0... nrows-1, 
       set it to the start of the ith row in linear storage */
    for ( i = 0; i < nrows; i++ ) {
        /* matrix_row_start is the address of (*matrix)[i] and
           ptr_to_row_in_storage is the address of the start of the 
           ith row in linear storage.
           Therefore, the following assignment changes the contents of 
           (*matrix)[i]  to store the start of the ith row in linear storage 
        */
        *matrix_row_start = (void*) ptr_to_row_in_storage;

        /* advance both pointers */
        matrix_row_start++;     /* next pointer in 2d array */
        ptr_to_row_in_storage +=  ncols * element_size; /* next row */
    }
    *errvalue = SUCCESS;
        

}

/******************************************************************************/
void print_matrix (
        void  **matrix,         /* matrix to be printed        */
        int     nrows,          /* number of rows in matrix    */
        int     ncols,          /* number of columns in matrix */
        MPI_Datatype dtype,     /* MPI type                    */
        FILE  *stream)          /* stream on which to print    */
{
    int i, j;
    int etype = get_type(dtype);
   
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            switch (etype) {
            case DOUBLE:
                fprintf (stream, "%6.3f ", ((double **)matrix)[i][j]); break;
            case FLOAT:
                fprintf (stream,"%6.3f ", ((float **)matrix)[i][j]);   break;
            case INT:
                fprintf (stream,"%6d ", ((int **)matrix)[i][j]);       break;
            case CHAR:
                fprintf (stream,"%6c ", ((char **)matrix)[i][j]);      break;
            case LONG:
                fprintf (stream,"%6ld ", ((long int **)matrix)[i][j]); break;
            case LLONG:
                fprintf (stream,"%6lld ", ((long long int **)matrix)[i][j]);
            }
        }
        fprintf(stream,"\n");
   }
}

/******************************************************************************/
void print_vector (
        void   *vector,         /* vector to be printed        */
        int     n,              /* number of elements in vector*/
        MPI_Datatype dtype,     /* MPI type                    */
        FILE  *stream)          /* stream on which to print    */
{
    int i;
    int etype = get_type(dtype);
   
    for (i = 0; i < n; i++) {
            switch (etype) {
            case DOUBLE:
                fprintf (stream, "%6.4f ", ((double *)vector)[i]); break;
            case FLOAT:
                fprintf (stream,"%6.4f ", ((float *)vector)[i]);   break;
            case INT:
                fprintf (stream,"%6d ", ((int *)vector)[i]);       break;
            case CHAR:
                fprintf (stream,"%6c ", ((char *)vector)[i]);      break;
            case LONG:
                fprintf (stream,"%6ld ", ((long int *)vector)[i]); break;
            case LLONG:
                fprintf (stream,"%6lld ", ((long long int *)vector)[i]);
            }
    }
}

void print_full_vector (
    void        *v,      /* IN - Address of vector */
    int          n,      /* IN - Elements in vector */
    MPI_Datatype dtype,  /* IN - Vector element type */
    MPI_Comm     comm)   /* IN - Communicator */
{
    int id;              /* Process rank */

    MPI_Comm_rank (comm, &id);
   
    if (0 == id) {
        print_vector (v, n, dtype, stdout);
        printf("\n");
    }
}

void print_block_vector (
    void        *v,       /* IN - Address of vector */
    MPI_Datatype dtype,   /* IN - Vector element type */
    int          n,       /* IN - Elements in vector */
    MPI_Comm     comm)    /* IN - Communicator */
{
    int        datum_size; /* Bytes per vector element */
    int        i;
    int        prompt;     /* Dummy variable */
    MPI_Status status;     /* Result of receive */
    void       *tmp;       /* Other process's subvector */
    int        id;         /* Process rank */
    int        p;          /* Number of processes */

    MPI_Comm_size (comm, &p);
    MPI_Comm_rank (comm, &id);
    datum_size = get_size (dtype);

    if ( 0 == id ) {
        print_vector (v,  size_of_block(id,n,p), dtype, stdout);
        if (p > 1) {
            tmp = malloc (size_of_block(p-1,n,p)*datum_size);
            for (i = 1; i < p; i++) {
                MPI_Send (&prompt, 1, MPI_INT, i, PROMPT_MSG, comm);
                MPI_Recv (tmp, size_of_block(i,n,p), dtype, i,
                          RESPONSE_MSG, comm, &status);
                print_vector (tmp,  size_of_block(i,n,p), dtype, stdout);
            }
            free (tmp);
        }
        printf ("\n\n");
    } 
    else {
        MPI_Recv (&prompt, 1, MPI_INT, 0, PROMPT_MSG, comm, &status);
        MPI_Send (v, size_of_block(id,n,p), dtype, 0,
                  RESPONSE_MSG, comm);
    }
}


/******************************************************************************/
/** init_communication_arrays() initialize arrays to pass to MPI gather/scatterv
 *  @param  int p            [IN]   Number of processes 
 *  @param  int n            [IN]   Total number of elements
 *  @param  int *count       [OUT]  Array of counts
 *  @return int *offset      [OUT]  Array of displacements 
 */
void init_communication_arrays (
    int p,          /* IN - Number of processes */
    int n,          /* IN - Total number of elements */
    int *count,     /* OUT - Array of counts */
    int *offset)    /* OUT - Array of displacements */
{
    int i;

    count[0]  = size_of_block(0,n,p);
    offset[0] = 0;
    for (i = 1; i < p; i++) {
        offset[i] = offset[i-1] + count[i-1];
        count[i]  = size_of_block(i,n,p);
    }
}

/******************************************************************************/
/** replicate_block_vector() copies a distributed vector into every process
 *  @param  void        *invec   [IN]   Block-distributed vector
 *  @param  int          n       [IN]   Total number of elements in vector
 *  @param  MPI_Datatype dtype   [IN]   MPI element type
 *  @param  void        *outvec  [OUT]  Replicated vector
 *  @param  MPI_Comm     comm    [IN]  Communicator
 */
void replicate_block_vector (
        void        *invec,  
        int          n,      
        MPI_Datatype dtype,  
        void        *outvec,   
        MPI_Comm     comm   
        )
{
    int *recv_count;  /* Elements contributed by each process */
    int *recv_offset; /* Displacement in concatenated array */
    int id;           /* Process id */
    int p;            /* Processes in communicator */

    MPI_Comm_size (comm, &p);
    MPI_Comm_rank (comm, &id);

    /* Allocate count and offset arrays and bail out if either fails. */
    recv_count  = malloc ( p * sizeof(int));
    recv_offset = malloc ( p * sizeof(int));
    if ( NULL == recv_offset || NULL == recv_count ) {
        printf ("malloc failed for process %d\n", id);
        MPI_Abort (MPI_COMM_WORLD, MALLOC_ERROR);
    }

    /* Fill the count and offset arrays to pass to MPI_Allgatherv
       so that the blocks are concatenated by process rank in the output
       vector */
    init_communication_arrays (p, n, recv_count, recv_offset);

    /* Use MPI_Allgatherv to copy the distributed blocks from invec in
       each process into a replicated outvec in each process. */
    MPI_Allgatherv (invec, recv_count[id], dtype, outvec, recv_count,
                    recv_offset, dtype, comm);

    /* Release the storage for the count and offset arrays. */
    free (recv_count);
    free (recv_offset);
}


/*******************************************************************************
                           Random Number Routines
*******************************************************************************/


/******************************************************************************/
/** init_random()  initializes the state for the C random() function
 *  @param  int    state_size  [IN]  Size of state array for random to use
 *  @return char*  a pointer to the state array allocated for random()
 *  @post          After this call, an array of size state_size*sizeof(char) has
 *                 been allocated and initialized by C initstate(). It must be
 *                 freed by calling free()
 */
char*  init_random( int state_size )
{
    char * state;
    state  = (char*) malloc ( state_size * sizeof(char));
    if ( NULL != state )
        initstate(time(NULL), state, state_size);
    return state;
}

void   finalize ( char* state )
{
    free (state );
}

/** uniform_random()  returns a uniformly distributed random number in [0,1]
 *  @return double  a pointer to the state array allocated for random()
 *  @pre           Either init_random() should have been called or srandom()
 */
double uniform_random()
{
    return (double) (random()) / RAND_MAX; 
}




