
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "utilities.h" 


#define PROMPT_MSG         1
#define RESPONSE_MSG       2



void read_and_distribute_matrix_byrows (
          char        *filename,       /* [IN]  name of file to read          */
          void      ***matrix,         /* [OUT] matrix to fill with data      */
          void       **matrix_storage, /* [OUT] linear storage for the matrix */
          MPI_Datatype dtype,          /* [IN]  matrix element type           */
          int         *nrows,          /* [OUT] number of rows in matrix      */
          int         *ncols,          /* [OUT] number of columns in matrix   */
          int         *errval,         /* [OUT] success/error code on return  */
          MPI_Comm     comm)           /* [IN] communicator handle            */
{

    int    id;              /* process rank process */
    int    p;               /* number of processes  in communicator group */
    size_t element_size;    /* number of bytes in matrix element type */
    int    mpi_initialized; /* flag to check if MPI_Init was called already */
    FILE   *file;           /* input file stream pointer */
    int    nlocal_rows;     /* number of rows calling process "owns" */
    MPI_Status   status;    /* result of MPI_Recv */
    const int    MSG_TAG=1;
    /* Make sure we are being called by a program that init-ed MPI */
    MPI_Initialized(&mpi_initialized);
    if ( !mpi_initialized ) {
       *errval = -1;
       return;
    }
  
    /* Get process rank and the number of processes in group */
    MPI_Comm_size (comm, &p);
    MPI_Comm_rank (comm, &id);
   
    /* Get the number of bytes in a matrix element */
    element_size = get_size (dtype);
    if ( element_size <= 0 ) {
       *errval = -1;
       return;
    }
       
    if ( p-1 == id ) { 
        /* Process p-1 opens the binary file containing the matrix and 
           reads the first two numbers, which are the number of rows and 
           columns respectively. */
        file = fopen (filename, "rb");
        if ( NULL == file ) {
            *nrows = 0;
            *ncols = 0;
        }
        else {
            fread (nrows, sizeof(int), 1, file);
            fread (ncols, sizeof(int), 1, file);
        }      
    }

    /* Process p-1  broadcasts the numbers of rows to all other processes. */
    MPI_Bcast (nrows, 1, MPI_INT, p-1, comm);

    if ( 0 == *nrows  ) {
       *errval = -1;
       return;
    }

    /* Process p-1  broadcasts the numbers of columns to all other processes. */
    MPI_Bcast (ncols, 1, MPI_INT, p-1, comm);

    /* Each process sets local_rows = the number of rows the process owns.
       The number of rows depends on id, *nrows, and p. It is the difference
       between the high address and the low address. */
    nlocal_rows = size_of_block( id, *nrows, p );

    /* Each process creates its linear storage and 2D matrix for accessing 
       the elements of its assigned rows. */
    alloc_matrix( nlocal_rows, *ncols, element_size,
                  matrix_storage, 
                  matrix,         
                  errval);
    if ( SUCCESS != *errval ) {
         MPI_Abort (comm, *errval);
    }

    if ( p-1 == id  ) {
        int nrows_to_send;     /* number of rows that p-1 sends to a process */
        int num_elements;      /* total number of matrix elements to send */
        size_t nelements_read; /* result of read operation */
        int i;                 /* loop index */

        /* For each process i, for i = 0 up to p-2,
           Process p-1 reads a consecutive chunk of bytes from the file 
           that contains the rows for process i and then sends that chunk 
           of bytes to process i. */
         for ( i = 0; i < p-1; i++) {
             nrows_to_send  = size_of_block( i, *nrows, p );
             num_elements   = nrows_to_send * (*ncols);
             nelements_read = fread (*matrix_storage, element_size,
                                         num_elements, file );

             /* Check that the number of items read matches the number 
                requested. If not, abort. */
             if ( nelements_read != num_elements )
                 MPI_Abort(comm, FILE_READ_ERROR);
             MPI_Send (*matrix_storage, num_elements, dtype,
                        i, MSG_TAG, comm);
        }
        /* Process p-1 reads the remainder of the file into its own 
           linear storage. */
        nelements_read = fread (*matrix_storage, element_size, 
                                    nlocal_rows * (*ncols), file);
        /* Check that the number of items read matches the number requested */
        if ( nelements_read != nlocal_rows * (*ncols) )
            MPI_Abort(comm, FILE_READ_ERROR);
        fclose (file);
    } 
    else /* what all other processes do */
         /* store the data sent by process p-1 into the linear storage array 
            for *matrix, which is *matrix_storage, not **matrix_storage!  
            The number of values expected is the number of rows for this 
            process times matrix width. 
        */
        MPI_Recv (*matrix_storage, nlocal_rows * (*ncols), 
                  dtype, p-1, MSG_TAG, comm, &status);
    errval = 0;
}

void read_and_distribute_square_matrix_byrows (
          char        *filename,       /* [IN]  name of file to read          */
          void      ***matrix,         /* [OUT] matrix to fill with data      */
          void       **matrix_storage, /* [OUT] linear storage for the matrix */
          MPI_Datatype dtype,          /* [IN]  matrix element type           */
          int         *nrows,          /* [OUT] size of matrix                */
          int         *errval,         /* [OUT] success/error code on return  */
          MPI_Comm     comm)           /* [IN] communicator handle            */
{

    int    id;              /* process rank process */
    int    p;               /* number of processes  in communicator group */
    size_t element_size;    /* number of bytes in matrix element type */
    int    mpi_initialized; /* flag to check if MPI_Init was called already */
    FILE   *file;           /* input file stream pointer */
    int    nlocal_rows;     /* number of rows calling process "owns" */
    MPI_Status   status;    /* result of MPI_Recv */
    const int    MSG_TAG=1;
    /* Make sure we are being called by a program that init-ed MPI */
    MPI_Initialized(&mpi_initialized);
    if ( !mpi_initialized ) {
       *errval = -1;
       return;
    }
  
    /* Get process rank and the number of processes in group */
    MPI_Comm_size (comm, &p);
    MPI_Comm_rank (comm, &id);
   
    /* Get the number of bytes in a matrix element */
    element_size = get_size (dtype);
    if ( element_size <= 0 ) {
       *errval = -1;
       return;
    }
       
    if ( p-1 == id ) { 
        /* Process p-1 opens the binary file containing the matrix and 
           reads the first two numbers, which are the number of rows and 
           columns respectively. */
        file = fopen (filename, "r");
        if ( NULL == file ) {
            *nrows = 0;
        }
        else {
            fread (nrows, sizeof(int), 1, file);
        }      
    }

    /* Process p-1  broadcasts the matrix size to all other processes. */
    MPI_Bcast (nrows, 1, MPI_INT, p-1, comm);

    if ( 0 == *nrows  ) {
       *errval = -1;
       return;
    }

    /* Each process sets local_rows = the number of rows the process owns.
       The number of rows depends on id, *nrows, and p. It is the difference
       between the high address and the low address. */
    nlocal_rows = size_of_block( id, *nrows, p );

    /* Each process creates its linear storage and 2D matrix for accessing 
       the elements of its assigned rows. */
    alloc_matrix( nlocal_rows, *nrows, element_size,
                  matrix_storage, 
                  matrix,         
                  errval);
    if ( SUCCESS != *errval ) {
         MPI_Abort (comm, *errval);
    }

    if ( p-1 == id  ) {
        int nrows_to_send;     /* number of rows that p-1 sends to a process */
        int num_elements;      /* total number of matrix elements to send */
        size_t nelements_read; /* result of read operation */
        int i;                 /* loop index */

        /* For each process i, for i = 0 up to p-2,
           Process p-1 reads a consecutive chunk of bytes from the file 
           that contains the rows for process i and then sends that chunk 
           of bytes to process i. */
         for ( i = 0; i < p-1; i++) {
             nrows_to_send  = size_of_block( i, *nrows, p );
             num_elements   = nrows_to_send * (*nrows);
             nelements_read = fread (*matrix_storage, element_size,
                                         num_elements, file );

             /* Check that the number of items read matches the number 
                requested. If not, abort. */
             if ( nelements_read != num_elements )
                 MPI_Abort(comm, FILE_READ_ERROR);
             MPI_Send (*matrix_storage, num_elements, dtype,
                        i, MSG_TAG, comm);
        }
        /* Process p-1 reads the remainder of the file into its own 
           linear storage. */
        nelements_read = fread (*matrix_storage, element_size, 
                                    nlocal_rows * (*nrows), file);
        /* Check that the number of items read matches the number requested */
        if ( nelements_read != nlocal_rows * (*nrows) )
            MPI_Abort(comm, FILE_READ_ERROR);
        fclose (file);
    } 
    else /* what all other processes do */
         /* store the data sent by process p-1 into the linear storage array 
            for *matrix, which is *matrix_storage, not **matrix_storage!  
            The number of values expected is the number of rows for this 
            process times matrix width. 
        */
        MPI_Recv (*matrix_storage, nlocal_rows * (*nrows), 
                  dtype, p-1, MSG_TAG, comm, &status);
    errval = 0;
}


void collect_and_print_matrix_byrows (
          void       **matrix,   /* [IN] matrix to print              */
          MPI_Datatype dtype,    /* [IN] matrix element type          */
          int          nrows,    /* [IN] number of rows in matrix     */
          int          ncols,    /* [IN] number of columns in matrix  */
          MPI_Comm     comm)     /* [IN] communicator handle          */
{
    int    id;              /* process rank process                         */
    int    p;               /* number of processes  in communicator group   */
    size_t element_size;    /* number of bytes in matrix element type       */
    int    nlocal_rows;     /* number of rows calling process "owns"        */
    MPI_Status   status;    /* result of MPI_Recv                           */
    void **submatrix_buffer;/* matrix to hold submatrices sent by processes */
    void  *buffer_storage;  /* linear storage for submatrix_buffer          */
    int    max_num_rows;    /* largest number of rows of any process        */
    int    prompt;          /* synchronizing variable                       */
    int    errval;          /* to hold error values                         */

    MPI_Comm_rank (comm, &id);
    MPI_Comm_size (comm, &p);
    
    nlocal_rows = size_of_block( id, nrows, p );

    if ( 0 == id ) {
        int i;

        /* Process 0 prints its rows first. */
        print_matrix (matrix, nlocal_rows, ncols, dtype, stdout);
        if (p > 1) {
            /* Get the number of bytes in a matrix element */
            element_size = get_size (dtype);
 
            /* Get the maximum number of rows used by any process,
               which is the number process p-1 uses. */
            max_num_rows = size_of_block( p-1, nrows, p ); 

            /* Allocate the 2D matrix and backing linear storage to hold
               arrays received from remaining processes */
            alloc_matrix( max_num_rows, ncols, element_size,
                  &buffer_storage, 
                  &submatrix_buffer,         
                  &errval);   
            if ( SUCCESS != errval ) {
                 MPI_Abort (comm, errval);
            }
         
            /* Request each other process to send its rows. Rather than just
               printing what it receives, which might flood the processor on
               which process 0 is running, process 0 prompts each other process
               to send its data and then waits for it. This is a form of lock
               step synchronization */ 
               
            for (i = 1; i < p; i++) {
                /* Calculate the number of elements to be received from
                   process i */
                int num_rows     = size_of_block( i, nrows, p );
                int num_elements = num_rows * ncols;

                /* Send a message to process i telling it to send data */
                MPI_Send (&prompt, 1, MPI_INT, i, PROMPT_MSG, comm);

                /* Wait for data to arrive from proess i */
                MPI_Recv (buffer_storage, num_elements, dtype,
                          i, RESPONSE_MSG, comm, &status);

                /* Print the matrix just received */
                print_matrix (submatrix_buffer,num_rows, ncols, dtype, stdout);
            }
            /* Free the allocated memory */
            free (submatrix_buffer);
            free (buffer_storage);
        }
        fprintf(stdout, "\n");
    } 
    else {
        /* Wait for prompt message from process 0 */
        MPI_Recv (&prompt, 1, MPI_INT, 0, PROMPT_MSG, comm, &status);

        /* On receiving it, send the matrix received by the call, which is
           the set of rows belonging to the process that called this 
           function. */
        MPI_Send (*matrix, nlocal_rows * ncols, dtype, 0, RESPONSE_MSG, comm);
   }
}


void read_and_distribute_2dblock_matrix (
          char        *filename,       /* [IN]  name of file to read          */
          void      ***matrix,         /* [OUT] matrix to fill with data      */
          void       **matrix_storage, /* [OUT] linear storage for the matrix */
          MPI_Datatype dtype,          /* [IN]  matrix element type           */
          int         *nrows,          /* [OUT] number of rows in matrix      */
          int         *ncols,          /* [OUT] number of columns in matrix   */
          int         *errval,         /* [OUT] sucess/error code on return   */
          MPI_Comm     cart_comm)      /* [IN]  communicator handle           */
{
    int    i,j,k;           /* various loop index variables                  */
    int    grid_id;         /* process rank in the cartesian grid            */
    int    p;               /* number of processes in the cartesian grid     */
    size_t element_size;    /* number of bytes in matrix element type        */
    int    mpi_initialized; /* flag to check if MPI_Init was called already  */
    FILE   *file;           /* input file stream pointer                     */
    int    nlocal_rows;     /* number of rows that calling process "owns"    */
    int    nlocal_cols;     /* number of columns that calling process "owns" */
    MPI_Status   status;    /* result of MPI_Recv call                       */
    int    dest_id;         /* rank of receiving process in cartesian grid   */
    int    grid_coord[2];   /* process coordinates in the grid               */
    int    grid_periodic[2];/* flags indicating if grid wraps around         */
    int    grid_size[2];    /* dimensions of grid                            */
    void*  buffer;          /* address of temp location to store rows        */
    int    block_coord[2];  /* coordinates in grid of current block          */
    void*  source_address;  /* address of block to be sent                   */
    void*  dest_address;    /* location where block is to be received        */
    
    /* Make sure we are being called by a program that init-ed MPI */
    MPI_Initialized(&mpi_initialized);
    if ( !mpi_initialized ) {
       *errval = -1;
       return;
    }
  
    /* Get process rank in grid and the number of processes in group */
    MPI_Comm_rank (cart_comm, &grid_id);  
    MPI_Comm_size (cart_comm, &p);            
   
    /* Get the number of bytes in a matrix element */
    element_size = get_size (dtype);
    if ( element_size <= 0 ) {
       *errval = -2;
       return;
    }
       
    /* Process 0 opens the file and reads the number of rows and columns */
    if ( 0 == grid_id ) { 
        /* Process 0 opens the binary file containing the matrix and 
           reads the first two numbers, which are the number of rows and 
           columns respectively. */
        file = fopen (filename, "r");
        if ( NULL == file ) {
            *nrows = 0;
            *ncols = 0;
        }
        else { /* successful open */
            fread (nrows, sizeof(int), 1, file);
            fread (ncols, sizeof(int), 1, file);
        }      
    }

    /* Process 0 broadcasts the numbers of rows to all other processes. */
    MPI_Bcast (nrows, 1, MPI_INT, 0, cart_comm);

    /* All processes check value of *nrows; if 0 it indicates failed open */
    if ( 0 == *nrows  ) {
       *errval = -3;
       return;
    }

    /* Process 0 broadcasts the numbers of columns to all other processes. 
       No need to check whether *ncols is zero. */
    MPI_Bcast (ncols, 1, MPI_INT, 0, cart_comm);

    /* All processes obtain the grid's topology so they can determine
       their block sizes. */
    MPI_Cart_get (cart_comm, 2, grid_size, grid_periodic, grid_coord);

    /* Each process sets nlocal_rows = the number of rows the process owns and
       local_cols to the number of columns it owns.
       The number of rows depends on the process's row coordinate, *nrows, 
       and the number of grid rows in total. This implements the formula
             blocks = floor((i+1)*n/p) - floor(i*n/p)
       where i is grid coordinate in given dimension, n is either total
       number of rows, or total number of columns, and p is the number of
       processes in grid in the given dimension.
    */
    nlocal_rows = size_of_block( grid_coord[0], *nrows, grid_size[0] );
    nlocal_cols = size_of_block( grid_coord[1], *ncols, grid_size[1] );

    /* Each process creates its linear storage and 2D matrix for accessing 
       the elements of its assigned rows. It needs storage for a 2D matrix
       of nlocal_rows by nlocal_cols elements. */
    alloc_matrix( nlocal_rows, nlocal_cols, element_size,
                  matrix_storage, 
                  matrix,         
                  errval);
    if ( SUCCESS != *errval ) {
         MPI_Abort (cart_comm, *errval);
    }

   /* Grid process 0 reads in the matrix one row at a time
      and distributes each row among the MPI processes. The first step is
      to allocate storage for one row of the matrix. */
    if ( 0 == grid_id ) {
        buffer = malloc (*ncols * element_size);
        if ( buffer == NULL ) { 
            MPI_Abort (cart_comm, *errval);
        }
    }

    /* This is the read and distribute loop. Process 0 will read a row
       and send it to the processes that are supposed to have it. It needs to
       break it into blocks, with successive blocks going to the processes
       in the same grid row but successive grid columns. */

    for (i = 0; i < grid_size[0]; i++) { /* for each grid row */
        /* Set block_coord[0] to the current grid row index */
        block_coord[0] = i;

        /* For every matrix row that is part of this grid row */
        for (j = 0; j < size_of_block(i, *nrows, grid_size[0] ); j++) {

            /* Process 0 reads  a row of the matrix */
            if ( 0 == grid_id ) {
                fread (buffer, element_size, *ncols, file);
            }

            /* Every process executes this loop. For each grid column within
               the current grid row ... */
            for (k = 0; k < grid_size[1]; k++) {
                block_coord[1] = k;

                /* Determine the grid id of the process in grid position
                   [i,k].  This is the destination process. Its id is returned 
                   in dest_id. */
                MPI_Cart_rank (cart_comm, block_coord, &dest_id);

                /* Process 0 needs to determine the start of the block to
                   be sent to process dest_id. This is the start in a row
                   with *ncols elements assigned to kth process out of
                   grid_size[1] many processes in that row. */
                if ( 0 == grid_id ) {
                    source_address = buffer + 
                                 ( (k*(*ncols))/grid_size[1] ) * element_size;
                    /* The process has to make sure it does not try to send to
                       itself. If so it does a memory copy instead. */
                    if (0 == dest_id ) {
                        /* It is sending to itself */
                        dest_address = (*matrix)[j];
                        memcpy (dest_address, source_address, 
                               nlocal_cols * element_size);                  
                    } 
                    else {
                        /* It is sending to another process */
                        int blocksize = size_of_block(k,*ncols, grid_size[1]);
                        MPI_Send (source_address,blocksize, dtype,
                                  dest_id, 0, cart_comm);
                    }
                }
                else if (grid_id == dest_id) {
                         MPI_Recv ((*matrix)[j], nlocal_cols, dtype, 0,
                          0, cart_comm, &status);
                }
            } /* end for k */
        } /* end for j */
    } /* for i */

    if (grid_id == 0) 
        free (buffer);
    *errval = 0;
}

void collect_and_print_2dblock_matrix (
   void       **a,            /* IN -2D matrix */
   MPI_Datatype dtype,        /* IN -Matrix element type */
   int          m,            /* IN -Matrix rows */
   int          n,            /* IN -Matrix columns */
   MPI_Comm     grid_comm,    /* IN - Communicator */
   FILE         *file)        /* IN - File stream*/
{
   void      *buffer;         /* Room to hold 1 matrix row */
   int        coords[2];      /* Grid coords of process
                                 sending elements */
   int        element_size;    /* Bytes per matrix element */
   int        els;            /* Elements received */
   int        grid_coords[2]; /* Coords of this process */
   int        grid_id;        /* Process rank in grid */
   int        grid_period[2]; /* Wraparound */
   int        grid_size[2];   /* Dims of process grid */
   int        i, j, k;
   void      *laddr;          /* Where to put subrow */
   int        local_cols;     /* Matrix cols on this proc */
   int        p;              /* Number of processes */
   int        src;            /* ID of proc with subrow */
   MPI_Status status;         /* Result of receive */

   MPI_Comm_rank (grid_comm, &grid_id);
   MPI_Comm_size (grid_comm, &p);
   element_size = get_size (dtype);

   MPI_Cart_get (grid_comm, 2, grid_size, grid_period,
      grid_coords);
   local_cols = size_of_block(grid_coords[1], n, grid_size[1]);
   if (0 == grid_id)
      buffer = malloc ( n * element_size);

   /* For each row of the process grid */
   for (i = 0; i < grid_size[0]; i++) {
      coords[0] = i;

      /* For each matrix row controlled by the process row */
      for (j = 0; j < size_of_block(i,m, grid_size[0] ); j++) {

         /* Collect the matrix row on grid process 0 and
            print it */
         if (0 == grid_id) {
            for (k = 0; k < grid_size[1]; k++) {
               coords[1] = k;
               MPI_Cart_rank (grid_comm, coords, &src);
               els = size_of_block(k,n, grid_size[1]);
               laddr = buffer +
                  ((k*n)/grid_size[1]) * element_size;
               if (src == 0) {
                  memcpy (laddr, a[j], els * element_size);
               } else {
                  MPI_Recv(laddr, els, dtype, src, 0,
                     grid_comm, &status);
               }
            }
            print_vector (buffer, n, dtype, file);
            fprintf (file,"\n");
         } 
         else if (grid_coords[0] == i) {
            MPI_Send (a[j], local_cols, dtype, 0, 0,
               grid_comm);
         }
      }
   }
   if (0 == grid_id) {
      free (buffer);
      printf ("\n");
   }
}

