#ifndef __MATRIX_IO_H__
#define __MATRIX_IO_H__

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/*******************************************************************************
                     I/O Routines for Matrices Decomposed by Rows

*******************************************************************************/

void read_and_distribute_square_matrix_byrows (
          char        *filename,       /* [IN]  name of file to read          */
          void      ***matrix,         /* [OUT] matrix to fill with data      */
          void       **matrix_storage, /* [OUT] linear storage for the matrix */
          MPI_Datatype dtype,          /* [IN]  matrix element type           */
          int         *nrows,          /* [OUT] size of matrix                */
          int         *errval,         /* [OUT] success/error code on return  */
          MPI_Comm     comm);           /* [IN] communicator handle            */

void read_and_distribute_matrix_byrows (
          char        *filename,       /* [IN]  name of file to read          */
          void      ***matrix,         /* [OUT] matrix to fill with data      */
          void       **matrix_storage, /* [OUT] linear storage for the matrix */
          MPI_Datatype dtype,          /* [IN]  matrix element type           */
          int         *nrows,          /* [OUT] number of rows in matrix      */
          int         *ncols,          /* [OUT] number of columns in matrix   */
          int         *errval,         /* [OUT] success/error code on return  */
          MPI_Comm     comm);          /* [IN] communicator handle            */



void collect_and_print_matrix_byrows (
          void       **matrix,         /* [IN] matrix to print                */
          MPI_Datatype dtype,          /* [IN] matrix element type            */
          int          nrows,          /* [IN] number of rows in matrix       */
          int          ncols,          /* [IN] number of columns in matrix    */
          MPI_Comm     comm);          /* [IN] communicator handle            */


/*******************************************************************************
                 I/O Routines for Matrices Decomposed by 2D Blocks

*******************************************************************************/


void read_and_distribute_2dblock_matrix (
          char        *filename,       /* [IN]  name of file to read          */
          void      ***matrix,         /* [OUT] matrix to fill with data      */
          void       **matrix_storage, /* [OUT] linear storage for the matrix */
          MPI_Datatype dtype,          /* [IN]  matrix element type           */
          int         *nrows,          /* [OUT] number of rows in matrix      */
          int         *ncols,          /* [OUT] number of columns in matrix   */
          int         *errval,         /* [OUT] success/error code on return  */
          MPI_Comm     cart_comm);     /* [IN]  communicator handle           */


void collect_and_print_2dblock_matrix  (
				 void       **a,            /* IN -2D matrix */
				 MPI_Datatype dtype,        /* IN -Matrix element type */
				 int          m,            /* IN -Matrix rows */
				 int          n,            /* IN -Matrix columns */
				 MPI_Comm     grid_comm);    /* IN - Communicator */


#endif /* __MATRIX_IO_H__ */

