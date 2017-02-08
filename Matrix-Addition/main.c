/***************************************************************
Title: AQ_HW6.c (main.c)
Author: Andres Quinones
Created on: Noveber 24th, 2016
Description: Adds two matrices of doubles in parallel
Purpose: Second MPI programming project
Usage: addmatrices <A> <B> [optional_output_file]
            where A and B are in binary format"
            sum_matrix is the name of an optional file to output the sum" 
           [Default behavior is to print output to stdout] \n"

Build With: make or make all. or the following.

mpicc -Wall -o addmatrices AQ_HW6.c matrix_io.c utilities.c -lm
  or if file is renamed 
    mpicc -Wall -o [new_program_name] AQ_HW6.c matrix_io.c utilities.c -lm

Executable name is addmatrices unless desired otherwise.



This program utilizes MPI routines implimented by Professor Stewart Weiss
at Hunter College with his permission.
These functions are contained in the included utilities.h and matrix_io.h
files.
Note: Some of these functions were slightly modified to better suit the task.
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

////////////////////////////////////////////////////////////
// These two files were provided by Professor Stewart Weiss.
//
// Contains commonly desired functions for matricies and vectors.
#include "utilities.h"
//
// Contains Matrix IO routines
#include "matrix_io.h"
////////////////////////////////////////////////////////////

#define PI 3.14159265358979323846


//Declaration for handle_Error. errorcode is code from read and distribute.
void handle_Error(int errcode);
void printUsage(char * args[]);



int main (int argc, char * argv[])
{	
	  double** A;         /* matrix to be analyzed       */
    double*  A_storage; /* backing storage for matrix    */
    double** B;         /* matrix to be analyzed       */
    double*  B_storage; /* backing storage for matrix    */
    int     nrowsA;            /* number of rows in matrix      */
    int     ncolsA;            /* number of columns in matrix   */
    int     nrowsB;            /* number of rows in matrix      */
    int     ncolsB;            /* number of columns in matrix   */
    int     id;               /* process rank                  */
    int     numProcesses;     /* number of processes           */
    int     i, j;             /* loop indices                  */
    int     error;            /* error exit value of calls     */
    int     rows;             /* number of rows on this process*/
    int     cols;             /* number of rows on this process*/

    int       grid_id;        /* rank of process within Cartesian grid        */
    int       grid_size[2];   /* grid dimensions                              */
    int       grid_coords[2]; /* process's coordinates in Cartesian grid      */
    MPI_Comm  grid_comm;      /* Cartesian communicator for grid of procs     */
    MPI_Comm  row_comm;       /* communicator for all processes in a grid row */
    MPI_Comm  col_comm;       /* communicator for all processes in a grid col */
    int       periodic[2];    /* array of wraparound flags when creating grid */

    double    max_seconds;    /* for timing program              */
    double    seconds;        /* Elapsed time for matrix-vector multiply */



	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &numProcesses);


	//Arg Checking on root
  if (id == 0)
  {

    if (argc == 3 || argc == 4)
      error = 0;

    else
    {
      error = 1;
      handle_Error(error);
      printUsage(argv);
      if (argc == 2)
        fprintf(stderr, "-- Missing B_matrix_file. \n");
    }
  }

  MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);

// If arguments good then continue procedure
  if (error == 0)
  {
    // Start timing
    seconds = - MPI_Wtime();
    
    /* Set up the communicators */
    grid_size[0] = 0;
    grid_size[1] = 0;
    MPI_Dims_create (numProcesses, 2, grid_size);
    
    /* The grid will not be periodic */
    periodic[0]  = 0;
    periodic[1]  = 0;

    /* Create the Cartesian commincator from all processes in program. */
    MPI_Cart_create (MPI_COMM_WORLD, 2, grid_size, 
                     periodic, 1, &grid_comm);

    /* Get the rank within the Cartesian communicator of this process */
    MPI_Comm_rank (grid_comm, &grid_id);

    /* Get the coordinates within the grid of this process */
    MPI_Cart_coords (grid_comm, grid_id, 2, grid_coords);

    /* Make comm for rows in grid comm */
    MPI_Comm_split (grid_comm, grid_coords[0], grid_coords[1], &row_comm);
   
    /* Make comm for cols in grid comm */
    MPI_Comm_split (grid_comm, grid_coords[1], grid_coords[0], &col_comm);

    // Read first matrix argument.
    read_and_distribute_2dblock_matrix (argv[1], (void *) &A,
      (void *) &A_storage, MPI_DOUBLE, &nrowsA, &ncolsA, &error, grid_comm);

    // If we have an error value at this point. Root Handles the code, otherwise continue procedure
    if (error != 0)
    {
      if (id == 0)
      {
        fprintf(stderr, "ERROR WITH %s: ", argv[1]);
        handle_Error(error);
        printUsage(argv);
      }
    }
    // printf(" READ SUCCESS \n");

    else 
      {
        // Read second matrix argument
        read_and_distribute_2dblock_matrix (argv[2], (void *) &B, (void *) &B_storage, MPI_DOUBLE, 
                                              &nrowsB, &ncolsB, &error, grid_comm);
        // Another error wall for second matrix success.
        if (error != 0)
        {
         if (id == 0)
         {
          fprintf(stderr, "ERROR WITH %s : ", argv[2]);
          handle_Error(error);
          printUsage(argv);
         }
        }

        // BOTH Read and distributes successful
        else 
        {
          // Make sure matrices are the same size. Print error if they arent
          if (id == 0)
            {
              if (ncolsA != ncolsB || nrowsA != nrowsB)
              {
                handle_Error(-4);
                error = -4;
                printUsage(argv);
              }
            }

            MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // Final Error check
            if (!error)
           {
              //No more possibilities for error. Completely safe to do work. 
              //Calculate local rows and columns.
              rows = size_of_block(grid_coords[0],nrowsA,grid_size[0]); 
              cols = size_of_block(grid_coords[1],ncolsA,grid_size[1]); 

              //Each processor does their portion of the matrix sum.
              for (i = 0; i < rows; i++)
                for (j = 0; j < cols; j++)
                  A[i][j] = A[i][j] + B[i][j];
                  
              // Stop timing.    
              seconds += MPI_Wtime();
              // Reduce max elapsed time to root
              MPI_Reduce (&seconds, &max_seconds, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

              //Depending on the existence of a third argument, each processor takes part in a collect_and_print operation but with different streams.
              if (argc == 4)
              {   
                  //Print to output file. Creates file if does not already exist
                  FILE * output_file;
                  output_file = fopen(argv[3], "w+");
                  collect_and_print_2dblock_matrix ( (void**) A, MPI_DOUBLE, nrowsA, ncolsA, grid_comm, output_file);
                  fclose(output_file);
              }
              else
                  //Print to standard output
                  collect_and_print_2dblock_matrix ( (void**) A, MPI_DOUBLE, nrowsA, ncolsA, grid_comm, stdout);

              // If root process then
              if (0 == id) {
              // also print elapsed time and work
                  fprintf (stderr, "Elapsed Time = %f sec ", max_seconds);
                  fprintf (stderr, " Work = %6.4f mflops \n", 2*ncolsA*ncolsA/(1000000.0*max_seconds));
              }
            }
          }
      }
  }
    MPI_Finalize();
}

// Original Read and distribute mainly returned -1. Slightly modified to code to return different errorcodes depending on the problem.
//All errors print out to stderr.
void handle_Error(int errcode)
{
  switch (errcode)
  {
  case (1) :
    {
      fprintf(stderr, "INVALID_ARGUMENTS\n");
      break;
    }
  case (-1) :
    {
      fprintf(stderr, "MPI_UNINIIALIZED ERROR\n");
      break;
    }
  case (-2) :
    {
      fprintf(stderr, "GET_SIZE Error\n");
      break;
    }
  case (-3) :
    {
      fprintf(stderr, "FILE_OPEN_ERROR\n");
      break;
    }
  case (-4) :
    {
      fprintf(stderr, "ARRAYS_NOT_SAME_SIZE\n");
      break;
    }
  default :
    {
      fprintf(stderr, "OTHER ERROR\n");
    }
  }
}

void printUsage(char * argv[])
{
  fprintf(stderr, "Usage: %s <A> <B> [sum_matrix]"
                  "\n    where A and B are in binary format"
                  "\n    sum_matrix is the name of an optional file to output the sum" 
                  "\n  [Default behavior is to print output to stdout] \n", argv[0]);
}

























