/***************************************************************
Title: AQ_HW5.c (main.c)
Author: Andres Quinones
Created on: November 6, 2016
Description: Computes what regions of a landscape (represented by
a matrix in binary format) are in the sun or shade for a given angle
that the suns rays make with the ground.
Purpose: Second MPI programming project
Usage: sunShine <theta> <matrix_filename>
theta is the angle the sun shines. 0 < theta < PI
matrix_filename is a matrix of heights in binary format
Build With: make or make all. or the following.

mpicc -Wall -o sunShine AQ_HW5.c matrix_io.c utilities.c -lm
  or if file is renamed 
    mpicc -Wall -o sunShine AQ_HW5.c matrix_io.c utilities.c -lm

Executable name is sunShine unless desired otherwise.



This program utilizes MPI routines implimented by Professor Stewart Weiss
at Hunter College with his permission.
These functions are contained in the included utilities.h and matrix_io.h
files.
***************************************************************/

#include <math.h>
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
//Contains Matrix IO routines
#include "matrix_io.h"
////////////////////////////////////////////////////////////
#define PI 3.14159265358979323846


//Declaration for handle_Error. errorcode is code from read and distribute.
void handle_Error(int errcode);


int main (int argc, char * argv[])
{	
	int** A;         /* matrix to be analyzed       */
    int*  A_storage; /* backing storage for matrix    */
    int     nrows;            /* number of rows in matrix      */
    int     ncols;            /* number of columns in matrix   */
    int     id;               /* process rank                  */
    int     numProcesses;     /* number of processes           */
    int     i, j;             /* loop indices                  */
    int     error;            /* error exit value of calls     */
    int     rows;             /* number of rows on this process*/
    double shadow_height_s = 0.0; /* height of shadow at start*/
    double shadow_height_m = 0.0; /* height of shadow in middle of block*/
    double theta = 0.0;    /* angle of sun*/
    double slope = 0.0;    /* slope of sun's light. */

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &numProcesses);

	//Arg Checking on root
  if (id == 0)
  {
    //If enough arguments
      if (argc == 3)
      {
        theta = atof(argv[1]);
        // Theta is 0.0 here -> arguments was in wrong order or user input 0.0. Invalid.
        // Theta < 0.0. Argument invalid. Set theta to invalid value
        // Theta > PI. Argument too big. Set theta to invalid value.
        if (theta <= 0.0 || theta > PI)
        {
          theta = 0.0;
          printf("Usage: %s <theta> <matrix_filename>\n theta is the angle the sun shines. 0 < theta < PI \n matrix_filename is a matrix of heights in binary format. Heights in matrix are expected to be positive integers. \n", argv[0]);
        }
        //Theta is valid -> get slope
        else
          slope = tan(theta);
      }
    //If wrong number of arguments
      else 
      {
        theta = 0.0;
        printf("Usage: %s <theta> <matrix_filename>\n theta is the angle the sun shines. 0 < theta < PI \n matrix_filename is a matrix of heights in binary format. Heights in matrix are expected to be positive integers. \n", argv[0]);
      }
  }

  MPI_Bcast(&theta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&slope, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

// If theta is valid do work. Otherwise do nothing and all processes simply finalize.
  if (theta != 0.0)
  {
    read_and_distribute_matrix_byrows (argv[2], 
                       (void *) &A,
                       (void *) &A_storage, 
                       MPI_INT, &nrows, &ncols, &error, 
                       MPI_COMM_WORLD);
    //If read and distribute fails. Do nothing. Have root display error.
    if (error != 0)
    {
      if (id == 0)
        handle_Error(error);
    }
    
    //Read and distribute successful -> Do work.
    else 
    {
      rows = size_of_block(id,nrows,numProcesses); 
     
     //If sun is shining from east
     if (theta > PI/2)
   		{
        //If shining from east. Invert angle. Shadow gets cast on 180degrees - theta.
        theta = PI - theta;
        double slope = tan(theta);
        //Process columns east to west.
   			for (i = 0; i < rows; i++) 
      	 for (j = ncols-1; j > -1; j--)
      		{
      		  if (j==ncols-1)
      		    { //NOTHING CASTING SHADOW ONTO THIS BLOCK. IN SUN CAST SHADOW TO OTHERS.
                shadow_height_s = (double) A[i][j];
                A[i][j] = 0;
      		  	}

      		  else
      		  	{ //Might be in shade. Caclulate height of shadow in middle of block.
      		  		shadow_height_m = shadow_height_s - (slope/2);
      		  		// If height of block is greater than or equal to height of shadow at 50%
                // Then shadow is casted on less than half of the block. 
                // This means that the block is in the sun and casts a new shadow to further blocks.
      		  		if ((double) A[i][j] > shadow_height_m)  
                  {
      		  			  shadow_height_s = (double) A[i][j];
                    A[i][j] = 0;
      		  		  }
      		  		// Otherwise. Block is in shade. but shadow might persist to future blocks or block might cast a bigger shadow
                else 
                  {
                    // Make shadow_height_s where the shadow ends at the end of the block.
                    // Also sets start of shadow for future blocks.
                    // Shadow might persist to future blocks.
                    shadow_height_s -= slope;
                    // If the height of the block is taller than the height of shadow at end of block
                    // Then the block is shaded between [50,100) percent. 
                    // Therefore this block casts a new taller shadow.
                    if (A[i][j] > shadow_height_s)
                      shadow_height_s = A[i][j];  
                    //Block in shade.
                    A[i][j] = 1; 
                  }
      		  	}
    		  }
    	}
      //If sun is shining from west
     else if (theta < PI/2)
  	   {
        //Process columns west to east
        for (i = 0; i < rows; i++) 
         for (j = 0; j < ncols; j++)
          {
            if (j==0)
              { //NOTHING CASTING SHADOW ONTO THIS BLOCK. IN SUN CAST SHADOW TO OTHERS.
                shadow_height_s = (double) A[i][j];
                A[i][j] = 0;
              }

            else
              { //Might be in shade. Caclulate height of shadow in middle of block.
                shadow_height_m = shadow_height_s - (slope/2);
                // If height of block is greater than height of shadow at 50%
                // Then shadow is casted on less than half of the block. 
                // This means that the block is in the sun and casts a new shadow to further blocks.
                if ((double) A[i][j] > shadow_height_m)  
                  {
                    shadow_height_s = (double) A[i][j];
                    A[i][j] = 0;
                  }
                // Otherwise. Block is in shade. but shadow might persist to future blocks or block might cast a bigger shadow
                else 
                  {
                    // Make shadow_height_s where the shadow ends at the end of the block.
                    // Also sets start of shadow for future blocks.
                    // Shadow might persist to future blocks.
                    shadow_height_s -= slope;
                    // If the height of the block is taller than the height of shadow at end of block
                    // Then the block is shaded between [50,100) percent. 
                    // Therefore this block casts a new taller shadow.
                    if (A[i][j] > shadow_height_s)
                      shadow_height_s = A[i][j];  
                    //Block in shade.
                    A[i][j] = 1; 
                  }
              }
          }
       }
      //Otherwise sun shining head on at 90 degrees. Everything is in the sun.
     else 
     	{
  	   		for (i = 0; i < rows; i++) 
      			for (j = 0; j < ncols; j++)
      				A[i][j] = 0;
      } 
      collect_and_print_matrix_byrows ((void **) A, MPI_INT, nrows, ncols,
             MPI_COMM_WORLD);  
    }
  }
    MPI_Finalize();
}

//Read and distribute mostly returns -1. But handle_error accomodates for
//Some codes that alloc_matrix might return (which is used by read and distribute.)
void handle_Error(int errcode)
{
  switch (errcode)
  {
  case (1) :
    {
      printf("MALLOC_ERROR\n");
      break;
    }
  case (2) :
    {
      printf("OPEN_FILE_ERROR\n");
      break;
    }
  case (3) :
    {
      printf("TYPE_ERROR\n");
      break;
    }
  case (4) :
    {
      printf("FILE_READ_ERROR\n");
      break;
    }
  default :
    printf("READ_AND_DISTRIBUTE FAILED\n");
  }
}

























