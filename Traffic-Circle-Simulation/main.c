/***************************************************************
Title: AQ_HW7.c (main)
Author: Andres Quinones
Created on: December 8th, 2016
Description: Simulates flow of traffic
Purpose: Fourth MPI programming project
Usage: 
Build With: make or make all. or the following.

mpicc -Wall -o traffic_circle AQ_HW7.c -lm

Executable name is traffic_circle unless desired otherwise.



This program utilizes some functions mostly used for monte carlo methods 
implimented by Professor Stewart Weiss
at Hunter College with his permission.
These functions are noted below.
***************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

//////////////////////////////////////////////////////////////////////////////
// These functions were implimented by Professor Stewart Weiss in his gen_exponential.c program
// 

/** init_sequence()  initializes a custom lag table for C random()
   Use the initstate() library function to create a custom lag table for 
   the random() function. random() uses a Lagged Fibonacci Generator which
   can be customized to have a bigger lag table. The argument to init_sequence
   is the number of entries in the lag table. 
   The memory must be released by calling finalize().
*/
char*  init_sequence( int state_size );

// Clear memory allocated for lag-table
void  finalize ( char* state );

/* Returns an exponentially distributed random value with parameter
   lambda using the random() library function using the method of
   inversion. 
*/

double gen_exponential_number( double lambda );
//////////////////////////////////////////////////////////////////////////////



/* Returns an uniform random value */
double gen_uniform_rand_var();

/* Modifies value of toggle. Toggle == 0 -> 1. Otherwise toggle -> 0 */
void flip_toggle(int * toggle);

/* Moves cars inside of array_current (old_traffic_circle). 
	If they reach their destination they leave the circle
	This updated circle is stored in array_next (updated_traffic_cricle)  */
void move_internal_cars(int * array_current, int * array_next, int size);

/* Simple function to set a position pos, in traffic_circle to a destination value*/
void add_car (int * traffic_circle, int pos, int destination);

/* Generates a destination for a car using a uniform random varible
	and comparing to the row that represents the probabilities for each exit in probabilities*/
int generate_destination(double * probabilities_2d, int * arry_exits, int row);

void handle_Error(int errcode);

char index_to_entrance(int index);

void printUsage(char * argv[]);


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//        START OF MAIN
int main (int argc, char * argv[])
{
	double f_time_btwn_arrivals[4];
	double d_trans_prob[16];

	int circle_current[16] = {-1,-1,-1,-1,			// Variables to represent current_timestep and next_timestep traffic circles
							  -1,-1,-1,-1,
							  -1,-1,-1,-1,
							  -1,-1,-1,-1};

	int circle_next[16] = {-1,-1,-1,-1,
						   -1,-1,-1,-1,  
						   -1,-1,-1,-1,  
						   -1,-1,-1,-1};

	int entrance_offsets[4] = {0,12,8,4};			// N,E,S,W offsets in traffic circle
	int arrival_counts[4] = {0,0,0,0};				// Counts arrivals at entrance
	int wait_counts[4] = {0,0,0,0};					// Counts number of cars that had to wait
	int queue[4] = {0,0,0,0};						// Queues for each entrance
	double avg_queue_accum[4] = {0.0,0.0,0.0,0.0};	// Local average queue sizes over time
	double global_avg_queue_accum[4] = {0.0,0.0,0.0,0.0};	// Global average queue sizes over time

	int i,j;										// Loop Indexes

	double next_arriv_time[4] = {0.0,0.0,0.0,0.0};	// Array to store next arrival times for each entrances
	int convergence = 0;							// Flag for convergence
	int global_max_iteration;						// Global Max_Iterations across all processes
	int local_max_iteration;						// Local Max_Iterations for each process
	int current_iteration = 1;						// Current iteration
	int curr_next_toggle = 0;						// Flag to help track which between circle_current and next is the last updated.
	double epsilon;									// Epsilon for determining convergence

	double prob_of_waiting[4] = {0.0,0.0,0.0,0.0};	// Local average probabibilities to wait
	double global_prob_of_waiting[4] = {0.0,0.0,0.0,0.0};  // Global average probabibilities to wait
    int     random_statesize = 1000;   /* size of lag table            */
    char    *rand_initstate;   /* address of lag table allocated       */

	int     id;               /* process rank                  */
    int     numProcesses;     /* number of processes           */
    int     error = 0;

    rand_initstate = init_sequence( random_statesize );  // Initialize lag table of size 1000

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &numProcesses);


	//Must do this on a per processor basis

    // Root process reads file and processes any errors
    if (id == 0)
    {
    	//Handle Params
    	if (argc == 2)
    	{
    		global_max_iteration = 1000000;
    		epsilon = 0.000000001;
    	}
    	else if (argc == 3)
    	{
    		global_max_iteration = atof(argv[2]);
    		epsilon = 0.000000001;
    	}
    	else if (argc == 4)
    	{
    		global_max_iteration = atof(argv[2]);
    		epsilon = atof(argv[3]);
    	}
    	else
    	{
    		error = -3;
    	}

    	if(error >= 0)
    	{
	    	FILE * input_file;
	    	input_file = fopen(argv[1], "r");
	    	if (input_file == NULL)
	    		error = -4;
	  		else
	  		{
	  			// Input f
		  		for (i=0; i < 4; i++)
		  		{
		  			int x = fscanf(input_file, "%lf", &(f_time_btwn_arrivals[i]));
		  			if (x != 1)
		  			{ 	
		  				error = -1;
		  				break;
		  			}
		  		}
		  		// If no error
		  		if(error >= 0)
		  		{
		  			// Proceed to input d
			  		for (i=0; i < 4; i++)
			  		{
			  			double sum = 0.0;
			  			for (j=0; j < 4; j++)
			  			{
			  				int it = i*4+j;
			  				int x = fscanf(input_file, "%lf", &(d_trans_prob[it]));
			  				if (x != 1)
			  				{ 
			  					error = -1;
			  					break;
			  				}
			  				else 
			  					sum += d_trans_prob[it];
			  			}
			  			if (error == -1)
			  				break;
			  			if (sum != 1.0)
			  			{
			  				error = -2;
			  				break;
			  			}
			  		}
			  		int x = fclose(input_file);
		  			if (x != 0)
		  				error = -5;
		  		}
		  	}		
	  	}
    }
    // Broadcast error value
    MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // If there is an error. Root handles error and all processes don't do work
    if (error < 0)
    {
    	if (id == 0)
    	{
    		handle_Error(error);
    		printUsage(argv);
    	}
    }
    // If no error -> Safe to do work.
    else 
	{
	    // Broadcast all needed data
	  	MPI_Bcast(f_time_btwn_arrivals, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  	MPI_Bcast(d_trans_prob, 16, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  	MPI_Bcast(&global_max_iteration, 1, MPI_INT, 0, MPI_COMM_WORLD);

	    // Each process computes local_max_iteration
	    local_max_iteration = global_max_iteration / numProcesses;
	    if (id < global_max_iteration % numProcesses)
	    	local_max_iteration++;

		// Loop until convergence or max iterations
		while (!convergence && current_iteration <= local_max_iteration)
		{
			int convergence_count = 0;
			int entrace_arrivals[4] = {0,0,0,0};	// Keep track of arrivals per timestep
			int * start, * next;

			// Logic to make: Previous iteration's next becomes this iteration's current and vise versa
			if (curr_next_toggle){
				start = circle_next;
				next = circle_current;
			}
			else {
				start = circle_current;
				next = circle_next;
			}


			// Cars arrive at entrances
			for (i=0; i<4; i++)
			{
				int c = 0;
				// Invert mean. Avoid calculating multiple times.
				double l = 1.0 / f_time_btwn_arrivals[i];
				// As long as next arrival time is < this current time step
				// Note: low f values might lead to many cars arriving at once 
				while ( (double) current_iteration > next_arriv_time[i])
				{
					c++;
					// When more than 1 car arrives at one time step. 
					// Because only 1 car can enter per timestep
					// the other cars will for sure be waiting before entering
					if(c > 1)
						wait_counts[i]++;
					arrival_counts[i]++;  	// Increment count for that entrance
					queue[i]++;  // Car enters queue
					entrace_arrivals[i] = 1;	// Flag car arrived for the entrance
					// Compute next arrival time with exponentially distributed random variable.
					next_arriv_time[i] += gen_exponential_number(l);	//
				}
			}
			// Cars inside circle move forward;
			move_internal_cars(start, next, 16);

			// For each entrance
			for (i=0; i<4;i++)
			{
				// If cars are waiting
				if(queue[i] > 0)
				{
					int slot_to_check = 0; 

					// We must check adjacent slot. If entrance is at 0. We check 15.
					if(entrance_offsets[i] == 0)
						slot_to_check = 15;
					// Otherwise we check entrance-1
					else
						slot_to_check = entrance_offsets[i] - 1;

					//If there is no incoming car, or incoming car is leaving. Then safe to enter circle.
					if(next[slot_to_check] < 0 || next[slot_to_check] == entrance_offsets[i])
					{
						// The car in front of queue enters the circle.
						// Note: A destination is generated from d_trans_prob
						add_car(next, entrance_offsets[i], generate_destination(d_trans_prob, entrance_offsets, i));
						queue[i]--;
					}
					// Otherwise car must wait. 
					// Any cars that just arrived. Had to wait before entering.
					else
					{
						if(entrace_arrivals[i] == 1)
							wait_counts[i]++;
					}
				}
			}

			// Compute new average from all new queue sizes. 
			// If this new average only differed by less than epsilon then it converged
			for (i=0; i<4; i++)
			{
				double old_avg = avg_queue_accum[i];
				avg_queue_accum[i] = (avg_queue_accum[i] * (current_iteration-1) + queue[i]) / current_iteration;
				double difference_in_avg = fabs(old_avg - avg_queue_accum[i]);
				// If difference in average is very tiny it converged. 
				// Difference > 0 ignores the first few steps when avg might be 0.
				if ( difference_in_avg < epsilon && difference_in_avg > 0.0) 
					convergence_count++;
			}


			// If all entrances converged. Signal to break
			if (convergence_count == 4)
				convergence = 1;

			current_iteration ++;
			// Flag to swap current and next 
			flip_toggle(&curr_next_toggle);
		}

		// Compute prob_of_waiting by #carswaited / #carsarrived
		for (i=0; i<4; i++)
		{
			prob_of_waiting[i] = (double) wait_counts[i] / (double) arrival_counts[i];
		}

		// Reduce prob_of_waiting, average
			 MPI_Reduce ( prob_of_waiting, global_prob_of_waiting, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	   		 MPI_Reduce ( avg_queue_accum, global_avg_queue_accum, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	   
	   	// Root displays output after computing average from sum of all proccesses local averages
	    if (id == 0)
	    {
			for (i = 0; i < 4; i++)
			{
				char entrance = index_to_entrance(i);
				global_prob_of_waiting[i] /= numProcesses;
	    		global_avg_queue_accum[i] /= numProcesses;
	    		fprintf(stderr, " %c: Avg Wait Probability: %0.2f \n", entrance, global_prob_of_waiting[i]);
	    		fprintf(stderr, " %c: Avg Queue Size: %0.2f \n", entrance, global_avg_queue_accum[i]);
			}
		}

	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	finalize(rand_initstate);
	return 0;
}
//        END OF MAIN
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////



char*  init_sequence( int state_size )
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

double gen_exponential_number( double lambda )
{
    double u = (double) (random()) / RAND_MAX; 
    return ( - log ( u ) / lambda);
}

/* Returns an uniform random value */
double gen_uniform_rand_var()
{
    return (double) (random()) / RAND_MAX; 
}

void flip_toggle(int * toggle)
{
	if ( *toggle == 0)
		*toggle = 1;
	else
		*toggle = 0;
}

void move_internal_cars(int * array_current, int * array_next, int size)
{
	for (int i = 0; i < size-1; i++)
	{
		if (array_current[i] > -1)
		{
			if (array_current[i] == i+1)
				array_next[i+1] = -1;
			else
				array_next[i+1] = array_current[i];
	    }
	}
	if (array_current[size-1] > -1)
	{
		if (array_current[size-1] == 0)
				array_next[0] = -1;
		else 
			array_next[0] = array_current[size-1];
	}
}

void add_car (int * traffic_circle, int pos, int destination)
{
	traffic_circle[pos] = destination;
}

int generate_destination(double * probabilities_2d, int * arry_exits, int row)
{
	double u = gen_uniform_rand_var();
	double sum = 0;
	int i;
	for (i=0; i < 4; i++)
	{
		sum += probabilities_2d[row * 4 + i];
		if (u < sum)
			break;
	}
	return arry_exits[i];
}

void handle_Error(int errcode)
{
  switch (errcode)
  {
  case (-3) :
    {
      fprintf(stderr, "INVALID_ARGUMENTS\n");
      break;
    }
  case (-1) :
    {
      fprintf(stderr, "READ_ERROR\n");
      break;
    }
  case (-2) :
    {
      fprintf(stderr, "INPUT_FILE_ROWS_DONT_SUM_TO_1\n");
      break;
    }
  case (-4) :
    {
      fprintf(stderr, "FILE_OPEN_ERROR\n");
      break;
    }
  case (-5) :
    {
      fprintf(stderr, "FILE_CLOSE_ERROR\n");
      break;
    }
  default :
    {
      fprintf(stderr, "OTHER ERROR\n");
    }
  }
}

char index_to_entrance(int index)
{
  switch (index)
  {
  case (0) :
    {
    	return 'N';
    }
  case (1) :
    {
    	return 'E';
    }
  case (2) :
    {
    	return 'S';
    }
  case (3) :
    {
    	return 'W';
    }
  default :
    {
		return 'O';
    }
  }
}

void printUsage(char * argv[])
{
  fprintf(stderr, "Usage: %s <input> [max_iterations] [epsilon_value]"
                  "\n    where input is a file containing parameters for D and F."
                  "\n                      for more info see project description)"
                  "\n    Optional: max_iterations is the max time steps to run in simulation in total over all processes" 
                  "\n    Optional: epsilon_value is the threshold for convergence." 
                  "\n  [Default: max_iterations = 1,000,000 Epsilon = 0.000000001  ] \n", argv[0]);
}


