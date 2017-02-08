/***************************************************************
Title: AQ_countPasswords_p.c
Author: Andres Quinones
Created on: October 18, 2016
Description: Counts the number of passwords that satisfy the pre-
determined constraints defined in Homework 4 for Parallel Compu-
ting with Professor Weiss at Hunter College
Purpose: First MPI programming assignment
Usage: combinations n [d] where 4 <= n <= 18 and d is a decimal digit
Build With: mpicc -Wall -o combinations AQ_countPasswords_p.c
	or if file is renamed 
		mpicc -Wall -o combinations <name_of_file>
***************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

// ************************************************************************************************************************************************************************
// Undesired Password Properties:
// Password should NOT have the properties that the function below checks for.

// numOfDigits is the number of elements in digit_array
// integer is the numeric value to check against numofdigits
// digitSum is the sum of the digits (presumably) in the password
int digitSum_is_Mod7111317(int digitSum);
int numericValue_is_MultipleOfLength(long long int integer, int numOfDigits);
//EDIT: Rest of properties written in validAgainstConstraints
// End Undesired Password Properties
// ************************************************************************************************************************************************************************

// Validator function.
// Takes needed values for the above constraints and checks if q_integer is a valid password.
// If password has undesired property then unvalid password (return 0)
// If new Undesired Property added above then it must be added into this function definition.
int validAgainstConstraints (long long int q_integer, int numOfDigits, int prohibited_Digit);


int main (int argc, char * argv[])
{
	//Using long long to support up to 18 digit passwords
	int prohibited = -1;
	long long int i;
	int id;
	int numProcesses;
	long long int maxNumber = 0;
	long long int validCount_Subtotal = 0;
	long long int validCount_GrandTotal = 0;
	double time = 0;
	double max_time = 0;
	int numberOfDigits;
	int invArgs=0;


	

	//Calculate largest number. argv[1] many 9s (ex. argv[1] == 4 then  9999)
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &numProcesses);

	//If root node process arguments. If error with args print usage.
	if (id == 0)
	{
		// if too many arguments or no argument for number of digits, print error and quit
		if (argc > 3 || argc < 2)
		{
			invArgs = -1;
			printf ("usage: %s n [d] where 4 <= n <= 18 and d is a decimal digit \n", argv[0]);
		}
		//If there is a prohibited digit, set prohibitied to that
		if (argc == 3)
			prohibited = atoi(argv[2]);
		//If there wasnt a supplied prohibited digit, prohibited remains -1 here. This is our flag that there was no supplied prohib digit.
		numberOfDigits = atoi(argv[1]);

		//If number of digits is too small or too big or prohibited is not a digit. Print error and quit
		if (numberOfDigits < 4 || numberOfDigits > 18 || prohibited > 9 || prohibited < -1)
		{
			invArgs = -1;
			printf ("usage: %s n [d] where 4 <= n <= 18 and d is a decimal digit. \n Note: d= -1 or no d argument allows all digits. \n", argv[0]);
		}

		//If args are valid, compute max number
		if (invArgs > -1)
		{
			long long int ten_exp=1;
			for (int i=0; i<numberOfDigits; i++, ten_exp *= 10)
				maxNumber += (9 * ten_exp);
		}
	}
	//Send these values \/\/ to all of the processes.
		MPI_Bcast(&numberOfDigits, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&prohibited, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&invArgs, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&maxNumber, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//If the args were valid perform parallel part. Otherwise do nothing.
	if (invArgs > -1)
	{
		
		//Start timing
		MPI_Barrier(MPI_COMM_WORLD);
		time = - MPI_Wtime();

		//Start Computation
		//Each mpi instance starts at it's id and jumps through all possibilities by a factor of p
		for (i =  id; i <= maxNumber; i += numProcesses)
		{	
			if (validAgainstConstraints(i, numberOfDigits, prohibited))
			{
				validCount_Subtotal++;
			}
		}
		//Reduce all subtotals into total in root process
		MPI_Reduce (&validCount_Subtotal, &validCount_GrandTotal, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		//Finish Computation

		//After reduction we have our answer in root so stop timing.
		time += MPI_Wtime();
		
		//Get max time of ANY process
		MPI_Reduce (&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		//If root process, print out answer 
		MPI_Barrier(MPI_COMM_WORLD);
		if(id == 0)
		{
			printf("%lld        %f seconds \n", validCount_GrandTotal, max_time);
		}
	}

	MPI_Finalize();
}


//Bunch of checks for digit sum == multiple of 7,11,13, or 17
int digitSum_is_Mod7111317(int digitSum)
{
	if ( (digitSum % 7) == 0 )
		return 1;
	if ( (digitSum % 11) == 0 )
		return 1;
	if ( (digitSum % 13) == 0 )
		return 1;
	if ( (digitSum % 17) == 0 )
		return 1;
	return 0;
}

//Check for number being a multiple of its length.
int numericValue_is_MultipleOfLength(long long int integer, int numOfDigits)
{
	if ((integer % (long long int) numOfDigits) == 0)
		return 1;
	return 0;
}

//If it satisfies all constraints return true >1
int validAgainstConstraints (long long int q_integer, int da_size, int prohibited_Digit)
{
	//Extract digits from integer and store them in an array for easier digit logic. 
	//If number is less than number of digits requested. Put 0s in front of the most significant bit.
	long long int count = q_integer;
	int i = da_size;
	int current_digit;
	int previous_digit = -1;
	int last_digit = count % 10;
	int digitSum = 0;
	int * appear_count = malloc (10 * sizeof(int));

	for (int x = 0; x < 10; x++)
		appear_count[x] = 0;

	//For each digit (Before returns within this loop, must free appear_count)
	while (i>0)
	{
		//Extract digits
		if (count > 0)
		{	
			current_digit = count % 10;
			count = count / 10;
		}
		//current_digit is a 0 that was in front of MSB
		else
			current_digit = 0;
	
		//If current digit is a prohibited digit, return false
		if (current_digit == prohibited_Digit)
		{
			free(appear_count);
			return 0;
		}
		//If current digit is the same as the previous digit return false
		//Note: for the first digit, this check doesnt pass because previous_digit
		//.      was initialized to -1
		if (current_digit == previous_digit)
		{
			free(appear_count);
			return 0;
		}
		// accumulate a sum of the digits
		digitSum += current_digit;
		//If we have seen 2 of current_digit so far then return false
		if (appear_count[current_digit] == 2)
		{
			free(appear_count);
			return 0;
		}
		//Otherwise, increment count for current_digit
		appear_count[current_digit]++;
		i--;
		//Update previous_digit to prepare for next call
		previous_digit = current_digit;
	}
	// appearcount not needed, free memory
	free(appear_count);
	//If last_digit == first digit return false
	if (last_digit == current_digit)
		return 0;
	//If sum of digits is a mod of the digits below return false
	if (digitSum_is_Mod7111317(digitSum))
		return 0;
	// If numeric value is a multiple of the length, return false
	if (numericValue_is_MultipleOfLength(q_integer, da_size))
		return 0;

	//If password doesnt have any of the above properties, then it is a valid password. 
	return 1;
}