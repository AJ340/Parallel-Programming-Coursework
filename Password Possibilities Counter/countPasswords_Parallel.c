/***************************************************************
Title: AQ_countpasswords_p.c
Author: Andres Quinones
Created on: October 18, 2016
Description: Counts the number of passwords that satisfy the pre-
determined constraints defined in Homework 4 for Parallel Compu-
ting with Professor Weiss at Hunter College
Purpose: First MPI programming assignment
Usage: combinations n [d] where 4 <= n <= 18 and d is a decimal digit
Build With: mpicc -Wall -o combinations AQ_countpasswords_p.c
	or if file is renamed 
		mpicc -Wall -o combinations <name_of_file>
***************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

// digit_array is an array of digits for some number i. 
// da_size is the number of elements in digit_array
int leftmost_is_SameAsRightmost(int * digit_array, int da_size);
int twoConsecutiveDigits_are_TheSame(int * digit_array, int da_size);
int a_digitAppearMoreThanTwice(int * digit_array, int da_size);
int digitSum_is_Mod7111317(int * digit_array, int da_size);
//prohib_D is a digit that isnt allowed in the password (represented in digit_array)
int passwordContainsProhibitedDigit(int * digit_array, int da_size, int prohib_D);
//integer is the numeric value to check against numofdigits
int numericValue_is_MultipleOfLength(long long int integer, int numOfDigits);

//Takes needed values for the above constraints and checks if q_integer is a valid password.
int validAgainstConstraints (long long int q_integer, int da_size, int prohibited_Digit);


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
	double startime = 0;
	double endtime = 0;

	// if too many arguments or no argument for number of digits, print error and quit
	if (argc > 3 || argc < 2)
	{
		printf ("usage: %s n [d] where 4 <= n <= 18 and d is a decimal digit \n", argv[0]);
		return 0;
	}
	//If there is a prohibited digit, set prohibitied to that
	if (argc == 3)
		prohibited = atoi(argv[2]);
	//If there wasnt a supplied prohibited digit, prohibited remains -1 here. This is our flag that there was no supplied prohib digit.
	int numberOfDigits = atoi(argv[1]);

	//If number of digits is too small or too big or prohibited is not a digit. Print error and quit
	if (numberOfDigits < 4 || numberOfDigits > 18 || prohibited > 9 || prohibited < 0)
	{
		printf ("usage: %s n [d] where 4 <= n <= 18 and d is a decimal digit \n", argv[0]);
		return 0;
	}

	//Calculate largest number. argv[1] many 9s (ex. argv[1] == 4 then  9999)
	long long int ten_exp=1;
	for (int i=0; i<numberOfDigits; i++, ten_exp *= 10)
		maxNumber += (9 * ten_exp);

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &numProcesses);
	//Start timing
	startime = MPI_Wtime();
	//Each mpi instance starts at it's id and jumps through all possibilities by a factor of p
	for (i = (long long int) id; i < maxNumber; i += numProcesses)
	{	
		if (validAgainstConstraints(i, numberOfDigits, prohibited))
			validCount_Subtotal++;
	}
	//Reduce all subtotals into total in root process
	MPI_Reduce (&validCount_Subtotal, &validCount_GrandTotal, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	//After reduction we have our answer so stop timing.
	endtime = MPI_Wtime();
	//If root process, print out answer
	if(id == 0)
		printf("%lld        %f seconds \n", validCount_GrandTotal, endtime-startime);
	MPI_Finalize();
}

//Check if leftmost == rightmost
int leftmost_is_SameAsRightmost(int * digit_array, int da_size)
{
	if (digit_array[0] == digit_array[da_size-1])
		return 1;
	return 0;
}

//Check for consecutive repeating digits
int twoConsecutiveDigits_are_TheSame(int * digit_array, int da_size)
{
	for (int x=0; x<da_size-1; x++)
		if(digit_array[x] == digit_array[x+1])
			return 1;
	return 0;
}

//Check for digit appearing more than 3 times
int a_digitAppearMoreThanTwice(int * digit_array, int da_size)
{
	int * appear_count = malloc (10 * sizeof(int));
	int ret_Val = 0;
	for (int x = 0; x < 10; x++)
		appear_count[x] = 0;

	//O(n). Store seen digits in appearcount. Whenever appearcount is > 2 for a digit we pass, return 1.
	for (int x=0; x<da_size; x++)
	{
		if (appear_count[digit_array[x]] == 2)
		{
			ret_Val = 1;
			break;
		}
		appear_count[digit_array[x]]++;
	}
	free(appear_count);
	return ret_Val;
}

//Bunch of checks for digit sum == multiple of 7,11,13, or 17
int digitSum_is_Mod7111317(int * digit_array, int da_size)
{
	int digitSum = 0;
	for (int x=0; x<da_size; x++)
		digitSum += digit_array[x];

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

//Check for password containing a prohibited digit
int passwordContainsProhibitedDigit(int * digit_array, int da_size, int prohib_D)
{
	//If prohib is negative then no provided prohibited digit so return false
	if (prohib_D < 0)
		return 0;
	//Check for prohibited digit accross all digits
	for (int x=0; x<da_size; x++)
	{
		if (digit_array[x] == prohib_D)
			return 1;
	}
	return 0;
}

//If it satisfies all constraints return true >1
int validAgainstConstraints (long long int q_integer, int da_size, int prohibited_Digit)
{
	int ret_Val = 1;
	int * digit_array = malloc (da_size * (sizeof(int)));
	//Extract digits from integer and store them in an array for easier digit logic. 
	//If number is less than number of digits requested. Put 0s in front of the most significant bit.
	for (int i=da_size-1, count = q_integer; i >= 0; i--,  count /= 10)
	{
		if (count == 0)
		{
			digit_array[i] = 0;
			continue;
		}
		digit_array[i] = count % 10;
	}
	//If it has a property it shouldnt have, return false. Otherwise true
	if (passwordContainsProhibitedDigit(digit_array, da_size, prohibited_Digit) ||
		leftmost_is_SameAsRightmost(digit_array, da_size) ||
		twoConsecutiveDigits_are_TheSame(digit_array, da_size) ||
		a_digitAppearMoreThanTwice(digit_array, da_size) ||
		digitSum_is_Mod7111317(digit_array, da_size) ||
		numericValue_is_MultipleOfLength(q_integer, da_size))
		ret_Val = 0;
	free(digit_array);
	return ret_Val;
}