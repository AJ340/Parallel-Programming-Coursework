#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int leftmost_is_SameAsRightmost(int * digit_array, int da_size)
{
	if (digit_array[0] == digit_array[da_size-1])
		return 1;
	return 0;
}

int twoConsecutiveDigits_are_TheSame(int * digit_array, int da_size)
{
	for (int x=0; x<da_size-1; x++)
		if(digit_array[x] == digit_array[x+1])
			return 1;
	return 0;
}

int a_digitAppearMoreThanTwice(int * digit_array, int da_size)
{
	int * appear_count = malloc (10 * sizeof(int));
	int ret_Val = 0;
	for (int x = 0; x < 10; x++)
		appear_count[x] = 0;

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

int digitSum_is_Mod7111317(int * digit_array, int da_size)
{
	int digitSum = 0;
	for (int x=0; x<da_size; x++)
	{
		digitSum += digit_array[x];
	}

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

int numericValue_is_MultipleOfLength(int integer, int numOfDigits)
{
	if ((integer % numOfDigits) == 0)
		return 1;
	return 0;
}

int password_LessThan_4Digits(int * digit_array, int numOfDigits)
{
	if (numOfDigits < 4)
		return 1;
	return 0;
}

int passwordContainsProhibitedDigit(int * digit_array, int da_size, int prohib_D)
{
	if (prohib_D < 0)
		return 0;

	for (int x=0; x<da_size; x++)
	{
		if (digit_array[x] == prohib_D)
			return 1;
	}
	return 0;
}


int validAgainstConstraints (int q_integer, int da_size, int prohibited_Digit)
{
	int ret_Val = 1;
	int * digit_array = malloc (da_size * (sizeof(int)));
	int current = da_size-1;
	int count = q_integer;
	while (count)
	{
		digit_array[current] = count % 10;
		current --;
    	count /= 10;
	}

	if (password_LessThan_4Digits(digit_array, da_size) ||
		passwordContainsProhibitedDigit(digit_array, da_size, prohibited_Digit) ||
		leftmost_is_SameAsRightmost(digit_array, da_size) ||
		twoConsecutiveDigits_are_TheSame(digit_array, da_size) ||
		a_digitAppearMoreThanTwice(digit_array, da_size) ||
		digitSum_is_Mod7111317(digit_array, da_size) ||
		numericValue_is_MultipleOfLength(q_integer, da_size))
		ret_Val = 0;
	free(digit_array);
	return ret_Val;
}

int main (int argc, char * argv[])
{
	int prohibited = -1;

	if (argc > 3)
	{
		printf ("Too many arguments");
		return 1;
	}
	if (argc == 3)
		prohibited = atoi(argv[2]);

	//printf("prohibited: %d \n", prohibited);
	int numberOfDigits = atoi(argv[1]);
	//printf("numberOfDigits: %d \n", numberOfDigits);
	int maxNumber = 0;
	int minNumber = 0;
	int ten_exp = 1;

	for (int i=0; i<numberOfDigits; i++)
	{
		maxNumber += 9 * ten_exp;
		minNumber = 1 * ten_exp;
		ten_exp *= 10;
	}

	//printf("maxNumber: %d \n", maxNumber);
	//printf("minNumber: %d \n", minNumber);



	int validCount = 0;
	for (int i= minNumber; i < maxNumber; i++)
	{	
		if (validAgainstConstraints(i, numberOfDigits, prohibited))
			validCount++;
	}

	printf("Total Number of Valid Passwords: %d \n", validCount);

}