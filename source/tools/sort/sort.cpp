/**
 * @file	sort.cpp
 * @date	13.01.2012
 * @author	Johannes Neitsch
 * @brief	short_discription
 *
 * some sort algorithm
 *
 */

/** short funtion discription
 * @param	paramter_name	parameter_discription
 * @return	return_discription
 *
 * extendet discription
 *
 */

#include "sort.h"


//check if array is sorted
bool check(double* array, int size)
{
	int j;

	for (j=0;j<size-1;j++)
		if (array[j]>array[j+1])
			return false;
	return true;
}

//quicksort of an array
void quickSort(double array[], int start, int end)
{
        int i = start;                          // index of left-to-right scan
        int k = end;                            // index of right-to-left scan

        if (end - start >= 1)                   // check that there are at least two elements to sort
        {
                double pivot = array[start];       // set the pivot as the first element in the partition

                while (k > i)                   // while the scan indices from left and right have not met,
                {
                        while (array[i] <= pivot && i <= end && k > i)  // from the left, look for the first
                                i++;                                    // element greater than the pivot
                        while (array[k] > pivot && k >= start && k >= i) // from the right, look for the first
                            k--;                                        // element not greater than the pivot
                        if (k > i)                                       // if the left seekindex is still smaller than
                                swap(array, i, k);                      // the right index, swap the corresponding elements
                }
                swap(array, start, k);          // after the indices have crossed, swap the last element in
                                                // the left partition with the pivot
                quickSort(array, start, k - 1); // quicksort the left partition
                quickSort(array, k + 1, end);   // quicksort the right partition
        }
        else    // if there is only one element in the partition, do not do any sorting
        {
                return;                     // the array is sorted, so exit
        }
}


//change two values in one array
void swap(double array[], int index1, int index2)
// pre: array is full and index1, index2 < array.length
// post: the values at indices 1 and 2 have been swapped
{
	double temp = array[index1];           // store the first value in a temp
	array[index1] = array[index2];      // copy the value of the second into the first
	array[index2] = temp;               // copy the value of the temp into the second
}


//change two values in one array with refference array
void swap(double array[], int index1, int index2, int array_reverenz[])
// pre: array is full and index1, index2 < array.length
// post: the values at indices 1 and 2 have been swapped
{
	double temp = array[index1];           // store the first value in a temp
	array[index1] = array[index2];      // copy the value of the second into the first
	array[index2] = temp;               // copy the value of the temp into the second

	int tmp = array_reverenz[index1];
	array_reverenz[index1] = array_reverenz[index2];
	array_reverenz[index2] = tmp;


}



//quicksort with reference
void quickSort(double array[], int start, int end, int array_reverenz[])
{
        int i = start;                          // index of left-to-right scan
        int k = end;                            // index of right-to-left scan

        if (end - start >= 1)                   // check that there are at least two elements to sort
        {
                double pivot = array[start];       // set the pivot as the first element in the partition

                while (k > i)                   // while the scan indices from left and right have not met,
                {
                        while (array[i] <= pivot && i <= end && k > i)  // from the left, look for the first
                                i++;                                    // element greater than the pivot
                        while (array[k] > pivot && k >= start && k >= i) // from the right, look for the first
                            k--;                                        // element not greater than the pivot
                        if (k > i)                                       // if the left seekindex is still smaller than
                                swap(array, i, k, array_reverenz);                      // the right index, swap the corresponding elements
                }
                swap(array, start, k, array_reverenz);          // after the indices have crossed, swap the last element in
                                                // the left partition with the pivot
                quickSort(array, start, k - 1, array_reverenz); // quicksort the left partition
                quickSort(array, k + 1, end, array_reverenz);   // quicksort the right partition
        }
        else    // if there is only one element in the partition, do not do any sorting
        {
	
                return;                     // the array is sorted, so exit
        }
}
