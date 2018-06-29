/**
 * @file	sort.h
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

#ifndef SORT_H_
#define SORT_H_

	bool check(double* array, int size);


	void swap(double array[], int index1, int index2);
	void quickSort(double array[], int start, int end);


	void quickSort(double array[], int start, int end, int array_reverenz[]);
	void swap(double array[], int index1, int index2, int array_reverenz[]);

#endif /* SORT_H_ */
