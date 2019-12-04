//
// Alexander Maksimenko, 29.11.19
//

#ifndef COMB_TYPE_H
#define COMB_TYPE_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>  
#include <math.h>

// The maximum of the number of columns in a matrix is not greater than 32 for unsigned long
#define MAX_COLUMNS 32
#define TRow uint32_t
//#define MAX_COLUMNS 64
//#define TRow uint64_t

class EvalCombType{
	// Main area for intermediate computing
	static TRow **matrices;  // The size is [MAX_COLUMNS+1][MAX_ROWS]; 
	// Additional area for a rank - the number of ones in an appropriate row of the main area
	static int **matrices_rangs; // The size is [MAX_COLUMNS+1][MAX_ROWS]; 
	static int max_columns; // For the allocation of memory
	// The main recursive procedure for finding the optimal permutation of rows and columns
	static int permutation_col_row(int columns, int rest_rows, int sorted_cols, TRow *cur_optimum);
public:
	static int init_error;
	// Before using the lexicographic_ordering() we have to create one object of this class for allocating a memory
	// Do not create two or more objects in one program!
	// max_columns -- the maximal number of columns
	// max_rows -- the maximal number of rows
	EvalCombType(int max_cols, int max_rows);
	~EvalCombType();
	// The entries of the array "matrix" are rows; bits are columns
	// rows is the number of rows
	// columns is the number of columns
	// For the given matrix the procedure lexicographic_ordering() looks for the minimal lexicographical ordering
	// and saves it in the optimum
	// For the correct processing, before call of the procedure, there must be created exactly one object of EvalCombType
	static int lexicographic_ordering(TRow *matrix, int rows, int columns, TRow *optimum);

	// Read the matrix from the file "filename"
	static int read_matrix (char *filename, int &rows, int &columns, TRow *matrix);

	// Write the matrix in the file
	static int write_matrix (char *filename, int rows, int columns, TRow *matrix);
};

#endif //COMB_TYPE_H
