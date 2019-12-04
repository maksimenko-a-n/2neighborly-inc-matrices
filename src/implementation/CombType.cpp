//
// Alexander Maksimenko, 29.11.19
//

#include <CombType.h>

// Init the static variables by zeroes
TRow **EvalCombType::matrices = NULL;
int **EvalCombType::matrices_rangs = NULL;
int EvalCombType::max_columns = 0;
int EvalCombType::init_error = 0;


// Before using the lexicographic_ordering() we have to create one object of this class for allocating a memory
// Do not create two or more objects in one program!
// max_columns -- the maximal number of columns
// max_rows -- the maximal number of rows
EvalCombType::EvalCombType(int max_cols, int max_rows)
{
	init_error = 0;
	if (matrices != NULL || matrices_rangs != NULL){
		printf ("Error: You are trying to create a second object of EvalCombType!\n");
		init_error = 1;
		return;
	}
	if (max_cols >= MAX_COLUMNS){
		init_error = 2;
		return;
	}
	max_columns = max_cols;
	matrices = (TRow**)calloc((max_columns+1), sizeof(TRow*));
	matrices_rangs = (int**)calloc((max_columns+1), sizeof(int*));
	if (matrices == NULL || matrices_rangs == NULL){
		init_error = 3;
		return;
	}
	long int i;
	for (i = 0; i <= max_columns; i++)
	{
		matrices[i] = (TRow*)malloc(sizeof(TRow)*max_rows);
		matrices_rangs[i] = (int*)malloc(sizeof(int)*max_rows);
		if (matrices[i] == NULL || matrices_rangs[i] == NULL){
			init_error = 4;
			return;
		}
	}
    //printf ("Allocate %dx%d memory\n", max_columns, max_rows);
};

// Free the memory
EvalCombType::~EvalCombType()
{
	long int i;
	if (matrices != NULL){
		for (i = 0; i <= max_columns; i++)
			free(matrices[i]);
		free (matrices);
	}
	if (matrices_rangs != NULL){
		for (i = 0; i <= max_columns; i++)
			free(matrices_rangs[i]);
		free (matrices_rangs);
	}
};
	

// The main recursive procedure for finding the optimal permutation of rows and columns
// columns is the number of columns
// rest_rows is the number of usorted rows
// sorted_cols is the number of sorted columns
// cur_optimum is an array (of size rest_rows) for the lexicographically minimal representative
// The data for sorting is in recursive_sort_array[sorted_cols]
int EvalCombType::permutation_col_row(int columns, int rest_rows, int sorted_cols, TRow *cur_optimum)
{
	// Important!!!
	// We suppose, that rang > 0 ALWAYS!!!
	// AND rows are already sorted by rangs (number of ones)
    // The procedure is designed in such a way that on each stage the rows with identical rang
	// are lexicographically sorted in the set of sorted columns
	TRow *old_matrix = matrices[sorted_cols];
	int *old_matrix_rangs = matrices_rangs[sorted_cols]; // The rangs of rows
	
	// row0_rang is the number of unsorted 1's in the 0th row
	TRow row0_rang = old_matrix_rangs[0]; 

	TRow sorted_bits; // sorted_bits contains 1's for the first sorted_cols bits
	sorted_bits = ((TRow)1 << sorted_cols) - 1;
	TRow row0_bits = old_matrix[0] & sorted_bits;

	// Look over rows with the same rang and collect all 1's in unsorted_bits
	TRow unsorted_bits = old_matrix[0]; // unsorted_bits are locations of unsorted 1's
	int j;
	for (j = 1; j < rest_rows; j++) {
		if (old_matrix_rangs[j] != row0_rang || ((old_matrix[j]) & sorted_bits) != row0_bits)
			break;
		unsorted_bits |= old_matrix[j];
	}
	// j is the number of rows with the minimal rang (is not used in the following)

	// new_matrix is an array for the next step of sorting
	TRow *new_matrix = matrices[sorted_cols+1];
	int *new_matrix_rangs = matrices_rangs[sorted_cols+1];

	int shift;
	TRow q1 = 1 << sorted_cols; // We will exchange this column with others
	TRow q2;
	// The cycle over all unsorted 1's
	// For all columns with 1's in rows with minimal rang
	for (q2 = q1, shift = 0; shift < columns - sorted_cols; shift++, q2 <<= 1) {
		if (!(unsorted_bits & q2)) // The column without 1
			continue;

		int new_rest_rows = rest_rows; // The number of rest rows for the next stage
		TRow q3 = ~((TRow)0) - (q1|q2); // Two zeroes among 1's
		int cur_rang = 0, beg_block = 0, end_block = 0, new_rang;
		TRow old_row, new_row;
		int isNotOptimum = 0;
		// For all rows
		// Exchange the columns q1 and q2. Recalculate rangs etc.
		for (j = 0; j < rest_rows; j++) {
			old_row = old_matrix[j]; // the current row
		    // Exchange columns q1 and q2
			new_row = (old_row & q3) | ((old_row & q2) >> shift) | ((old_row & q1) << shift);
			
			// Calculate the new rang of the row (rang = the number of unsorted 1's)
			// The new rang may be left the same or decreased by 1
			new_rang = old_matrix_rangs[j] - ((new_row >> sorted_cols) & 1);

			if (new_rang == 0) { 
			    // We have found the row with minimal rang
				// If the matrix has such rows, 
				// then the first found is lexicographically minimal
				int ind = rest_rows - new_rest_rows;
				if (cur_optimum[ind] < new_row){ // The current value is worse than optimum
					//goto End_of_big_cycle; 
					isNotOptimum = 1;
					break;
				}
				if (cur_optimum[ind] > new_row) {
				    // If the current value can be better, then update the optimum
					cur_optimum[ind] = new_row;
					if (new_rest_rows > 1)
						cur_optimum[ind+1] = ~((TRow)0); // No comparisons are required in this branch
				}

				new_rest_rows--;
			}
			else {
				if (new_rang >= cur_rang) {
					// If the rang has increased (w.r.t. previous rows), then we fix the beginning of new block
					if (new_rang > cur_rang) 
					    beg_block = end_block;
					// Add a record to the block of rows with identical rang
					new_matrix[end_block] = new_row;
					new_matrix_rangs[end_block] = new_rang;
					cur_rang = new_rang; 
				}
				else {	// if (new_rang < cur_rang) => We need to insert the record before the last block
					int l;
				    // Shift the last block by 1
					for (l = end_block; l > beg_block; l--) 
						new_matrix[l] = new_matrix[l-1];
					new_matrix_rangs[end_block] = cur_rang;
					// Insert the new record before the beginning of block
					new_matrix[beg_block] = new_row;
					new_matrix_rangs[beg_block] = new_rang;
					beg_block++;
					// cur_rang is not changed
				}
				end_block++;
			}
		}
		if(isNotOptimum)
			continue;
		sorted_bits = (1 << (sorted_cols+1)) - 1;
		// Shift (press) all new_matrix_rangs[0] of 1's to the right
		TRow shifted_bits = ((1 << new_matrix_rangs[0]) - 1) << (sorted_cols+1);
		TRow next_row = ((new_matrix[0]) & sorted_bits) | shifted_bits;
		if (new_rest_rows > 1) { // Continue the search
			if (cur_optimum[rest_rows-new_rest_rows] >= next_row )  // The new round of recursion
			    permutation_col_row(columns, new_rest_rows, sorted_cols + 1, cur_optimum + rest_rows - new_rest_rows);
		}
		else {
			if (new_rest_rows == 1) { // There are no rows for sorting
				if (cur_optimum[rest_rows-1] > next_row)
					cur_optimum[rest_rows-1] = next_row;
			}
		}
		
		//End_of_big_cycle: ;
	}
	return 0;
}

// Для массива matrix (ячейка - строка, биты - столбцы) находим 
// лексикографически минимальную перестановку строк и столбцов
// rows - число строк
// columns - число столбцов
// Для корректной работы перед вызовом этой процедуры
// требуется выделить память под промежуточные вычисления 
// с помощью array_alloc(int max_columns, int max_rows)
// а после работы процедуры нужно освободить память с помощью array_free(int max_columns)
// Результат работы процедуры сохраняется в глобальный массив optimum
// The entries of the array "matrix" are rows; bits are columns
// rows is the number of rows
// columns is the number of columns
// For the given matrix the procedure lexicographic_ordering() looks for the minimal lexicographical ordering
// and saves it in the optimum
// For the correct processing, before call of the procedure, there must be created exactly one object of EvalCombType<TRow>
int EvalCombType::lexicographic_ordering(TRow *matrix, int rows, int columns, TRow *optimum)
{
	int i, j;
	for (j = 0; j <= columns; j++) // Zero the array for sorting the rows by the number of 1's
		matrices_rangs[j][0] = 0;
	
	int rang;
	TRow mask;
	for (i = 0; i < rows; i++) {
		mask = matrix[i];
		// Count the number of 1's in matrix[i]
		for ( rang = 0; mask > 0; mask >>= 1)
			rang += (mask & 1);
		// Add the result to the row "rang" of the table of rangs
		matrices[rang][matrices_rangs[rang][0]] = matrix[i];
		matrices_rangs[rang][0] += 1;
	}

    // Write to optimum the rows with only zeroes
	int num_of_zeros = matrices_rangs[0][0]; // The number of rows containing only 0's
    for (j = 0; j < num_of_zeros; j++)
	    optimum[j] = 0;
	if (num_of_zeros == rows)
		return 0;

	// Collect the rows of the rang table into matrices[0] (sort them by rangs)
	int last_row = 0;
	for (j = 1; j <= columns; j++) {
		for (i = 0; i < matrices_rangs[j][0]; i++)	{
			matrices[0][last_row] = matrices[j][i]; // Save rows
			matrices_rangs[0][last_row] = j; // Save rangs
			last_row++;
		}
	}
		
	// At the beginning, we have no the optimal solution
	memset (optimum + num_of_zeros, ~((unsigned char)0), sizeof(TRow)*(rows - num_of_zeros));
	// Call the recursive procedure
	permutation_col_row (columns, last_row, 0, optimum + num_of_zeros);
	return 0;
}

// Read the matrix from the file "filename"
// Variables rows and columns will get the appropriate values
int EvalCombType::read_matrix (char *filename, int &rows, int &columns, TRow *matrix){
	FILE *incfile = fopen(filename, "r");
	if (incfile == NULL){
		printf ("Input error: Cann't open file %s\n", filename);
		return 1;
	}
	char buffer[4096], *endptr;
	long line = 0, ret;
	do{
		if (fgets ( buffer, 4095, incfile ) == NULL){
			printf ("Input error: Cann't find \"begin\" in file %s\n", filename);
			return 2;
		}	
		line++;
		endptr = strchr(buffer,'b');
		ret = memcmp(endptr, "begin", 5);
		//printf ("cmp = %d; %s", ret, buffer);
	}while(ret != 0);
	
	if (fgets ( buffer, 4095, incfile ) == NULL){
		printf ("Input error: Unexpected end of file %s\n", filename);
		return 3;
	}	
	line++;
	
	int N = strtol (buffer, &endptr, 10);
	int M = strtol (endptr, NULL, 10);
	
	if (N < 2){
		printf ("Input error: The number %d of rows is less than 2\n", N);
		return 4;
	}
	
	if (M < 2 || M > MAX_COLUMNS){
		printf ("Input error: The number of columns = %d, but must be in [2,%d]\n", M, MAX_COLUMNS);
		return 5;
	}
	printf ("rows = %d, columns = %d\n", N, M);

	TRow x, y;
	int i, j, read;
	for (i = 0, x = 1; i < N; i++, x <<= 1){
		if (fgets ( buffer, 4095, incfile ) == NULL){
			printf ("Input error: Unexpected end of file %s\n", filename);
			return 6;
		}	
		line++;
		matrix[i] = 0;
		endptr = buffer;
		for (j = 0, y = 1; j < M; j++, y <<= 1){
			read = strtol (endptr, &endptr, 10);
			if (read != 0 && read != 1){
				printf ("Input error: Wrong format in line %d\n", line);
				return 7;
			}	
			if (read)
				matrix[i] += y;
		}
	}

	fclose (incfile);
	columns = M;
	rows = N;
	return 0;
}

// Write the matrix to the file "filename"
int EvalCombType::write_matrix (char *filename, int rows, int columns, TRow *matrix){
	FILE *incfile = fopen(filename, "w");
	if (incfile == NULL){
		printf ("Output error: Cann't open file %s for writing\n", filename);
		return 1;
	}
	int i, j, x;
	//fprintf (incfile, "Ordered matrix\n");
	fprintf (incfile, "begin\n%d %d\n", rows, columns);
	for (i = 0; i < rows; i++){
		for (j = 0, x = 1 << (columns-1); j < columns; j++, x >>= 1){
			fprintf (incfile, " %d", (matrix[i]&x)/x);
		}
		fprintf (incfile, "\n");
	}
	fprintf (incfile, "end\n");
	fclose (incfile);
	return 0;
}

