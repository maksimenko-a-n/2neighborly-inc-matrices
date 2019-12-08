//
// Alexander Maksimenko, 29.11.19
//

#include <Explorer.h>

void printM(const TRow *matrix, int rows, int cols){
	for (int i = 0; i < rows; i++) {
		TRow row = matrix[i];
		for (int j = 0; j < cols; j++, row >>= 1){
			std::cout << (row & 1u);
		}
		std::cout << std::endl;
	}
}

void evalSignature(uint8_t dimension, uint8_t rows, uint8_t cols, TRow *matrix, uint8_t *signature){
	memset(signature, 0, siglen);
	signature[0] = cols;
	signature[1] = rows;
	for (int i = 0; i < rows; i++)
		signature[__builtin_popcount(matrix[i]) - dimension + 2] += 1;
}

int evalColSignature(uint8_t dimension, uint8_t rows, uint8_t cols, TRow *matrix, uint8_t *signature){
	uint8_t columns[MAX_COLUMNS];
	memset(columns, 0, cols);
	for (int i = 0; i < rows; i++){
		auto row = matrix[i];
		for (int j = 0; row; j++, row >>= 1)
			columns[j] += (row & 1u);
	}
	memset(signature, 0, siglen);
	for (int j = 0; j < cols; j++){
		int i = columns[j] - dimension;
		if (i < 0) return 1;
		signature[i] += 1;
	}
	return 0;
}

CombType::CombType(Polytope &polytope, FILE *logfile){
	log_file = logfile;
	dimension = polytope.getDimension();
	auto rows = polytope.getCountFacets();
	auto cols = polytope.getCountVertex();
	auto matrix = polytope.getData(); 
	evalSignature(dimension, rows, cols, matrix, signature);
	evalColSignature(dimension, rows, cols, matrix, colSignature);
	TRow M[MAX_FACETS];
	memcpy(M, matrix, rows * sizeof(TRow));
	memset(M+rows, 0, (MAX_FACETS - rows) * sizeof(TRow));	
	// Recursively generate all matrices
	genMatrices(M, rows, cols, 0);
}

// Recursively generate all permutations of columns, 
// for every permutation sort the rows and append to matrices
void CombType::genMatrices(TRow *in_matrix, int rows, int cols, int curcol){
	TRow matrix[MAX_FACETS];
	memcpy(matrix, in_matrix, MAX_FACETS * sizeof(TRow));
	if (curcol == cols - 1){
		IncidenceMatrix incM(matrix, rows);
		incM.sort();
		matrices.insert(incM);
		if (matrices.size() % 1000000 == 0){
			fprintf (log_file, " %ld", matrices.size());
			fflush (log_file);
		}
		return;
	}
	for (int c = curcol; c < cols; c++){
		for (int r = 0; r < rows; r++){
			TRow source = matrix[r];
			TRow xored = ((source >> curcol) ^ (source >> c)) & 1u; 
			// To swap two sets, we need again to XOR the xored with source
			matrix[r] ^= ((xored << curcol) | (xored << c)); 
		}
		genMatrices(matrix, rows, cols, curcol+1);
	}
}

Explorer::Explorer(int dim, int nVert, int nFact, int minfv, int src_index, bool is_log) : nVertices(nVert), nFacets(nFact), min_facets_in_vert(minfv){	
	// Open log-file
    std::string log_file_name = std::to_string(dim) + "d" + 
			std::to_string(nVertices) + "v" + std::to_string(nFacets) + "f-from(" +
			std::to_string(src_index) + ").log";
	if (is_log)
		log_file = fopen(log_file_name.c_str(), "w");
	else
		log_file = stdout;
};


// Test if x is a proper subset of a row of the matrix
bool Explorer::is_proper_subset(TRow x, TRow *matrix, int rows){
	for (int i = 0; i < rows; i++){
		if (x != matrix[i] && x == (x & matrix[i]))
				return true;
	}
	return false;
}

// Test if x is a subset of a row of the matrix or is equal to a row
bool Explorer::is_subseteq(TRow x, TRow *matrix, int rows){
	for (int i = 0; i < rows; i++){
		if (x == (x & matrix[i]))
				return true;
	}
	return false;
}

// Suppose rows of the matrix is a vertices and all of them are adjacent
// Test if x can be appended to the matrix and the new matrix will remain 2-neighborly
// edges is the array of edges (intersections for every pair of rows)
bool Explorer::test2n(TRow x, TRow *matrix, int rows, std::vector<TRow> *edges){
	for (int i = 0; i < rows; i++){
		TRow e = x & matrix[i]; // edge
		// Test the edge
		for (int j = 0; j < rows; j++)
			if (i != j && (e == (e & matrix[j]))){
				return false;
			}
	}
	if (edges != NULL){
		// Test the old edges
		int N = edges->size();
		TRow * M = edges->data();
		for (int i = 0; i < N; i++){
			if (M[i] == (x & M[i]))
				return false;
		}
	}	
	return true;
}

// Read facets from the file
int Explorer::read_facets(const char* fname){
	// Read 2-neighborly polytopes with the given parameters
	// They will be used as facets of desirable polytopes
    fprintf(log_file, "Read from %s", fname);
    Polytope::readFromFile(fname, facets_list, nFacets-1, nVertices-1, false);
    fprintf(log_file, " %ld polytopes\n", facets_list.size());
	if (facets_list.size() == 0){
		fprintf(log_file, "There are no polytopes!\n");
		return 1;
	}
	dimension = facets_list[0].getDimension();
	fprintf (log_file, "Source dimension: %d\n", dimension);
	// Preprocessing of combinatorial types
    fprintf(log_file, "Generating permutations for facets:");
	for (auto &facet : facets_list){
		types_list.push_back(CombType(facet,log_file));
		fprintf(log_file, " %ld;", types_list.back().matrices.size());
	}
	fprintf(log_file, " %ld facets.\n", types_list.size());
	dimension++;
	return 0;
}

// Generate feasible rows
void Explorer::find_new_facets(TRow *matrix){
	// Generate new facets (rows)
	newfacets.clear();
	TRow maxf = 1u << nVerticesSrc;
	for (unsigned int newfacet = 0; newfacet < maxf; newfacet++){
		// The number of ones must be not less than the dimension
		if (__builtin_popcount(newfacet) < dimension - addcols)
			continue; 		
		// The new row must be a subset of some of the source rows
		if(is_proper_subset(newfacet, matrix, nFacetsSrc))
			newfacets.push_back(newfacet);
	}
	//std::cout << "Feasible rows: " << newfacets.size() << "\n";
}


// Generate feasible columns
void Explorer::find_new_vertices(TRow *matrix){
	// Find new vertices (columns)
	newvertices.clear();
	TRow maxv = 1 << (nFacets - 1);
	for (unsigned int newcol = 0; newcol < maxv; newcol++){
		int nOnes = __builtin_popcount(newcol);
		// The number of ones must be not less than min_facets_in_vert
		if (nOnes < min_facets_in_vert)	continue; 		
		// The number of zeroes must be not less than 3 (else this is a pyramid)
		if (nFacets - nOnes < 3) continue; 		
		// Test 2-neighborliness
		//if (test2n(newcol, &tmatrix, &edges))
		if (test2n(newcol, matrix, nVerticesSrc))
			newvertices.push_back(newcol);
	}
	//std::cout << "Feasible vertices: " << newvertices.size() << "\n";
}

// Final test for the number of vertices in every facet
bool Explorer::test_vert_in_facets(TRow *matrix, int min_ones){
	// Test the vertices in the source facets
	int r = 0;
	for ( ; r < nFacetsSrc; r++){
		int nOnes = __builtin_popcount(matrix[r]);
		if (nOnes > nVerticesSrc)
			return false;
	}
	// Test added facets
	int size = nFacets - 1;
	//for (int r = nFacetsSrc; r < size; r++){
	for (; r < size; r++){
		int nOnes = __builtin_popcount(matrix[r]);
		if(nOnes < min_ones || nOnes > nVerticesSrc)
			return false;
	}
	return true;
}

// Append current column (vertex) to facets
void Explorer::appendColInFacets(int current){
	TRow bit = 1u;
    for (int i = 0; i < nFacets; i++, bit <<= 1) {
		if (tmatrix[current] & bit)
			curFacets[i].appendColumn(tmatrix[current]);
	}
}

// Remove current column (vertex) from facets
void Explorer::removeColFromFacets(int current){
	TRow bit = 1u;
    for (int i = 0; i < nFacets; i++, bit <<= 1) {
		if (tmatrix[current] & bit)
			curFacets[i].removeLastCol();
	}
}

// Returns 0/1-vector of length nFacets. It characterizes if the i-th row generates a good facet (1) or not (0).
TRow Explorer::good_facets(TRow mask, int current){
	//cnt1++;
	TRow bit = 1u;
    for (int i = 0; i < nFacets; i++, bit <<= 1) {
		//if (((mask & bit) == 0) && (tmatrix[current] & bit)){
		if (tmatrix[current] & bit){
			//cnt2++;
			// These two lines take almost all of the computational time
			curFacet = curFacets[i];
			curFacet.matrix[i] = 0;
			curFacet.removeSubrows();
			// ATTENTION!!! May be used in comparing
			//memset(curFacet.matrix + nRows, 0, (MAX_FACETS - nRows) * sizeof(TRow));
			//curFacet.getFacet(i);
			evalSignature(dimension-1, curFacet.getCountRow(), curFacet.getCountColumn(), curFacet.matrix, signature);
			// Check facets
			for (auto &type : types_list){
				if (!type.isEqualSignature(signature))
					continue;
				//if (evalColSignature(dimension-1, curFacet.getCountRow(), curFacet.getCountColumn(), curFacet.matrix, colSignature))
				//	continue;
				//if(!type.isEqualColSignature(colSignature))
				//	continue;
				//cnt3++;
				if (type == curFacet){
					mask |= bit;
					break;
				}
			}
		}
    }
	return mask;
}

// Add a column (vertex) to the matrix
void Explorer::add_vertex(int current, int index, TRow mask){
	int iterations_index = nFacets - 1 - nFacetsSrc + current - nVerticesSrc;
	print_status(iterations_index);
	if (current < nVertices - 1){
		//if (nFacets <= nFacetsSrc+1 && current < nVertices-1)
		/*
		if (current < nVertices-3){
			std::cout << " " << cnt;
			std::cout.flush();
		}*/
		int Nv = newvertices.size();
		for ( ; index < Nv; index++){
			iterations[iterations_index] = index;
			auto vertex = newvertices[index];
			// A vertex (column) cann't have 1 for the row that is already a facet
			if (vertex & mask)
				continue;
			if (!test2n(vertex, tmatrix.data(), current, &edges))
				continue;
			// Add the vertex to the polytope and to the tmatrix
			poly.setLastColumn(current, vertex);
			tmatrix[current] = vertex;
			// Test the number of ones in added rows
			if (!test_vert_in_facets(poly.getData(), dimension - nVertices + current + 1))
				continue;
			// Add new edges
			int edges_size = edges.size();
			for (int i = 0; i < current; i++)
				edges.push_back(tmatrix[i] & vertex);
			// Before running good_facets() we have to update the array of current facets
			appendColInFacets(current);
			add_vertex(current+1, index+1, good_facets(mask, current));
			// Restore the array of current facets in the previous state
			removeColFromFacets(current);
			edges.resize(edges_size);
		}
		// Clear columns
		poly.removeColumnsFrom(current);
	}	
	else{
		TRow vertex = (1u << nFacets) - 1 - mask;
		// Test the number of ones (facets) in the vertex
		if (__builtin_popcount(vertex) < min_facets_in_vert) return; 		
		if (!test2n(vertex, tmatrix.data(), current, &edges)) return;
		// Add the vertex to the polytope and to the tmatrix
		poly.setColumn(current, vertex);
		tmatrix[current] = vertex;
		TRow *matrix_data = poly.getData();
		// Test sum of 1s in every row. It must be <= nVerticesSrc
		if(!test_vert_in_facets(matrix_data, dimension)) return;
		// Every added row cann't be a subset of others
		for (int r = nFacetsSrc; r < nFacets - 1; r++){
			if(is_subseteq(matrix_data[r], matrix_data, r))
				return;
		}
		cnt++;
		// Before running good_facets() we have to update the array of current facets
		appendColInFacets(current);
		if (good_facets(mask, current) == ((1u << nFacets) - 1)){
		//if(Checker::isFacetsPolytopeInVector(poly, facets_list) ){
			poly.printToStream(file_result);
			file_result.flush();
			// Add the polytope to the array of polytopes
			polytopes.push_back(poly);
			//polytopes.push_back(std::make_shared<Polytope>(newpolytope));
			fprintf (log_file, "  good! ");
		}
		// Restore the array of current facets in the previous state
		removeColFromFacets(current);
	}
}

// Add columns (vertices) to the matrix
void Explorer::add_vertices(){
	// Generate feasible columns (vertices)
	find_new_vertices(tmatrix.data());
	//if (nFacets <= nFacetsSrc+1)
	//	std::cout << "Feasible cols: " << newvertices.size() << std::endl;
	// Add the space for vertices
	for (int i = 0; i < addcols; i++)
		tmatrix.push_back(0);
	// Generate the array of current facets
	curFacets.clear();
	auto matrix = poly.getData();
	for (int f = 0; f < nFacets; f++){
		IncidenceMatrix facet(nFacets);
		auto row = matrix[f];
		for (int i = 0; i < nVerticesSrc; i++, row >>= 1){
			if (row & 1u) // Append only columns with ones
				facet.appendColumn(tmatrix[i]);
		}	
		curFacets.push_back(facet);
	}
	// Edges for testing 2-neighborliness
	edges.clear();
	// The last row in the matrix is already a facet
	// Hence mask = 1u << (nFacets - 1)
	add_vertex(nVerticesSrc, 0, 1u << (nFacets - 1));
}

// Add new rows (facets) to the matrix
void Explorer::add_facets(int currentF, int indexFrom){
	if (currentF < nFacets - 1){
		int Nf = newfacets.size();
		//bool isPrintStatus = (0 < nFacets - currentF + addcols - 4);
		//bool isPrintStatus = (0 < nFacets - currentF + addcols - 6);
		//bool isPrintStatus = (0 < nFacets - currentF + addcols - 8);
		for ( ; indexFrom < Nf; indexFrom++){
			iterations[currentF - nFacetsSrc] = indexFrom;
			/*
			if (isPrintStatus){
				if(currentF == nFacetsSrc)
					std::cout << " (" << indexFrom << " / " << Nf << ")";
				else
					std::cout << " " << cnt;
				std::cout.flush();
			}*/
			// Init the facet
			poly.setRow(currentF, newfacets[indexFrom]);
			//add_facets(currentF+1, indexFrom+1);
			add_facets(currentF+1, indexFrom); // ATTENTION!!! indexFrom+1
		}
	}
	else{
		// Transpose matrix
		tmatrix = poly.getTransposedMatrix();
		// Test the number of facets in a vertex
		for (int v = 0; v < tmatrix.size(); v++){
			int nOnes = __builtin_popcount(tmatrix[v]);
			if (nOnes < min_facets_in_vert)	return;
			// The number of zeroes must be not less than 3 (else this is a pyramid)
			if (nFacets - nOnes < 3) return; 		
		}
		// Recursively add vertices
		add_vertices();
	}
}

// Remove duplicates from the polytopes and write them to the file
void Explorer::writeToFile(const std::string &fname){
	Polytope::removeDuplicates(polytopes);
    Polytope::printToFile(polytopes, fname);
}

// The main process of enumerating matrices with the given submatrix (a facet from the given list)
int Explorer::evaluate(int src_facet_index){
	// Parameters for status printing
	print_time = CLOCKS_PER_SEC * 5;
	max_print_time = CLOCKS_PER_SEC * 1200;
	// Time is running
	start = clock();
	next_time = start + print_time;
	
	// Init the polytope
	poly = facets_list[src_facet_index];
	poly.setDimension(dimension);
	/*
	//ATTENTION!!!
    IncidenceMatrix source_matrix;
	for (int i = 0; i < 15; i++)
		source_matrix.appendRow(specialM[i]);
	poly = std::make_shared<Polytope>(dimension, source_matrix);
	*/
	//std::cout << "Source:" << std::endl;
	//poly.printToStream(std::cout);
	
	nFacetsSrc = poly.getCountFacets();
	nVerticesSrc = poly.getCountVertex();
    addcols = nVertices - nVerticesSrc;

	// Open files for output
    std::string file_name_result = std::to_string(dimension) + "d" + 
			std::to_string(nVertices) + "v" + std::to_string(nFacets) + "f-from-" +
			std::to_string(nVerticesSrc) + "v" + std::to_string(nFacetsSrc) + "f-" + std::to_string(min_facets_in_vert) + "mfv";
    //std::string log_file_name = file_name_result + ".log";
    //log_file = fopen(log_file_name.c_str(), "w");
	//log_file = stdout;
	fprintf (log_file, "V=%d, F=%d, minFinV=%d, Index=%d: vsrc=%d, fsrc=%d\n", nVertices, nFacets, min_facets_in_vert, src_facet_index, nVerticesSrc, nFacetsSrc);
	fflush(log_file);

    std::string file_name_tmp = file_name_result + ".tmp";
    file_result = std::ofstream(file_name_tmp);

    //std::string file_name_time = file_name_result + ".time";
    //Timer timer(file_name_time);
	file_name_result += ".res";

	//auto pmatrix = poly.getMatrix();
	//std::cout << "Source:" << std::endl;
	//printM(poly.getData(), poly.getCountRow(), poly.getCountColumn());
	// Generate the list of new rows
	if(nFacetsSrc < nFacets - 1)
		find_new_facets(poly.getData());

	// Append facets to the polytope
    int addrows = nFacets - nFacetsSrc;
	for (int i = 1; i < addrows; i++)
		poly.appendRow(0);
	// Append the current (source) facet (row)
	poly.appendRow((1u << nVerticesSrc) - 1);
	// The matrix for the second stage (May be not needed)
	//matrix2 = *pmatrix;

	// Run main process
    fprintf(log_file, "Check polytopes:");
	cnt = 0;
	add_facets(nFacetsSrc, 0);
	fprintf(log_file, "\n");
	fprintf (log_file, "Total checked: %ld\n", cnt);
	
	// Close the temporary file
    file_result.close();
	fprintf(log_file, "Found: %ld\n", polytopes.size());
	// Sort and write results
	writeToFile(file_name_result);
	fprintf(log_file, "Good: %ld\n", polytopes.size());
	//remove (file_name_tmp.c_str());
	// Print the elapsed time
	clock_t time = clock() - start;
	float sec = (float)time / CLOCKS_PER_SEC;
    fprintf (log_file, "Time: %4.3f sec\n", sec);
	return 0;
}

inline void Explorer::print_status(int iter_index){
	clock_t t = clock();
    if (t >= next_time){
		//print_cur_stat();
		fprintf (log_file, " [");
		if (nFacets - nFacetsSrc > 1)
			fprintf (log_file, "%d/%ld", iterations[0], newfacets.size());
		int i = 1;
		for ( ; i < nFacets - nFacetsSrc - 1; i++)
			fprintf (log_file, " %d", iterations[i]);
		if (i < iter_index)
			fprintf (log_file, "; %d/%ld", iterations[i], newvertices.size());
		i++;
		for ( ; i < iter_index; i++)
			fprintf (log_file, " %d", iterations[i]);
		fprintf (log_file, "]");
		double elapsed_time = (double)(t - start) / CLOCKS_PER_SEC;
		//double left_time = elapsed_time * (inputlen - sent) / sent;
		//long sec = round(left_time);
		long sec = round(elapsed_time);
		int min = sec / 60;
		int hour = min / 60;
		fprintf (log_file, " (%d:%02d:%02ld)", hour, min%60, sec%60);
        fflush (log_file);
        //fflush (logf);
        //fflush (outf);
		if (print_time < max_print_time){
			print_time *= 2;
			if (print_time > max_print_time)
				print_time = max_print_time;
		}	
		next_time = t + print_time;
    }    
}
