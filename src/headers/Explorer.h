//
// Alexander Maksimenko, 29.11.19
//

#ifndef EXPLORER_H
#define EXPLORER_H

#include <CombType.h>
#include <IncidenceMatrix.h>
#include <Polytope.h>
#include <set>

#define siglen 8

// For debugging and testing the number of functions calls
//unsigned long long cnt1 = 0, cnt2 = 0, cnt3 = 0, cnt4 = 0, cnt5 = 0;


// Evaluate the signature of an incidence matrix.
// This is needed for fast checking of identity of matrices
void evalSignature(uint8_t dimension, uint8_t rows, uint8_t cols, TRow *matrix, uint8_t *signature);

// Like a signature it evaluates signature of columns.
// The test shows that it has little new information w.r.t signature
int evalColSignature(uint8_t dimension, uint8_t rows, uint8_t cols, TRow *matrix, uint8_t *signature);

// For comparing matrices in the set CombType::matrices
struct CmpMatrices {
	bool operator() (const IncidenceMatrix& left, const IncidenceMatrix& right) const {
		//return memcmp(left.matrix, right.matrix, MAX_FACETS * sizeof(TRow)) < 0;
		return memcmp(left.matrix, right.matrix, right.nRows * sizeof(TRow)) < 0;
	}
};

// Special class for speeding up the test of combinatorial equivalence of facets
class CombType {
	public:
	int dimension; // the dimension of the facet
	uint8_t signature[siglen]; // The signature of the matrix. Evaluated by evalSignature()
	uint8_t colSignature[siglen]; // The signature of the matrix columns. Is not used at the current time
	std::set<IncidenceMatrix, CmpMatrices> matrices; // The set of matrices for all permutations of columns
	// When creating the object the procedure genMatrices() generates a large set of matrices
	CombType(Polytope &polytope);
    inline bool isEqualSignature(uint8_t *sig){return (memcmp(signature, sig, siglen) == 0);};	
    inline bool isEqualColSignature(uint8_t *sig){return (memcmp(colSignature, sig, siglen) == 0);};	
	// Recursively generate all permutations of columns, 
	// for every permutation sort the rows and append to matrices
	void genMatrices(TRow *in_matrix, int rows, int cols, int curcol);
	bool operator==(IncidenceMatrix &compared_matrix){
		//if (isEqualSignature(sig))
		return (matrices.find(compared_matrix) != matrices.end());
		//return false;
	}
};

class Explorer {
	uint8_t dimension; // the dimension of the new polytope
	int nVerticesSrc;  // the num of vertices of the source facet
	int nFacetsSrc;    // the num of facets of the source facet
	int nVertices;     // the num of vertices of the new polytope
	int nFacets;       // the num of facets of the new polytope
	int min_facets_in_vert; // the minimal num of 1's in every column
    int addcols; // addcols = nVertices - nVerticesSrc;
	Polytope poly; // The new polytope
	std::vector<TRow> newfacets; // The set of feasible rows
	std::vector<TRow> newvertices; // The set of feasible columns
	std::vector<TRow> tmatrix; // The transposed incidence matrix of the new polytope
	std::vector<TRow> edges;   // The current edges. Every edge is an intersection of two columns (vertices)
	std::vector<Polytope> polytopes; // The set of all new polytopes
	std::ofstream file_result; // Where to write results
	unsigned long cnt; // For counting the number of final testings
public:
	std::vector<Polytope> facets_list; // The given list of facets (polytopes in dimension - 1)
	std::vector<CombType> types_list;  // Combinatorial types of facets for speeding the identity testing
	std::vector<IncidenceMatrix> curFacets; // The list of current facets (matrices) of the new polytope
	IncidenceMatrix curFacet; // The current tested facet of the new polytope
	uint8_t signature[siglen]; // The signature of the matrix. Evaluated by evalSignature()
	uint8_t colSignature[siglen]; // The signature of the matrix columns. Is not used at the current time

    Explorer(int nVert, int nFact, int minfv) : nVertices(nVert), nFacets(nFact), min_facets_in_vert(minfv){};
	// Read facets from file
	int read_facets(const char* fname);
	// The main process of enumerating matrices with the given submatrix (a facet from the given list)
	int evaluate(int src_facet_index = -1);
	// Generate all feasible rows for the new matrix
	void find_new_facets(TRow *matrix);
	// Generate all feasible columns for the new matrix
	void find_new_vertices(TRow *matrix);
	// Final test for the number of vertices in every facet
	bool test_vert_in_facets(TRow *matrix, int min_ones);
	// Add a column (vertex) to the matrix
	void add_vertex(int current, int index, TRow mask);
	// Add columns (vertices) to the matrix
	void add_vertices();
	// Add new rows (facets) to the matrix
	void add_facets(int current, int index);
	// Append current column (vertex) to facets
	void appendColInFacets(int current);
	// Remove current column (vertex) from facets
	void removeColFromFacets(int current);
	// Returns 0/1-vector of length nFacets. 
	// It characterizes if the i-th row generates a good facet (1) or not (0).
	TRow good_facets(TRow mask, int current);
	// Test if x is a proper subset of a row of the matrix
	bool is_proper_subset(TRow x, TRow *matrix, int rows);
	// Test if x is a subset of a row of the matrix or is equal to a row
	bool is_subseteq(TRow x, TRow *matrix, int rows);
	// Suppose rows of the matrix is a vertices and all of them are adjacent
	// Test if x can be appended to the matrix and the new matrix will remain 2-neighborly
	// edges is the array of edges (intersections for every pair of rows)
	bool test2n(TRow x, TRow *matrix, int rows, std::vector<TRow> *edges = NULL);
	// Remove duplicates from the polytopes and write them to the file
	void writeToFile(const std::string &fname);
	unsigned long src_list_size(){return facets_list.size();};
};

#endif //EXPLORER_H
