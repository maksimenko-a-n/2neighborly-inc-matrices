//
// Alexander Maksimenko, 29.11.19
//

#ifndef INCIDENCE_MATRIX_H
#define INCIDENCE_MATRIX_H

#include <stdint.h>  // uint8_t
#include <vector>
#include <fstream>
#include <iostream>  // std::cout()
#include <string>
#include <cstring> // memcmp(), memset(), memcpy()
#include <assert.h>
//#include <bliss/graph.hh>

#include <CombType.h>

//#define TRow unsigned long

#define MAX_FACETS 20


class IncidenceMatrix {
public:

    //IncidenceMatrix() = default;
    ~IncidenceMatrix() = default;

	IncidenceMatrix(int rows = 0);
	IncidenceMatrix(TRow *m, int rows = MAX_FACETS);
	void init(TRow *m, int rows = MAX_FACETS);

    uint8_t nColumns = 0;
	uint8_t nRows = 0;
	TRow matrix[MAX_FACETS];
    bool isCanon = false;

    TRow *getData(){return matrix;};
    unsigned int getCountColumn();
    unsigned int getCountRow();

    void appendRow( TRow new_row );
    void removeRow(unsigned int index);
    TRow getRow(unsigned int index);
    void setRow(unsigned int index, TRow new_value);

    void appendColumn(TRow new_column );
    void removeColumn(unsigned int index);
    void removeColumnsFrom(unsigned int index);
	void removeLastCol(){removeColumnsFrom(nColumns-1);};
    TRow getColumn(unsigned int index);
    void setColumn(unsigned int index, TRow new_value);
    void setLastColumn(unsigned int index, TRow new_value);
    void updateCountColumn();

    void readMatrix(std::istream &i_stream);
    void printToStream(std::ostream &o_stream);

    void sort();
    void removeSubrows();
    void removeSubcols();
	//void getFacet(int index_row);
    std::vector<TRow> getTransposedMatrix();

	void amax_canon(){
		// optimum will be the lexicographically minimal representation of a matrix
		// (we can permute rows and columns)
		TRow optimum[MAX_FACETS]; 
		// Find the optimum
		EvalCombType::lexicographic_ordering(matrix, nRows, nColumns, optimum);
		// Save the optimum to matrix
		memcpy(matrix, optimum, nRows * sizeof(TRow));
	}
	// You can include <bliss/graph.hh> and use the following function
	// for getting the canonical representation of a matrix
    //void bliss_canon(); // You need also uncomment the description of bliss_canon() in the IncidenceMatrix.cpp
    void canonize(){
		if( isCanon ) return;
		//bliss_canon();
		amax_canon();
		isCanon = true;
	};

    bool operator==(IncidenceMatrix &other);
    bool operator<(IncidenceMatrix &other);
    bool operator>(IncidenceMatrix &other);
};

#endif //INCIDENCE_MATRIX_H
