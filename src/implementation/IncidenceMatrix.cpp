//
// Alexander Maksimenko, 29.11.19
//

#include <IncidenceMatrix.h>
#include <algorithm> // std::sort()


inline int compInt (const void *a, const void *b){
	return (long)(*(TRow *)a) - *(TRow *)b;
}

void IncidenceMatrix::sort() {
	qsort(getData(), nRows, sizeof (TRow), compInt);
}

IncidenceMatrix::IncidenceMatrix(int rows){		
	memset(matrix, 0, MAX_FACETS * sizeof(TRow));
	nRows = rows;
	nColumns = 0;
}

IncidenceMatrix::IncidenceMatrix(TRow *m, int rows){
	init(m, rows);
}

void IncidenceMatrix::init(TRow *m, int rows){
	memcpy(matrix, m, rows * sizeof(TRow));
	memset(matrix + rows, 0, (MAX_FACETS - rows) * sizeof(TRow));
	nRows = rows;
}
/*
void IncidenceMatrix::getFacet(int index_row) {
    matrix[index_row] = 0;
	removeSubrows();
	// ATTENTION!!! May be used in comparing
	//memset(matrix + nRows, 0, (MAX_FACETS - nRows) * sizeof(TRow));
}
*/

void IncidenceMatrix::removeSubrows(){
    sort();
    for (unsigned int i = 0; i < nRows; i++) {
        auto i_row = matrix[i];
        for (unsigned int j = nRows-1; j > i; j--) {
            if( ( i_row & matrix[j] ) == i_row ) {
                matrix[i] = 0;
                break;
            }
        }
    }
	int current = 0;
    for (unsigned int from = 0; from < nRows; from++) {
		auto row = matrix[from];
		if (row){
			matrix[current++] = row;
		}
    }
	nRows = current;
    //isCanon = false;
}

// ATTENTION!!! This procedure was not tested!
void IncidenceMatrix::removeSubcols(){
    TRow transposed_matrix[MAX_FACETS];
	// Transpose matrix
    for (unsigned int i = 0; i < nColumns; ++i)
        transposed_matrix[i] = getColumn(i);
	memcpy(getData(), transposed_matrix, nColumns * sizeof(TRow));
	std::swap(nRows, nColumns);
	removeSubrows();
	// Again transpose matrix
    for (unsigned int i = 0; i < nColumns; ++i)
        transposed_matrix[i] = getColumn(i);
	memcpy(getData(), transposed_matrix, nColumns * sizeof(TRow));
	std::swap(nRows, nColumns);
    //updateCountColumn();
}

void IncidenceMatrix::appendRow(TRow new_row) {
    //isCanon = false;
	matrix[nRows] = new_row;
	nRows++;
}

TRow IncidenceMatrix::getRow(unsigned int index) {
    if(nRows > index)
        return matrix[index];
    return 0;
}

unsigned int IncidenceMatrix::getCountRow() {
    return nRows;
}

unsigned int IncidenceMatrix::getCountColumn() {
    return nColumns;
}

void IncidenceMatrix::appendColumn(TRow new_column ){
    //isCanon = false;
    for(int i = 0; i < nRows; i++, new_column >>= 1) {
        matrix[i] |= ((new_column & 1u) << nColumns);
    }
    nColumns++;
}

void IncidenceMatrix::removeColumn(unsigned int index) {
    //isCanon = false;
    auto right = (1u << index) - 1, left = ~((1u << (index + 1)) - 1);
    for(int i = 0; i < nRows; i++){
		TRow row = matrix[i];
        matrix[i] = ((left & row) >> 1u) + (right & row);
    }
    nColumns--;
}

void IncidenceMatrix::removeColumnsFrom(unsigned int index) {
    //isCanon = false;
    //auto right = (1u << index) - 1, left = ((1u << nColumns) - 1) - ((1u << (index + 1)) - 1);
    auto right = (1u << index) - 1;
    for(int i = 0; i < nRows; i++){
        matrix[i] &= right;
    }
    nColumns = index;
}

void IncidenceMatrix::printToStream(std::ostream &o_stream) {
    if ( nColumns == 0 || nRows == 0){
        o_stream << "IncidenceMatrix empty" << std::endl << std::endl;
        return;
    }
    for(int i = 0; i < nRows; i++){
		TRow row = matrix[i];
        //std::string buf;
        for (int j = 0; j < nColumns; j++){
            //buf.push_back((row & 1u) ? '1' : '0');
			o_stream << ((row & 1u) ? '1' : '0');
            row >>= 1u;
        }
        //std::reverse(buf.begin(), buf.end());
        o_stream << std::endl;
    }
    o_stream << std::endl;
}

TRow IncidenceMatrix::getColumn(unsigned int index) {
    TRow result = 0;
    //TRow bit = static_cast<TRow>(1) << index;
	//TRow c = 1;
    for(int i = 0; i < nRows; i++){ //, c <<= 1){
        //result |= ((matrix[i] & bit) ? c : 0);
        result |= (((matrix[i] >> index) & 1u) << i);
    }
    return result;
}

void IncidenceMatrix::removeRow(unsigned int index) {
    //isCanon = false;
	memcpy(matrix + index, matrix + index + 1, (nRows - index - 1) * sizeof(TRow));
	nRows--;
}

void IncidenceMatrix::setColumn(unsigned int index, TRow new_value){
    //isCanon = false;
    TRow mask = ~(1u << index);
    for(int i = 0; i < nRows; i++, new_value >>= 1) {
        matrix[i] &= mask; // Clear the column
        matrix[i] |= ((new_value & 1u) << index);
    }
	if (nColumns <= index)
		nColumns = index+1;
}

void IncidenceMatrix::setLastColumn(unsigned int index, TRow new_value){
    //isCanon = false;
    //TRow mask = ~(1u << index);
    auto right = (1u << index) - 1;
    for(int i = 0; i < nRows; i++, new_value >>= 1) {
        matrix[i] &= right; // Clear columns
        matrix[i] |= ((new_value & 1u) << index);
    }
    nColumns = index+1;
}


void IncidenceMatrix::setRow(unsigned int index, TRow new_value) {
    //isCanon = false;
    matrix[index] = new_value;
}


void IncidenceMatrix::updateCountColumn() {
	TRow mask = 0;
	for (int r = 0; r < nRows; r++)
		mask |= matrix[r];
	for (nColumns = 0; mask != 0; mask >>= 1u, nColumns++);
}

void IncidenceMatrix::readMatrix(std::istream &i_stream) {
	nRows = 0;
	nColumns = 0;
    if (!i_stream.eof()) {
        do {
            std::string buff;
            std::getline(i_stream, buff);
			//std::cout << buff << std::endl;
            if( buff.size() < 3 )
                break;
            try {
				std::reverse(buff.begin(), buff.end());
                appendRow(std::stoull(buff, nullptr, 2));
            }
            catch (const std::invalid_argument &e) {
                std::cout << "Catch error in IncidenceMatrix::readMatrix()" << std::endl
                          << "(" << e.what() << ") line: " << buff << std::endl;
                break;
            }

        } while (!i_stream.eof());
		updateCountColumn();
    } 
	else
        std::cout << "[IncidenceMatrix] Конец файла" << std::endl;
}

std::vector<TRow> IncidenceMatrix::getTransposedMatrix() {
    std::vector<TRow> result;
	int cols = getCountColumn();
    result.reserve(cols);
    for (unsigned int i = 0; i < cols; ++i) {
        result.push_back(getColumn(i));
    }
    return result;
}

/*
// The straightforward way of canonizing a matrix is very slow 
void IncidenceMatrix::amax_canon(TRow *inmat, int curcol){
    unsigned int nRows = getCountRow();
    unsigned int cols = getCountColumn();

	if (curcol == cols - 1){
		qsort(inmat, nRows, sizeof (TRow), compInt);
		if (memcmp (matrix.data(), inmat, nRows * sizeof(TRow)) > 0)
			memcpy(matrix.data(), inmat, nRows * sizeof(TRow));
		return;
	}
	TRow outmat[MAX_FACETS];
	memcpy(outmat, inmat, nRows * sizeof(TRow));
	for (int b = curcol; b < cols; b++){
		for (int r = 0; r < nRows; r++){ // Swap two bits curcol and b
			TRow x = outmat[r];
			TRow xored = ((x >> curcol) ^ (x >> b)) & 1u; 
			// To swap two sets, we need to again XOR the xored with original
			outmat[r] = x ^ ((xored << curcol) | (xored << b)); 
		}
		amax_canon(outmat, curcol+1);
	}
}
*/
/*
void IncidenceMatrix::bliss_canon() {
	//cnt_bliss++;
    bliss::Digraph g(nColumns + nRows);
    bliss::Digraph::SplittingHeuristic shs_directed = bliss::Digraph::shs_fsm;
    g.set_splitting_heuristic(shs_directed);
    g.set_verbose_level(0);
    g.set_verbose_file(stdout);
    g.set_failure_recording(true);
    g.set_component_recursion(true);

// Copy data from matrix to the graph g
    unsigned int one = 1u << (nColumns - 1u);
    for (int f = 0; f < nRows; f++) {
        unsigned int row = matrix[f];
        for (int v = 0; v < nColumns; v++, row <<= 1u) {
            if (row & one)
                g.add_edge(nColumns + f, v);
        }
    }

    bliss::Stats stats;
// Canonical labeling
    const unsigned int *cl = g.canonical_form(stats, nullptr, stdout);
// Permute to canonical labeling
    bliss::Digraph *cf = g.permute(cl);

    one = 1;
    for (int f = 0; f < nRows; f++) {
        std::vector<unsigned int> *edges_out = &(cf->vertices[nColumns + f].edges_out);
        int len = edges_out->size();
        unsigned int *vs = edges_out->data();
        matrix[f] = 0;
        for (int v = 0; v < len; v++) {
            matrix[f] |= one << vs[v];
        }
    }
    delete cf;
}
*/

bool IncidenceMatrix::operator==(IncidenceMatrix &other) {
    if( nColumns != other.nColumns || nRows != other.nRows)
        return false;
	canonize();
	other.canonize();
	return memcmp(getData(), other.getData(), nRows * sizeof(TRow)) == 0;
}

bool IncidenceMatrix::operator<(IncidenceMatrix &other) {
    if( nColumns != other.nColumns || nRows != other.nRows)
        return ((nColumns != other.nColumns) ? (nColumns < other.nColumns) : (nRows < other.nRows));
	canonize();
	other.canonize();
	return memcmp(getData(), other.getData(), nRows * sizeof(TRow)) < 0;
}

bool IncidenceMatrix::operator>(IncidenceMatrix &other) {
    if( nColumns != other.nColumns || nRows != other.nRows)
        return ((nColumns != other.nColumns) ? (nColumns > other.nColumns) : (nRows > other.nRows));
	canonize();
	other.canonize();
	return memcmp(getData(), other.getData(), nRows * sizeof(TRow)) > 0;
}
