//
// Alexander Maksimenko, 29.11.19
//

#include <Polytope.h>
#include <IncidenceMatrix.h>
#include <algorithm> // std::sort()


// Read polytopes from file_name and return the array.
// max_facets is the upper bound for the number of facets
// if vertices_fix == true, then the number of vertices is fixed
// else vertices is the upper bound for the number of polytope's vertices
int Polytope::readFromFile(const std::string& file_name, std::vector<Polytope> &result, 
                unsigned int max_facets, unsigned int vertices, bool vertices_fix) {

    std::ifstream file_in(file_name);
    if(!file_in.is_open()) {
        std::cout << "File [" << file_name << "] not found" << std::endl;
        return 1;
    }

    while (!file_in.eof()) {
        Polytope polyhedron;
        if( polyhedron.readFromStream(file_in) == 0) {
            if (polyhedron.getCountFacets() <= max_facets && polyhedron.getCountVertex() <= vertices) {
				//printf (" (%d,%d)", polyhedron.getCountVertex(), polyhedron.getCountFacets());
                if (!vertices_fix || polyhedron.getCountVertex() == vertices)
                    result.push_back(polyhedron);
            }
        } else {
            std::cout << "Error: File: [" << file_name << "] position: " << file_in.tellg() << std::endl;
        }
    }
    file_in.close();
	return 0;
}

int Polytope::readFromFile(const std::string &file_name, std::vector<Polytope> &result) {
    std::ifstream file_in(file_name);
    if(!file_in.is_open()) {
        std::cout << "File [" << file_name << "] not found" << std::endl;
        return 1;
    }

    while (!file_in.eof()) {
        Polytope polyhedron; 
        if( polyhedron.readFromStream(file_in) == 0) {
            result.push_back(polyhedron);
        } else {
            std::cout << "Error: File: [" << file_name << "] position: " << file_in.tellg() << std::endl;
        }
    }
    file_in.close();
    return 0;
}

int Polytope::readFromStream(std::istream &i_stream) {
    if (!i_stream.eof()) {
		dimension = 0;
        std::string buff;
        do{
			std::getline(i_stream, buff);
			//std::cout << buff << std::endl;
        }while(buff.empty() && !i_stream.eof());
        try {
			dimension = std::stoul(buff, nullptr, 10);
        }
        catch (const std::invalid_argument &e) {
            std::cout << "Non-correct FILE" << std::endl
                      << "(" << e.what() << ") line: " << buff << std::endl;
            return 2;
        }
		readMatrix(i_stream);
		//printf ("Read matrix: %d, %d, %d\n", dimension, getCountColumn(), getCountRow());
		if (dimension == 0 || getCountRow() == 0)
			return 3;
        return 0;
    } else {
        std::cout << "[Polytope] End file" << std::endl;
    }
    return 1;
}


void Polytope::printToStream(std::ostream& out_stream) {
	unsigned int d = dimension;
	out_stream << d << " " << getCountVertex() << " " << getCountFacets() << std::endl;
    IncidenceMatrix::printToStream(out_stream);
}


void Polytope::removeDuplicates(std::vector<Polytope> &polytopes){
	std::sort(polytopes.begin(), polytopes.end());
    for(int i = 1; i < polytopes.size(); i++ ){
        if (polytopes[i-1] == polytopes[i]) {
            polytopes.erase(polytopes.begin() + i);
            i--;
        }
    }
}

void Polytope::printToFile(std::vector<Polytope> &polyhedrons, const std::string &file_name) {
    std::ofstream out_file(file_name);
    for (auto &p : polyhedrons) {
        p.printToStream(out_file);
    }
    out_file.close();
}

/*

void Polytope::getVertexFigure(unsigned int index_column, Polytope *polytope) {
    //if( !matrix ) return nullptr;
    std::shared_ptr<IncidenceMatrix> new_matrix = std::make_shared<IncidenceMatrix>();
    unsigned int new_dimension = dimension - 1;
	unsigned int selected_column = matrix->getColumn(index_column);
    for(unsigned int i = 0; selected_column != 0; i++, selected_column >>= 1) {
        if (selected_column & 1u) {
            new_matrix->appendRow(matrix->getRow(i));
        }
    }
    new_matrix->removeColumn(index_column);
	new_matrix->removeSubcols();
	polytope->setDimension(new_dimension);
	polytope->setMatrix(new_matrix);
    //return std::make_shared<Polytope>(new_dimension, new_matrix);
}

void Polytope::getPolytopeFacet(unsigned int index_row, Polytope *polytope) {
    //if( !matrix )  return nullptr;
    std::shared_ptr<IncidenceMatrix> new_matrix = std::make_shared<IncidenceMatrix>();
    unsigned int new_dimension = dimension - 1;
    TRow select_row = matrix->getRow(index_row);
    for(unsigned int i = 0; select_row != 0; select_row >>= 1u, i++){
        if( select_row & 1u ) {
            new_matrix->appendColumn(matrix->getColumn(i));
        }
    }
    new_matrix->removeRow(index_row);
	new_matrix->removeSubrows();
	polytope->setDimension(new_dimension);
	polytope->setMatrix(new_matrix);
    //return std::make_shared<Polytope>(new_dimension, new_matrix);
}

std::shared_ptr<IncidenceMatrix> Polytope::getMatrix() {
    //if( !matrix ) return nullptr;
    return matrix;
}

void Polytope::setMatrix(std::shared_ptr<IncidenceMatrix> new_matrix) {
    //if( !matrix )  return;
    //if( checkIncidenceMatrix( new_matrix, dimension ) )
        matrix = std::move(new_matrix);
    //else
    //    std::cout << "new matrix no check" << std::endl;
}

bool Polytope::checkIncidenceMatrix(const std::shared_ptr<IncidenceMatrix>& incidenceMatrix, unsigned int dimension) {
    //if( !incidenceMatrix )  return false;
    for (unsigned int i = 0; i < incidenceMatrix->getCountRow(); ++i) {
        auto i_row = incidenceMatrix->getRow(i);
        for (unsigned int j = i + 1; j < incidenceMatrix->getCountRow(); j++) {
            auto j_row = incidenceMatrix->getRow(j);
			auto intersection = i_row & j_row;
            //if( i_row != j_row &&  (( i_row & j_row ) == j_row || ( i_row & j_row ) == i_row )) {
            if( intersection == i_row || intersection == j_row) {
                return false;
            }
        }
    }
    return true;
}
*/
