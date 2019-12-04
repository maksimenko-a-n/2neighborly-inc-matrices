//
// Alexander Maksimenko, 29.11.19
//

#ifndef POLYTOPE_H
#define POLYTOPE_H

#include <IncidenceMatrix.h>

class Polytope: public IncidenceMatrix {
    uint8_t dimension = 0;

public:

	Polytope() = default;

    uint8_t getDimension(){return dimension;};
    void setDimension(uint8_t dim){dimension = dim;};
    unsigned int getCountFacets(){return getCountRow();};
    unsigned int getCountVertex(){return getCountColumn();};

    int readFromStream(std::istream &i_stream);
    static int readFromFile(const std::string& file_name, std::vector<Polytope> &result);
    static int readFromFile(const std::string& file_name, std::vector<Polytope> &result, 
						unsigned int max_row, unsigned int max_column, bool count_column_fix);

    static void removeDuplicates(std::vector<Polytope> &polytopes);
    void printToStream(std::ostream& out_stream);
    static void printToFile(std::vector<Polytope> &polyhedrons, const std::string &file_name);

    //void getVertexFigure(unsigned int index_column, Polytope *polytope);
    //void getPolytopeFacet(unsigned int index_row, Polytope *polytope);
};


#endif //POLYTOPE_H
