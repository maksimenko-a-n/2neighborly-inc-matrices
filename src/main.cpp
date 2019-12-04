#include <CombType.h>
#include <IncidenceMatrix.h>
#include <Polytope.h>
#include <Explorer.h>

#include <chrono>
using namespace std::chrono;


// The program has 5 parameters (the last parameter is optional):
//  1) the input file with incidence matrices of facets (2-neighborly polytopes);
//  2) the number of vertices V of desirable polytope;
//  3) the number of facets F of desirable polytope;
//  4) the minimal number of facets in a vertex;
//  5) the index of the source facet (by default the last facet of the input list will be choosen).
// From the input file, only facets with vertices < V and facets < F will be choosen
// The output is a file with desirable polytopes
int main(int argc, char *argv[])
{
    setlocale(LC_ALL, "rus");

	// Read parameters
    if(argc < 5)
    {
        std::cout << "Usage: " << argv[0] << " [file name] [vertices] [facets] [min f in vertices]" << std::endl;
        return 1;
    }
    unsigned int vertices = 0;
    unsigned int facets = 0;
    unsigned int min_facets_in_vert = 0;
    int src_index = -1; // Index of the source facet for the construction
    try{
        vertices = std::stoul(argv[2]);
        facets = std::stoul(argv[3]);
        min_facets_in_vert = std::stoul(argv[4]);
    }
    catch(...){
        std::cout << "Usage: " << argv[0] << " [file name] [vertices] [facets] [min f in vertices]" << std::endl;
        return 1;
    }
	if (vertices < 3 || vertices > 32){
		printf ("Wrong number of vertices: %d!\n", vertices);
		return 2;
	}
	if (facets < 3 || facets > 32){
		printf ("Wrong number of facets: %d!\n", facets);
		return 2;
	}
	if (min_facets_in_vert < 3 || min_facets_in_vert > 32){
		printf ("Wrong number of min_facets_in_vert: %d!\n", min_facets_in_vert);
		return 2;
	}
    if(argc > 5)
    {
		try{
			src_index = std::stoi(argv[5]);
        }
		catch(...){
			std::cout << "The fifth argument must be an index of the facet-source" << std::endl;
			return 1;
		}
    }
	// Create an object for evaluation
	Explorer exp(vertices, facets, min_facets_in_vert);
	// Read 2-neighborly polytopes from the given file
	// They will be used as facets of desirable polytopes
	if (exp.read_facets(argv[1]))
		return 1;
	// Correct the index of a source facet
	unsigned long len = exp.src_list_size();
	if (src_index < 0)
		src_index += len;
	if (src_index < 0 || src_index >= len)
		src_index = len-1;

    // Allocate the memory for evaluation of combinatorial type (lexicographically minimal permutation of rows and columns)
	EvalCombType evalCT(vertices, facets);
	if(EvalCombType::init_error){
		printf ("Error %d when initializing EvalCombType!\n", EvalCombType::init_error);
		return 2;
	}
		
	// The main evaluation
	steady_clock::time_point timeBegin = steady_clock::now();
	exp.evaluate(src_index); 
    steady_clock::time_point timeEnd = steady_clock::now();


    duration<double> time_span = duration_cast<duration<double>>(timeEnd - timeBegin);
    std::cout << "Time: " << time_span.count() << " seconds\n";
//	printf ("cnt1 = %lld, cnt2 = %lld, cnt3 = %lld, cnt4 = %lld, cnt5 = %lld", cnt1, cnt2, cnt3, cnt4, cnt5);

    return 0;
}
