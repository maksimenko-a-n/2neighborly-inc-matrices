# 2neighborly-inc-matrices
The program generates incidence matrices of 2-neighborly d-polytopes for the given list of incidence matrices of facets (2-neighborly (d-1)-polytopes).
Input -- 4 parameters: the name of file with facets, the number of vertices of a new polytope, the number of facets for it, and the mimimal number of facets incident to a vertex (the absolute minimum is d).
The fifth parameter is optional -- it is an index of a source facet for constructing a polytope. If the 5th parameter is absent, then the last facet from the given list is taken.

For compilation use the "CMakeLists.txt". For Windows it can be done with CMake. For Linux: cmake CMakeLists.txt

The file "2neighborlyPolytopes.txt" consists main information (V-representations and incidence matrices) about all 2-neighborly d-polytopes with at most d+9 facets
