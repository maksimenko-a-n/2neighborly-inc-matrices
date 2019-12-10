# 2neighborly-inc-matrices
The program generates incidence matrices of 2-neighborly d-polytopes for the given list of incidence matrices of facets (2-neighborly (d-1)-polytopes).
The theory described in https://arxiv.org/abs/1912.03900
Input -- 4 parameters: the name of file with facets, the number of vertices of a new polytope, the number of facets for it, and the mimimal number of facets incident to a vertex (the absolute minimum is d).
The fifth parameter is optional -- it is an index of a source facet for constructing a polytope. If the 5th parameter is absent, then the last facet from the given list is taken.

Examples of tested commands:
incmat.exe 4d2n.txt 8 12 8
incmat.exe 5d2n.txt 10 14 9 -l
incmat.exe 5d2n.txt 10 15 9 3
incmat.exe 6d2n.txt 14 16 11 7 -l

For compilation use the "CMakeLists.txt". For Windows the makefile can be done with CMake. For Linux: cmake CMakeLists.txt
After that you can "make" the executable file in Linux or, for example, "mingw32-make" (if the compiler is MinGW) for Windows.

The file "2neighborlyPolytopes.txt" consists main information (V-representations and incidence matrices) about all 2-neighborly d-polytopes with at most d+9 facets.

Files [k]d2n.txt consist incidence matrices for 2-neighborly k-polytopes

