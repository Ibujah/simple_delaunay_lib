# Delaunay triangulation and tetrahedralization in Rust

This software is an implementation of Delaunay 2D and 3D algorithms, for coding exercice. It is not fully optimized, but should at least work correctly.


## Delaunay 2D

A list of 2D points is ordered along an Hilbert curve.

A first (non flat) triangle is inserted, then next points are inserted one by one.

For each point, a walk inside the Delaunay graph is done until a triangle containing the point is found. The point in inserted within the triangle, and a sequence of edge flips is done until the graph is fully Delaunay.


## Delaunay 3D

A list of 3D points is ordered along an Hilbert curve.

A first (non flat) tetrahedron is inserted, then next points are inserted one by one.

For each point, a walk inside the Delaunay graph is done until a triangle containing the point is found. Then, every neighbor tetrahedron which sphere contains the point is removed. A new set of tetrahedra is then inserted to replace the hole in the graph (Bowyer-Watson algorithm).


## Sources

David F. Watson, *Computing the n-dimensional Delaunay tessellation with application to Voronoi polytopes*, The Computer Journal, vol. 24, no 2, 1981, p. 167-172

Adrian Bowyer, *Computing Dirichlet tessellations, The Computer Journal*, vol. 24, no 2, 1981, p. 162-166

Olivier Devilliers, *Delaunay Triangulation, Theory vs practice*, EuroCG12 [link here](https://inria.hal.science/hal-00850561/PDF/EuroCG12-devillers-slides.pdf) 

SI, Hang, *TetGen, a Delaunay-based quality tetrahedral mesh generator*, ACM Transactions on Mathematical Software, 2015, vol. 41, no 2, p. 11.

Gavrilova, Marina, Ratschek, Helmut, et Rokne, Jon G. *Exact computation of Delaunay and power triangulations*, Reliable Computing, 2000, vol. 6, p. 39-60.

[Hilbert curve](https://en.wikipedia.org/wiki/Hilbert_curve)
