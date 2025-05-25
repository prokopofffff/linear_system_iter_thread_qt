#ifndef COMMON_H
#define COMMON_H

// Structure to hold problem parameters and context
struct ApproximationContext {
    double a, b, c, d;       // Domain boundaries [a,b]×[c,d]
    int n_x, n_y;            // Number of interpolation points
    double h_x, h_y;         // Grid step sizes
    int k;                   // Function identifier
    double eps;              // Required accuracy
    int m_i;                 // Maximum number of iterations
    int p;                   // Number of computational threads
};

// Structure to represent a point
struct Point {
    double x, y;
};

// Structure to represent a triangle
struct Triangle {
    int vertices[3];         // Indices of vertices in the grid
    Point centroid;          // Center of gravity of the triangle
};

// Get index of a grid point (i,j) in the solution array
inline int getIndex(int i, int j, int n_y) {
    return i * (n_y + 1) + j;
}

#endif // COMMON_H
