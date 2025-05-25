#include "functions.h"
#include <cmath>
#include <iostream>
#include <algorithm>

// Evaluate the function with the given k at point (x,y)
double evaluateFunction(int k, double x, double y) {
    switch (k) {
        case 0:
            return 1.0;                         // f(x,y) = 1
        case 1:
            return x;                           // f(x,y) = x
        case 2:
            return y;                           // f(x,y) = y
        case 3:
            return x + y;                       // f(x,y) = x + y
        case 4:
            return std::sqrt(x*x + y*y);        // f(x,y) = sqrt(x² + y²)
        case 5:
            return x*x + y*y;                   // f(x,y) = x² + y²
        case 6:
            return std::exp(x*x - y*y);         // f(x,y) = e^(x² - y²)
        case 7:
            return 1.0 / (25.0 * (x*x + y*y) + 1.0);  // f(x,y) = 1/(25(x² + y²) + 1)
        default:
            std::cerr << "Invalid function index k=" << k << std::endl;
            return 0.0;
    }
}

// Generate all triangles in the grid
void generateTriangles(Triangle* triangles, const ApproximationContext& context) {
    int n_x = context.n_x;
    int n_y = context.n_y;
    double a = context.a;
    double c = context.c;
    double h_x = context.h_x;
    double h_y = context.h_y;

    int triangleIdx = 0;

    // Generate two triangles for each grid cell
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
            // Bottom-left, bottom-right, top-left, top-right
            int bl = getIndex(i, j, n_y);
            int br = getIndex(i+1, j, n_y);
            int tl = getIndex(i, j+1, n_y);
            int tr = getIndex(i+1, j+1, n_y);

            // First triangle (bottom-left, bottom-right, top-left)
            triangles[triangleIdx].vertices[0] = bl;
            triangles[triangleIdx].vertices[1] = br;
            triangles[triangleIdx].vertices[2] = tl;

            // Calculate centroid (center of gravity)
            double x1 = a + i * h_x;
            double y1 = c + j * h_y;
            double x2 = a + (i+1) * h_x;
            double y2 = c + j * h_y;
            double x3 = a + i * h_x;
            double y3 = c + (j+1) * h_y;

            triangles[triangleIdx].centroid.x = (x1 + x2 + x3) / 3.0;
            triangles[triangleIdx].centroid.y = (y1 + y2 + y3) / 3.0;

            triangleIdx++;

            // Second triangle (bottom-right, top-right, top-left)
            triangles[triangleIdx].vertices[0] = br;
            triangles[triangleIdx].vertices[1] = tr;
            triangles[triangleIdx].vertices[2] = tl;

            // Calculate centroid (center of gravity)
            x1 = a + (i+1) * h_x;
            y1 = c + j * h_y;
            x2 = a + (i+1) * h_x;
            y2 = c + (j+1) * h_y;
            x3 = a + i * h_x;
            y3 = c + (j+1) * h_y;

            triangles[triangleIdx].centroid.x = (x1 + x2 + x3) / 3.0;
            triangles[triangleIdx].centroid.y = (y1 + y2 + y3) / 3.0;

            triangleIdx++;
        }
    }
}

// Barycentric coordinates for a point (x,y) in a triangle
void barycentricCoords(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3,
                       double& lambda1, double& lambda2, double& lambda3) {
    double det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
    lambda1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det;
    lambda2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det;
    lambda3 = 1.0 - lambda1 - lambda2;
}

// Find triangle containing the point (x,y) and return its index
int findTriangle(double x, double y, const ApproximationContext& context) {
    double a = context.a;
    double c = context.c;
    double h_x = context.h_x;
    double h_y = context.h_y;

    // For efficiency, we can directly calculate which cell contains the point
    int i = std::min(std::max(static_cast<int>((x - a) / h_x), 0), context.n_x - 1);
    int j = std::min(std::max(static_cast<int>((y - c) / h_y), 0), context.n_y - 1);

    // Check which of the two triangles in the cell contains the point
    int firstTriangleIdx = 2 * (i * context.n_y + j);
    int secondTriangleIdx = firstTriangleIdx + 1;

    // Calculate vertices of the first triangle
    double x1 = a + i * h_x;
    double y1 = c + j * h_y;
    double x2 = a + (i+1) * h_x;
    double y2 = c + j * h_y;
    double x3 = a + i * h_x;
    double y3 = c + (j+1) * h_y;

    // Calculate barycentric coordinates
    double lambda1, lambda2, lambda3;
    barycentricCoords(x, y, x1, y1, x2, y2, x3, y3, lambda1, lambda2, lambda3);

    // If all barycentric coordinates are non-negative, the point is in the first triangle
    if (lambda1 >= -1e-10 && lambda2 >= -1e-10 && lambda3 >= -1e-10) {
        return firstTriangleIdx;
    } else {
        // Must be in the second triangle
        return secondTriangleIdx;
    }
}

// Evaluate the approximating function at point (x,y) using solution coefficients
double evaluateApproximation(const double* solution, double x, double y, const ApproximationContext& context) {
    // Get number of triangles
    int numTriangles = 2 * context.n_x * context.n_y;

    // Allocate memory for triangles if needed
    Triangle* triangles = new Triangle[numTriangles];
    generateTriangles(triangles, context);

    // Find which triangle contains the point
    int triangleIdx = findTriangle(x, y, context);

    // Get vertices of the triangle
    int v1 = triangles[triangleIdx].vertices[0];
    int v2 = triangles[triangleIdx].vertices[1];
    int v3 = triangles[triangleIdx].vertices[2];

    // Get coordinates of vertices
    int n_y = context.n_y;
    double a = context.a;
    double c = context.c;
    double h_x = context.h_x;
    double h_y = context.h_y;

    // Convert vertex indices to grid indices
    int i1 = v1 / (n_y + 1);
    int j1 = v1 % (n_y + 1);
    int i2 = v2 / (n_y + 1);
    int j2 = v2 % (n_y + 1);
    int i3 = v3 / (n_y + 1);
    int j3 = v3 % (n_y + 1);

    // Get actual coordinates
    double x1 = a + i1 * h_x;
    double y1 = c + j1 * h_y;
    double x2 = a + i2 * h_x;
    double y2 = c + j2 * h_y;
    double x3 = a + i3 * h_x;
    double y3 = c + j3 * h_y;

    // Calculate barycentric coordinates
    double lambda1, lambda2, lambda3;
    barycentricCoords(x, y, x1, y1, x2, y2, x3, y3, lambda1, lambda2, lambda3);

    // Interpolate function value using barycentric coordinates
    double result = lambda1 * solution[v1] + lambda2 * solution[v2] + lambda3 * solution[v3];

    // Clean up memory
    delete[] triangles;

    return result;
}
