#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "common.h"

// Evaluate the function with the given k at point (x,y)
double evaluateFunction(int k, double x, double y);

// Evaluate the approximating function at point (x,y) using solution coefficients
double evaluateApproximation(const double* solution, double x, double y, const ApproximationContext& context);

// Get all triangles in the grid for the approximation
void generateTriangles(Triangle* triangles, const ApproximationContext& context);

#endif // FUNCTIONS_H
