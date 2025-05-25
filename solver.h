#ifndef SOLVER_H
#define SOLVER_H

#include "common.h"
#include "matrix.h"

// Solve the system Ax = b using the minimum error method with Jacobi preconditioner
int solveSystem(const SparseMatrix& A, const double* b, double* x, const ApproximationContext& context);

// Calculate C1 error rate (max |f(P_l)-P_f(P_l)|) at triangle centroids
double calculateC1Error(const double* solution, const ApproximationContext& context);

// Calculate L1 error rate (sum |f(P_l)-P_f(P_l)|*area/2) at triangle centroids
double calculateL1Error(const double* solution, const ApproximationContext& context);

// Calculate C2 error rate (max |f(x_i,y_j)-P_f(x_i,y_j)|) at grid points
double calculateC2Error(const double* solution, const ApproximationContext& context);

// Calculate L2 error rate (sum |f(x_i,y_j)-P_f(x_i,y_j)|*h_x*h_y) at grid points
double calculateL2Error(const double* solution, const ApproximationContext& context);

#endif // SOLVER_H
