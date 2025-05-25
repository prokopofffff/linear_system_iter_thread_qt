#ifndef MATRIX_H
#define MATRIX_H

#include "common.h"

// Modified Sparse Row (MSR) format for sparse matrix
struct SparseMatrix {
    int n;              // Matrix size
    int* ia;            // Row pointers (size n+2)
    int* ja;            // Column indices
    double* a;          // Non-zero values
    int nnz;            // Number of non-zero elements
};

// Function to build the structure of the sparse matrix (MSR format)
void buildMatrixStructure(SparseMatrix& matrix, const ApproximationContext& context);

// Function to calculate the Gram matrix of the basis from Courant functions
void calculateGramMatrix(SparseMatrix& matrix, const ApproximationContext& context);

// Function to calculate the right-hand side of the system
void calculateRightHandSide(double* b, const ApproximationContext& context);

// Free memory allocated for the sparse matrix
void freeMatrix(SparseMatrix& matrix);

// Sparse matrix-vector multiplication y = A*x
void sparseMV(const SparseMatrix& A, const double* x, double* y);

// Calculate the Jacobi preconditioner (diagonal elements of A)
void calculateJacobiPreconditioner(const SparseMatrix& A, double* M);

// Solve the system M*x = b where M is the Jacobi preconditioner (diagonal matrix)
void solvePreconditioner(const double* M, const double* b, double* x, int n);

#endif // MATRIX_H
