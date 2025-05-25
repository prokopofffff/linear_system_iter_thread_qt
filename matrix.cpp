#include "matrix.h"
#include "functions.h"
#include <cstring>
#include <algorithm>
#include <iostream>
#include <vector>
#include <pthread.h>

// Thread argument structure for parallel computation
struct ThreadArg {
    SparseMatrix* matrix;
    ApproximationContext context;
    int thread_id;
    double* b;  // For right-hand side calculation
};

// Build the structure of the sparse matrix (MSR format)
void buildMatrixStructure(SparseMatrix& matrix, const ApproximationContext& context) {
    int n_x = context.n_x;
    int n_y = context.n_y;
    int n = (n_x + 1) * (n_y + 1);  // Total number of nodes/unknowns

    // Allocate memory for the MSR format
    matrix.n = n;
    matrix.ia = new int[n + 2]; // Row pointers (size n+2 for MSR format)

    // First, we need to determine the structure (count non-zeros)
    int* nnz_per_row = new int[n];
    std::memset(nnz_per_row, 0, n * sizeof(int));

    // Count non-zeros per row
    for (int i = 0; i <= n_x; i++) {
        for (int j = 0; j <= n_y; j++) {
            int row = getIndex(i, j, n_y);

            // Count non-zeros in this row (including diagonal)
            // For each node, we consider connections to neighboring nodes through triangles

            // Add diagonal element
            nnz_per_row[row]++;

            // Add off-diagonal elements (neighbors)
            // Left neighbor
            if (i > 0) {
                nnz_per_row[row]++;
            }
            // Right neighbor
            if (i < n_x) {
                nnz_per_row[row]++;
            }
            // Bottom neighbor
            if (j > 0) {
                nnz_per_row[row]++;
            }
            // Top neighbor
            if (j < n_y) {
                nnz_per_row[row]++;
            }
            // Diagonal neighbors
            if (i > 0 && j > 0) {
                nnz_per_row[row]++;
            }
            if (i < n_x && j > 0) {
                nnz_per_row[row]++;
            }
            if (i > 0 && j < n_y) {
                nnz_per_row[row]++;
            }
            if (i < n_x && j < n_y) {
                nnz_per_row[row]++;
            }
        }
    }

    // Set up row pointers
    matrix.ia[0] = 0;
    matrix.ia[1] = n; // First n entries are diagonal elements

    for (int i = 0; i < n; i++) {
        matrix.ia[i + 2] = matrix.ia[i + 1] + nnz_per_row[i] - 1; // -1 because diagonal is stored separately
    }

    // Calculate total number of non-zero elements
    matrix.nnz = matrix.ia[n + 1];

    // Allocate memory for column indices and values
    matrix.ja = new int[matrix.nnz];
    matrix.a = new double[matrix.nnz + n]; // +n for diagonal elements

    // Initialize diagonal elements (first n entries in a)
    for (int i = 0; i < n; i++) {
        matrix.a[i] = 0.0;
    }

    // Set up column indices (structure only, values will be filled later)
    int* pos = new int[n];
    for (int i = 0; i < n; i++) {
        pos[i] = matrix.ia[i + 1];
    }

    for (int i = 0; i <= n_x; i++) {
        for (int j = 0; j <= n_y; j++) {
            int row = getIndex(i, j, n_y);

            // Add off-diagonal elements (neighbors)
            // Left neighbor
            if (i > 0) {
                int col = getIndex(i-1, j, n_y);
                matrix.ja[pos[row]++] = col;
            }
            // Right neighbor
            if (i < n_x) {
                int col = getIndex(i+1, j, n_y);
                matrix.ja[pos[row]++] = col;
            }
            // Bottom neighbor
            if (j > 0) {
                int col = getIndex(i, j-1, n_y);
                matrix.ja[pos[row]++] = col;
            }
            // Top neighbor
            if (j < n_y) {
                int col = getIndex(i, j+1, n_y);
                matrix.ja[pos[row]++] = col;
            }
            // Diagonal neighbors
            if (i > 0 && j > 0) {
                int col = getIndex(i-1, j-1, n_y);
                matrix.ja[pos[row]++] = col;
            }
            if (i < n_x && j > 0) {
                int col = getIndex(i+1, j-1, n_y);
                matrix.ja[pos[row]++] = col;
            }
            if (i > 0 && j < n_y) {
                int col = getIndex(i-1, j+1, n_y);
                matrix.ja[pos[row]++] = col;
            }
            if (i < n_x && j < n_y) {
                int col = getIndex(i+1, j+1, n_y);
                matrix.ja[pos[row]++] = col;
            }
        }
    }

    // Sort column indices for each row
    for (int i = 0; i < n; i++) {
        std::sort(&matrix.ja[matrix.ia[i + 1]], &matrix.ja[matrix.ia[i + 2]]);
    }

    // Initialize all values to zero
    std::memset(matrix.a + n, 0, matrix.nnz * sizeof(double));

    // Clean up
    delete[] nnz_per_row;
    delete[] pos;
}

// Parallel function to calculate part of the Gram matrix
void* calculateGramMatrixThread(void* arg) {
    ThreadArg* targ = static_cast<ThreadArg*>(arg);
    SparseMatrix* matrix = targ->matrix;
    ApproximationContext context = targ->context;
    int thread_id = targ->thread_id;
    int p = context.p;
    int n_x = context.n_x;
    int n_y = context.n_y;
    double h_x = context.h_x;
    double h_y = context.h_y;
    int n = matrix->n;

    // Calculate which rows this thread is responsible for
    int rows_per_thread = n / p;
    int start_row = thread_id * rows_per_thread;
    int end_row = (thread_id == p - 1) ? n : (thread_id + 1) * rows_per_thread;

    // Generate triangles
    int numTriangles = 2 * n_x * n_y;
    Triangle* triangles = new Triangle[numTriangles];
    generateTriangles(triangles, context);

    // For each triangle
    for (int t = 0; t < numTriangles; t++) {
        int v1 = triangles[t].vertices[0];
        int v2 = triangles[t].vertices[1];
        int v3 = triangles[t].vertices[2];

        // Check if any vertex is in our range
        if ((v1 >= start_row && v1 < end_row) ||
            (v2 >= start_row && v2 < end_row) ||
            (v3 >= start_row && v3 < end_row)) {

            // Calculate the area of the triangle
            int i1 = v1 / (n_y + 1);
            int j1 = v1 % (n_y + 1);
            int i2 = v2 / (n_y + 1);
            int j2 = v2 % (n_y + 1);
            int i3 = v3 / (n_y + 1);
            int j3 = v3 % (n_y + 1);

            double x1 = context.a + i1 * h_x;
            double y1 = context.c + j1 * h_y;
            double x2 = context.a + i2 * h_x;
            double y2 = context.c + j2 * h_y;
            double x3 = context.a + i3 * h_x;
            double y3 = context.c + j3 * h_y;

            double area = 0.5 * std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

            // Update Gram matrix elements
            // For simplicity, we use a constant value for each triangle
            // In a real implementation, you would integrate over the triangle

            if (v1 >= start_row && v1 < end_row) {
                // Diagonal element
                matrix->a[v1] += area / 6.0;

                // Off-diagonal elements
                int pos_v2 = -1, pos_v3 = -1;

                // Find position of v2 in row v1
                for (int pos = matrix->ia[v1 + 1]; pos < matrix->ia[v1 + 2]; pos++) {
                    if (matrix->ja[pos] == v2) {
                        pos_v2 = pos;
                        break;
                    }
                }

                // Find position of v3 in row v1
                for (int pos = matrix->ia[v1 + 1]; pos < matrix->ia[v1 + 2]; pos++) {
                    if (matrix->ja[pos] == v3) {
                        pos_v3 = pos;
                        break;
                    }
                }

                if (pos_v2 >= 0) {
                    matrix->a[n + pos_v2] += area / 12.0;
                }
                if (pos_v3 >= 0) {
                    matrix->a[n + pos_v3] += area / 12.0;
                }
            }

            if (v2 >= start_row && v2 < end_row) {
                // Diagonal element
                matrix->a[v2] += area / 6.0;

                // Off-diagonal elements
                int pos_v1 = -1, pos_v3 = -1;

                // Find position of v1 in row v2
                for (int pos = matrix->ia[v2 + 1]; pos < matrix->ia[v2 + 2]; pos++) {
                    if (matrix->ja[pos] == v1) {
                        pos_v1 = pos;
                        break;
                    }
                }

                // Find position of v3 in row v2
                for (int pos = matrix->ia[v2 + 1]; pos < matrix->ia[v2 + 2]; pos++) {
                    if (matrix->ja[pos] == v3) {
                        pos_v3 = pos;
                        break;
                    }
                }

                if (pos_v1 >= 0) {
                    matrix->a[n + pos_v1] += area / 12.0;
                }
                if (pos_v3 >= 0) {
                    matrix->a[n + pos_v3] += area / 12.0;
                }
            }

            if (v3 >= start_row && v3 < end_row) {
                // Diagonal element
                matrix->a[v3] += area / 6.0;

                // Off-diagonal elements
                int pos_v1 = -1, pos_v2 = -1;

                // Find position of v1 in row v3
                for (int pos = matrix->ia[v3 + 1]; pos < matrix->ia[v3 + 2]; pos++) {
                    if (matrix->ja[pos] == v1) {
                        pos_v1 = pos;
                        break;
                    }
                }

                // Find position of v2 in row v3
                for (int pos = matrix->ia[v3 + 1]; pos < matrix->ia[v3 + 2]; pos++) {
                    if (matrix->ja[pos] == v2) {
                        pos_v2 = pos;
                        break;
                    }
                }

                if (pos_v1 >= 0) {
                    matrix->a[n + pos_v1] += area / 12.0;
                }
                if (pos_v2 >= 0) {
                    matrix->a[n + pos_v2] += area / 12.0;
                }
            }
        }
    }

    // Clean up
    delete[] triangles;

    return nullptr;
}

// Calculate the Gram matrix in parallel
void calculateGramMatrix(SparseMatrix& matrix, const ApproximationContext& context) {
    int p = context.p;
    pthread_t* threads = new pthread_t[p];
    ThreadArg* thread_args = new ThreadArg[p];

    for (int i = 0; i < p; i++) {
        thread_args[i].matrix = &matrix;
        thread_args[i].context = context;
        thread_args[i].thread_id = i;

        pthread_create(&threads[i], nullptr, calculateGramMatrixThread, &thread_args[i]);
    }

    // Wait for all threads to complete
    for (int i = 0; i < p; i++) {
        pthread_join(threads[i], nullptr);
    }

    // Clean up
    delete[] threads;
    delete[] thread_args;
}

// Thread function for calculating part of the right-hand side
void* calculateRHSThread(void* arg) {
    ThreadArg* targ = static_cast<ThreadArg*>(arg);
    double* b = targ->b;
    ApproximationContext context = targ->context;
    int thread_id = targ->thread_id;
    int p = context.p;
    int n_x = context.n_x;
    int n_y = context.n_y;
    double h_x = context.h_x;
    double h_y = context.h_y;
    int n = (n_x + 1) * (n_y + 1);

    // Calculate which triangles this thread is responsible for
    int numTriangles = 2 * n_x * n_y;
    int triangles_per_thread = numTriangles / p;
    int start_triangle = thread_id * triangles_per_thread;
    int end_triangle = (thread_id == p - 1) ? numTriangles : (thread_id + 1) * triangles_per_thread;

    // Generate all triangles
    Triangle* triangles = new Triangle[numTriangles];
    generateTriangles(triangles, context);

    // Thread-local copy of the right-hand side
    double* local_b = new double[n];
    std::memset(local_b, 0, n * sizeof(double));

    // For each triangle in this thread's range
    for (int t = start_triangle; t < end_triangle; t++) {
        Triangle& triangle = triangles[t];
        int v1 = triangle.vertices[0];
        int v2 = triangle.vertices[1];
        int v3 = triangle.vertices[2];

        // Get coordinates of vertices
        int i1 = v1 / (n_y + 1);
        int j1 = v1 % (n_y + 1);
        int i2 = v2 / (n_y + 1);
        int j2 = v2 % (n_y + 1);
        int i3 = v3 / (n_y + 1);
        int j3 = v3 % (n_y + 1);

        double x1 = context.a + i1 * h_x;
        double y1 = context.c + j1 * h_y;
        double x2 = context.a + i2 * h_x;
        double y2 = context.c + j2 * h_y;
        double x3 = context.a + i3 * h_x;
        double y3 = context.c + j3 * h_y;

        // Calculate area of the triangle
        double area = 0.5 * std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

        // Get function value at the centroid
        double f_val = evaluateFunction(context.k, triangle.centroid.x, triangle.centroid.y);

        // Distribute function value to vertices
        local_b[v1] += f_val * area / 3.0;
        local_b[v2] += f_val * area / 3.0;
        local_b[v3] += f_val * area / 3.0;
    }

    // Critical section: add local contributions to global b
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&mutex);
    for (int i = 0; i < n; i++) {
        b[i] += local_b[i];
    }
    pthread_mutex_unlock(&mutex);
    pthread_mutex_destroy(&mutex);

    // Clean up
    delete[] triangles;
    delete[] local_b;

    return nullptr;
}

// Calculate the right-hand side of the system in parallel
void calculateRightHandSide(double* b, const ApproximationContext& context) {
    int n_x = context.n_x;
    int n_y = context.n_y;
    int n = (n_x + 1) * (n_y + 1);
    int p = context.p;

    // Initialize right-hand side to zero
    std::memset(b, 0, n * sizeof(double));

    // Create threads
    pthread_t* threads = new pthread_t[p];
    ThreadArg* thread_args = new ThreadArg[p];

    for (int i = 0; i < p; i++) {
        thread_args[i].b = b;
        thread_args[i].context = context;
        thread_args[i].thread_id = i;

        pthread_create(&threads[i], nullptr, calculateRHSThread, &thread_args[i]);
    }

    // Wait for all threads to complete
    for (int i = 0; i < p; i++) {
        pthread_join(threads[i], nullptr);
    }

    // Clean up
    delete[] threads;
    delete[] thread_args;
}

// Free memory allocated for the sparse matrix
void freeMatrix(SparseMatrix& matrix) {
    delete[] matrix.ia;
    delete[] matrix.ja;
    delete[] matrix.a;
}

// Sparse matrix-vector multiplication y = A*x
void sparseMV(const SparseMatrix& A, const double* x, double* y) {
    int n = A.n;

    // First handle diagonal elements
    for (int i = 0; i < n; i++) {
        y[i] = A.a[i] * x[i];
    }

    // Then handle off-diagonal elements
    for (int i = 0; i < n; i++) {
        for (int j = A.ia[i + 1]; j < A.ia[i + 2]; j++) {
            y[i] += A.a[n + j] * x[A.ja[j]];
        }
    }
}

// Calculate the Jacobi preconditioner (diagonal elements of A)
void calculateJacobiPreconditioner(const SparseMatrix& A, double* M) {
    int n = A.n;

    for (int i = 0; i < n; i++) {
        M[i] = A.a[i];  // Diagonal elements are stored separately
    }
}

// Solve the system M*x = b where M is the Jacobi preconditioner (diagonal matrix)
void solvePreconditioner(const double* M, const double* b, double* x, int n) {
    for (int i = 0; i < n; i++) {
        x[i] = (std::abs(M[i]) > 1e-10) ? b[i] / M[i] : b[i];
    }
}
