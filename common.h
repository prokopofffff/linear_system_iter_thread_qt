#ifndef COMMON_H
#define COMMON_H

#include <pthread.h>

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

struct SparseMatrix;

// Thread argument structure for parallel computation
struct ThreadArg {
    SparseMatrix* matrix;
    ApproximationContext context;
    int thread_id;
    double* b;  // For right-hand side calculation
};

// Structure for parallel error calculation
struct ErrorArg {
    const double* solution;
    ApproximationContext context;
    int thread_id;
    double result;  // Local result for the thread
};

// Global thread manager for persistent thread handling
struct ThreadManager {
    pthread_t* threads;
    ThreadArg* matrix_thread_args;
    ErrorArg* error_thread_args;
    int num_threads;
    bool initialized;

    ThreadManager() : threads(nullptr), matrix_thread_args(nullptr), error_thread_args(nullptr), num_threads(0), initialized(false) {}

    // Initialize threads with given number
    void initialize(int p);

    // Cleanup threads and free memory
    void cleanup();

    // Setup thread arguments for matrix operations
    void setupMatrixThreads(SparseMatrix* matrix, const ApproximationContext& context, double* b = nullptr);

    // Setup thread arguments for error calculations
    void setupErrorThreads(const double* solution, const ApproximationContext& context);

    // Execute matrix calculation threads (Gram matrix)
    void executeMatrixThreads(void* (*thread_func)(void*));

    // Execute error calculation threads
    void executeErrorThreads(void* (*thread_func)(void*));

    // Wait for all threads to complete
    void waitForCompletion();
};

// Global thread manager instance
extern ThreadManager g_thread_manager;

#endif // COMMON_H
