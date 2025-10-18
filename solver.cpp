#include "solver.h"
#include "functions.h"
#include <cstring>
#include <cmath>
#include <iostream>
#include <pthread.h>
#include <algorithm>


// Solve the system Ax = b using the minimum error method with Jacobi preconditioner
int solveSystem(const SparseMatrix& A, const double* b, double* x, const ApproximationContext& context) {
    int n = A.n;
    int m_i = context.m_i;
    double eps = context.eps;

    // Initialize solution to zero
    std::memset(x, 0, n * sizeof(double));

    // Allocate temporary vectors
    double* r = new double[n];     // Residual
    double* z = new double[n];     // Preconditioned residual
    double* p = new double[n];     // Search direction
    double* q = new double[n];     // A*p
    double* M = new double[n];     // Jacobi preconditioner

    // Calculate Jacobi preconditioner
    calculateJacobiPreconditioner(A, M);

    // Initialize residual: r = b - A*x
    sparseMV(A, x, r);
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - r[i];
    }

    // Solve M*z = r
    solvePreconditioner(M, r, z, n);

    // Initialize search direction: p = z
    std::memcpy(p, z, n * sizeof(double));

    // Calculate initial residual norm
    double r_norm_squared = 0.0;
    for (int i = 0; i < n; i++) {
        r_norm_squared += r[i] * r[i];
    }
    double r0_norm = std::sqrt(r_norm_squared);
    double r_norm = r0_norm;

    // Iteration counter
    int it = 0;

    // Iterate until convergence or maximum iterations
    while (it < m_i && r_norm / r0_norm > eps) {
        // q = A*p
        sparseMV(A, p, q);

        // Calculate step size alpha = (r,z) / (p,q)
        double rz = 0.0;
        double pq = 0.0;
        for (int i = 0; i < n; i++) {
            rz += r[i] * z[i];
            pq += p[i] * q[i];
        }
        double alpha = rz / std::max(pq, 1e-14);

        // Update solution: x = x + alpha*p
        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
        }

        // Update residual: r = r - alpha*q
        for (int i = 0; i < n; i++) {
            r[i] -= alpha * q[i];
        }

        // Calculate new residual norm
        r_norm_squared = 0.0;
        for (int i = 0; i < n; i++) {
            r_norm_squared += r[i] * r[i];
        }
        r_norm = std::sqrt(r_norm_squared);

        // Solve M*z = r
        solvePreconditioner(M, r, z, n);

        // Calculate beta = (r_new,z_new) / (r_old,z_old)
        double rz_new = 0.0;
        for (int i = 0; i < n; i++) {
            rz_new += r[i] * z[i];
        }
        double beta = rz_new / std::max(rz, 1e-14);

        // Update search direction: p = z + beta*p
        for (int i = 0; i < n; i++) {
            p[i] = z[i] + beta * p[i];
        }

        // Update rz for next iteration
        rz = rz_new;

        it++;
    }

    // Clean up
    delete[] r;
    delete[] z;
    delete[] p;
    delete[] q;
    delete[] M;

    return it;
}

// Thread function for calculating C1 error (max error at triangle centroids)
void* calculateC1ErrorThread(void* arg) {
    ErrorArg* earg = static_cast<ErrorArg*>(arg);
    const double* solution = earg->solution;
    ApproximationContext context = earg->context;
    int thread_id = earg->thread_id;
    int p = context.p;
    int n_x = context.n_x;
    int n_y = context.n_y;

    // Number of triangles
    int numTriangles = 2 * n_x * n_y;
    int triangles_per_thread = numTriangles / p;
    int start_triangle = thread_id * triangles_per_thread;
    int end_triangle = (thread_id == p - 1) ? numTriangles : (thread_id + 1) * triangles_per_thread;

    // Generate triangles
    Triangle* triangles = new Triangle[numTriangles];
    generateTriangles(triangles, context);

    // Find maximum error in this thread's range
    double max_error = 0.0;
    for (int t = start_triangle; t < end_triangle; t++) {
        Triangle& triangle = triangles[t];
        double x = triangle.centroid.x;
        double y = triangle.centroid.y;

        // Exact function value
        double f_val = evaluateFunction(context.k, x, y);

        // Approximated function value
        double p_f_val = evaluateApproximation(solution, x, y, context);

        // Calculate error
        double error = std::abs(f_val - p_f_val);

        // Update maximum error
        max_error = std::max(max_error, error);
    }

    // Save thread's result
    earg->result = max_error;

    // Clean up
    delete[] triangles;

    return nullptr;
}

// Calculate C1 error rate (max |f(P_l)-P_f(P_l)|) at triangle centroids
double calculateC1Error(const double* solution, const ApproximationContext& context) {
    // Setup thread arguments using global thread manager
    g_thread_manager.setupErrorThreads(solution, context);

    // Execute threads for C1 error calculation
    g_thread_manager.executeErrorThreads(calculateC1ErrorThread);

    // Wait for all threads to complete
    g_thread_manager.waitForCompletion();

    // Find global maximum
    double max_error = 0.0;
    for (int i = 0; i < context.p; i++) {
        max_error = std::max(max_error, g_thread_manager.error_thread_args[i].result);
    }

    return max_error;
}

// Thread function for calculating L1 error (sum of errors at triangle centroids)
void* calculateL1ErrorThread(void* arg) {
    ErrorArg* earg = static_cast<ErrorArg*>(arg);
    const double* solution = earg->solution;
    ApproximationContext context = earg->context;
    int thread_id = earg->thread_id;
    int p = context.p;
    int n_x = context.n_x;
    int n_y = context.n_y;
    double h_x = context.h_x;
    double h_y = context.h_y;

    // Number of triangles
    int numTriangles = 2 * n_x * n_y;
    int triangles_per_thread = numTriangles / p;
    int start_triangle = thread_id * triangles_per_thread;
    int end_triangle = (thread_id == p - 1) ? numTriangles : (thread_id + 1) * triangles_per_thread;

    // Generate triangles
    Triangle* triangles = new Triangle[numTriangles];
    generateTriangles(triangles, context);

    // Calculate sum of errors in this thread's range
    double sum_error = 0.0;
    for (int t = start_triangle; t < end_triangle; t++) {
        Triangle& triangle = triangles[t];
        double x = triangle.centroid.x;
        double y = triangle.centroid.y;

        // Exact function value
        double f_val = evaluateFunction(context.k, x, y);

        // Approximated function value
        double p_f_val = evaluateApproximation(solution, x, y, context);

        // Calculate error
        double error = std::abs(f_val - p_f_val);

        // Add weighted error to sum (area = h_x*h_y/2 for each triangle)
        sum_error += error * h_x * h_y / 2.0;
    }

    // Save thread's result
    earg->result = sum_error;

    // Clean up
    delete[] triangles;

    return nullptr;
}

// Calculate L1 error rate (sum |f(P_l)-P_f(P_l)|*area/2) at triangle centroids
double calculateL1Error(const double* solution, const ApproximationContext& context) {
    // Setup thread arguments using global thread manager
    g_thread_manager.setupErrorThreads(solution, context);

    // Execute threads for L1 error calculation
    g_thread_manager.executeErrorThreads(calculateL1ErrorThread);

    // Wait for all threads to complete
    g_thread_manager.waitForCompletion();

    // Sum errors from all threads
    double total_error = 0.0;
    for (int i = 0; i < context.p; i++) {
        total_error += g_thread_manager.error_thread_args[i].result;
    }

    return total_error;
}

// Thread function for calculating C2 error (max error at grid points)
void* calculateC2ErrorThread(void* arg) {
    ErrorArg* earg = static_cast<ErrorArg*>(arg);
    const double* solution = earg->solution;
    ApproximationContext context = earg->context;
    int thread_id = earg->thread_id;
    int p = context.p;
    int n_x = context.n_x;
    int n_y = context.n_y;

    // Calculate which grid points this thread is responsible for
    int total_points = (n_x + 1) * (n_y + 1);
    int points_per_thread = total_points / p;
    int start_point = thread_id * points_per_thread;
    int end_point = (thread_id == p - 1) ? total_points : (thread_id + 1) * points_per_thread;

    // Find maximum error in this thread's range
    double max_error = 0.0;
    for (int point = start_point; point < end_point; point++) {
        int i = point / (n_y + 1);
        int j = point % (n_y + 1);

        double x = context.a + i * context.h_x;
        double y = context.c + j * context.h_y;

        // Exact function value
        double f_val = evaluateFunction(context.k, x, y);

        // Approximated function value
        double p_f_val = solution[getIndex(i, j, n_y)];

        // Calculate error
        double error = std::abs(f_val - p_f_val);

        // Update maximum error
        max_error = std::max(max_error, error);
    }

    // Save thread's result
    earg->result = max_error;

    return nullptr;
}

// Calculate C2 error rate (max |f(x_i,y_j)-P_f(x_i,y_j)|) at grid points
double calculateC2Error(const double* solution, const ApproximationContext& context) {

    // Setup thread arguments using global thread manager
    g_thread_manager.setupErrorThreads(solution, context);

    // Execute threads for C2 error calculation
    g_thread_manager.executeErrorThreads(calculateC2ErrorThread);

    // Wait for all threads to complete
    g_thread_manager.waitForCompletion();

    // Find global maximum
    double max_error = 0.0;
    for (int i = 0; i < context.p; i++) {
        max_error = std::max(max_error, g_thread_manager.error_thread_args[i].result);
    }

    return max_error;
}

// Thread function for calculating L2 error (sum of errors at grid points)
void* calculateL2ErrorThread(void* arg) {
    ErrorArg* earg = static_cast<ErrorArg*>(arg);
    const double* solution = earg->solution;
    ApproximationContext context = earg->context;
    int thread_id = earg->thread_id;
    int p = context.p;
    int n_x = context.n_x;
    int n_y = context.n_y;
    double h_x = context.h_x;
    double h_y = context.h_y;

    // Calculate which grid points this thread is responsible for
    int total_points = (n_x + 1) * (n_y + 1);
    int points_per_thread = total_points / p;
    int start_point = thread_id * points_per_thread;
    int end_point = (thread_id == p - 1) ? total_points : (thread_id + 1) * points_per_thread;

    // Calculate sum of errors in this thread's range
    double sum_error = 0.0;
    for (int point = start_point; point < end_point; point++) {
        int i = point / (n_y + 1);
        int j = point % (n_y + 1);

        double x = context.a + i * h_x;
        double y = context.c + j * h_y;

        // Exact function value
        double f_val = evaluateFunction(context.k, x, y);

        // Approximated function value
        double p_f_val = solution[getIndex(i, j, n_y)];

        // Calculate error
        double error = std::abs(f_val - p_f_val);

        // Add weighted error to sum
        sum_error += error * h_x * h_y;
    }

    // Save thread's result
    earg->result = sum_error;

    return nullptr;
}

// Calculate L2 error rate (sum |f(x_i,y_j)-P_f(x_i,y_j)|*h_x*h_y) at grid points
double calculateL2Error(const double* solution, const ApproximationContext& context) {
    // Setup thread arguments using global thread manager
    g_thread_manager.setupErrorThreads(solution, context);

    // Execute threads for L2 error calculation
    g_thread_manager.executeErrorThreads(calculateL2ErrorThread);

    // Wait for all threads to complete
    g_thread_manager.waitForCompletion();

    // Sum errors from all threads
    double total_error = 0.0;
    for (int i = 0; i < context.p; i++) {
        total_error += g_thread_manager.error_thread_args[i].result;
    }

    return total_error;
}
