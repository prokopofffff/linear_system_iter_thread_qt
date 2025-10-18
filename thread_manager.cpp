#include "common.h"
#include "matrix.h"
#include "solver.h"
#include <pthread.h>
#include <iostream>

// Global thread manager instance
ThreadManager g_thread_manager;

// Initialize threads with given number
void ThreadManager::initialize(int p) {
    if (initialized) {
        cleanup();
    }

    num_threads = p;
    threads = new pthread_t[p];
    matrix_thread_args = new ThreadArg[p];
    error_thread_args = new ErrorArg[p];
    initialized = true;

    std::cout << "Thread manager initialized with " << p << " threads" << std::endl;
}

// Cleanup threads and free memory
void ThreadManager::cleanup() {
    if (!initialized) {
        return;
    }

    waitForCompletion();

    delete[] threads;
    delete[] matrix_thread_args;
    delete[] error_thread_args;

    threads = nullptr;
    matrix_thread_args = nullptr;
    error_thread_args = nullptr;
    num_threads = 0;
    initialized = false;

    std::cout << "Thread manager cleaned up" << std::endl;
}

// Setup thread arguments for matrix operations
void ThreadManager::setupMatrixThreads(SparseMatrix* matrix, const ApproximationContext& context, double* b) {
    if (!initialized) {
        std::cerr << "Thread manager not initialized!" << std::endl;
        return;
    }

    for (int i = 0; i < num_threads; i++) {
        matrix_thread_args[i].matrix = matrix;
        matrix_thread_args[i].context = context;
        matrix_thread_args[i].thread_id = i;
        matrix_thread_args[i].b = b;
    }
}

// Setup thread arguments for error calculations
void ThreadManager::setupErrorThreads(const double* solution, const ApproximationContext& context) {
    if (!initialized) {
        std::cerr << "Thread manager not initialized!" << std::endl;
        return;
    }

    for (int i = 0; i < num_threads; i++) {
        error_thread_args[i].solution = solution;
        error_thread_args[i].context = context;
        error_thread_args[i].thread_id = i;
        error_thread_args[i].result = 0.0;
    }
}

// Execute matrix calculation threads (Gram matrix)
void ThreadManager::executeMatrixThreads(void* (*thread_func)(void*)) {
    if (!initialized) {
        std::cerr << "Thread manager not initialized!" << std::endl;
        return;
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_create(&threads[i], nullptr, thread_func, &matrix_thread_args[i]);
    }
}

// Execute error calculation threads
void ThreadManager::executeErrorThreads(void* (*thread_func)(void*)) {
    if (!initialized) {
        std::cerr << "Thread manager not initialized!" << std::endl;
        return;
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_create(&threads[i], nullptr, thread_func, &error_thread_args[i]);
    }
}

// Wait for all threads to complete
void ThreadManager::waitForCompletion() {
    if (!initialized) {
        return;
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], nullptr);
    }
}
