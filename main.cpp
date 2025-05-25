#include <QApplication>
#include "window.h"
#include <iostream>
#include <cmath>
#include <vector>
#include "matrix.h"
#include "approximation.h"
#include "solver.h"
#include "functions.h"

void printUsage(const char* programName) {
    std::cerr << "Usage: " << programName << " a b c d nx ny mx my k eps mi p" << std::endl;
    std::cerr << "  a, b, c, d: boundaries of area [a,b]×[c,d] (double)" << std::endl;
    std::cerr << "  nx, ny: number of interpolation points on X and Y axes (int)" << std::endl;
    std::cerr << "  mx, my: number of visualization points on X and Y axes (int)" << std::endl;
    std::cerr << "  k: function to approximate (int, 0-7)" << std::endl;
    std::cerr << "  eps: accuracy of the solution (double)" << std::endl;
    std::cerr << "  mi: maximum number of iterations (int)" << std::endl;
    std::cerr << "  p: number of computational threads (int)" << std::endl;
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);

    if (argc != 13) {
        printUsage(argv[0]);
        return 1;
    }

    // Parse command line arguments
    double a = std::atof(argv[1]);
    double b = std::atof(argv[2]);
    double c = std::atof(argv[3]);
    double d = std::atof(argv[4]);
    int nx = std::atoi(argv[5]);
    int ny = std::atoi(argv[6]);
    int mx = std::atoi(argv[7]);
    int my = std::atoi(argv[8]);
    int k = std::atoi(argv[9]);
    double eps = std::atof(argv[10]);
    int mi = std::atoi(argv[11]);
    int p = std::atoi(argv[12]);

    if (a >= b || c >= d || nx <= 0 || ny <= 0 || mx <= 0 || my <= 0 || k < 0 || k > 7 || eps <= 0 || mi <= 0 || p <= 0) {
        std::cerr << "Invalid argument values!" << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    // Calculate grid step sizes
    double h_x = (b - a) / nx;
    double h_y = (d - c) / ny;

    // Create context for the problem
    ApproximationContext context;
    context.a = a;
    context.b = b;
    context.c = c;
    context.d = d;
    context.n_x = nx;
    context.n_y = ny;
    context.h_x = h_x;
    context.h_y = h_y;
    context.k = k;
    context.eps = eps;
    context.m_i = mi;
    context.p = p;

    // Build the sparse matrix and solve the system
    SparseMatrix matrix;
    buildMatrixStructure(matrix, context);
    calculateGramMatrix(matrix, context);

    double* b_vector = new double[(nx + 1) * (ny + 1)];
    calculateRightHandSide(b_vector, context);

    double* solution = new double[(nx + 1) * (ny + 1)];
    solveSystem(matrix, b_vector, solution, context);

    // Generate visualization data for function, approximation, and error
    std::vector<double> functionData((mx + 1) * (my + 1));
    std::vector<double> approximationData((mx + 1) * (my + 1));
    std::vector<double> errorData((mx + 1) * (my + 1));
    double fmax = 0, amax = 0, emax = 0;
    for (int i = 0; i <= mx; ++i) {
        double x = a + i * (b - a) / mx;
        for (int j = 0; j <= my; ++j) {
            double y = c + j * (d - c) / my;
            int idx = i * (my + 1) + j;
            double fval = evaluateFunction(k, x, y);
            double aval = evaluateApproximation(solution, x, y, context);
            double eval = fval - aval;
            functionData[idx] = fval;
            approximationData[idx] = aval;
            errorData[idx] = eval;
            fmax = std::max(fmax, std::abs(fval));
            amax = std::max(amax, std::abs(aval));
            emax = std::max(emax, std::abs(eval));
        }
    }
    std::cout << "Max |f|: " << fmax << ", Max |approx|: " << amax << ", Max |error|: " << emax << std::endl;

    // Create and show the visualization window
    VisualizationWindow window;
    window.setParameters(a, b, c, d, nx, ny, mx, my, k, eps, mi, p);
    window.setFunctionData(functionData, mx, my);
    window.setApproximationData(approximationData, mx, my);
    window.setErrorData(errorData, mx, my);
    window.show();

    // Cleanup
    delete[] b_vector;
    delete[] solution;
    freeMatrix(matrix);

    return app.exec();
}
