#include "window.h"
#include "graphics.h"
#include <QPainter>
#include <QDebug>
#include <algorithm>
#include <QString>
#include <cmath>
#include "matrix.h"
#include "approximation.h"
#include "solver.h"
#include "functions.h"

VisualizationWindow::VisualizationWindow(QWidget* parent)
    : QMainWindow(parent), a_(0), b_(1), c_(0), d_(1), nx_(10), ny_(10), mx_(10), my_(10), k_(0), mi_(100), p_(1), eps_(1e-6), currentMode_(0), scaleLevel_(0), perturbation_(0)
{
    setWindowTitle("Sparse Matrix Visualization");
    resize(800, 600);
}

void VisualizationWindow::setParameters(double a, double b, double c, double d, int nx, int ny, int mx, int my, int k, double eps, int mi, int p) {
    a_ = a; b_ = b; c_ = c; d_ = d;
    nx_ = nx; ny_ = ny; mx_ = mx; my_ = my;
    k_ = k; eps_ = eps; mi_ = mi; p_ = p;
    scaleLevel_ = 0;
    perturbation_ = 0;
    recomputeData();
    updateVisualization();
}

void VisualizationWindow::setFunctionData(const std::vector<double>& data, int mx, int my) {
    functionData_ = data;
    mx_ = mx; my_ = my;
    updateVisualization();
}

void VisualizationWindow::setApproximationData(const std::vector<double>& data, int mx, int my) {
    approximationData_ = data;
    mx_ = mx; my_ = my;
    updateVisualization();
}

void VisualizationWindow::setErrorData(const std::vector<double>& data, int mx, int my) {
    errorData_ = data;
    mx_ = mx; my_ = my;
    updateVisualization();
}

void VisualizationWindow::paintEvent(QPaintEvent*) {
    QPainter painter(this);
    painter.fillRect(rect(), Qt::white);
    if (!image_.isNull()) {
        painter.drawImage(0, 0, image_);
    }
    displayState(painter);
    drawPaletteLegend(painter);
}

void VisualizationWindow::keyPressEvent(QKeyEvent* event) {
    bool needRecompute = false;
    switch (event->key()) {
        case Qt::Key_0:
            k_ = (k_ + 1) % 8; // Cycle function index 0-7
            needRecompute = true;
            break;
        case Qt::Key_1:
            currentMode_ = (currentMode_ + 1) % 3; // Cycle mode: 0=function, 1=approx, 2=error
            updateVisualization();
            break;
        case Qt::Key_2:
            scaleLevel_++;
            needRecompute = true;
            break;
        case Qt::Key_3:
            if (scaleLevel_ > 0) scaleLevel_--;
            needRecompute = true;
            break;
        case Qt::Key_4:
            nx_ *= 2; ny_ *= 2;
            needRecompute = true;
            break;
        case Qt::Key_5:
            if (nx_ > 2) nx_ /= 2;
            if (ny_ > 2) ny_ /= 2;
            needRecompute = true;
            break;
        case Qt::Key_6:
            perturbation_++;
            needRecompute = true;
            break;
        case Qt::Key_7:
            if (perturbation_ > 0) perturbation_--;
            needRecompute = true;
            break;
        case Qt::Key_8:
            mx_ *= 2; my_ *= 2;
            needRecompute = true;
            break;
        case Qt::Key_9:
            if (mx_ > 2) mx_ /= 2;
            if (my_ > 2) my_ /= 2;
            needRecompute = true;
            break;
        default:
            QMainWindow::keyPressEvent(event);
    }
    if (needRecompute) {
        recomputeData();
        updateVisualization();
    }
}

void VisualizationWindow::recomputeData() {
    // Set up context
    double h_x = (b_ - a_) / nx_;
    double h_y = (d_ - c_) / ny_;
    context_.a = a_;
    context_.b = b_;
    context_.c = c_;
    context_.d = d_;
    context_.n_x = nx_;
    context_.n_y = ny_;
    context_.h_x = h_x;
    context_.h_y = h_y;
    context_.k = k_;
    context_.eps = eps_;
    context_.m_i = mi_;
    context_.p = p_;

    // Build and solve system
    SparseMatrix matrix;
    buildMatrixStructure(matrix, context_);
    calculateGramMatrix(matrix, context_);
    std::vector<double> b_vector((nx_ + 1) * (ny_ + 1));
    calculateRightHandSide(b_vector.data(), context_);

    // Apply perturbation at center
    int cx = nx_ / 2;
    int cy = ny_ / 2;
    int centerIdx = cx * (ny_ + 1) + cy;
    // Find max|f| in the initial area
    double fmax = 0.0;
    for (int i = 0; i <= nx_; ++i) {
        double x = a_ + i * h_x;
        for (int j = 0; j <= ny_; ++j) {
            double y = c_ + j * h_y;
            double val = evaluateFunction(k_, x, y);
            fmax = std::max(fmax, std::abs(val));
        }
    }
    b_vector[centerIdx] += perturbation_ * 0.1 * fmax;

    solution_.resize((nx_ + 1) * (ny_ + 1));
    solveSystem(matrix, b_vector.data(), solution_.data(), context_);

    // Generate visualization data for function, approximation, and error
    functionData_.resize((mx_ + 1) * (my_ + 1));
    approximationData_.resize((mx_ + 1) * (my_ + 1));
    errorData_.resize((mx_ + 1) * (my_ + 1));
    double a1 = a_, b1 = b_, c1 = c_, d1 = d_;
    double zoom = std::pow(0.5, scaleLevel_);
    double cx_area = (a_ + b_) / 2.0;
    double cy_area = (c_ + d_) / 2.0;
    double w = (b_ - a_) * zoom / 2.0;
    double h = (d_ - c_) * zoom / 2.0;
    a1 = cx_area - w; b1 = cx_area + w; c1 = cy_area - h; d1 = cy_area + h;
    for (int i = 0; i <= mx_; ++i) {
        double x = a1 + i * (b1 - a1) / mx_;
        for (int j = 0; j <= my_; ++j) {
            double y = c1 + j * (d1 - c1) / my_;
            int idx = i * (my_ + 1) + j;
            double fval = evaluateFunction(k_, x, y);
            double aval = evaluateApproximation(solution_.data(), x, y, context_);
            double eval = fval - aval;
            functionData_[idx] = fval;
            approximationData_[idx] = aval;
            errorData_[idx] = eval;
        }
    }
    freeMatrix(matrix);
}

void VisualizationWindow::updateVisualization() {
    // Choose data and compute min/max
    const std::vector<double>* data = nullptr;
    QString modeStr;
    switch (currentMode_) {
        case 0: data = &functionData_; modeStr = "Function"; break;
        case 1: data = &approximationData_; modeStr = "Approximation"; break;
        case 2: data = &errorData_; modeStr = "Error"; break;
    }
    if (!data || data->empty()) return;
    double minV = (*std::min_element(data->begin(), data->end()));
    double maxV = (*std::max_element(data->begin(), data->end()));
    // Prepare image
    image_ = QImage(width(), height(), QImage::Format_RGB32);
    image_.fill(Qt::white);
    Graphics::drawGrid(image_, *data, mx_, my_, minV, maxV);
    update();
    // Print max/min to console
    double maxAbs = std::max(std::abs(minV), std::abs(maxV));
    qDebug() << modeStr << "max|F|:" << maxAbs;
}

void VisualizationWindow::drawPaletteLegend(QPainter& painter) {
    int W = width();
    int H = height();
    int legendH = 20;
    QRect legendRect(10, H - legendH - 10, W - 20, legendH);
    for (int i = 0; i < legendRect.width(); ++i) {
        double alpha = (double)i / (legendRect.width() - 1);
        QColor color = Graphics::palette(alpha, 0.0, 1.0);
        painter.setPen(color);
        painter.drawLine(legendRect.left() + i, legendRect.top(), legendRect.left() + i, legendRect.bottom());
    }
    painter.setPen(Qt::black);
    painter.drawRect(legendRect);
}

void VisualizationWindow::displayState(QPainter& painter) {
    painter.setPen(Qt::red);
    int y = 20;
    painter.drawText(10, y, QString("Function k=%1").arg(k_)); y += 20;
    painter.drawText(10, y, QString("Mode: %1").arg(currentMode_ == 0 ? "Function" : currentMode_ == 1 ? "Approximation" : "Error")); y += 20;
    painter.drawText(10, y, QString("Scale level: %1").arg(scaleLevel_)); y += 20;
    painter.drawText(10, y, QString("nx=%1 ny=%2").arg(nx_).arg(ny_)); y += 20;
    painter.drawText(10, y, QString("mx=%1 my=%2").arg(mx_).arg(my_)); y += 20;
    painter.drawText(10, y, QString("Perturbation: %1").arg(perturbation_));
}
