#ifndef WINDOW_H
#define WINDOW_H

#include <QMainWindow>
#include <QImage>
#include <QKeyEvent>
#include <vector>
#include "common.h"

class VisualizationWindow : public QMainWindow {
    Q_OBJECT
public:
    VisualizationWindow(QWidget* parent = nullptr);
    void setParameters(double a, double b, double c, double d, int nx, int ny, int mx, int my, int k, double eps, int mi, int p);
    void setFunctionData(const std::vector<double>& data, int mx, int my);
    void setApproximationData(const std::vector<double>& data, int mx, int my);
    void setErrorData(const std::vector<double>& data, int mx, int my);

protected:
    void paintEvent(QPaintEvent* event) override;
    void keyPressEvent(QKeyEvent* event) override;

private:
    void updateVisualization();
    void drawPaletteLegend(QPainter& painter);
    void displayState(QPainter& painter);
    void recomputeData();
    // Parameters
    double a_, b_, c_, d_;
    int nx_, ny_, mx_, my_, k_, mi_, p_;
    double eps_;
    // Data
    std::vector<double> functionData_;
    std::vector<double> approximationData_;
    std::vector<double> errorData_;
    // Visualization
    QImage image_;
    int currentMode_; // 0: function, 1: approximation, 2: error
    int scaleLevel_;
    int perturbation_;
    // For recomputation
    ApproximationContext context_;
    std::vector<double> solution_;
};

#endif // WINDOW_H
