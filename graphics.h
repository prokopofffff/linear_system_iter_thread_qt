#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <QColor>
#include <QImage>
#include <QPointF>
#include <vector>

namespace Graphics {
    QColor palette(double value, double min, double max);
    void drawTriangle(QImage& image, const QPointF& p1, const QPointF& p2, const QPointF& p3, const QColor& color);
    void drawGrid(QImage& image, const std::vector<double>& data, int mx, int my, double min, double max);
}

#endif // GRAPHICS_H
