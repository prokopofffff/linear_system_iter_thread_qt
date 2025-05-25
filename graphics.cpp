#include "graphics.h"
#include <QPainter>
#include <QPolygonF>
#include <algorithm>

namespace Graphics {
QColor palette(double value, double min, double max) {
    // Simple grayscale mapping
    double alpha = (max > min) ? (value - min) / (max - min) : 0.5;
    alpha = std::clamp(alpha, 0.0, 1.0);
    int v = static_cast<int>(alpha * 255);
    return QColor(v, v, v);
}

void drawTriangle(QImage& image, const QPointF& p1, const QPointF& p2, const QPointF& p3, const QColor& color) {
    QPainter painter(&image);
    QPolygonF triangle;
    triangle << p1 << p2 << p3;
    painter.setBrush(color);
    painter.setPen(Qt::NoPen);
    painter.drawPolygon(triangle);
}

void drawGrid(QImage& image, const std::vector<double>& data, int mx, int my, double min, double max) {
    int W = image.width();
    int H = image.height();
    // For now, map [0,mx]x[0,my] to [0,W-1]x[0,H-1]
    for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {
            // Grid points in data: (i,j), (i+1,j), (i,j+1), (i+1,j+1)
            int idx_bl = i * (my + 1) + j;
            int idx_br = (i + 1) * (my + 1) + j;
            int idx_tl = i * (my + 1) + (j + 1);
            int idx_tr = (i + 1) * (my + 1) + (j + 1);

            // Pixel coordinates
            QPointF bl((double)i / mx * (W - 1), (double)j / my * (H - 1));
            QPointF br((double)(i + 1) / mx * (W - 1), (double)j / my * (H - 1));
            QPointF tl((double)i / mx * (W - 1), (double)(j + 1) / my * (H - 1));
            QPointF tr((double)(i + 1) / mx * (W - 1), (double)(j + 1) / my * (H - 1));

            // First triangle: (bl, br, tl)
            double v1 = data[idx_bl];
            double v2 = data[idx_br];
            double v3 = data[idx_tl];
            double centroid1 = (v1 + v2 + v3) / 3.0;
            QColor color1 = palette(centroid1, min, max);
            drawTriangle(image, bl, br, tl, color1);

            // Second triangle: (br, tr, tl)
            v1 = data[idx_br];
            v2 = data[idx_tr];
            v3 = data[idx_tl];
            double centroid2 = (v1 + v2 + v3) / 3.0;
            QColor color2 = palette(centroid2, min, max);
            drawTriangle(image, br, tr, tl, color2);
        }
    }
}

} // namespace Graphics
