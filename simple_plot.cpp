#include "implot/implot.h"
#include "implot/implot_internal.h"
#include <iostream>


int main() {
    double x[5] = {1, 2, 3, 4, 5};
    double y[5] = {2, 4, 6, 8, 10};

    //ImPlot::SetPlotLimits(0, 10, 0, 10);
    ImPlot::BeginPlot("My Scatter Plot", "X", "Y");
    ImPlot::PlotScatter("My Points", x, y, 5);
    ImPlot::EndPlot();

    return 0;
}