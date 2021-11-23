#include "kdtree2.hpp"
#include <stdio.h>
#include <cmath>
#include <random>

int main()
{
    int n = 500000;
    int small = 2;
    kdtree::dim2::Points _points(n + small);
    int x,y;
    for(int i = 0; i < small; ++i)
    {
        double theta = 2.0 * M_PI / small * i;
        double r = 1.0;
        _points[i] = {r * std::cos(theta), r * std::sin(theta)};
    }
    for(int i = 0; i < n; ++i)
    {
        double theta = 2.0 * M_PI / n * i;
        double r = 7.0;
        _points[i + small] = {r * std::cos(theta), r * std::sin(theta)};
    }


    double px = 0.0;
    double py = 0.0;
    double radius = 2.0;

    auto access = [&](int i)
    {
        return _points[i];
    };


    kdtree::dim2::Tree tree{_points.size(), access};
    tree.make2DTree(0, _points.size(), 0);

    kdtree::dim2::Points res;
    tree.find_radius(0, px, py , radius, 0, res);

    printf("%d\n", res.size());
    for(auto val : res)
    {
        printf("(%f, %f)\n", val[0], val[1]);
    }

    res.clear();
    tree.find_radius(0, px, py , radius, 0, res);
    printf("%d\n", res.size());

    res.clear();
    tree.find_radius(0, px, py , radius, 0, res);
    printf("%d\n", res.size());
}