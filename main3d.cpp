#include "kdtree3.hpp"
#include <stdio.h>
#include <cmath>
#include <random>

int main()
{
    int n = 1000000;
    int small = 10;
    kdtree::dim3::Points _points(n + small);
    int x,y;

    std::random_device seed;
    std::mt19937 mt(seed());

    // 一様実数分布
    // [0.0, 1.0)の値の範囲で、等確率に実数を生成する
    std::uniform_real_distribution<> dist(0.0, 1.0);
    auto rand = [&](){return dist(mt);};


    for(int i = 0; i < small; ++i)
    {
        double theta = M_PI * rand();
        double phi = 2 * M_PI * rand();
        double r = 1.0 * rand();
        _points[i] = {r * std::sin(theta)*std::cos(phi), 
                      r * std::sin(theta)*std::sin(phi),
                      r * std::cos(theta)};
    }
    for(int i = 0; i < n; ++i)
    {
        double theta = M_PI * rand();
        double phi = 2 * M_PI * rand();
        double r = 6.0 + rand();
        _points[i + small] = {r * std::sin(theta)*std::cos(phi), 
                              r * std::sin(theta)*std::sin(phi),
                              r * std::cos(theta)};
    }


    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    double radius = 2.0;

    auto access = [&](int i)
    {
        return _points[i];
    };


    kdtree::dim3::Tree tree{_points.size(), access};
    tree.make3DTree(0, _points.size(), 0);

    kdtree::dim3::Points res;
    tree.find_radius(0, px, py, pz , radius, 0, res);

    printf("%d\n", res.size());
    for(auto val : res)
    {
        printf("(%f, %f)\n", val[0], val[1]);
    }

    res.clear();
    tree.find_radius(0, px, py, pz , radius, 0, res);
    printf("%d\n", res.size());

    res.clear();
    tree.find_radius(0, px, py, pz , radius, 0, res);
    printf("%d\n", res.size());
}