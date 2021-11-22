#include "kdtree1.hpp"
#include "kdtree2.hpp"
#include "kdtree3.hpp"
#include <iostream>

int main()
{
    std::vector<double> points 
    = {1, 3, 5 ,6, 10, 13, 14, 16, 19, 21};


    kdtree::dim1::Tree tree{int(points.size()), points};
    tree.make1DTree(0, points.size());


    std::vector<double> res;
    tree.find_radius(0, 19, 8, res);

    for(auto val : res)
    {
        printf("%f\n", val);
    }

}