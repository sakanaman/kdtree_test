#ifndef __KDTREE_DIM1__
#define __KDTREE_DIM1__

#include <vector>
#include <algorithm>
#include <functional>

// 1d-tree

namespace kdtree
{

namespace dim1
{

const int NIL = -1;

class Node
{
public:
    int location;
    int left_info = 0;
    int right_info = 0;
};


class Tree
{
public:
    Tree(int point_num, const std::vector<double>& points);
    int make1DTree(int l, int r);
    void find_radius(const int v, const double p, 
                     const double radius, std::vector<double>& res) const;
private:
    int np = 0;
    std::vector<Node> _nodes;
    std::vector<double> _points;
};

Tree::Tree(int point_num, const std::vector<double>& points)
{
    _nodes.resize(point_num * 3);

    _points = points;

    std::sort(_points.begin(), _points.end());
}

int Tree::make1DTree(int l, int r)
{
    if(r <= l) return NIL;

    int mid = (l + r) / 2;


    int t = np++;

    _nodes[t].location = mid;
    _nodes[t].left_info = make1DTree(l, mid);
    _nodes[t].right_info = make1DTree(mid + 1, r);

    return t;
}

void Tree::find_radius(const int v,
                       const double p,
                       const double radius, 
                       std::vector<double>& res) const
{
    double x = _points[_nodes[v].location];

    if(p - radius < x && x < p + radius)
    {
        res.push_back(x);
    }

    if(_nodes[v].left_info != NIL && p - radius <= x)
    {
        find_radius(_nodes[v].left_info, p, radius, res);
    }

    if(_nodes[v].right_info != NIL && x <= p + radius)
    {
        find_radius(_nodes[v].right_info, p, radius, res);
    }
}
    
} // namespace dim1

    
} // namespace kdtree







#endif
