#ifndef __KDTREE_DIM2__
#define __KDTREE_DIM2__

#include <vector>
#include <array>
#include <functional>

namespace kdtree
{

namespace dim2
{

const int NIL = -1;

using Point = std::array<double, 2>;
using Points = std::vector<std::array<double,2>>;
using iPoints = std::vector<int>;

class Node
{
public:
    int location;
    int left;
    int right;
};


using Nodes = std::vector<Node>;
using Accessor = std::function<Point(int)>;


class Tree
{
public:
    Tree(int points_num, const Accessor& access);

    int make2DTree(const int l,
                   const int r,
                   const int depth);

    
    void find_radius(const int v,
                     const double px, double py, const double r,
                     int depth,
                     Points& res);

private:

    bool is_inside_radius(double x, double y, double px, double py, double r)
    {
        double dx = x - px;
        double dy = y - py;

        return dx*dx + dy*dy < r*r;
    }

    int np = 0;
    Accessor _access;
    iPoints _points;
    Nodes _nodes;
};


Tree::Tree(int points_num, const Accessor& access)
    :_access(access), _nodes(points_num)
{
    _points.resize(points_num);

    for(int i = 0; i < points_num; ++i)
    {
        _points[i] = i;
    }
}


int Tree::make2DTree(const int l, const int r,
                      const int depth)
{
    int axis = depth % 2;
    auto criteria = [&](int a, int b)
    {
        Point one = _access(a);
        Point other = _access(b);

        return one[axis] < other[axis];
    };

    if(r <= l)
    {
        return NIL;
    }

    auto begin_iter = _points.begin();

    int mid = (l + r) / 2;

    // [l,r)の中央値で分割する
    std::nth_element(begin_iter + l, begin_iter + mid, begin_iter + r, criteria);

    int t = np++;
    _nodes[t].location = mid;
    _nodes[t].left = make2DTree(l, mid, depth+1);
    _nodes[t].right = make2DTree(mid+1, r, depth+1);

    return t;
}

void Tree::find_radius(const int v,
                       const double px, double py, const double r,
                       int depth,
                       Points& res)
{
    Point coord = _access(_points[_nodes[v].location]);

    if(is_inside_radius(coord[0], coord[1], px, py, r))
    {
        res.push_back(coord);
    }


    if(depth % 2 == 0)
    {
        if(_nodes[v].left != NIL && px - r <= coord[0])
        {
            find_radius(_nodes[v].left, px, py, r, depth+1, res);
        }
        if(_nodes[v].right != NIL && coord[0] <= px + r)
        {
            find_radius(_nodes[v].right, px, py, r, depth+1, res);
        }
    }
    else
    {
        if(_nodes[v].left != NIL && py - r <= coord[1])
        {
            find_radius(_nodes[v].left, px, py, r, depth+1, res);
        }
        if(_nodes[v].right != NIL && coord[1] <= py + r)
        {
            find_radius(_nodes[v].right, px, py, r, depth+1, res);
        }
    }
}                    

} // namespace dim2


} // namespace kdtree


#endif