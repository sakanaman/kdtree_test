#ifndef __KDTREE_DIM3__
#define __KDTREE_DIM3__

#include <vector>
#include <array>
#include <functional>

namespace kdtree
{

namespace dim3
{

const int NIL = -1;

using Point = std::array<double, 3>;
using Points = std::vector<std::array<double,3>>;
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

    int make3DTree(const int l,
                   const int r,
                   const int depth);

    
    void find_radius(const int v,
                     const double px, double py, double pz, const double r,
                     int depth,
                     Points& res);

private:

    bool is_inside_radius(double x, double y, double z, double px, double py, double pz, double r)
    {
        double dx = x - px;
        double dy = y - py;
        double dz = z - pz;

        return dx*dx + dy*dy + dz*dz < r*r;
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


int Tree::make3DTree(const int l, const int r,
                      const int depth)
{
    auto criteria_x = [&](int a, int b)
    {
        Point one = _access(a);
        Point other = _access(b);

        return one[0] < other[0];
    };
    auto criteria_y = [&](int a, int b)
    {
        Point one = _access(a);
        Point other = _access(b);

        return one[1] < other[1];
    };

    auto criteria_z = [&](int a, int b)
    {
        Point one = _access(a);
        Point other = _access(b);

        return one[2] < other[2];
    };

    if(r <= l)
    {
        return NIL;
    }

    auto begin_iter = _points.begin();
    int mid = (l + r) / 2;

    // [l,r)の中央値で分割する
    if(depth % 3 == 0)
    {
        std::nth_element(begin_iter + l, begin_iter + mid, begin_iter + r, criteria_x);
    }   
    else if(depth % 3 == 1)
    {
        std::nth_element(begin_iter + l, begin_iter + mid, begin_iter + r, criteria_y);
    }
    else
    {
        std::nth_element(begin_iter + l, begin_iter + mid, begin_iter + r, criteria_z);
    }

    int t = np++;
    _nodes[t].location = mid;
    _nodes[t].left = make3DTree(l, mid, depth+1);
    _nodes[t].right = make3DTree(mid+1, r, depth+1);

    return t;
}

void Tree::find_radius(const int v,
                       const double px, double py, double pz, const double r,
                       int depth,
                       Points& res)
{
    Point coord = _access(_points[_nodes[v].location]);

    if(is_inside_radius(coord[0], coord[1], coord[2], px, py, pz, r))
    {
        res.push_back(coord);
    }


    if(depth % 3 == 0)
    {
        if(_nodes[v].left != NIL && px - r <= coord[0])
        {
            find_radius(_nodes[v].left, px, py, pz, r, depth+1, res);
        }
        if(_nodes[v].right != NIL && coord[0] <= px + r)
        {
            find_radius(_nodes[v].right, px, py, pz, r, depth+1, res);
        }
    }
    else if (depth % 3 == 1)
    {
        if(_nodes[v].left != NIL && py - r <= coord[1])
        {
            find_radius(_nodes[v].left, px, py, pz, r, depth+1, res);
        }
        if(_nodes[v].right != NIL && coord[1] <= py + r)
        {
            find_radius(_nodes[v].right, px, py, pz, r, depth+1, res);
        }
    }
    else
    {
        if(_nodes[v].left != NIL && pz - r <= coord[2])
        {
            find_radius(_nodes[v].left, px, py, pz, r, depth+1, res);
        }
        if(_nodes[v].right != NIL && coord[2] <= pz + r)
        {
            find_radius(_nodes[v].right, px, py, pz, r, depth+1, res);
        }
    }
}                    

} // namespace dim2


} // namespace kdtree


#endif