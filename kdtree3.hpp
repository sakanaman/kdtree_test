#ifndef __KDTREE_DIM3__
#define __KDTREE_DIM3__

#include <vector>
#include <array>
#include <functional>
#include <cmath>
#include <queue>

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

    void find_knn(const int k, const Point& p, Points& res);
    
    void find_knn_imp(const int v, const int k, const Point& p, int depth, 
                      std::priority_queue<int, std::vector<int>, std::function<bool(int, int)>>& prio_queue);

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
    int axis = depth % 3;
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


double dist(const Point& v1, const Point& v2)
{
    double dx = v1[0] - v2[0];
    double dy = v1[1] - v2[1];
    double dz = v1[2] - v2[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}


void Tree::find_knn(const int k, const Point& p, Points& res)
{
    auto compare = [&](int a, int b)
    {
        Point one = _access(a);
        Point other = _access(b);
        return dist(one, p) < dist(other, p);
    };
    std::priority_queue<int, std::vector<int>, std::function<bool(int, int)>> prio_queue{compare};

    find_knn_imp(0, k, p, 0, prio_queue);

    res.resize(k);
    for(int i = 0; i < k; ++i)
    {
        res[i] = _access(prio_queue.top());
        prio_queue.pop();
    }
}

void Tree::find_knn_imp(const int v, const int k, const Point& p, int depth, 
                        std::priority_queue<int, std::vector<int>, std::function<bool(int, int)>>& prio_queue 
                       )
{
    int now_index = _points[_nodes[v].location];
    Point now_coord = _access(now_index);
    double now_dist = dist(now_coord, p);

    if(prio_queue.size() < k)
    {
        prio_queue.push(now_index);
        if(_nodes[v].left != NIL)
        {
            find_knn_imp(_nodes[v].left, k, p, depth+1, prio_queue);
        }
        if(_nodes[v].right != NIL)
        {
            find_knn_imp(_nodes[v].right, k, p, depth+1, prio_queue);
        }
    }
    else
    {
        int top_index = prio_queue.top();
        Point top_coord = _access(top_index);
        double top_dist = dist(p, top_coord);
        if(now_dist < top_dist)
        {
            prio_queue.pop();
            prio_queue.push(now_index);
            // renew top_*
            top_index = prio_queue.top();
            top_coord = _access(top_index);
            top_dist = dist(p, top_coord);
        }
        int axis = depth % 3;
        if(_nodes[v].left != NIL && p[axis] - top_dist < now_coord[axis])
        {
            find_knn_imp(_nodes[v].left, k, p, depth+1, prio_queue);
        }
        if(_nodes[v].right != NIL && now_coord[axis] < p[axis] + top_dist)
        {
            find_knn_imp(_nodes[v].right, k, p, depth+1, prio_queue);
        }
    }
}                        



} // namespace dim2


} // namespace kdtree


#endif