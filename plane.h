#ifndef PLANE
#define PLANE

#include <Dense>
#include <Eigenvalues>
#include <Core>
#include <map>
#include <StdVector>
#include <set>

class plane
{
    public:
        plane(){points.resize(0,3);}
        Eigen::Vector3d normal; // oriented contrary to origin
        double distance;
        Eigen::MatrixXd points;
        Eigen::VectorXi points_boundary; // indices of points which belong to boundary
        std::set<std::pair<int,int>> pixels;
        std::set<int> connected;
        std::set<int> obstructing; // list of planes indices obstructing the current plane
        std::set<int> obstructed; // list of planes indices obstructed by the current plane
        void computeNormal();
        void appendPoint(Eigen::Vector3d pt);
        void appendPixel(std::pair<int,int> p);
        std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator seed;
        int index;

    private:
};

#include "plane.inl"

#endif // plane
