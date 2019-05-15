#ifndef CORNER
#define CORNER

#include <Dense>
#include <Eigenvalues>
#include <Core>
#include <map>
#include <StdVector>
#include <set>
#include "plane.h"
#include "intersection.h"

class corner
{
    public:
        corner(){planes.resize(3); lines.resize(2);}
        std::vector<plane*> planes;
        std::vector<intersection*> lines;
        bool isPlaneIntersection;
        bool isLineIntersection;
        Eigen::Vector3d pt;
        void computePointFromPlanes();
        void computePointFromLines(int plane_idx);
        void setPlanes(plane* p1, plane* p2, plane* p3);
        void setLines(intersection* l1, intersection* l2);
        std::vector<int> corners_connected;

    private:

};

#include "corner.inl"

#endif // corner
