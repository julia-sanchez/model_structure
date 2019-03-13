#ifndef CORNER
#define CORNER

#include <Dense>
#include <Eigenvalues>
#include <Core>
#include <map>
#include <StdVector>
#include <set>
#include "plane.h"

class corner
{
    public:
        corner(){planes.resize(3);}
        std::vector<plane*> planes;
        Eigen::Vector3d pt;
        void computePoint();
        void setPlanes(plane* p1, plane* p2, plane* p3);

    private:

};

#include "corner.inl"

#endif // corner
