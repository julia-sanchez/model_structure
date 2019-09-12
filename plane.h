#ifndef PLANE
#define PLANE

#include <Dense>
#include <Eigenvalues>
#include <Core>
#include <map>
#include <StdVector>
#include <set>

// PCL INCLUDES
#include <pcl/point_types.h>
#include <pcl/common/common.h>
#include <pcl/search/kdtree.h>
#include <pcl/filters/uniform_sampling.h>
#include <pcl/io/pcd_io.h>

class plane
{
    public:
        plane(){points.resize(0,3);}
        Eigen::Vector3d normal; // oriented contrary to origin
        double distance;
        std::vector<Eigen::Vector3d> pts;
        Eigen::MatrixXd points;
        Eigen::VectorXi points_boundary; // indices of points which belong to boundary
        std::vector<std::pair<int,int>> pixels;
        std::set<int> connected;
        std::set<int> obstructing; // list of planes indices obstructing the current plane
        std::set<int> obstructed; // list of planes indices obstructed by the current plane
        void computeNormal();
        void appendPoint(Eigen::Vector3d pt);
        void appendPixel(std::pair<int,int> p);
        std::pair<std::pair<int,int>, Eigen::Vector3d> seed;
        int index;
        void mean();
        pcl::PointCloud<pcl::PointXYZ> cloud;
        std::set<int> intersections_indices;
        Eigen::Vector3d mean_point_;
        std::vector<Eigen::Vector3d> ordered_corners;

    private:
};

#include "plane.inl"

#endif // plane
