#ifndef INTERSECTION
#define INTERSECTION

#include <Dense>
#include <Core>
#include "plane.h"
#include "corner.h"
#include <set>
#include <Geometry>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/common/common.h>
#include <time.h>

const double eps = 0.00001;

class intersection
{
    public:
        intersection()
        {
            delta_phi = 0.01;
            delta_theta = 0.01;
            isConnection = false;
            isObstruction = false;
            isOpening = false;
            isLine = true;
            Nrow = (int)((2*M_PI+2*eps)/delta_phi);
            Ncol = (int)((M_PI+2*eps)/delta_theta);
        }
        intersection(plane* p1, plane* p2, float dp, float dt)
        {
            isConnection = false;
            isObstruction = false;
            isOpening = false;
            isLine = false;
            plane_ref = p1;
            plane_neigh = p2;
            delta_phi = dp;
            delta_theta = dt;
            Nrow = (int)((2*M_PI+2*eps)/delta_phi);
            Ncol = (int)((M_PI+2*eps)/delta_theta);
        }
        void setPlanes(plane* p1, plane* p2){plane_ref = p1; plane_neigh = p2;}
        Eigen::Vector3d tangente;
        Eigen::Vector3d normal;
        Eigen::Vector2d tangente2D;
        Eigen::Vector2d normal2D;
        double distance;
        double distance2D;
        plane* plane_ref;
        plane* plane_neigh;
        std::vector<std::pair<float,float>> theoritical_phi_theta;
        std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> points;
        std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> other_points;
        std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> points2D;
        std::set<std::pair<int,int>> theoritical_pixels;
        std::vector<std::pair<int,int>> pixels;
        std::vector<std::pair<int,int>> other_pixels;
        std::map<std::pair<int,int>, int> pixel2idx_point;
        std::map<int, std::pair<int,int>> idx_point2pixel;
        float start;
        float end;
        std::set<float> lim;
        void computeTheoriticalLineFeatures();
        void computeTheoriticalPhiTheta();
        void computeLim(std::vector<corner>& corners, std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> possible_inter3D);
        void setDeltas(float dp, float dt);
        void definePlaneConnection();
        bool isConnection;
        bool isObstruction;
        bool isOpening;
        bool isLine;
        Eigen::Vector3d pt;
        Eigen::Affine3d rot;
        Eigen::Affine3d trans_z;
        Eigen::Affine3d trans_xy;
        Eigen::Affine3d rot_inv;
        Eigen::Affine3d trans_z_inv;
        Eigen::Affine3d trans_xy_inv;
        Eigen::Vector2d pt_mean;
        std::vector<int> theoriticalLineIntersection;

    private:
        int Nrow;
        int Ncol;
        float delta_phi;
        float delta_theta;
        void ProjectBoundaryOnPlane();
        void Hough(double error_angle, double error_distance);
        void RANSAC(double error, int tests);

};

#include "intersection.inl"

#endif // intersection
