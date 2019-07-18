#ifndef INTERSECTION
#define INTERSECTION

#include <Dense>
#include <Core>
#include "plane.h"
#include <set>
#include <Geometry>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/common/common.h>
#include <time.h>

const double eps = 0.00001;
const double max_plane_distance = 0.05;
const double max_line_distance = 0.02;
const int ransac_iterations = 100000;
const int min_number_points_on_line = 5;
const int min_number_pixels = 5;
const double min_line_length = 0.1;
const double perc_pixels_belonging_to_theoretical = 0.1;
const double line_margin = 0.02;
const int max_dist_between_pixels_in_line = 20;

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
            theoreticalLim.resize(0);
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
            theoreticalLim.resize(0);
        }
        int index;
        void setPlanes(plane* p1, plane* p2){plane_ref = p1; plane_neigh = p2;}
        Eigen::Vector3d tangente;
        Eigen::Vector3d normal;
        Eigen::Vector3d pt_mean;
        Eigen::Vector3d tangente_sister;
        Eigen::Vector3d normal_sister;
        Eigen::Vector3d pt_mean_sister;
        double distance_sister;
        Eigen::Vector2d tangente2D;
        Eigen::Vector2d normal2D;
        double distance;
        double distance2D;
        plane* plane_ref;
        plane* plane_neigh;
        std::vector<std::pair<float,float>> theoritical_phi_theta;
        std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> points;
        std::set<int> indices_self;
        std::set<int> indices_sister;
        std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> other_points;
        std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> points2D;
        std::set<std::pair<int,int>> theoritical_pixels;
        std::vector<std::pair<int,int>> pixels;
        std::vector<std::pair<int,int>> other_pixels;
        std::map<std::pair<int,int>, int> pixel2idx_point;
        std::map<int, std::pair<int,int>> idx_point2pixel;
        float start;
        float end;
        Eigen::Vector3d start_pt;
        Eigen::Vector3d end_pt;
        std::pair<int,int> start_pixel;
        std::pair<int,int> end_pixel;
        std::set<float> lim;
        void computeTheoriticalLineFeaturesObstruction(std::multimap<std::pair<int,int>,std::pair<int,int>> neighborPix2currentPix);
        void computeTheoriticalLineFeaturesOpening();
        void computeTheoriticalLineFeaturesConnection();
        void computeTheoriticalPhiTheta();
        void computeLim();
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
        std::vector<Eigen::Vector3d> theoreticalLim;
        std::vector<int> possible_corner_indices;
        std::vector<int> corner_indices;
        bool replaced_start;
        bool replaced_end;
        void ProjectLineFeaturesOnPlane(plane* p);
        std::set<int> repeated;
        std::set<int> not_repeated;
        void selectCorners(int rad, float delta_phi, float delta_theta, Eigen::Vector3d rot_axis, Eigen::Vector3d axis_init_phi);
        float length;
        void correctObstructions(intersection& sister);
        intersection export_sister();
        int idx_sister;
    private:
        int Nrow;
        int Ncol;
        float delta_phi;
        float delta_theta;
        void ProjectBoundaryOnPlane(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> pts_to_projet, plane* plane_to_project);
//        void Hough(double error_angle, double error_distance);
        void RANSAC(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& points_init, std::vector<std::pair<int,int>>& pixels_init, double error, int tests, std::set<int>& indices_line, std::set<int>& repeated );
        void SeparatePoints(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& points_of_plane_ref, std::vector<std::pair<int,int>>& pixels_of_plane_ref, std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& points_of_plane_neigh, std::vector<std::pair<int,int>>& pixels_of_plane_neigh);

};

#include "intersection.inl"

#endif // intersection
