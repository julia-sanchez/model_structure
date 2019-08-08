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
const double max_line_distance = 0.03;
const int ransac_iterations = 50000;
const int min_number_points_on_line = 4;
const int min_number_pixels = 4;
const double min_line_length = 0.03;
const double perc_pixels_belonging_to_theoretical = 0.1;
const double line_margin = 0.01;
const int max_dist_between_pixels_in_line = 5;
const float max_dist_between_points_in_line = 0.1;
const float max_spat_dist_for_line_connection = 0.15;
const float max_pix_dist_for_line_connection = 10;

class intersection
{
    public:
        intersection()
        {
            delta_phi = 0.01;
            delta_theta = 0.01;
            isConnection = false;
            isObstruction = false;
            isObject = false;
            isOpening = false;
            isLine = true;
            Nrow = (int)((2*M_PI+2*eps)/delta_phi);
            Ncol = (int)((M_PI+2*eps)/delta_theta);
            theoreticalLim.resize(0);
            max_pixel_diff_start = max_pix_dist_for_line_connection;
            max_spatial_diff_start = max_spat_dist_for_line_connection;
            max_pixel_diff_end = max_pix_dist_for_line_connection;
            max_spatial_diff_end = max_spat_dist_for_line_connection;
            has_sister = false;
            has_points_after_ls = false;
            new_start_pt = Eigen::Vector3d::Zero();
            new_end_pt = Eigen::Vector3d::Zero();
            isConnected = false;
        }
        intersection(plane* p1, plane* p2, float dp, float dt)
        {
            isConnection = false;
            isObstruction = false;
            isObject = false;
            isOpening = false;
            isLine = false;
            plane_ref = p1;
            plane_neigh = p2;
            delta_phi = dp;
            delta_theta = dt;
            Nrow = (int)((2*M_PI+2*eps)/delta_phi);
            Ncol = (int)((M_PI+2*eps)/delta_theta);
            theoreticalLim.resize(0);
            max_pixel_diff_start = max_pix_dist_for_line_connection;
            max_spatial_diff_start = max_spat_dist_for_line_connection;
            max_pixel_diff_end = max_pix_dist_for_line_connection;
            max_spatial_diff_end = max_spat_dist_for_line_connection;
            has_sister = false;
            has_points_after_ls = false;
            new_start_pt = Eigen::Vector3d::Zero();
            new_end_pt = Eigen::Vector3d::Zero();
            isConnected = false;
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
        std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> remaining_points;
        std::vector<std::pair<int,int>> remaining_pixels;
        std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> points2D;
        std::set<std::pair<int,int>> theoritical_pixels;
        std::vector<std::pair<int,int>> pixels;
        std::vector<std::pair<int,int>> other_pixels;
        std::map<std::pair<int,int>, int> pixel2idx_point;
        std::map<int, std::pair<int,int>> idx_point2pixel;
        float start;
        float end;
        Eigen::Vector3d start_pt;
        Eigen::Vector3d new_start_pt;
        Eigen::Vector3d end_pt;
        Eigen::Vector3d new_end_pt;
        std::pair<int,int> start_pixel;
        std::pair<int,int> end_pixel;
        std::set<float> lim;
        void computeTheoriticalLineFeaturesObstruction(std::multimap<std::pair<int,int>,std::pair<int,int>> neighborPix2currentPix);
        void computeTheoriticalLineFeaturesObject(std::multimap<std::pair<int,int>,std::pair<int,int>> neighborPix2currentPix);
        void computeTheoriticalLineFeaturesOpening();
        void computeTheoriticalLineFeaturesConnection();
        void computeTheoriticalPhiTheta();
        bool computeLim();
        void setDeltas(float dp, float dt);
        void definePlaneConnection();
        bool isConnection;
        bool isObstruction;
        bool isObject;
        bool isOpening;
        bool isLine;
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
        float max_pixel_diff_start;
        float max_spatial_diff_start;
        float max_pixel_diff_end;
        float max_spatial_diff_end;
        void ProjectLineFeaturesOnPlane();
        std::set<int> repeated;
        std::set<int> not_repeated;
        void selectCorners(int rad, float delta_phi, float delta_theta, Eigen::Vector3d rot_axis, Eigen::Vector3d axis_init_phi);
        float length;
        void correctObstructions(intersection& sister);
        intersection export_sister();
        int idx_sister;
        std::vector<intersection> export_sisters();
        bool isDoubled(std::vector<intersection> all_edges);
        bool has_sister;
        bool isNearCorner(std::vector<Eigen::Vector3d> corners_pt);
        bool isConnected;
    private:
        int Nrow;
        int Ncol;
        float delta_phi;
        float delta_theta;
        void ProjectBoundaryOnPlane(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> pts_to_projet, plane* plane_to_project);
//        void Hough(double error_angle, double error_distance);
        void RANSAC(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& points_init, std::vector<std::pair<int,int>>& pixels_init, int tests, std::set<int>& indices_line, std::set<int>& repeated );
        void SeparatePoints(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& points_of_plane_ref, std::vector<std::pair<int,int>>& pixels_of_plane_ref, std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& points_of_plane_neigh, std::vector<std::pair<int,int>>& pixels_of_plane_neigh, std::set<int>& repeated_ref, std::set<int>& repeated_neigh);
        void searchLine( std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& points_of_plane, std::vector<std::pair<int,int>>& pixels_of_plane, plane* p, std::set<int>& indices_line_plane, std::set<int>& repeated_plane);
        bool has_points_after_ls;
};

#include "intersection.inl"

#endif // intersection
