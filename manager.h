#ifndef MANAGER
#define MANAGER

#include <Dense>
#include <Core>
#include <Eigen>
#include <map>
#include <Eigen/StdVector>
#include "save_image_pgm.h"
#include "plane.h"
#include "intersection.h"
#include "corner.h"
#include <stdio.h>

#include <chrono>
#include "manager3D.h"
#include "manager2D.h"

const int min_number_of_pixels = 50;                                           // min number of pixel for a cluster to represent a plane
const int pixel_radius_for_growing = 2;                                        // in case of soft obstruction over all the length  of the plane OR if black pixels (no return of laser bim)
const double normals_similarity_threshold_for_cleaning_when_one_cluster = 0.9;   // when one point belongs to one unique cluster if normal too different erase
const double normals_similarity_threshold_to_select_seed = 0.99;                 // comparison of seed normal to neighbors normals (seed must be on a precise zone of a plane)
const double normals_similarity_to_add_neighbors = 0.9;  // (when the region growing goes on another wall on a band) check the neighborhs normals to be sure not to go too far
const double not_parallel_dot_threshold = 0.9961;
const double angle_parallel_threshold = 15;
const int incertainty_pixels = 6;
const double min_dist_planes =0.1; //min distance between two parallel planes (if closer gathering)
const double min_polygone_angle = 10 * M_PI/180;

class manager
{
public:
    manager()
    {
        inter_remaining.isLine = false;
        inter_remaining.isObstruction = true;
        lim_theta.first = -1;
        lim_theta.second = -1;
        seeds_pixels = Eigen::MatrixXi::Zero(Nrow, Ncol);
    }
    manager(typename pcl::PointCloud<pcl::PointNormal>::Ptr c, int Nr, int Nc, double pa, double ta, int mc, Eigen::Vector3d ra, Eigen::Vector3d aip);

//    ----------------------------------------------------------------------------------------------------------------------------------------

    void setMap(std::map<std::pair<double, double>, std::pair<Eigen::Vector3d, Eigen::Vector3d>> p){PhiTheta2PtN = p;}
    int getMaxCol(){return max_col;}
    std::vector<Eigen::MatrixXi, Eigen::aligned_allocator<Eigen::MatrixXi>> getClusterizedImage(){return image_clusterized;}
    Eigen::MatrixXi getClusterizedIndicesImage(){return image_clusterized_indices;}
    std::vector<Eigen::MatrixXi, Eigen::aligned_allocator<Eigen::MatrixXi>> getInitImage(){return image_init;}
    void setMAxCol(int mc){max_col = mc;}
    void searchClusters(double thresh_plane_belonging, int radius, double thresh_neigh_for_seed);
    void init2Image();
    void clusterized2Image();
    void cleanClusters();
    std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator> getNeighbors(std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator it, int rad);
    void initSearchCluster(double rad);
    pcl::PointCloud<pcl::PointXYZ> extractBoundCloud();
    void extractBoundImage();
    void searchClusters2(double thresh_plane_belonging);
    void computeLines();
    void DefineIntersectionsPixPoints(int idx_plane);
    void computeTheoriticalPlanesIntersections();
    void setLimTheta(std::pair<int,int> lt){lim_theta = lt;}
    std::multimap<std::pair<int,int>,std::pair<int,int>> neighborPix2currentPix;
    void computeParallelLinesIntersections();
    Eigen::MatrixXi seeds_pixels;
    void order_polygone_corners();
    void export_mesh();

private:
    pcl::PointCloud<pcl::PointNormal>::Ptr cloud;
    int max_col;
    Eigen::Vector3d rot_axis;
    void setRotAxis(Eigen::Vector3d axis){rot_axis = axis;}
    Eigen::Vector3d axis_init_phi;
    void setAxisInitPhi(Eigen::Vector3d axis){axis_init_phi = axis;}
    int Nrow;
    int Ncol;
    double phi_app;
    double theta_app;
    float delta_phi;
    float delta_theta;
    //maps which link (phi, theta) or (X,Y) to each individual (xyz)(nx,ny,nz)
    std::map<std::pair<double, double>, std::pair<Eigen::Vector3d, Eigen::Vector3d>> PhiTheta2PtN; // phi theta
    std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>> XY2PtN;
    //we associate one mode plane to each point
    std::multimap<std::pair<int,int>, std::pair<int,int>> XY2Plane_idx; // to keep double values
    std::map<std::pair<int,int>, Eigen::Vector3d> XY2Plane_clusterized_clean; // to process points with various modes and associate one unique mode to them
    std::vector<std::pair<std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator>, Eigen::Vector3d>> regions; // vector of ( (vector of iterators) , (plane))
    std::vector<Eigen::MatrixXi, Eigen::aligned_allocator<Eigen::MatrixXi>> image_clusterized;
    Eigen::MatrixXi image_clusterized_indices;
    std::vector<Eigen::MatrixXi, Eigen::aligned_allocator<Eigen::MatrixXi>> image_init;
    std::map<std::pair<int,int>, int> boundary;
    std::vector<plane> planes;
    std::vector<intersection> intersections;
    std::vector<std::set<int>> edges;
    std::vector<intersection> all_edges;
    std::vector<corner> possible_corners;
    std::vector<corner> corners;
    bool isSeed(std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator it, int radius, double thresh_neigh_for_seed);
    void quantifyPhiTheta();
    Eigen::MatrixXi doubt;
    void computeConnections(int idx_plane);
    manager3D man3D;
    manager2D man2D;
    Eigen::MatrixXi processed_planes;
    void remove_intersection(int* k);
    intersection inter_remaining;
    std::vector<std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>> possible_inter3D;
    void computeTheoriticalLinesIntersections();
    std::pair<int,int> lim_theta;
    void detect_margin();
    Eigen::MatrixXi all_boundaries_image;
    bool replaceLim2(intersection& inter1, intersection& inter2, Eigen::Vector3d potential_corner_pt);
    bool replaceLim3(intersection& inter1, intersection& inter2, intersection& inter3, Eigen::Vector3d potential_corner_pt);
    void fill_edges();
    bool arePlanesClose(plane pi, plane pj);
    void actualizeChanged();
    bool DoesCrossPoly(std::vector<Eigen::Vector3d> jonction, std::vector<std::vector<Eigen::Vector3d>> polys, Eigen::Vector3d normal);
    bool isLineConnectedInPlane(int idx_line, int idx_plane, int end);
    void correctLinesCrossing();
    void clean_edges();
    std::vector<intersection> edges_removed;
    void projectOnPlane(intersection& inter, Eigen::Affine3d rot, Eigen::Vector2d& pinit, Eigen::Vector2d& pfin, Eigen::Vector2d& tangente2D, Eigen::Vector2d& normal2D, double* distance);
    void fill_edges_in_plane();
    void fusionCorners();
    void fill_edges_in_planes();
    bool DoesCross(Eigen::Affine3d rot, std::vector<Eigen::Vector3d> jonction, std::vector<Eigen::Vector3d> vec_tested);
    bool DoesCrossLines(std::vector<Eigen::Vector3d> jonction, plane& p);
    void recoverEqualStartEnd();

};

#include "manager.inl"

#endif // MANAGER
