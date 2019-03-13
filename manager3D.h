#ifndef MANAGER3D
#define MANAGER3D

#include <vector>
#include <string>
#include <iostream>

// PCL INCLUDES

#include <pcl/io/pcd_io.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/point_types.h>
#include <pcl/common/common.h>
#include <pcl/common/transforms.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/filters/conditional_removal.h>
#include <pcl/filters/random_sample.h>
#include <pcl/search/kdtree.h>

#define USE_OPENMP_FOR_NORMEST

const double epsilon = 0.00001;
const double theta_margin = 5*M_PI/180;

class manager3D
{
public:
    manager3D();
    manager3D(typename pcl::PointCloud<pcl::PointNormal>::Ptr c){cloud = c;}

    //----------------------------------------------------------------------------------------------------------------------------------------
    void setInputCloud(typename pcl::PointCloud<pcl::PointNormal>::Ptr c){cloud = c; tree.setInputCloud(cloud);}
    pcl::PointCloud<pcl::PointNormal>::Ptr getInputCloud(){return cloud;}
    std::map<std::pair<double, double>, std::pair<Eigen::Vector3d, Eigen::Vector3d>> getMap(){return planes;}
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> getImageBound(){return image_boundaries;}
    void setRotAxis(Eigen::Vector3d axis){rot_axis = axis;}
    int getSize() const {return cloud->size();}
    void getNormals(double);
    void sample(double);
    void createMap();
    bool isSeed(Eigen::Vector3d& pt, Eigen::Vector3d& normal, double radius, double threshold);
    Eigen::Vector3d rot_axis;
    Eigen::Vector3d axis_init_phi;

private:
    pcl::PointCloud<pcl::PointNormal>::Ptr cloud;
    pcl::search::KdTree<pcl::PointNormal> tree;
    std::map<std::pair<double, double>, std::pair<Eigen::Vector3d,Eigen::Vector3d>> planes;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_boundaries;
    void orient();
    void pcn2pc(pcl::PointCloud<pcl::PointXYZ>::Ptr);
    void insertInMap(double phi, double theta, int i);
};

#include "manager3D.inl"

#endif // manager3D
