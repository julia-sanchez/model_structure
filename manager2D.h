#ifndef MANAGER2D
#define MANAGER2D

#include <Dense>
#include <Core>
#include <Eigen>
#include <map>
#include "save_image_pgm.h"
#include <chrono>
#include "save_image_ppm.h"
#include <Eigen/StdVector>
#include "coeffWiseMulti.h"
#include "plane.h"

//const double epsilon = 0.00001;
const double margin_factor = 1.1;

class manager2D
{
public:
    manager2D(Eigen::MatrixXi p, int mc){image_clusterized = p; max_col = mc; Nrow = image_clusterized.rows(); Ncol = image_clusterized.cols();}
    manager2D(int mc);
    manager2D();

    //----------------------------------------------------------------------------------------------------------------------------------------

    void setImageClusterized(Eigen::MatrixXi p, std::pair<int,int> lt);
    int getMaxCol(){return max_col;}
    void setMaxCol(int mc){max_col = mc;}
    Eigen::MatrixXi getBoundariesImage();
    std::map<std::pair<int,int>, int> getBoundaries(){return XY2planeIdx_boundary;}
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> getImBinarized(){return image_bool;}
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> getImBinarizedMorpho(){return image_morpho;}
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> getOtherImBinarized(){return image_other_bool;}
    void closure(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> SE);
    void computeBoundaries(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> SE, Eigen::MatrixXi all_boundaries);
    void binarize(int p);
    std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator> getNeighbors(std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator it, int rad);
    int image_clusterized_idx;
    std::multimap<std::pair<int,int>,std::pair<int,int>> neighborPix2currentPix;
    std::multimap<std::pair<int,int>,std::pair<int,int>> getNeighborPix2currentPix(){return neighborPix2currentPix;};

private:
    int max_col;
    int Nrow;
    int Ncol;
    Eigen::MatrixXi image_clusterized; //(intensity = plane index)
    Eigen::MatrixXi image_grad; //(intensity = plane index)
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_bool;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_morpho;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_other_bool;
    std::map<std::pair<int,int>, int> XY2planeIdx_boundary;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> mynot(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_in);
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> morpho(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_in, Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> SE);
    void detect_margin();
    std::pair<int,int> lim_theta;
    std::multimap<std::pair<int,int>,std::pair<int,int>> reference;
};

#include "manager2D.inl"

#endif // manager2D
