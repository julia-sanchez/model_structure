#include "manager.h"

//#include "save_image_pgm.h"
#include "save_image_ppm.h"
#include <chrono>
#include <Eigen/StdVector>

int vecMed(std::vector<int> vec);

int main(int argc, char *argv[])
{
    if(argc != 10)
    {
        std::cout<< "file_pcd radius_for_normals radius_for_density_filter Ncol Nrow max_color luz depth_map_type hokuyo_ou_autre"<<std::endl<<std::endl;
        std::cout<< "on doit préciser si c'est hokuyo parce que le 0 n'est pas à l'origine du scanner (y)"<<std::endl<<std::endl;
    }

    // LOAD AND PREPROCESS CLOUD-------------------------------------------------------------------------------------------------------------------------------------------------

    std::string file_name = argv[1];
    pcl::PointCloud<pcl::PointNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointNormal>);

    if( pcl::io::loadPCDFile<pcl::PointNormal>( file_name, *cloud ) == -1 )
    {
        PCL_ERROR ("Couldn't read file test_pcd.pcd \n");
    }

    int hokuyo = atoi(argv[11]);
    std::cout<< "Pretransfrom to put sensor in (0,0)"<<std::endl;
    Eigen::Matrix4f rotation_transform = Eigen::Matrix4f::Identity();
    Eigen::Vector3d axis = {0,0,1};

    if(hokuyo)
    {
        axis = {0,1,0};
        std::cout<< "HOKUYO DETECTED"<<std::endl<<std::endl;
        rotation_transform (0,3) = 0.027;
        rotation_transform (1,3) = 0;
        rotation_transform (2,3) = -0.18;
//        pcl::transformPointCloudWithNormals(*cloud, *cloud, rotation_transform);

        for(int i = 0 ; i<cloud->points.size(); ++i)
        {
            double distance = sqrt(cloud->points[i].x*cloud->points[i].x + cloud->points[i].y*cloud->points[i].y + cloud->points[i].z*cloud->points[i].z);
            if(distance<0.3)
            {
                cloud->points.erase(cloud->points.begin()+i);
                --i;
            }
        }

        cloud->width = cloud->points.size();
    }

    pcl::io::savePCDFileASCII ("cloud_init.pcd", *cloud);
    system("bash pcd2csv.sh cloud_init.pcd");

    std::cout<< "Stop initial transformation"<<std::endl;

    // CONVERTIR EN CARTE DE PROFONDEUR-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    double normal_radius = atof(argv[2]);
    int max_col =atoi(argv[7]);
    int Nrow=atoi(argv[3]); //phi
    int Ncol=atoi(argv[4]); //theta
    double phi_app = atof(argv[5]);
    phi_app *= M_PI/180;
    double theta_app = atof(argv[6]);
    theta_app *= M_PI/180;
    double thresh_plane_belonging = atof(argv[8]);
    int radius_for_seed = atoi(argv[9]);
    double thresh_neigh_for_seed = atof(argv[10]);

    // FILTRE MORPHO POUR BOUCHER LES PIXELS NOIRS (SANS VALEUR)------------------------------------------------------------------------------------------------------------------------
    manager man(cloud, Nrow, Ncol, phi_app, theta_app, max_col);
    std::cout<< "Creating map"<<std::endl;
    man.initSearchCluster(normal_radius, axis);
    std::cout<<std::endl<< "Saving init"<<std::endl;
    man.init2Image();
    std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>, Eigen::aligned_allocator<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>> init = man.getInitImage();
    save_image_ppm("init", "", init, max_col);
    std::cout<< "Stop creating map"<<std::endl<<std::endl;

    std::cout<<std::endl<< "Start searching clusters"<<std::endl;
    man.searchClusters(thresh_plane_belonging, radius_for_seed, thresh_neigh_for_seed);
    std::cout<<std::endl<< "Stop searching clusters"<<std::endl<<std::endl;

    std::cout<<std::endl<< "Start cleaning clusters"<<std::endl;
    man.cleanClusters();
    std::cout<<std::endl<< "Stop cleaning clusters"<<std::endl<<std::endl;

    //do it twice with the updated normal to be sure one piece of plane is not avoided because of the high length of the plane and the uncertainty of the normal

//    std::cout<<std::endl<< "Research clusters"<<std::endl;
//    man.searchClusters2(0.06);
//    std::cout<<std::endl<< "Stop researching clusters"<<std::endl<<std::endl;

//    std::cout<<std::endl<< "Start recleaning clusters"<<std::endl;
//    man.cleanClusters();
//    std::cout<<std::endl<< "Stop recleaning clusters"<<std::endl<<std::endl;

    //save and display initial image
//    std::cout<<std::endl<< "Saving init"<<std::endl;
//    man.init2Image();
//    std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>, Eigen::aligned_allocator<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>> init = man.getInitImage();
//    save_image_ppm("init", "", init, max_col);

    //save and display clusterized image
    std::cout<<std::endl<< "Saving clusterized"<<std::endl;
    man.clusterized2Image(); // fill image_clusterized with colors of planes features and image_clusterized_indices with indices
    std::vector<Eigen::MatrixXi, Eigen::aligned_allocator<Eigen::MatrixXi>> clusterized = man.getClusterizedImage();
    save_image_ppm("clusterized", "", clusterized, max_col);
    Eigen::MatrixXi clusterized_indices = man.getClusterizedIndicesImage();
    save_image_pgm("clusterized_indices", "", clusterized_indices, clusterized_indices.maxCoeff());

    std::cout<<std::endl<< "Extract boundaries"<<std::endl;
    man.extractBoundImage();
//    Eigen::Vector3f normal = {1,0,0};
//    man.loadBoundary(("results/cloud_boundary_5.pcd"), normal);
//    man.searchLinesInBoundary(error_angle, error_distance);
}
