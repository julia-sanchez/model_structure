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
        std::cout<< "file_pcd radius_for_normals Nrow Ncol phi_app theta_app max_col thresh_plane_belonging radius_for_seed thresh_neigh_for_seed lim_theta hokuyo_ou_autre"<<std::endl<<std::endl;
        //                           0.15        1200 1200   360      180      500              0.15               5               0.03                0            0
        std::cout<< "on doit préciser si c'est hokuyo parce que le 0 n'est pas à l'origine du scanner (y)"<<std::endl<<std::endl;
    }

    // LOAD AND PREPROCESS CLOUD-------------------------------------------------------------------------------------------------------------------------------------------------

    std::string file_name = argv[1];
    pcl::PointCloud<pcl::PointNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointNormal>);

    if( pcl::io::loadPCDFile<pcl::PointNormal>( file_name, *cloud ) == -1 )
    {
        PCL_ERROR ("Couldn't read file test_pcd.pcd \n");
    }

    int hokuyo = atoi(argv[12]);
    std::cout<< "Pretransfrom to put sensor in (0,0)"<<std::endl;
    Eigen::Matrix4f rotation_transform = Eigen::Matrix4f::Identity();
    Eigen::Vector3d rot_axis = {0,0,1};
    Eigen::Vector3d axis_init_phi = {1,0,0};

    if(hokuyo)
    {
        rot_axis = {0,1,0};
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
    std::pair<int,int> lim_theta;
    if(atoi(argv[11])>epsilon)
    {
        lim_theta.first = atoi(argv[11]);
        lim_theta.second = Ncol - atoi(argv[11]);
    }
    else
    {
        lim_theta.first = -1;
        lim_theta.second = -1;
    }

    manager man(cloud, Nrow, Ncol, phi_app, theta_app, max_col, rot_axis, axis_init_phi);
    man.setLimTheta(lim_theta);
    std::cout<< "Creating map"<<std::endl;
    man.initSearchCluster(normal_radius);
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

    //do it a second time with the updated normal and a new seed to be sure one piece of plane is not avoided because of the high length of the plane and the uncertainty of the normal

//    std::cout<<std::endl<< "Research clusters"<<std::endl;
//    man.searchClusters2(0.03);
//    std::cout<<std::endl<< "Stop researching clusters"<<std::endl<<std::endl;

//    std::cout<<std::endl<< "Start recleaning clusters"<<std::endl;
//    man.cleanClusters();
//    std::cout<<std::endl<< "Stop recleaning clusters"<<std::endl<<std::endl;

    //save and display clusterized image
    std::cout<<std::endl<< "Saving clusterized"<<std::endl;
    man.clusterized2Image(); // fill image_clusterized with colors of planes features and image_clusterized_indices with indices
    std::vector<Eigen::MatrixXi, Eigen::aligned_allocator<Eigen::MatrixXi>> clusterized = man.getClusterizedImage();
    save_image_ppm("clusterized", "", clusterized, max_col);
    Eigen::MatrixXi clusterized_indices = man.getClusterizedIndicesImage();
    save_image_pgm("clusterized_indices", "", clusterized_indices, clusterized_indices.maxCoeff());

//    std::vector<Eigen::MatrixXi, Eigen::aligned_allocator<Eigen::MatrixXi>> clusterized_with_seeds(3);
//    for(int i = 0; i<3; ++i )
//        clusterized_with_seeds[i] = Eigen::MatrixXi::Zero(Nrow, Ncol);

//    for(int i = 0; i < Nrow; ++i)
//    {
//        for(int j = 0; j < Ncol; ++j)
//        {
//            if(man.seeds_pixels(i,j))
//            {
//                clusterized_with_seeds[0](i,j) = max_col;
//                clusterized_with_seeds[1](i,j) = max_col;
//                clusterized_with_seeds[2](i,j) = max_col;
//            }
//            else
//            {
//                clusterized_with_seeds[0](i,j) = clusterized[0](i,j);
//                clusterized_with_seeds[1](i,j) = clusterized[1](i,j);
//                clusterized_with_seeds[2](i,j) = clusterized[2](i,j);
//            }
//        }
//    }
//    save_image_ppm("seeds", "", clusterized_with_seeds, max_col);

    std::cout<<std::endl<< "Extract boundaries"<<std::endl;
    man.extractBoundImage();
}
