manager::manager(typename pcl::PointCloud<pcl::PointNormal>::Ptr c, int Nr, int Nc, double pa, double ta, int mc, Eigen::Vector3d ra, Eigen::Vector3d aip)
{
    cloud = c;
    man3D.setInputCloud(c);
    Nrow = Nr;
    Ncol = Nc;
    phi_app = pa;
    theta_app = ta;
    max_col = mc;
    image_clusterized_indices = Eigen::MatrixXi::Zero(Nrow, Ncol);
    image_clusterized.resize(3);
    image_init.resize(3);
    for(int i = 0; i<3; ++i )
    {
        image_clusterized[i] = Eigen::MatrixXi::Zero(Nrow, Ncol);
        image_init[i] = Eigen::MatrixXi::Zero(Nrow, Ncol);
    }
    inter_remaining.isLine = false;
    inter_remaining.isObstruction = true;
    lim_theta.first = -1;
    lim_theta.second = -1;
    rot_axis = ra;
    axis_init_phi = aip;
}

void manager::quantifyPhiTheta()
{
    int n = 0;
    XY2PtN.clear();
    delta_phi = (phi_app+epsilon)/Nrow;
    delta_theta = (theta_app+epsilon)/Ncol;
    for(auto it = PhiTheta2PtN.begin(); it!=PhiTheta2PtN.end(); ++it)
    {
        int X = (int)(it->first.first/delta_phi);
        int Y = (int)(it->first.second/delta_theta);
        auto idx_is_inserted = XY2PtN.insert(std::make_pair(std::make_pair(X,Y), it->second));
        if(!idx_is_inserted.second)
            ++n;
    }
    std::cout<<"Number of points not inserted into XY2PtN because same X and same Y (due to 3d noise affecting small theta): "<<n<<std::endl<<std::endl;
}

void manager::initSearchCluster(double rad)
{
    man3D.getNormals(rad);
    man3D.setRotAxis(rot_axis);
    man3D.setAxisInitPhi(axis_init_phi);
    man3D.createMap();
    PhiTheta2PtN = man3D.getMap();
    quantifyPhiTheta();
    detect_margin();
}

void manager::init2Image()
{
    double maxi_dist = epsilon;
    for(auto it = XY2PtN.begin(); it!=XY2PtN.end(); ++it)
    {
        double dist = abs(it->second.second.dot(it->second.first));
        if(maxi_dist<dist)
            maxi_dist = dist;
    }

    maxi_dist *= margin_factor;

    for(auto it = XY2PtN.begin(); it!=XY2PtN.end(); ++it)
    {
        Eigen::Vector3i rgb;
        Eigen::Vector3d temp;
        temp = (it->second.second + Eigen::Vector3d::Ones())/2;
        double dist = abs(it->second.second.dot(it->second.first));
        rgb = ( max_col*(1 - dist/maxi_dist)*temp ) .cast<int>();
        rgb = (rgb.cast<double>() * (1-dist/maxi_dist)).cast<int>();
        image_init[0](it->first.first, it->first.second) = rgb(0);
        image_init[1](it->first.first, it->first.second) = rgb(1);
        image_init[2](it->first.first, it->first.second) = rgb(2);
    }
}

void manager::clusterized2Image()
{
    double maxi_dist = epsilon;
    for(int i = 0; i<planes.size(); ++i)
    {
        double dist = planes[i].distance;
        if(maxi_dist<dist)
            maxi_dist = dist;
    }

    maxi_dist *= margin_factor;

    //first fill image_clusterized_indices with detected clusters
    for(int i = 0; i<planes.size(); ++i)
    {
        Eigen::Vector3i rgb;
        Eigen::Vector3d temp;
        temp = (planes[i].normal + Eigen::Vector3d::Ones())/2;
        double dist = planes[i].distance;
        temp = max_col*(1 - dist/maxi_dist)*temp;
        rgb = temp.cast<int>();

        for(auto it_pix = planes[i].pixels.begin(); it_pix != planes[i].pixels.end(); ++it_pix)
        {
            image_clusterized_indices(it_pix->first, it_pix->second) = planes[i].index+1;
            image_clusterized[0](it_pix->first, it_pix->second) = rgb(0);
            image_clusterized[1](it_pix->first, it_pix->second) = rgb(1);
            image_clusterized[2](it_pix->first, it_pix->second) = rgb(2);
        }
    }

    //then fill remaining points which do not belong to any cluster (other objects)
    for(int i = 0; i<Nrow; ++i)
    {
        for(int j = 0; j<Ncol; ++j)
        {
            if(image_clusterized_indices(i,j)==0) //pixel not filled with existing cluster....
            {
                auto idx_found_XY2PtN = XY2PtN.find(std::make_pair(i,j));
                if(idx_found_XY2PtN!=XY2PtN.end()) // ... but existing in first place
                    image_clusterized_indices(i,j) = planes.size()+1; // create new cluster (index : planes.size()+1) and put indice of this new cluster
            }
        }
    }
}

bool manager::isSeed(std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator it, int radius, double thresh_neigh_for_seed)
{
    //les voisins sont choisis sur l'image donc ils peuvent etre trèès proches en points 3D
    std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator> neigh = getNeighbors(it, radius);
    Eigen::Vector3d mean_pt;
    mean_pt = it->second.first;

    for(int i = 0; i< neigh.size(); ++i)
    {
        double dist_neigh = abs(neigh[i]->second.first.dot(neigh[i]->second.second));
        double dist = abs(mean_pt.dot(it->second.second));
        if( abs(dist- dist_neigh) > thresh_neigh_for_seed)
            return false;
//        double dist = abs((mean_pt-neigh[i]->second.first).dot(it->second.second));
//        if( dist > thresh_neigh_for_seed || abs(neigh[i]->second.second.dot(it->second.second))<normals_similarity_threshold_to_select_seed)
//            return false;
    }
    return true;
}

void manager::searchClusters(double thresh_plane_belonging, int radius, double thresh_neigh_for_seed) // create the vector of "regions" in which there is a vector of points with their associated mode
{
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> processed_seed = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> processed_growing = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();
    pcl::PointCloud<pcl::PointXYZ> pc_clus;
    regions.resize(0);

    for(auto it = XY2PtN.begin(); it != XY2PtN.end(); ++it)
    {
        if(!processed_seed(it->first.first, it->first.second) && it->first.second > (lim_theta.first+10) && it->first.second < (lim_theta.second-10)) // in order not to process one point which has already been selected in a cluster
        {
            pc_clus.points.resize(0);
            if(isSeed(it, radius, thresh_neigh_for_seed))
            {
                processed_growing = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();
                std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator> region; // vector of iterator
                region.push_back(it);

                for(int i = 0; i<region.size(); ++i)
                {
                    std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator> neighbors = getNeighbors(region[i], pixel_radius_for_growing);
                    bool at_least_one_inlier = false;
                    for(int k = 0; k<neighbors.size(); ++k)
                    {
                        if(!processed_growing(neighbors[k]->first.first, neighbors[k]->first.second)) //in order not to process twice one neighbor pixel
                        {
                            if(abs(neighbors[k]->second.second.dot(it->second.second))>normals_similarity_to_add_neighbors)
                            {
                                at_least_one_inlier = true;
                                break;
                            }
                        }
                    }
                    if(at_least_one_inlier)
                    {
                        for(int k = 0; k<neighbors.size(); ++k)
                        {
                            if(!processed_growing(neighbors[k]->first.first, neighbors[k]->first.second)) //in order not to process twice one neighbor pixel
                            {
                                if(abs((neighbors[k]->second.first-it->second.first).dot(it->second.second)) < thresh_plane_belonging)// && acos(neighbors[k]->second.dot(it->second)/(neighbors[k]->second.norm()*it->second.norm())) < thresh_angle )
                                {
                                    if(neighbors[k]->first.first == 607 && neighbors[k]->first.second == 147)
                                    {
                                        std::cout<<"pixel bizarre"<<std::endl;
                                        std::cout<<abs((neighbors[k]->second.first-it->second.first).dot(it->second.second))<<std::endl;
                                        std::cout<<it->second.second.transpose()<<std::endl;
                                    }
                                    processed_growing(neighbors[k]->first.first, neighbors[k]->first.second) = true;

                                    region.push_back(neighbors[k]);

//                                    pcl::PointXYZ pt_clus;
//                                    pt_clus.x = neighbors[k]->second.first(0);
//                                    pt_clus.y = neighbors[k]->second.first(1);
//                                    pt_clus.z = neighbors[k]->second.first(2);
//                                    pc_clus.points.push_back(pt_clus);
                                }
                            }
                        }
                    }
                }

                if(region.size()>min_number_of_pixels)
                {
                    processed_seed = processed_seed || processed_growing;
                    Eigen::Vector3d cluster_normal = it->second.second;
                    double cluster_distance = abs(it->second.first.dot(it->second.second));
                    cluster_normal /= cluster_normal.norm();
                    cluster_normal *= cluster_distance;
                    regions.push_back(std::make_pair(region, cluster_normal));
                    std::cout<<"New region number: "<< regions.size()-1 << std::endl<< "seed : "<<it->first.first<<" "<<it->first.second<<"-->"<<it->second.first.transpose()<<std::endl<<"normal : "<<it->second.second.transpose()<<std::endl<<"distance : "<<cluster_distance<<std::endl<<"number of points : "<<region.size()<<std::endl<<std::endl;
//                    std::vector<Eigen::MatrixXi, Eigen::aligned_allocator<Eigen::MatrixXi>> test(3);
//                    test[0] = coeffWiseMulti(image_init[0], processed_seed.cast<int>());
//                    test[1] = coeffWiseMulti(image_init[1], processed_seed.cast<int>());
//                    test[2] = coeffWiseMulti(image_init[2], processed_seed.cast<int>());
//                    save_image_ppm("processed", "", test,max_col);
//                    save_image_pgm("processed", "", processed_growing.cast<int>(), 1);

////                    pc_clus.width = pc_clus.points.size();
////                    pc_clus.height = 1;
////                    pc_clus.is_dense=false;
////                    pcl::io::savePCDFileASCII ("cloud_clus.pcd", pc_clus);

////                    system("bash pcd2csv.sh cloud_clus.pcd");

//                    std::cout<<"searchClusters 1 "<<std::endl<<std::endl;
//                    getchar();
//                    getchar();
                }
            }
        }
    }
}

void manager::cleanClusters() // from regions, keep all points in XY2Plane_idx and then clean points having various modes
{
    //regions = vector of ( (vector of iterators) , (plane)) (with iterator = XY2PtN (points on regions))
    XY2Plane_idx.clear();
    for(int k = 0; k < regions.size(); ++k) // for each cluster
    {
        for(int i = 0; i<regions[k].first.size(); ++i)
            XY2Plane_idx.insert(std::make_pair(regions[k].first[i]->first, std::make_pair(k,i))); // Plane_idx = pair(idx_region, idx_XY2PtN)
    }

    std::cout<<"Number of points in XY2Plane_idx "<<XY2Plane_idx.size()<<std::endl<<std::endl;

    std::cout<<"Start cleaning points having various clusters "<<std::endl<<std::endl;
    int n =0;
    planes.clear();
    planes.resize(regions.size());
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> share_points(planes.size(),planes.size());
    std::pair<std::pair<int,int>, std::pair<int,int>> inserted;
    for(auto it_XY2Plane_idx = XY2Plane_idx.begin(); it_XY2Plane_idx!=XY2Plane_idx.end(); ++it_XY2Plane_idx) // for each cluster
    {
        auto it_XY2PtN = regions[it_XY2Plane_idx->second.first].first[it_XY2Plane_idx->second.second];

        std::vector<std::multimap<std::pair<int,int>, std::pair<int,int>>::iterator> sameXY; // points XY with all associated clusters
        int Xinit = it_XY2Plane_idx->first.first;
        int Yinit = it_XY2Plane_idx->first.second;
        while(it_XY2Plane_idx->first.first == Xinit && it_XY2Plane_idx->first.second == Yinit)
        {
            sameXY.push_back(it_XY2Plane_idx);
            ++it_XY2Plane_idx;
        }

        --it_XY2Plane_idx;

        if(sameXY.size()>1) // points which belong to various clusters
        {

            //Count points which have multiple modes
            ++n;

            //Select the best cluster (mode most similar)

            //order and select the most likely knowing its own normal
            std::vector<int> planes_indices;
            std::map<double, std::pair<int,int>> dotty; // pour mapper
            for(int k = 0; k<sameXY.size(); ++k)
            {
                planes_indices.push_back(sameXY[k]->second.first);
                double dot = (regions[sameXY[k]->second.first].second/regions[sameXY[k]->second.first].second.norm()).dot(it_XY2PtN->second.second); // compute scalar product between real normal and region normal
                dotty.insert(std::make_pair(dot, sameXY[k]->second));                     // keep result in map to order scalar product and select the region normal which is more similar
            }

            for(int i = 0; i<planes_indices.size()-1; ++i)
            {
                for(int j = i+1; j<planes_indices.size(); ++j)
                {
                    share_points(planes_indices[i],planes_indices[j]) = true;
                    share_points(planes_indices[j],planes_indices[i]) = true;
                }
            }


            auto it_map = dotty.end();
            --it_map;                                                               // it_map is the ptr of [ scalar_product, iter] with iter corresponding to the more similar region normal
            if(abs(regions[it_map->second.first].second.norm()-abs(it_XY2PtN->second.first.dot((regions[it_map->second.first].second/regions[it_map->second.first].second.norm())))) > 0.1)
                --it_map;

            if( abs( ((regions[it_map->second.first].second/regions[it_map->second.first].second.norm()).dot(regions[it_map->second.first].first[it_map->second.second]->second.second)) ) > normals_similarity_threshold_for_cleaning_when_one_cluster )
                inserted = std::make_pair(it_XY2Plane_idx->first, it_map->second);
            else
                inserted = std::make_pair(std::make_pair(-1,-1), std::make_pair(-1,-1));
        }
        else // points which belong to one unique cluster
        {
            if( abs( ((regions[it_XY2Plane_idx->second.first].second/regions[it_XY2Plane_idx->second.first].second.norm()).dot(it_XY2PtN->second.second)) ) > normals_similarity_threshold_for_cleaning_when_one_cluster )
                inserted = *it_XY2Plane_idx;
            else
                inserted = std::make_pair(std::make_pair(-1,-1), std::make_pair(-1,-1));
        }

        // To recompute planes features with cleaned points :
        if(inserted.first.first>=0) // if the point was inserted
        {
            planes[inserted.second.first].appendPoint(it_XY2PtN->second.first);
            planes[inserted.second.first].appendPixel(inserted.first);
        }
    }

    std::cout<<"Stop cleaning points having various clusters "<<std::endl<<std::endl;

    //------------------------------------------------------------------------------------------------------------------

    std::cout<<"Start recomputing normal and distance "<<std::endl<<std::endl;
    //Recompute features

    for(int i = 0; i<planes.size(); ++i)
            planes[i].computeNormal();                      //recompute normal with points detected in plane

    std::cout<<"Stop recomputing normal and distance "<<std::endl<<std::endl;

    //------------------------------------------------------------------------------------------------------------------
    //gather planes with similar features

    std::vector<int> to_erase;

    for(int i = 0; i<planes.size()-1; ++i)
    {
        for(int j = i+1; j<planes.size(); ++j)
        {
            if(acos(planes[i].normal.dot(planes[j].normal))<5*M_PI/180 && abs(planes[i].distance-planes[j].distance)<0.1 && share_points(i,j)) // if planes similar...
            {
                for(int k = 0; k< planes[j].pts.size(); ++k)
                {
                    planes[i].appendPoint(planes[j].pts[k]); //gather points
                    planes[i].appendPixel(planes[j].pixels[k]); // gather pixels
                }
                planes[i].computeNormal(); //recompute normal
                to_erase.push_back(j); //erase second plane
            }
        }
    }

    for(int k = 0; k<to_erase.size(); ++k)
    {
        planes.erase(planes.begin()+to_erase[k]);
        regions.erase(regions.begin()+to_erase[k]);
    }

    for(int i = 0; i<planes.size(); ++i)
        planes[i].index = i;

    //------------------------------------------------------------------------------------------------------------------

    processed_planes = Eigen::MatrixXi::Zero(planes.size()+1, planes.size()+1).cast<bool>();
    std::cout<<"Number of points which had multiple modes : "<<n<<std::endl<<std::endl;
}


void manager::searchClusters2(double thresh_plane_belonging) // create the vector of "regions" in which there is a vector of points with their associated mode
{
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> processed_seed = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> processed_growing = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();
    pcl::PointCloud<pcl::PointXYZ> pc_clus;

    regions.resize(0);

    for(int planes_it = 0; planes_it<planes.size(); ++planes_it)
    {
        planes[planes_it].mean();       //compute new seed (point closest from the mean position)
        pc_clus.points.resize(0);
        processed_growing = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();
        std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator> region; // vector of iterators of XY2PtN
        auto iter_seed = XY2PtN.find(planes[planes_it].seed.first);
        region.push_back(iter_seed);

        for(int i = 0; i<region.size(); ++i)
        {
            std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator> neighbors = getNeighbors(region[i], pixel_radius_for_growing);
            bool at_least_one_inlier = true;
            for(int k = 0; k<neighbors.size(); ++k)
            {
                if(abs(neighbors[k]->second.second.dot(planes[planes_it].normal))>normals_similarity_to_add_neighbors && abs((neighbors[k]->second.first-planes[planes_it].seed.second).dot(planes[planes_it].normal)) < thresh_plane_belonging)
                {
                    at_least_one_inlier = false;
                    break;
                }
            }
            if(!at_least_one_inlier)
            {
                for(int k = 0; k<neighbors.size(); ++k)
                {
                    if(!processed_growing(neighbors[k]->first.first, neighbors[k]->first.second)) //in order not to process twice one pixel
                    {
                            processed_growing(neighbors[k]->first.first, neighbors[k]->first.second) = true;
                            region.push_back(neighbors[k]);

//                            pcl::PointXYZ pt_clus;
//                            pt_clus.x = neighbors[k]->second.first(0);
//                            pt_clus.y = neighbors[k]->second.first(1);
//                            pt_clus.z = neighbors[k]->second.first(2);
//                            pc_clus.points.push_back(pt_clus);
                    }
                }
            }
        }

        if(region.size()>min_number_of_pixels)
        {
            Eigen::Vector3d cluster_normal = planes[planes_it].normal;
            double cluster_distance = abs(planes[planes_it].seed.second.dot(planes[planes_it].normal));
            cluster_normal /= cluster_normal.norm();
            cluster_normal *= cluster_distance;
            regions.push_back(std::make_pair(region, cluster_normal));

            processed_seed = processed_seed || processed_growing;
            std::cout<<"New region number: "<< regions.size()-1 << std::endl<< "seed : "<<planes[planes_it].seed.first.first<<" "<<planes[planes_it].seed.first.second<<"-->"<<planes[planes_it].seed.second.transpose()<<std::endl<<"normal : "<<planes[planes_it].normal.transpose()<<std::endl<<"distance : "<<cluster_distance<<std::endl<<"number of points : "<<region.size()<<std::endl<<std::endl;
//            std::vector<Eigen::MatrixXi, Eigen::aligned_allocator<Eigen::MatrixXi>> test(3);
//            test[0] = coeffWiseMulti(image_init[0], processed_seed.cast<int>());
//            test[1] = coeffWiseMulti(image_init[1], processed_seed.cast<int>());
//            test[2] = coeffWiseMulti(image_init[2], processed_seed.cast<int>());
//            save_image_ppm("processed", "", test,max_col);
//            save_image_pgm("processed", "", processed_growing.cast<int>(), 1);

//            pc_clus.width = pc_clus.points.size();
//            pc_clus.height = 1;
//            pc_clus.is_dense=false;
//            pcl::io::savePCDFileASCII ("cloud_clus.pcd", pc_clus);

//            system("bash pcd2csv.sh cloud_clus.pcd");
//            getchar();
//            getchar();
//            getchar();
        }
    }
}


std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator> manager::getNeighbors(std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator it, int rad)
{
    std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator> neighbors;
    int initx;
    int finitx;
    if(it->first.first-rad>0)
        initx = it->first.first-rad;
    else
        initx = Nrow - (rad - it->first.first);

    if(it->first.first+rad<Nrow)
        finitx = it->first.first+rad;
    else
        finitx = it->first.first+rad - Nrow;

    int inity = std::max(it->first.second-rad, 0);
    int finity = std::min(it->first.second+rad, Ncol-1);

    for(int i = initx; i != finitx; ++i)
    {
        if(i>=Nrow)
            i = -1;
        else
        {
            for(int j = inity; j <= finity; ++j)
            {
                auto idx = XY2PtN.find(std::make_pair(i,j));
                if(idx != XY2PtN.end())
                    neighbors.push_back(idx);
            }
        }
    }
    return neighbors;
}

bool cmp_vector3d(Eigen::Vector3d a, Eigen::Vector3d b)
{
    if(abs(a(0)-b(0))>epsilon)
        return (a(0)<b(0));
    else if (abs(a(1)-b(1))>epsilon)
        return (a(1)<b(1));
    else
        return (a(2)<b(2));
}

void manager::extractBoundImage()
{
    man2D.setImageClusterized(image_clusterized_indices, lim_theta); // image_clusterized_indices is a matrix which contains the indices of planes (+1 to let 0 be the non cluster)
    man2D.setMaxCol(max_col); //maximum color
    pcl::PointCloud<pcl::PointXYZ> all_boundaries;
    edges.resize(planes.size()+1);
    system("rm results/binarized*.pgm");
    system("rm results/morpho*.pgm");
    system("rm results/boundary*.pgm");
    system("rm results/cloud_boundary_*.pcd");
    system("rm results/all_boundaries.pcd");
    system("rm results/*.csv");
    system("rm theoritical_line*.csv");
    system("rm Vertices.csv");
    system("rm cloud_boundary_*.pcd");
    system("rm cloud_boundary_*.csv");
    system("rm real_*.csv");
    system("rm theoretical_lines.csv");
    system("rm real_lines.csv");
    system("rm not_lines.csv");
    system("rm not_line*.csv");

    std::stringstream stm;
    int n = 0;
    all_edges.resize(0);

    all_boundaries_image = Eigen::MatrixXi::Zero(image_clusterized_indices.rows(), image_clusterized_indices.cols());

    for(int k = 0; k<planes.size(); ++k)
    {
        std::cout<<std::endl<<"PLANE NUMBER : "<<k<<std::endl;
        //binarize
        man2D.binarize(k);
        Eigen::Matrix<bool,Eigen::Dynamic, Eigen::Dynamic> binarized = man2D.getImBinarized();
        save_image_pgm("binarized", std::to_string(k), binarized.cast<int>()*max_col, max_col);

        //apply closure
        Eigen::Matrix<bool, 5, 5> SE_morpho;
        SE_morpho<<true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true;
        man2D.closure(SE_morpho);
        Eigen::Matrix<bool,Eigen::Dynamic, Eigen::Dynamic> morpho = man2D.getImBinarizedMorpho();
        save_image_pgm("morpho", std::to_string(k), morpho.cast<int>()*max_col, max_col);

        //compute boundary on image
        Eigen::Matrix<bool, 3, 3> SE_grad;
        SE_grad<<true, true, true, true, true, true, true, true, true;
        man2D.computeBoundaries(SE_grad, all_boundaries_image);
        boundary = man2D.getBoundaries();

        Eigen::MatrixXi boundary_image = man2D.getBoundariesImage();
        all_boundaries_image += boundary_image;
        save_image_pgm("boundary", std::to_string(k), boundary_image, max_col);
//        save_image_pgm("all_boundary", std::to_string(k), all_boundaries_image, max_col);

        if(boundary.size()>0)
        {
            pcl::PointCloud<pcl::PointXYZ> cloud_boundary = extractBoundCloud();
            all_boundaries += cloud_boundary;

            std::cout<<"Start computing theoretical connexions"<<std::endl;
            neighborPix2currentPix = man2D.getNeighborPix2currentPix();
            DefineIntersectionsPixPoints(k);
            computeLines(); // compute theoritical line corresponding to each pixel of the boundary + fill pixels corresponding to this intersection
            std::cout<<"Stop computing theoretical connexions"<<std::endl<<std::endl;

            for(int i = 0; i<intersections.size(); ++i)
            {
//                if(intersections[i].plane_ref->index == 14 && intersections[i].plane_neigh->index == 14)
//                {
//                    Eigen::MatrixXi image_test = Eigen::MatrixXi::Zero(Nrow, Ncol);
//                    for(auto it = intersections[i].pixels.begin(); it != intersections[i].pixels.end(); ++it)
//                        image_test(it->first, it->second) = 1;
//                    save_image_pgm("pixels_of_lines", "", image_test, 1);
//                    getchar();
//                    getchar();
//                    getchar();
//                }
                if(intersections[i].isLine)
                {
                    std::cout<<"intersection number : "<<i<<std::endl<<std::endl;
                    //correct obstruction lines
                    if(intersections[i].isObstruction && intersections[i].plane_ref->index != intersections[i].plane_neigh->index) // obstruction by other plane
                    {
                        intersection sister = intersections[i].export_sister();
                        intersections[i].correctObstructions(sister);
                        all_edges.push_back(intersections[i]);
                        std::cout<<"Start computing real lim of line n°"<<n<<" between planes "<<intersections[i].plane_ref->index<<" and "<<intersections[i].plane_neigh->index<<std::endl<<std::endl;
                        all_edges[n].computeLim();
                        all_edges[n].index = n;
                        // gather intersections by planes in "edges". index = plane . edges is used in corner detection
                        edges[all_edges[n].plane_ref->index].push_back(n);


                        all_edges.push_back(sister);
                        ++n;
                        std::cout<<"Start computing real lim of line n°"<<n<<" between planes "<<intersections[i].plane_ref->index<<" and "<<intersections[i].plane_neigh->index<<std::endl<<std::endl;
                        all_edges[n].computeLim();
                        all_edges[n].index = n;
                        // gather intersections by planes in "edges". index = plane . edges is used in corner detection
                        edges[all_edges[n].plane_ref->index].push_back(n);
                    }
                    else if(intersections[i].isObstruction && intersections[i].plane_ref->index == intersections[i].plane_neigh->index) // obstruction by object
                    {
                        intersection sister = intersections[i].export_sister();
                        intersections[i].correctObstructions(sister);
                        all_edges.push_back(intersections[i]);
                        std::cout<<"Start computing real lim of line n°"<<n<<" between plane "<<intersections[i].plane_ref->index<<" and other object"<<std::endl<<std::endl;
                        all_edges[n].computeLim();
                        all_edges[n].index = n;
                        edges[all_edges[n].plane_ref->index].push_back(n);
                    }
                    else
                    {
                        all_edges.push_back(intersections[i]);
                        //correct obstructions
                        std::cout<<"Start computing real lim of line n°"<<n<<" between planes "<<intersections[i].plane_ref->index<<" and "<<intersections[i].plane_neigh->index<<std::endl<<std::endl;
                        all_edges[n].computeLim();
                        all_edges[n].index = n;

                        // gather intersections by planes in "edges". index = plane . edges is used in corner detection
                        edges[all_edges[n].plane_ref->index].push_back(n);
                        if(all_edges[n].plane_ref->index != all_edges[n].plane_neigh->index)
                            edges[all_edges[n].plane_neigh->index].push_back(n);
                    }
                }
                else                        //keep it in other indices of "edges"
                {
                    edges[planes.size()].push_back(n);
                    all_edges.push_back(intersections[i]);
                }

                all_edges[n].plane_ref->intersections_indices.push_back(n);
                all_edges[n].plane_neigh->intersections_indices.push_back(n);

                ++n;
            }

            std::cout<<"-------------------------------------------------------------------"<<std::endl;
            std::cout<<"Number of edges : "<<intersections.size()<<std::endl;
            std::cout<<"Stop computing connexions"<<std::endl<<std::endl;

            //visualization------------------------------------------------------------------------------------------------------------------
            stm.str("");
            stm<<"cloud_boundary_"<<k<<".pcd";
            pcl::io::savePCDFileASCII (stm.str(), cloud_boundary);

            std::stringstream stm2("");
            stm2<<"bash pcd2csv.sh "<<stm.str();
            std::string command = stm2.str();
            system(command.c_str());

    //        save_image_pgm("binarized", std::to_string(k), binarized.cast<int>()*max_col, max_col);
    //        save_image_pgm("morpho", std::to_string(k), morpho.cast<int>()*max_col, max_col);
    //        Eigen::MatrixXi boundary_image = man2D.getBoundariesImage();
    //        save_image_pgm("boundary", std::to_string(k), boundary_image, max_col);
            //--------------------------------------------------------------------------------------------------------------------------------
        }
    }

    std::cout<<"Start computing theoritical corners"<<std::endl;
    possible_corners.resize(0);
    computeTheoriticalPlanesIntersections(); // for 3 planes intersection
    computeTheoriticalLinesIntersections(); // for 2 edges of one plane intersection
    std::cout<<"Stop computing theoritical corners"<<std::endl<<std::endl;

    float delta = 0.01;

    stm.str("");
    stm<<"Vertices.csv";
    std::ofstream file_vertices(stm.str());
    corners.resize(0);
    std::set<int> possible_corners_indices;
    //select corners among theoretical ones + add limits of line to possible_corners + draw lines
    for(int k = 0; k<all_edges.size(); ++k)
    {
        if(all_edges[k].plane_ref->index == 5)
        {
            std::ofstream f_theo_corners("theoretical_corners.csv");
            for (int i_theo_corners = 0; i_theo_corners < all_edges[k].theoreticalLim.size(); ++i_theo_corners)
                f_theo_corners<<all_edges[k].theoreticalLim[i_theo_corners].transpose()<<"\n";

            f_theo_corners.close();
            getchar();
            getchar();
        }

        if(all_edges[k].isLine)
        {
            std::cout<<"Intersection number : "<<k<<";   planes concerned : "<<all_edges[k].plane_ref->index<<" and "<<all_edges[k].plane_neigh->index <<std::endl;
            all_edges[k].selectCorners(incertainty_pixels, delta_phi, delta_theta, rot_axis, axis_init_phi);
            if(!all_edges[k].replaced_start && !all_edges[k].replaced_end)
            {
                corner c1;
                c1.planes.push_back(all_edges[k].plane_ref);
                c1.planes_indices.insert(all_edges[k].plane_ref->index);
                c1.planes.push_back(all_edges[k].plane_neigh);
                c1.planes_indices.insert(all_edges[k].plane_neigh->index);
                c1.pt = all_edges[k].start_pt;
                possible_corners.push_back(c1);
                //-----
                corner c2;
                c2.planes.push_back(all_edges[k].plane_ref);
                c2.planes_indices.insert(all_edges[k].plane_ref->index);
                c2.planes.push_back(all_edges[k].plane_neigh);
                c2.planes_indices.insert(all_edges[k].plane_neigh->index);
                c2.pt = all_edges[k].end_pt;
                possible_corners.push_back(c2);

                possible_corners_indices.insert(possible_corners.size()-2);
                possible_corners_indices.insert(possible_corners.size()-1);
            }
            else if(all_edges[k].replaced_start && !all_edges[k].replaced_end)
            {
                possible_corners_indices.insert(all_edges[k].corner_indices[0]);
                corner c;
                c.planes.push_back(all_edges[k].plane_ref);
                c.planes_indices.insert(all_edges[k].plane_ref->index);
                c.planes.push_back(all_edges[k].plane_neigh);
                c.planes_indices.insert(all_edges[k].plane_neigh->index);
                c.pt = all_edges[k].end_pt;
                possible_corners.push_back(c);
                possible_corners_indices.insert(possible_corners.size()-1);
            }
            else if(all_edges[k].replaced_end && !all_edges[k].replaced_start)
            {
                possible_corners_indices.insert(all_edges[k].corner_indices[0]);
                corner c;
                c.planes.push_back(all_edges[k].plane_ref);
                c.planes_indices.insert(all_edges[k].plane_ref->index);
                c.planes.push_back(all_edges[k].plane_neigh);
                c.planes_indices.insert(all_edges[k].plane_neigh->index);
                c.pt = all_edges[k].start_pt;
                possible_corners.push_back(c);
                possible_corners_indices.insert(possible_corners.size()-1);
            }
            else if(all_edges[k].replaced_start && all_edges[k].replaced_end)
            {
                possible_corners_indices.insert(all_edges[k].corner_indices[0]);
                possible_corners_indices.insert(all_edges[k].corner_indices[1]);
            }

            stm.str("");
            stm<<"theoritical_line_"<<k<<".csv";
            std::ofstream file(stm.str());
            int n = 0;
            while (n*delta+all_edges[k].start<all_edges[k].end)
            {
                file << ((n*delta+all_edges[k].start) * all_edges[k].tangente + all_edges[k].distance*all_edges[k].normal).transpose()<<"\n";
                ++n;
            }
            file.close();

            file.open("theoretical_lines.csv",  std::ofstream::app);
            n = 0;
            while (n*delta+all_edges[k].start<all_edges[k].end)
            {
                file << ((n*delta+all_edges[k].start) * all_edges[k].tangente + all_edges[k].distance*all_edges[k].normal).transpose()<<"\n";
                ++n;
            }
            file.close();

            stm.str("");
            stm<<"real_line_"<<k<<".csv";
            std::ofstream file_real_intersection(stm.str());
            for(int pt_idx = 0; pt_idx<all_edges[k].points.size(); ++pt_idx)
                file_real_intersection << all_edges[k].points[pt_idx](0)<<","<<all_edges[k].points[pt_idx](1)<<","<<all_edges[k].points[pt_idx](2)<<"\n";
            file_real_intersection.close();

            file.open("real_lines.csv",  std::ofstream::app);
            for(int pt_idx = 0; pt_idx<all_edges[k].points.size(); ++pt_idx)
                file << all_edges[k].points[pt_idx](0)<<","<<all_edges[k].points[pt_idx](1)<<","<<all_edges[k].points[pt_idx](2)<<"\n";
            file.close();

        }
        else
        {
            //save all_edges[k] which is not a line
            stm.str("");
            stm<<"not_line_"<<k<<".csv";
            std::ofstream file_not_line(stm.str());
            for(int pt_idx = 0; pt_idx<all_edges[k].points.size(); ++pt_idx)
                file_not_line << all_edges[k].points[pt_idx](0)<<","<<all_edges[k].points[pt_idx](1)<<","<<all_edges[k].points[pt_idx](2)<<"\n";
            file_not_line.close();

            if(all_edges[k].plane_ref->index != planes.size())
            {
                file_not_line.open("not_lines.csv",  std::ofstream::app);
                for(int pt_idx = 0; pt_idx<all_edges[k].points.size(); ++pt_idx)
                    file_not_line << all_edges[k].points[pt_idx](0)<<","<<all_edges[k].points[pt_idx](1)<<","<<all_edges[k].points[pt_idx](2)<<"\n";
                file_not_line.close();
            }
        }
    }

    //save corners in "corners" vector + draw corners

    for(auto it_ci = possible_corners_indices.begin(); it_ci != possible_corners_indices.end(); ++it_ci)
    {
        corners.push_back(possible_corners[*it_ci]);
        file_vertices << possible_corners[*it_ci].pt(0)<<","<<possible_corners[*it_ci].pt(1)<<","<<possible_corners[*it_ci].pt(2)<<"\n";
    }
    file_vertices.close();

    pcl::io::savePCDFileASCII ("results/all_boundaries.pcd", all_boundaries);
    system("bash pcd2csv.sh results/all_boundaries.pcd");
}

pcl::PointCloud<pcl::PointXYZ> manager::extractBoundCloud()
{
    std::cout<<"Start converting boundaries from image to 3D pc"<<std::endl;

    pcl::PointCloud<pcl::PointXYZ> cloud_bound;
    for (int i = 0; i < Nrow; ++i)
    {
        for (int j = 0; j < Ncol; ++j)
        {
            auto idx_found_boundary = boundary.find(std::make_pair(i,j));
            if(idx_found_boundary != boundary.end())
            {
                auto idx_found_XY2PtN = XY2PtN.find(std::make_pair(i,j));
                if (idx_found_XY2PtN != XY2PtN.end())
                {
                    Eigen::Vector3d pt= idx_found_XY2PtN->second.first;
                    pcl::PointXYZ pcl_pt;
                    pcl_pt.x = pt(0);
                    pcl_pt.y = pt(1);
                    pcl_pt.z = pt(2);
                    cloud_bound.points.push_back(pcl_pt);
                }
                else
                    std::cout<<"Error extracting bound cloud from boudary image : Can't find this boundary point in 3D pointcloud : "<< i<<" "<<j<<std::endl<<std::endl;
            }
        }
    }

    cloud_bound.width    = cloud_bound.points.size();
    cloud_bound.height   = 1;
    cloud_bound.is_dense = false;
    cloud_bound.points.resize (cloud_bound.width * cloud_bound.height);

    std::cout<<"Number of points on boundary : "<<cloud_bound.size()<<std::endl;

    std::cout<<"Stop converting boundaries from image to 3D pc"<<std::endl<<std::endl;
    return cloud_bound;
}

//--------------------------------------------------------------------------------------------------------------------

void manager::DefineIntersectionsPixPoints(int idx_plane)
{
    ++idx_plane;
    intersections.resize(0);
    int intersections_number = 0;
    int idx_plane_neighbor;
    std::map<int, int> idx_plane2idx_intersection;

    std::cout<<"size boundary before :"<<boundary.size()<<std::endl<<std::endl;

    //Remember : If i am stydying plane P and a pixel px belongs to P and has a neighbor pixel belonging to a neighbor plane Pn -> the index of px is the index of Pn

    //1_ process neighboring planes ( connection or obstruction)
    for(auto it_boundary = boundary.begin(); it_boundary != boundary.end();)
    {
        idx_plane_neighbor = it_boundary->second;

        if(idx_plane_neighbor != idx_plane && idx_plane_neighbor != planes.size()+1)
        {
            if(!processed_planes(idx_plane_neighbor-1, idx_plane-1)) // if an intersection with this neighbor plane has already been processed when sutdying the neighbor plane
            {
                //Compute theoritical intersection with first pixel encountered
                if(!processed_planes(idx_plane-1, idx_plane_neighbor-1)) // to not add new intersection if it has not already been processed with other pixel
                {
                    processed_planes(idx_plane-1, idx_plane_neighbor-1) = true; // to process only the first pixel encountered of this color (from this neighbor plane)
                    intersection inter (&planes[idx_plane-1], &planes[idx_plane_neighbor-1], delta_phi, delta_theta);
                    intersections.push_back(inter);
                    std::cout<<" intersections_number : "<<intersections_number<<"  is a relation with other plane"<<std::endl;
                    idx_plane2idx_intersection.insert(std::make_pair(idx_plane_neighbor-1, intersections_number)); //save which neighboring plane correspond to which intersection number
                    ++intersections_number;
                }

                auto idx_found_XY2PtN = XY2PtN.find(it_boundary->first);
                // put current boundary point into the attribute "pixels" of the specific intersection

                intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].pixels.push_back(it_boundary->first);
                intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].points.push_back(idx_found_XY2PtN->second.first);

                if(image_clusterized_indices(it_boundary->first.first, it_boundary->first.second) == idx_plane)
                    intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].indices_self.insert(intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].pixels.size()-1);
                else
                    intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].indices_sister.insert(intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].pixels.size()-1);
            }

            //erase neighboring points of the same studied plane (opening) or of "other objects" to be sure not to detect false openings/obstruction by object
            int rad = 2;
            int min_ki = std::max(it_boundary->first.first-rad, 0);
            int max_ki = std::min(it_boundary->first.first+rad, Nrow-1);
            int min_kj = std::max(it_boundary->first.second-rad, 0);
            int max_kj = std::min(it_boundary->first.second+rad, Ncol-1);

            for(int ki = min_ki; ki<=max_ki; ++ki)
            {
                for(int kj = min_kj; kj<=max_kj; ++kj)
                {
                    auto it_found_boundary = boundary.find(std::make_pair(ki, kj));
                    if(it_found_boundary != boundary.end())
                    {
                        if(it_found_boundary->second == idx_plane || it_found_boundary->second-1 == planes.size()) //opening || other object
                            boundary.erase(it_found_boundary);
                    }
                }
            }

            //erase current studied point
            it_boundary = boundary.erase(it_boundary);
        }
        else
            ++it_boundary;
    }

    std::cout<<"size boundary after connection pixels :"<<boundary.size()<<std::endl<<std::endl;

    //2_ process obstruction by object
    for(auto it_boundary = boundary.begin(); it_boundary != boundary.end();)
    {
        idx_plane_neighbor = it_boundary->second;

        if(idx_plane_neighbor-1 == planes.size()) // if there is an obstructing OBJECT
        {
            if(!processed_planes(idx_plane_neighbor-1, idx_plane-1)) // to not process intersection if the neighbor plane has already been processed
            {
                //Compute theoritical intersection with first pixel encountered
                if(!processed_planes(idx_plane-1, idx_plane_neighbor-1)) // to not add new intersection if it has already been processed with other pixel
                {
                    processed_planes(idx_plane-1, idx_plane_neighbor-1) = true; // to process only the first pixel encountered of this color (from this neighbor plane)
                    intersection inter (&planes[idx_plane-1], &planes[idx_plane-1], delta_phi, delta_theta);
                    inter.isObstruction = true;
                    intersections.push_back(inter);
                    idx_plane2idx_intersection.insert(std::make_pair(idx_plane_neighbor-1, intersections_number)); //save which neighboring plane correspond to which intersection number
                    std::cout<<" intersections_number : "<<intersections_number<<"  is an obstruction by object"<<std::endl;
                    ++intersections_number;
                }

                auto idx_found_XY2PtN = XY2PtN.find(it_boundary->first);

                intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].other_pixels.push_back(it_boundary->first);
                intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].other_points.push_back(idx_found_XY2PtN->second.first);

                if(image_clusterized_indices(it_boundary->first.first, it_boundary->first.second) == idx_plane)
                    intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].indices_self.insert(intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].other_pixels.size()-1);
                else
                    intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].indices_sister.insert(intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].other_pixels.size()-1);
            }

            //erase neighboring points of the same studied plane to be sure not to detect false openings
            int rad = 2;
            int min_ki = std::max(it_boundary->first.first-rad, 0);
            int max_ki = std::min(it_boundary->first.first+rad, Nrow-1);
            int min_kj = std::max(it_boundary->first.second-rad, 0);
            int max_kj = std::min(it_boundary->first.second+rad, Ncol-1);

            for(int ki = min_ki; ki<=max_ki; ++ki)
            {
                for(int kj = min_kj; kj<=max_kj; ++kj)
                {
                    auto it_found_boundary = boundary.find(std::make_pair(ki, kj));
                    if(it_found_boundary != boundary.end())
                    {
                        if(it_found_boundary->second == idx_plane)
                            boundary.erase(it_found_boundary);
                    }
                }
            }

            //erase current studied point
            it_boundary = boundary.erase(it_boundary);
        }
        else
            ++it_boundary;
    }

    std::cout<<"size boundary after obstruction by object pixels :"<<boundary.size()<<std::endl<<std::endl;

    //3_ process remaining = openings
    bool first_time = true;
    for(auto it_boundary = boundary.begin(); it_boundary != boundary.end();)
    {
        idx_plane_neighbor = it_boundary->second;

        if(idx_plane_neighbor == idx_plane) // if it belongs to an opening
        {
            if(!processed_planes(idx_plane_neighbor-1, idx_plane-1)) // to not process intersection if the neighbor plane has already been processed
            {
                //Compute theoritical intersection with first pixel encountered
                if(first_time) // to not process the pixels if the intersection has already been processed with other pixel
                {
                    intersection inter (&planes[idx_plane-1], &planes[idx_plane_neighbor-1], delta_phi, delta_theta);
                    inter.isOpening = true;
                    intersections.push_back(inter);
                    idx_plane2idx_intersection.insert(std::make_pair(idx_plane_neighbor-1, intersections_number)); //save which neighboring plane correspond to which intersection number
                    std::cout<<" intersections_number : "<<intersections_number<<"  is an opening"<<std::endl;
                    ++intersections_number;
                    first_time = false;
                }

                // put current boundary point into the attribute "pixels" of the specific intersection
                intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].other_pixels.push_back(it_boundary->first);
                // put current boundary point into the attribute "points" of the specific intersection
                auto idx_found_XY2PtN = XY2PtN.find(it_boundary->first);
                if(idx_found_XY2PtN!=XY2PtN.end())
                    intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].other_points.push_back(idx_found_XY2PtN->second.first);
                else
                    std::cout<<"Error creating intersection from opening pixels : did not find this boundary point in XY2PtN : "<< it_boundary->first.first <<" "<<it_boundary->first.second<<std::endl<<std::endl;
            }
            it_boundary = boundary.erase(it_boundary);
        }
        else
            ++it_boundary;
    }

    std::cout<<"size boundary after opening pixels :"<<boundary.size()<<std::endl<<std::endl;

    //--------------------------------------------------------------------------------------------------------------------------------------
}


void manager::computeLines()
{
    //Define the type of connection and compute limits of lines
    for (int k = 0; k<intersections.size(); ++k)
    {
        std::cout<<"-------------------------------------------------------"<<std::endl<<std::endl;
        std::cout<<"Intersection number : "<<k<<std::endl<<std::endl;
        std::cout<<"plane_ref : "<<intersections[k].plane_ref->index<<"  plane_neigh : "<<intersections[k].plane_neigh->index<<std::endl<<std::endl;
        std::cout<<"Number of pixels of this intersection : "<<intersections[k].pixels.size()<<std::endl;
        std::cout<<"Number of other_pixels of this intersection : "<<intersections[k].other_pixels.size()<<std::endl<<std::endl;

        //-------------------------------------------------------------------------------------------------------------------------------------------------
        // if opening
        if(intersections[k].isOpening)
        {
            std::cout<<"It is an opening"<<std::endl<<std::endl;
            int N_other_pixels_before;
            int N_other_pixels_after;

            if(intersections[k].other_pixels.size()>min_number_pixels)
            {
                N_other_pixels_before = intersections[k].other_pixels.size();
                intersections[k].computeTheoriticalLineFeaturesOpening(); //RANSAC to detect the main line of the intersection
                N_other_pixels_after = intersections[k].other_pixels.size();

                if(intersections[k].pixels.size()>min_number_pixels) // if the line found contains sufficient pixels
                {
                    std::cout<<"length : "<<intersections[k].length<<std::endl;
                    std::cout<<"number of pixels : "<<intersections[k].pixels.size()<<std::endl;
                    std::cout<<"number of other_pixels : "<<intersections[k].other_pixels.size()<<std::endl;
                    std::cout<<"normal : "<<intersections[k].normal.transpose()<<"   tangente : "<<intersections[k].tangente.transpose()<<"  distance : "<<intersections[k].distance<<"  initial point : "<<intersections[k].pt.transpose()<<std::endl<<std::endl;
                    intersections[k].isLine = true;

                    //create new intersections to look for more lines in remaining points
                    if(intersections[k].other_pixels.size()>min_number_pixels && (N_other_pixels_before-N_other_pixels_after)>0)
                    {
                        intersection inter(intersections[k].plane_ref, intersections[k].plane_neigh, delta_phi, delta_theta);
                        inter.other_points = intersections[k].other_points;
                        inter.other_pixels = intersections[k].other_pixels;
                        inter.isObstruction = false;
                        inter.repeated = intersections[k].repeated;
                        inter.isOpening  = true;
                        inter.isLine = false;
                        intersections.push_back(inter);
                        std::cout<<"The baby is : "<<intersections.size()-1<<std::endl<<std::endl;
                    }
                    else
                    {
                        //make reappear points not used in all_boundaries_image for other walls to use it
                        //put remaining points to inter_remaining of this plane
                        for(auto it_not_repeated = intersections[k].not_repeated.begin(); it_not_repeated != intersections[k].not_repeated.end(); ++it_not_repeated)
                        {
                            std::pair<int,int> pix = intersections[k].other_pixels[*it_not_repeated];
                            all_boundaries_image(pix.first, pix.second) = 0;
                            inter_remaining.pixels.push_back(pix);
                            Eigen::Vector3d point = intersections[k].other_points[*it_not_repeated];
                            inter_remaining.points.push_back(point);
                        }
                    }
                }
                else
                {
                    std::cout<<"intersection n°"<<k<<" could not be modelized by a line, the best line found only contains = "<<intersections[k].points.size()<<" points"<<std::endl;
                    if(intersections[k].pixels.size() + intersections[k].not_repeated.size()>min_number_pixels) //if there remain enough points we consider this intersection is not modelisable by a line but it still exists
                    {
                        std::cout<<"However, it has many points so it must not be a line but an obstruction not geometric."<<std::endl;
                        //make reappear points not used in all_boundaries_image for other walls to use it
                        //put remaining points to inter_remaining of this plane
                        for(auto it_not_repeated = intersections[k].not_repeated.begin(); it_not_repeated != intersections[k].not_repeated.end(); ++it_not_repeated)
                        {
                            std::pair<int,int> pix = intersections[k].other_pixels[*it_not_repeated];
                            intersections[k].pixels.push_back(pix);
                            all_boundaries_image(pix.first, pix.second) = 0;
                            Eigen::Vector3d point = intersections[k].other_points[*it_not_repeated];
                            intersections[k].points.push_back(point);
                        }
                        intersections[k].isLine = false;
                    }
                    else // if there does not remain enough points we consider they are noise around the other intersection and we remove them
                        remove_intersection(&k);
                }
            }
            else
                remove_intersection(&k);
        }
        else if(intersections[k].isObstruction) //if detected obstruction (obstruction baby)
        {
            std::cout<<"It is an obstruction baby between plane ref : "<<intersections[k].plane_ref->index<<" and plane neigh : "<<intersections[k].plane_neigh->index<<std::endl<<std::endl;

            int N_other_pixels_before;
            int N_other_pixels_after;

            if(intersections[k].other_pixels.size()>min_number_pixels)
            {
                N_other_pixels_before = intersections[k].other_pixels.size();
                intersections[k].computeTheoriticalLineFeaturesObstruction(neighborPix2currentPix); //RANSAC to detect the main line of the intersection
                N_other_pixels_after = intersections[k].other_pixels.size();
            }

            if(intersections[k].pixels.size()>min_number_pixels) // if the line found contains sufficient pixels
            {
                std::cout<<"length : "<<intersections[k].length<<std::endl;
                std::cout<<"number of pixels : "<<intersections[k].pixels.size()<<std::endl;
                std::cout<<"number of other_pixels : "<<intersections[k].other_pixels.size()<<std::endl;
                std::cout<<"normal : "<<intersections[k].normal.transpose()<<"   tangente : "<<intersections[k].tangente.transpose()<<"  distance : "<<intersections[k].distance<<"  initial point : "<<intersections[k].pt.transpose()<<std::endl<<std::endl;
                processed_planes(intersections[k].plane_neigh->index, intersections[k].plane_ref->index) = true;
                intersections[k].isLine = true;

                //create new intersections to look for more lines in remaining points
                if(intersections[k].other_pixels.size()>min_number_pixels && (N_other_pixels_before-N_other_pixels_after)>0)
                {
                    intersection inter(intersections[k].plane_ref, intersections[k].plane_neigh, delta_phi, delta_theta);
                    inter.other_points = intersections[k].other_points;
                    inter.other_pixels = intersections[k].other_pixels;
                    inter.isObstruction = intersections[k].isObstruction;
                    inter.repeated = intersections[k].repeated;
                    inter.indices_self = intersections[k].indices_self;
                    inter.indices_sister = intersections[k].indices_sister;
                    inter.isOpening  = intersections[k].isOpening;
                    inter.isLine = false;
                    intersections.push_back(inter);
                    std::cout<<"The baby is : "<<intersections.size()-1<<std::endl<<std::endl;
                    std::cout<<"It contains "<<intersections[k].repeated.size()<<" repeated points"<<std::endl<<std::endl;
                }
                else
                {
                    //make reappear points not used in all_boundaries_image for other walls to use it
                    //put remaining points to inter_remaining of this plane
                    for(auto it_not_repeated = intersections[k].not_repeated.begin(); it_not_repeated != intersections[k].not_repeated.end(); ++it_not_repeated)
                    {
                        std::pair<int,int> pix = intersections[k].other_pixels[*it_not_repeated];
                        all_boundaries_image(pix.first, pix.second) = 0;
                        inter_remaining.pixels.push_back(pix);
                        Eigen::Vector3d point = intersections[k].other_points[*it_not_repeated];
                        inter_remaining.points.push_back(point);
                    }
                }
            }
            else
            {
                std::cout<<"intersection n°"<<k<<" does not contain enough points onto the plane studied"<<std::endl;
                if(intersections[k].pixels.size() + intersections[k].not_repeated.size()>min_number_pixels) //if there remain enough points we consider this intersection is not modelisable by a line but it still exists
                {
                    std::cout<<"intersection n°"<<k<<" can not be modelize by a line because the best line found only contains = "<<intersections[k].points.size()<<" points"<<std::endl;
                    //make reappear points not used in all_boundaries_image for other walls to use it
                    //put remaining points to inter_remaining of this plane
                    for(auto it_not_repeated = intersections[k].not_repeated.begin(); it_not_repeated != intersections[k].not_repeated.end(); ++it_not_repeated)
                    {
                        std::pair<int,int> pix = intersections[k].other_pixels[*it_not_repeated];
                        intersections[k].pixels.push_back(pix);
                        all_boundaries_image(pix.first, pix.second) = 0;
                        Eigen::Vector3d point = intersections[k].other_points[*it_not_repeated];
                        intersections[k].points.push_back(point);
                    }
                    intersections[k].isLine = false;
                }
                else // if there does not remain enough points we consider they are noise around the other intersection and we remove them
                    remove_intersection(&k);
            }
        }
        //-----------------------------------------------------------------------------------------------------------------------------------
        else             //if connection or not detected obstruction
        {
            std::cout<<"this intersection is a connection or an obstruction which has not been detected yet between planes n°"<<intersections[k].plane_ref->index<<" and plane n°"<<intersections[k].plane_neigh->index<<std::endl;
            if(intersections[k].pixels.size()>min_number_pixels)
            {
                //affichage des pixels neighbors qui ont une correspondance dans ref
//                if(intersections[k].plane_ref->index == 3 && intersections[k].plane_neigh->index == 5)
//                {
//                    Eigen::MatrixXi test_image = Eigen::MatrixXi::Zero(Ncol, Nrow);
//                    for(auto it_neigh_pixel = intersections[k].pixels.begin(); it_neigh_pixel != intersections[k].pixels.end(); ++it_neigh_pixel)
//                    {
//                        auto pair_of_it_neighborPix2currentPix = neighborPix2currentPix.equal_range(*it_neigh_pixel);
//                        if( pair_of_it_neighborPix2currentPix.first != pair_of_it_neighborPix2currentPix.second)
//                            test_image(it_neigh_pixel->first, it_neigh_pixel->second) = 1;
//                    }
//                    save_image_pgm("test_neighborPix2currentPix", "", test_image, 1);
//                    getchar();
//                    getchar();
//                }

                intersections[k].computeTheoriticalLineFeaturesConnection(); // does not touch points/pixels other_points/other_pixels pt_mean -> just modify normal/tangente
                std::cout<<"normal : "<<intersections[k].normal.transpose()<<"   tangente : "<<intersections[k].tangente.transpose()<<"  distance : "<<intersections[k].distance<<"  initial point : "<<intersections[k].pt.transpose()<<std::endl<<std::endl;
                processed_planes(intersections[k].plane_neigh->index, intersections[k].plane_ref->index) = true;

                intersections[k].definePlaneConnection();               // define if the intersection is a real connection between planes or an obstruction (obstructed or obstructing) and if it is a line
                                                                        //move all points which do not belong to a connection from points to other_points

                //create new intersections with other_points to search obstructions

                if(intersections[k].other_pixels.size()>min_number_pixels)
                {
                    std::cout<<"there is an obstruction baby"<<std::endl;
                    intersection inter (intersections[k].plane_ref, intersections[k].plane_neigh, delta_phi, delta_theta);
                    inter.other_points = intersections[k].other_points;
                    inter.other_pixels = intersections[k].other_pixels;
                    inter.indices_self = intersections[k].indices_self;
                    inter.indices_sister = intersections[k].indices_sister;
                    inter.isObstruction = true;
                    intersections.push_back(inter);

                    std::cout<<"The hypothetical baby obstruction is : "<<intersections.size()-1<<" with "<<intersections[intersections.size()-1].other_points.size()<< " points"<<std::endl<<std::endl;

                }
                else
                {
                    std::cout<<"there is not an obstruction baby"<<std::endl;
                    for(auto it_pixels = intersections[k].other_pixels.begin(); it_pixels != intersections[k].other_pixels.end(); ++it_pixels)
                    {
                        all_boundaries_image(it_pixels->first, it_pixels->second) = 0;
                        inter_remaining.pixels.push_back(*it_pixels);
                    }
                    for(auto it_points = intersections[k].other_points.begin(); it_points != intersections[k].other_points.end(); ++it_points)
                        inter_remaining.points.push_back(*it_points);
                }

                if(intersections[k].isConnection)
                {
                    std::cout<<"Conclusion : it is a connection with plane n°"<<intersections[k].plane_neigh->index<<std::endl<<std::endl;
                    intersections[k].isLine = true;
                    //remaining points may be obstructions
                }
                else
                {
                    std::cout<<"Conclusion : it is not a connection with plane n°"<<intersections[k].plane_neigh->index<<std::endl<<std::endl;
                    intersections.erase(intersections.begin()+k); // real erase because points are all splitted into ref_ref and neigh_neigh
                    --k;
                }
            }
        }
    }

//    inter_remaining.pixels.clear();
//    inter_remaining.points.clear();
//    inter_remaining.setPlanes(&planes[idx_plane-1], &planes[idx_plane-1]);
//    //3_ process remaining points (details + obstructing objects) and put it in a last intersection
//    for(auto it_boundary = boundary.begin(); it_boundary != boundary.end(); ++it_boundary)
//    {
//        inter_remaining.pixels.push_back(it_boundary->first);
//        // put current boundary point into the attribute "points" of the specific intersection
//        auto idx_found_XY2PtN = XY2PtN.find(it_boundary->first);
//        if(idx_found_XY2PtN!=XY2PtN.end())
//            inter_remaining.points.push_back(idx_found_XY2PtN->second.first);
//    }
//    intersections.push_back(inter_remaining);
}

void manager::remove_intersection(int* k)
{
    //make reappear points not used in all_boundaries_image for other walls to use it
    //put remaining points to inter_remaining of this plane
    for(auto it_not_repeated = intersections[*k].not_repeated.begin(); it_not_repeated != intersections[*k].not_repeated.end(); ++it_not_repeated)
    {
        std::pair<int,int> pix = intersections[*k].other_pixels[*it_not_repeated];
        all_boundaries_image(pix.first, pix.second) = 0;
        inter_remaining.pixels.push_back(pix);
        Eigen::Vector3d point = intersections[*k].other_points[*it_not_repeated];
        inter_remaining.points.push_back(point);
    }
    for(auto it_pixels = intersections[*k].pixels.begin(); it_pixels != intersections[*k].pixels.end(); ++it_pixels)
    {
        all_boundaries_image(it_pixels->first, it_pixels->second) = 0;
        inter_remaining.pixels.push_back(*it_pixels);
    }
    for(int i = 0; i<intersections[*k].points.size(); ++i)
        inter_remaining.points.push_back(intersections[*k].points[i]);

    intersections.erase(intersections.begin() + *k);
    --(*k);
}

void manager::computeTheoriticalLinesIntersections()
{
    for(int plane_idx = 0; plane_idx<planes.size(); ++plane_idx)
    {
        if(edges[plane_idx].size()>=2)
        {
            for(int edge_idx = 0; edge_idx<edges[plane_idx].size()-1; ++edge_idx)
            {
                for(int edge_idx_other = edge_idx+1; edge_idx_other<edges[plane_idx].size(); ++edge_idx_other)
                {
                    if(!all_edges[edges[plane_idx][edge_idx]].isConnection || !all_edges[edges[plane_idx][edge_idx_other]].isConnection)
                    {
                        if(abs(all_edges[edges[plane_idx][edge_idx]].tangente.dot(all_edges[edges[plane_idx][edge_idx_other]].tangente)) < parallel_dot_threshold)
                        {
                            corner c;
                            c.setLines(&all_edges[edges[plane_idx][edge_idx]], &all_edges[edges[plane_idx][edge_idx_other]]);
                            c.computePointFromLines(plane_idx); // compute intersection point between lines on studied plane
                            c.isPlaneIntersection = false;
                            c.isLineIntersection = true;
                            possible_corners.push_back(c);

                            all_edges[edges[plane_idx][edge_idx]].possible_corner_indices.push_back(possible_corners.size()-1);
                            all_edges[edges[plane_idx][edge_idx]].theoreticalLim.push_back(c.pt);
                            all_edges[edges[plane_idx][edge_idx_other]].possible_corner_indices.push_back(possible_corners.size()-1);
                            all_edges[edges[plane_idx][edge_idx_other]].theoreticalLim.push_back(c.pt);
                        }
                    }
                }
            }
        }
    }
}

void manager::computeTheoriticalPlanesIntersections()
{
    std::vector<std::vector<std::vector<bool>>> already_treated(planes.size(), std::vector<std::vector<bool>>(planes.size(), std::vector<bool>(planes.size(), false)));
    for(int plane_idx = 0; plane_idx<planes.size(); ++plane_idx)
    {
        if(edges[plane_idx].size()>=2)
        {
            for(int edge_idx = 0; edge_idx<edges[plane_idx].size()-1; ++edge_idx)
            {
                int i = edges[plane_idx][edge_idx];
                if(all_edges[i].isConnection)
                {
                    for(int edge_idx_other = edge_idx+1; edge_idx_other<edges[plane_idx].size(); ++edge_idx_other)
                    {
                        int j = edges[plane_idx][edge_idx_other];
                        if(all_edges[j].isConnection)
                        {
                            plane* first_plane;
                            plane* second_plane;
                            plane* third_plane;

                            if(all_edges[i].plane_ref->index == plane_idx)
                            {
                                first_plane = all_edges[i].plane_ref;
                                second_plane = all_edges[i].plane_neigh;
                            }
                            else
                            {
                                first_plane = all_edges[i].plane_neigh;
                                second_plane = all_edges[i].plane_ref;
                            }
                            if(all_edges[j].plane_ref->index == plane_idx)
                                third_plane = all_edges[j].plane_neigh;
                            else
                                third_plane = all_edges[j].plane_ref;

                            int k;
                            int third_edge_idx;
                            for(third_edge_idx = 0; third_edge_idx<edges[third_plane->index].size(); ++third_edge_idx)
                            {
                                k = edges[third_plane->index][third_edge_idx];
                                if(all_edges[k].isConnection && (all_edges[k].plane_ref->index == second_plane->index || all_edges[k].plane_neigh->index == second_plane->index))
                                    break;
                            }

                            if(!already_treated[first_plane->index][second_plane->index][third_plane->index] && third_edge_idx!=edges[third_plane->index].size())
                            {
                                corner c;
                                c.setPlanes(first_plane, second_plane, third_plane);
                                c.isPlaneIntersection = true;
                                c.isLineIntersection = false;
                                c.computePointFromPlanes();
                                possible_corners.push_back(c);

                                all_edges[i].possible_corner_indices.push_back(possible_corners.size()-1);
                                all_edges[i].theoreticalLim.push_back(c.pt);
                                all_edges[j].possible_corner_indices.push_back(possible_corners.size()-1);
                                all_edges[j].theoreticalLim.push_back(c.pt);
                                all_edges[k].possible_corner_indices.push_back(possible_corners.size()-1);
                                all_edges[k].theoreticalLim.push_back(c.pt);
                                already_treated[first_plane->index][second_plane->index][third_plane->index] = true;
                                already_treated[first_plane->index][third_plane->index][second_plane->index] = true;
                                already_treated[second_plane->index][first_plane->index][third_plane->index] = true;
                                already_treated[second_plane->index][third_plane->index][first_plane->index] = true;
                                already_treated[third_plane->index][second_plane->index][first_plane->index] = true;
                                already_treated[third_plane->index][first_plane->index][second_plane->index] = true;
                            }
                        }
                    }
                }
            }
        }
    }

//    for(int i = 0; i<all_edges.size()-1; ++i)
//    {
//        if(all_edges[i].isConnection)
//        {
//            for(int j = i+1; j<all_edges.size(); ++j)
//            {
//                if(all_edges[j].isConnection)
//                {
//                    bool common_plane = all_edges[i].plane_ref->index == all_edges[j].plane_ref->index || all_edges[i].plane_neigh->index == all_edges[j].plane_ref->index || all_edges[i].plane_ref->index == all_edges[j].plane_neigh->index || all_edges[i].plane_neigh->index == all_edges[j].plane_neigh->index;
//                    if(common_plane)
//                    {
//                        corner c;
//                        if(all_edges[i].plane_ref->index == all_edges[j].plane_ref->index || all_edges[i].plane_neigh->index == all_edges[j].plane_ref->index)
//                            c.setPlanes(all_edges[i].plane_ref, all_edges[i].plane_neigh, all_edges[j].plane_neigh);
//                        else if(all_edges[i].plane_ref->index == all_edges[j].plane_neigh->index || all_edges[i].plane_neigh->index == all_edges[j].plane_neigh->index)
//                            c.setPlanes(all_edges[i].plane_ref, all_edges[i].plane_neigh, all_edges[j].plane_ref);

//                        c.isPlaneIntersection = true;
//                        c.isLineIntersection = false;
//                        c.computePointFromPlanes();
//                        possible_corners.push_back(c);
//                        all_edges[i].possible_corner_indices.push_back(possible_corners.size()-1);
//                        all_edges[i].theoreticalLim.push_back(c.pt);
//                        all_edges[j].possible_corner_indices.push_back(possible_corners.size()-1);
//                        all_edges[j].theoreticalLim.push_back(c.pt);
//                    }
//                }
//            }
//        }
//    }
}


void manager::detect_margin()
{
    if(lim_theta.first == -1)
    {
        int temp;
        int min_theta = 0;
        int max_theta = Ncol;
        for (auto it_XY2PtN = XY2PtN.begin(); it_XY2PtN != XY2PtN.end(); ++it_XY2PtN)
        {
            if(it_XY2PtN->first.second>min_theta)
                min_theta = it_XY2PtN->first.second;

            temp = it_XY2PtN->first.first ;
            while(it_XY2PtN->first.first == temp)
                ++it_XY2PtN;
            --it_XY2PtN;

            if(it_XY2PtN->first.second<max_theta)
                max_theta = it_XY2PtN->first.second;
        }

        for (auto it_XY2PtN = XY2PtN.begin(); it_XY2PtN != XY2PtN.end();)
        {
            if(it_XY2PtN->first.second<min_theta || it_XY2PtN->first.second>max_theta)
                it_XY2PtN = XY2PtN.erase(it_XY2PtN);
            else
                ++it_XY2PtN;
        }

        std::cout<<"min_theta : "<<min_theta<<"   max_theta : "<<max_theta<<std::endl<<std::endl;
        lim_theta = std::make_pair(min_theta, max_theta);
    }
}
