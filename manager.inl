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
    inter_remaining.isObject = true;
    lim_theta.first = -1;
    lim_theta.second = -1;
    rot_axis = ra;
    axis_init_phi = aip;
    seeds_pixels = Eigen::MatrixXi::Zero(Nrow, Ncol);
}

std::pair<int,int> pt2XY(Eigen::Vector3d pt, float delta_phi, float delta_theta, Eigen::Vector3d rot_axis, Eigen::Vector3d axis_init_phi)
{
    float theta = acos(pt.dot(rot_axis)/pt.norm());
    float phi = atan2((rot_axis.cross(axis_init_phi)).dot(pt), axis_init_phi.dot(pt));
    if (phi<0)
        phi += 2*M_PI;

    return std::make_pair((int)(phi/delta_phi), (int)(theta/delta_theta));
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

    image_clusterized_indices = Eigen::MatrixXi::Zero(Nrow, Ncol);
    for(int i = 0; i<3; ++i )
        image_clusterized[i] = Eigen::MatrixXi::Zero(Nrow, Ncol);

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
        int test = 0;
        for(auto it_pix = planes[i].pixels.begin(); it_pix != planes[i].pixels.end(); ++it_pix)
        {
            image_clusterized_indices(it_pix->first, it_pix->second) = planes[i].index+1;
            image_clusterized[0](it_pix->first, it_pix->second) = rgb(0);
            image_clusterized[1](it_pix->first, it_pix->second) = rgb(1);
            image_clusterized[2](it_pix->first, it_pix->second) = rgb(2);
            ++test;
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

    for(int i = 0; i< neigh.size(); ++i)
    {
        double dist = abs((it->second.first-neigh[i]->second.first).dot(it->second.second));
        if( dist > thresh_neigh_for_seed || abs(neigh[i]->second.second.dot(it->second.second))<normals_similarity_threshold_to_select_seed)
            return false;
    }

    //-------------------
    Eigen::Vector3d mean_pt = Eigen::Vector3d::Zero();
    for(int i = 0; i <neigh.size(); ++i)
        mean_pt += neigh[i]->second.first;

    mean_pt /= neigh.size();
    it->second.first = mean_pt;

    Eigen::Vector3d mean_norm = Eigen::Vector3d::Zero();
    for(int i = 0; i <neigh.size(); ++i)
        mean_norm += neigh[i]->second.second;

    mean_norm /= neigh.size();
    mean_norm /= mean_norm.norm();
    it->second.second= mean_norm;
    //-------------------
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
                seeds_pixels(it->first.first, it->first.second) = 1;
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
                                neighbors[k]->second.second /= neighbors[k]->second.second.norm();
                                it->second.second /= it->second.second.norm();
//                                double local_distance_seed = abs(it->second.first.dot(it->second.second));
//                                double local_distance_neigh = abs(neighbors[k]->second.first.dot(neighbors[k]->second.second));
//                                double dotty = abs(it->second.second.dot(neighbors[k]->second.second));
//                                if(dotty>1)
//                                    dotty = 1;
//                                if(abs(local_distance_seed-local_distance_neigh)<0.1 && acos(dotty)<45*M_PI/180)
                                if(abs((neighbors[k]->second.first-it->second.first).dot(it->second.second)) < thresh_plane_belonging)// && acos(neighbors[k]->second.dot(it->second)/(neighbors[k]->second.norm()*it->second.norm())) < thresh_angle )
                                {
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
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> share_points =  Eigen::MatrixXi::Zero(planes.size(),planes.size()).cast<bool>();
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

    for(int i = 0; i<planes.size(); ++i)
        planes[i].index = i;

    clusterized2Image();
    save_image_ppm("clusterized_before_gathering", "", image_clusterized, max_col);
    save_image_pgm("clusterized_indices_before_gathering", "", image_clusterized_indices, planes.size()+1);

    //------------------------------------------------------------------------------------------------------------------
    //gather planes with similar features

    std::cout<<"Start gathering superimposed planes "<<std::endl<<std::endl;
    std::map<int,int> to_erase_map; // key -> erased // value->completed (gathered)
    std::set<int> to_erase_idx;

    for(int i = 0; i<planes.size()-1; ++i)
    {
        for(int j = i+1; j<planes.size(); ++j)
        {
            Eigen::Vector3d mean_normal = (planes[i].normal + planes[j].normal);
            mean_normal /= mean_normal.norm();
            if(acos(abs(planes[i].normal.dot(planes[j].normal)))<10*M_PI/180 && abs((planes[i].mean_point_-planes[j].mean_point_).dot(mean_normal))<min_dist_planes) // if planes similar...
            {
                if(share_points(i,j)) //if planes share points
                {
                    auto it_i_in_to_erase = to_erase_map.find(i);
                    if( it_i_in_to_erase == to_erase_map.end()) // to check if the reference plane i has not been already gathered with other
                    {
                        for(int k = 0; k< planes[j].pts.size(); ++k)
                        {
                            planes[i].appendPoint(planes[j].pts[k]); //gather points
                            planes[i].appendPixel(planes[j].pixels[k]); // gather pixels
                        }
                        planes[i].computeNormal(); //recompute normal
                        to_erase_map.insert(std::make_pair(j,i)); //erase second plane
                        to_erase_idx.insert(j);
                    }
                    else
                    {
                        for(int k = 0; k< planes[j].pts.size(); ++k)
                        {
                            planes[it_i_in_to_erase->second].appendPoint(planes[j].pts[k]); //gather points
                            planes[it_i_in_to_erase->second].appendPixel(planes[j].pixels[k]); // gather pixels
                        }
                        planes[it_i_in_to_erase->second].computeNormal(); //recompute normal
                        to_erase_map.insert(std::make_pair(j,it_i_in_to_erase->second)); //erase second plane
                        to_erase_idx.insert(j);
                    }
                }
                else if (arePlanesClose(planes[i], planes[j]))
                {
                    auto it_i_in_to_erase = to_erase_map.find(i);
                    if( it_i_in_to_erase == to_erase_map.end()) // to check if the reference plane i has not been already gathered with other
                    {
                        for(int k = 0; k< planes[j].pts.size(); ++k)
                        {
                            planes[i].appendPoint(planes[j].pts[k]); //gather points
                            planes[i].appendPixel(planes[j].pixels[k]); // gather pixels
                        }
                        planes[i].computeNormal(); //recompute normal
                        to_erase_map.insert(std::make_pair(j,i)); //erase second plane
                        to_erase_idx.insert(j);
                    }
                    else
                    {
                        for(int k = 0; k< planes[j].pts.size(); ++k)
                        {
                            planes[it_i_in_to_erase->second].appendPoint(planes[j].pts[k]); //gather points
                            planes[it_i_in_to_erase->second].appendPixel(planes[j].pixels[k]); // gather pixels
                        }
                        planes[it_i_in_to_erase->second].computeNormal(); //recompute normal
                        to_erase_map.insert(std::make_pair(j,it_i_in_to_erase->second)); //erase second plane
                        to_erase_idx.insert(j);
                    }
                }
            }
        }
    }

    std::cout<<"Stop gathering superimposed planes "<<std::endl<<std::endl;

    //------------------------------------------------------------------------------------------------------------------
    //ultimate cleaning now I have corrected plane features
    auto it_to_erase = to_erase_idx.begin();
    std::vector<std::pair<int,int>> pixels_temp;
    std::vector<Eigen::Vector3d> points_temp;
    for(int i = 0; i<planes.size(); ++i)
    {
        points_temp.clear();
        pixels_temp.clear();
        if(i != *it_to_erase)
        {
            for (int k = 0; k<planes[i].pts.size(); ++k)
            {
                if(abs(abs(planes[i].pts[k].dot(planes[i].normal))-planes[i].distance) < max_plane_distance)
                {
                    points_temp.push_back(planes[i].pts[k]);
                    pixels_temp.push_back(planes[i].pixels[k]);
                }
            }
            planes[i].pts = points_temp;
            planes[i].pixels = pixels_temp;
        }
        else
        {
            ++it_to_erase;
            if(it_to_erase == to_erase_idx.end())
                --it_to_erase;
        }
    }

    for(int i = 0; i<planes.size(); ++i)
    {
        if(planes[i].pts.size() < min_number_of_pixels)
            to_erase_idx.insert(i);
    }

    //------------------------------------------------------------------------------------------------------------------
    auto it_to_erase_start = to_erase_idx.end();
    --it_to_erase_start;
    auto it_to_erase_end = to_erase_idx.begin();
    --it_to_erase_end;
    for(it_to_erase = it_to_erase_start; it_to_erase != it_to_erase_end; --it_to_erase)
    {
        planes.erase(planes.begin() + *it_to_erase);
        regions.erase(regions.begin() + *it_to_erase);
    }

    for(int i = 0; i<planes.size(); ++i)
        planes[i].index = i;

    processed_planes = Eigen::MatrixXi::Zero(planes.size()+1, planes.size()+1).cast<bool>();
    //------------------------------------------------------------------------------------------------------------------

    std::cout<<"------------------------------------------------------------------------"<<std::endl;
    std::cout<<"planes features after cleaning"<<std::endl<<std::endl;
    std::cout<<"------------------------------------------------------------------------"<<std::endl;
    for (int i = 0; i<planes.size(); ++i)
    {
        std::cout<<"Plane "<<i<<"\n";
        std::cout<<"normal : "<<planes[i].normal.transpose()<<std::endl;
        std::cout<<"pixel : "<<planes[i].pixels[0].first<<" "<<planes[i].pixels[0].second<<std::endl;
        std::cout<<"\n";
    }

}


bool manager::arePlanesClose(plane pi, plane pj)
{
    std::cout<<"Are planes close?"<<std::endl<<std::endl;
    for (int i = 0; i<pi.pixels.size(); ++i)
    {
        for (int j = 0; j < pj.pixels.size(); ++j)
        {
            if(sqrt(pow(pi.pixels[i].first-pj.pixels[j].first,2) + pow(pi.pixels[i].second-pj.pixels[j].second,2)) <= max_dist_between_pixels_in_line)
                return true;
        }
    }
    return false;
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
    std::cout<<"removing little objects (normals noise)"<<std::endl<<std::endl;
    man2D.removeLittleObjects(50);
    image_clusterized_indices = man2D.getImageClusterized();
    save_image_pgm("filtered_clusterized", "", image_clusterized_indices, planes.size()+1);
//    image_clusterized = pgm2ppm(image_clusterized_indices);
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
        save_image_pgm("binarized", std::to_string(k), binarized.cast<int>(), 1);

        //apply closure
        int rad = 2;
        man2D.morpho(rad);
        Eigen::MatrixXi morpho = man2D.getImBinarizedMorpho();
        save_image_pgm("morpho", std::to_string(k), morpho, 1);

        //compute boundary on image
        man2D.computeBoundaries(all_boundaries_image);
        boundary = man2D.getBoundaries();

        Eigen::MatrixXi boundary_image = man2D.getBoundariesImage();

        save_image_pgm("boundary", std::to_string(k), boundary_image, planes.size()+1);

        if(boundary.size()>0)
        {
            pcl::PointCloud<pcl::PointXYZ> cloud_boundary = extractBoundCloud();
            all_boundaries += cloud_boundary;

            std::cout<<"Start computing theoretical line connexions between planes"<<std::endl;
            neighborPix2currentPix = man2D.getNeighborPix2currentPix();
            DefineIntersectionsPixPoints(k);
            computeLines(); // compute theoritical line corresponding to each pixel of the boundary + fill pixels corresponding to this intersection
            std::cout<<"Stop computing theoretical connexions"<<std::endl<<std::endl;
            for(int i = 0; i<intersections.size(); ++i)
            {
                if(intersections[i].isLine)
                {
                    std::cout<<"intersection number : "<<i<<std::endl<<std::endl;
                    //correct obstruction lines
                    if(intersections[i].isObstruction) // obstruction by other plane
                    {
                        intersection sister = intersections[i].export_sister();
                        intersections[i].correctObstructions(sister);
                        bool ok = intersections[i].computeLim();
                        if(ok)
                            all_edges.push_back(intersections[i]);

                        ok = sister.computeLim();
                        if(ok)
                            all_edges.push_back(sister);
                    }
                    else if(intersections[i].isConnection || intersections[i].isOpening) // connection or opening -> nothing to refine
                    {
                        bool ok = intersections[i].computeLim();
                        if(ok)
                            all_edges.push_back(intersections[i]);
                    }
                }
                else                        //keep it in other indices of "edges"
                    all_edges.push_back(intersections[i]);
            }

            for(int i = 0; i<intersections.size(); ++i)
            {
                if(intersections[i].isLine)
                {
                    if(intersections[i].isObject) // interaction with object
                    {
                        std::cout<<"intersection number : "<<i<<std::endl<<std::endl;
                        if(intersections[i].has_sister)
                        {
                            intersection sister = intersections[i].export_sister();
                            intersections[i].correctObstructions(sister);
                        }

                        std::cout<<"Start computing visible lim of line n°"<<n<<" between plane "<<intersections[i].plane_ref->index<<" and other object"<<std::endl<<std::endl;
                        bool ok = intersections[i].computeLim();

                        if(ok)
                        {
                            if(!intersections[i].isDoubled(all_edges))
                                all_edges.push_back(intersections[i]);
                        }
                    }
                }
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
        }
    }

    pcl::io::savePCDFileASCII ("results/all_boundaries.pcd", all_boundaries);
    system("bash pcd2csv.sh results/all_boundaries.pcd");

    //fill edge vector for corner computation (will be done 3 times)
    fill_edges();
    std::cout<<"Start computing theoritical corners (intersection of 3 planes)"<<std::endl;
    computeTheoriticalPlanesIntersections(); // for 3 planes intersection
    std::cout<<"Stop computing theoritical corners (intersection of 3 planes)"<<std::endl<<std::endl;

    //--------------------------------//--------------------------------//--------------------------------
    std::cout<<"Removing objects found near corners : "<<std::endl<<std::endl;

    std::vector<Eigen::Vector3d> corners_pt (possible_corners.size());
    for(int corners_idx = 0; corners_idx < possible_corners.size(); ++corners_idx)
        corners_pt[corners_idx] = possible_corners[corners_idx].pt;

    std::vector<intersection> all_edges_temp;
    for(int k = 0; k<all_edges.size(); ++k)
    {
        if(all_edges[k].isObject)
        {
            if(!all_edges[k].isNearCorner(corners_pt))
                all_edges_temp.push_back(all_edges[k]);
        }
        else
            all_edges_temp.push_back(all_edges[k]);
    }
    all_edges = all_edges_temp;

    //refill edge vector for corner computation (2nd time)
    possible_corners.clear();
    for(int k = 0; k<all_edges.size(); ++k)
    {
        if(all_edges[k].isLine)
        {
            all_edges[k].index = k;
            all_edges[k].new_start_pt = all_edges[k].start_pt;
            all_edges[k].new_end_pt = all_edges[k].end_pt;
            all_edges[k].start = all_edges[k].start_pt.dot(all_edges[k].tangente);
            all_edges[k].end = all_edges[k].end_pt.dot(all_edges[k].tangente);
            all_edges[k].start_changed = false;
            all_edges[k].end_changed = false;
            all_edges[k].max_spatial_diff_start = max_spat_dist_for_line_connection;
            all_edges[k].max_spatial_diff_end = max_spat_dist_for_line_connection;
            all_edges[k].max_pixel_diff_start = max_pix_dist_for_line_connection;
            all_edges[k].max_pixel_diff_end = max_pix_dist_for_line_connection;
        }
    }

    fill_edges();

    std::cout<<"Start computing theoritical corners (intersection of 3 planes)"<<std::endl;
    computeTheoriticalPlanesIntersections(); // for 3 planes intersection
    std::cout<<"Stop computing theoritical corners (intersection of 3 planes)"<<std::endl<<std::endl;

    std::cout<<"reStart computing theoritical corners (intersections between lines of a same plane)"<<std::endl;
    computeTheoriticalLinesIntersections(); // for 2 edges of one plane intersection
    std::cout<<"reStop computing theoritical corners (intersections between lines of a same plane)"<<std::endl<<std::endl;

    //fusion corners
    for(int c = 0 ; c < possible_corners.size()-1; ++c)
    {
        if(possible_corners[c].isLineIntersection) // select first corner possibly fusionable
        {
            int idx_connection = 100000;
            bool no_connection_has_corner = false;
            bool no_connection_is_connected_twice = false;
            //check which connection line is studied
            if( possible_corners[c].lines[0]->isConnection && !possible_corners[c].lines[1]->isConnection )
            {
                idx_connection = possible_corners[c].lines[0]->index;
                no_connection_has_corner = (possible_corners[c].lines[1]->new_start_pt - possible_corners[c].pt).norm() < epsilon || (possible_corners[c].lines[1]->new_end_pt - possible_corners[c].pt).norm() < epsilon;
                no_connection_is_connected_twice = possible_corners[c].lines[1]->start_changed && possible_corners[c].lines[1]->end_changed;
            }
            else if ( !possible_corners[c].lines[0]->isConnection && possible_corners[c].lines[1]->isConnection )
            {
                no_connection_has_corner = (possible_corners[c].lines[0]->new_start_pt - possible_corners[c].pt).norm() < epsilon || (possible_corners[c].lines[0]->new_end_pt - possible_corners[c].pt).norm() < epsilon;
                idx_connection = possible_corners[c].lines[1]->index;
                no_connection_is_connected_twice = possible_corners[c].lines[0]->start_changed && possible_corners[c].lines[0]->end_changed;
            }

            if(idx_connection != 100000 && no_connection_has_corner && no_connection_is_connected_twice) // one of the lines of the corner is a connection
            {
                // check whether the connection is connected to other lines by start or by end
                bool start_of_connection = false;

                if((all_edges[idx_connection].new_start_pt - possible_corners[c].pt).norm() < (all_edges[idx_connection].new_end_pt - possible_corners[c].pt).norm())
                    start_of_connection = true;
                else
                    start_of_connection = false;

                for(int c1 = c+1 ; c1 < possible_corners.size(); ++c1)
                {
                    //check if there is other corner with same connection line involved and which is a lineintersection
                    if(possible_corners[c1].isLineIntersection && (possible_corners[c1].pt - possible_corners[c].pt).norm() < 0.06 && (possible_corners[c1].lines[0]->index == idx_connection  || possible_corners[c1].lines[1]->index == idx_connection))
                    {
                        bool no_connection1_has_corner = false;
                        bool no_connection1_is_connected_twice = false;
                        if (possible_corners[c1].lines[0]->index == idx_connection)
                        {
                            no_connection1_has_corner = (possible_corners[c1].lines[1]->new_start_pt - possible_corners[c1].pt).norm() < epsilon || (possible_corners[c1].lines[1]->new_end_pt - possible_corners[c1].pt).norm() < epsilon;
                            no_connection1_is_connected_twice = possible_corners[c1].lines[1]->start_changed && possible_corners[c1].lines[1]->end_changed;
                        }
                        else
                        {
                            no_connection1_has_corner = (possible_corners[c1].lines[0]->new_start_pt - possible_corners[c1].pt).norm() < epsilon || (possible_corners[c1].lines[0]->new_end_pt - possible_corners[c1].pt).norm() < epsilon;
                            no_connection1_is_connected_twice = possible_corners[c1].lines[0]->start_changed && possible_corners[c1].lines[0]->end_changed;
                        }
                        bool start_of_connection1 = false;
                        if((all_edges[idx_connection].new_start_pt - possible_corners[c1].pt).norm() < (all_edges[idx_connection].new_end_pt - possible_corners[c1].pt).norm())
                            start_of_connection1 = true;
                        else
                            start_of_connection1 = false;

                        //check that the corner point of c1 is on the same side of the connection line
                        if( (start_of_connection1 && start_of_connection ) || ( !start_of_connection1 && !start_of_connection ) && no_connection1_has_corner && no_connection1_is_connected_twice)
                        {
                            Eigen::Vector3d mean_pt = (possible_corners[c1].pt + possible_corners[c].pt)/2;

                            //actualize lines new_points
                            if( (possible_corners[c].lines[0]->new_start_pt - possible_corners[c].pt).norm()< epsilon )
                                possible_corners[c].lines[0]->new_start_pt = mean_pt;
                            else if( (possible_corners[c].lines[0]->new_end_pt - possible_corners[c].pt).norm()< epsilon )
                                possible_corners[c].lines[0]->new_end_pt = mean_pt;

                            if( (possible_corners[c].lines[1]->new_start_pt - possible_corners[c].pt).norm()< epsilon )
                                possible_corners[c].lines[1]->new_start_pt = mean_pt;
                            else if( (possible_corners[c].lines[1]->new_end_pt - possible_corners[c].pt).norm()< epsilon )
                                possible_corners[c].lines[1]->new_end_pt = mean_pt;

                            if( (possible_corners[c1].lines[0]->new_start_pt - possible_corners[c1].pt).norm()< epsilon )
                                possible_corners[c1].lines[0]->new_start_pt = mean_pt;
                            else if( (possible_corners[c1].lines[0]->new_end_pt - possible_corners[c1].pt).norm()< epsilon )
                                possible_corners[c1].lines[0]->new_end_pt = mean_pt;

                            if( (possible_corners[c1].lines[1]->new_start_pt - possible_corners[c1].pt).norm()< epsilon )
                                possible_corners[c1].lines[1]->new_start_pt = mean_pt;
                            else if( (possible_corners[c1].lines[1]->new_end_pt - possible_corners[c1].pt).norm()< epsilon )
                                possible_corners[c1].lines[1]->new_end_pt = mean_pt;

                            //actualize corner
                            if(possible_corners[c1].lines[0]->index != possible_corners[c].lines[0]->index && possible_corners[c1].lines[0]->index != possible_corners[c].lines[1]->index)
                                possible_corners[c].lines.push_back(possible_corners[c1].lines[0]);
                            else if(possible_corners[c1].lines[1]->index != possible_corners[c].lines[0]->index && possible_corners[c1].lines[1]->index != possible_corners[c].lines[1]->index)
                                possible_corners[c].lines.push_back(possible_corners[c1].lines[1]);
                            possible_corners[c].pt = mean_pt;
                        }
                    }
                }
            }
        }
    }



    //--------------------------------//--------------------------------//--------------------------------

    std::cout<<"Actualize start_changed and end_changed : "<<std::endl<<std::endl;
    actualizeChanged();

    //--------------------------------//--------------------------------//--------------------------------

    std::cout<<"reStart computing theoritical corners (two parallel lines connected)"<<std::endl;
    computeParallelLinesIntersections();
    std::cout<<"reStop computing theoritical corners (two parallel lines connected)"<<std::endl<<std::endl;

    //--------------------------------//--------------------------------//--------------------------------

    std::cout<<"Removing 1 or 2 segments connected with nothing : "<<std::endl<<std::endl;

    all_edges_temp.clear();
    std::vector<intersection> edges_removed;
    for(int k = 0; k<all_edges.size(); ++k)
    {
        std::cout<<"line : "<<k<<std::endl;
        if(all_edges[k].isLine)
        {
            std::cout<<"line : "<<k<<"  has start changed ? :"<<all_edges[k].start_changed<<"  has end changed ? :"<<all_edges[k].end_changed<<std::endl;
            bool shared_corner = false;
            if((all_edges[k].start_changed && all_edges[k].end_changed))        //start and end have changed
                all_edges_temp.push_back(all_edges[k]);
            else if(all_edges[k].start_changed || all_edges[k].end_changed)     //start or end has not changed and the other has
            {
                std::cout<<"line : "<<k<<std::endl;
                Eigen::Vector3d uniq_connexion;
                if(all_edges[k].start_changed && !all_edges[k].end_changed)
                    uniq_connexion = all_edges[k].new_start_pt;
                else if(!all_edges[k].start_changed && all_edges[k].end_changed)
                    uniq_connexion = all_edges[k].new_end_pt;

                for(int l = 0; l<all_edges.size(); ++l)
                {
                    if(all_edges[l].isLine && k != l && (all_edges[l].start_changed && all_edges[l].end_changed))
                    {
                        std::cout<<"lines : "<<k<<" "<<l<<std::endl;
                        if( (all_edges[l].new_start_pt-uniq_connexion).norm()<epsilon || (all_edges[l].new_end_pt-uniq_connexion).norm()<epsilon)
                            shared_corner = true;
                    }
                    if(all_edges[l].isLine && k != l)
                    {
                        if( (all_edges[l].new_start_pt-uniq_connexion).norm()<epsilon || (all_edges[l].new_end_pt-uniq_connexion).norm()<epsilon)
                            std::cout<<"There is a case : "<<k<<" "<<l<<std::endl<<std::endl;
                    }
                    if(shared_corner)
                        break;
                }
                if(shared_corner)
                    all_edges_temp.push_back(all_edges[k]);
                else
                    edges_removed.push_back(all_edges[k]);
            }
            else
                edges_removed.push_back(all_edges[k]);
        }
        else
            all_edges_temp.push_back(all_edges[k]);
    }
//    all_edges = all_edges_temp;

    //--------------------------------//--------------------------------//--------------------------------

    actualizeChanged();
    for(int k = 0; k<all_edges.size(); ++k)
    {
        all_edges[k].index = k;
        all_edges[k].plane_ref->intersections_indices.insert(k);
        if(all_edges[k].isConnection)
            all_edges[k].plane_neigh->intersections_indices.insert(k);
    }

    fill_edges();

//    std::cout<<" check lines extremities and numerical problems"<<std::endl;
//    for(int k = 0; k<edges.size(); ++k)
//    {
//        std::cout<<"plane n°"<<k<<std::endl;
//        for(auto it = edges[k].begin(); it != edges[k].end(); ++it)
//            std::cout<<all_edges[*it].new_start_pt.transpose()<<"       "<<all_edges[*it].new_end_pt.transpose()<<"\n";
//        std::cout<<"\n";
//    }
    for(int k = 0; k<edges.size(); ++k)
    {
        std::cout<<"Intersections of plane n°"<<k<<std::endl;
        for(auto it = edges[k].begin(); it != edges[k].end(); ++it)
            std::cout<<*it<<" ";
        std::cout<<"\n"<<"\n";
    }

    //remove unused corners and create corner vector
    for(int it_corners = 0; it_corners < possible_corners.size(); ++it_corners)
    {
        bool added = false;
        for(int k = 0; k<all_edges.size(); ++k)
        {
            if(all_edges[k].isLine)
            {
                if( ((all_edges[k].new_start_pt -possible_corners[it_corners].pt).norm() < epsilon && all_edges[k].start_changed) || ((all_edges[k].new_end_pt -possible_corners[it_corners].pt).norm() < epsilon && all_edges[k].end_changed))
                {
                    if(!added)
                    {
                        corners.push_back(possible_corners[it_corners]);
                        added = true;
                    }
                    all_edges[k].indices_corners.push_back(corners.size()-1);
                    corners[corners.size()-1].indices_lines.insert(k);
                }
            }
        }
    }

    //create indices_connected_lines_start and indices_connected_lines_end for each edge

    for(int it_corners = 0; it_corners < corners.size(); ++it_corners)
    {
        auto it_end = corners[it_corners].indices_lines.end();
        --it_end;
        for(auto it_indices_lines = corners[it_corners].indices_lines.begin(); it_indices_lines != it_end; ++it_indices_lines)
        {
            auto it_indices_lines_other = it_indices_lines;
            ++it_indices_lines_other;
            for(;it_indices_lines_other != corners[it_corners].indices_lines.end(); ++it_indices_lines_other)
            {
                if( (corners[it_corners].pt - all_edges[*it_indices_lines].new_start_pt).norm() < epsilon)
                    all_edges[*it_indices_lines].indices_connected_lines_start.insert(*it_indices_lines_other);
                else
                    all_edges[*it_indices_lines].indices_connected_lines_end.insert(*it_indices_lines_other);

                if( (corners[it_corners].pt - all_edges[*it_indices_lines_other].new_start_pt).norm() < epsilon)
                    all_edges[*it_indices_lines_other].indices_connected_lines_start.insert(*it_indices_lines);
                else
                    all_edges[*it_indices_lines_other].indices_connected_lines_end.insert(*it_indices_lines);
            }
        }
    }

    for(int it_corners = 0; it_corners < corners.size(); ++it_corners)
    {
        if(corners[it_corners].isLineIntersection)
        {
            std::cout<<"corner n° : "<<it_corners<<" -> isLineintersection"<<std::endl;
            std::cout<<"\t lines n° : ";
            for(auto it_indices_lines = corners[it_corners].indices_lines.begin(); it_indices_lines != corners[it_corners].indices_lines.end(); ++it_indices_lines)
                std::cout<< *it_indices_lines<<" ";
            std::cout<<std::endl;
        }
        if(corners[it_corners].isLineContinuity)
        {
            std::cout<<"corner n° : "<<it_corners<<" -> isLineContinuity"<<std::endl;
            std::cout<<"\t lines n° : ";
            for(auto it_indices_lines = corners[it_corners].indices_lines.begin(); it_indices_lines != corners[it_corners].indices_lines.end(); ++it_indices_lines)
                std::cout<< *it_indices_lines<<" ";
            std::cout<<std::endl;
        }
    }


    std::cout<<"Start drawing corners"<<std::endl<<std::endl;
    float delta = 0.01;
    std::ofstream file_vertices("Vertices.csv");

    for(int it_corners = 0; it_corners < corners.size(); ++it_corners)
    {
        Eigen::Vector3d color;
        if(corners[it_corners].isPlaneIntersection)
            color = {255,0,0};
        else if (corners[it_corners].isLineIntersection)
            color = {150,70,0};
        else if (corners[it_corners].isLineContinuity)
            color = {0,150,0};
        file_vertices << corners[it_corners].pt.transpose()<<" "<<color.transpose()<<"\n";
    }
    file_vertices.close();

    std::ofstream file_pot_vertices("Potential_vertices.csv");
    //Draw potential corners
    for(int it_corners = 0; it_corners < possible_corners.size(); ++it_corners)
        file_pot_vertices << possible_corners[it_corners].pt(0)<<","<<possible_corners[it_corners].pt(1)<<","<<possible_corners[it_corners].pt(2)<<"\n";
    file_pot_vertices.close();


    std::cout<<"Start drawing lines"<<std::endl<<std::endl;
    //draw lines
    for(int k = 0; k<all_edges.size(); ++k)
    {
        all_edges[k].start = all_edges[k].new_start_pt.dot(all_edges[k].tangente);
        all_edges[k].end = all_edges[k].new_end_pt.dot(all_edges[k].tangente);

        if(all_edges[k].isLine)
        {
            Eigen::Vector3i color;
            if(all_edges[k].isConnection)
                color = {255, 0, 0};
            else if(all_edges[k].isObstruction)
                color = {0, 255, 0};
            else if(all_edges[k].isObject)
                color = {0, 0, 255};
            else if(all_edges[k].isOpening)
                color = {0, 255, 255};
            std::cout<<"Drawing intersection number : "<<k<<";   planes concerned : "<<all_edges[k].plane_ref->index<<" and "<<all_edges[k].plane_neigh->index <<std::endl;
            stm.str("");
            stm<<"theoritical_line_"<<k<<".csv";
            std::ofstream file(stm.str());
            int n = 0;
            while (n*delta+all_edges[k].start<=all_edges[k].end)
            {
                file << ((n*delta+all_edges[k].start) * all_edges[k].tangente + all_edges[k].distance*all_edges[k].normal).transpose()<<" "<<color.transpose()<<"\n";
                ++n;
            }
            file.close();

            file.open("theoretical_lines.csv",  std::ofstream::app);
            n = 0;
            while (n*delta+all_edges[k].start<=all_edges[k].end)
            {
                file << ((n*delta+all_edges[k].start) * all_edges[k].tangente + all_edges[k].distance*all_edges[k].normal).transpose()<<" "<<color.transpose()<<"\n";
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
        else if(!all_edges[k].isLine)
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
    std::ofstream file_remaining("remaining.csv");
    for (int i = 0; i<inter_remaining.points.size(); ++i)
        file_remaining<< inter_remaining.points[i](0)<<","<<inter_remaining.points[i](1)<<","<<inter_remaining.points[i](2)<<"\n";
    file_remaining.close();

    //draw removed lines

    system("rm removed_lines.csv");
    for(int k = 0; k<edges_removed.size(); ++k)
    {
        edges_removed[k].start = edges_removed[k].new_start_pt.dot(edges_removed[k].tangente);
        edges_removed[k].end = edges_removed[k].new_end_pt.dot(edges_removed[k].tangente);
            Eigen::Vector3i color;
            if(edges_removed[k].isConnection)
                color = {255, 0, 0};
            else if(edges_removed[k].isObstruction)
                color = {0, 255, 0};
            else if(edges_removed[k].isObject)
                color = {0, 0, 255};
            else if(edges_removed[k].isOpening)
                color = {0, 255, 255};
            std::cout<<"Drawing intersection number : "<<k<<";   planes concerned : "<<edges_removed[k].plane_ref->index<<" and "<<edges_removed[k].plane_neigh->index <<std::endl;
            std::ofstream file;
            file.open("removed_lines.csv",  std::ofstream::app);
            n = 0;
            while (n*delta+edges_removed[k].start<=edges_removed[k].end)
            {
                file << ((n*delta+edges_removed[k].start) * edges_removed[k].tangente + edges_removed[k].distance*edges_removed[k].normal).transpose()<<" "<<color.transpose()<<"\n";
                ++n;
            }
            file.close();
    }

    export_mesh();
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
    int rad = 1;

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

        if(idx_plane_neighbor-1 == planes.size()) // if there is an OBJECT not planar
        {
            if(!processed_planes(idx_plane_neighbor-1, idx_plane-1)) // to not process intersection if the neighbor plane has already been processed
            {
                //Compute theoritical intersection with first pixel encountered
                if(!processed_planes(idx_plane-1, idx_plane_neighbor-1)) // to not add new intersection if it has already been processed with other pixel
                {
                    processed_planes(idx_plane-1, idx_plane_neighbor-1) = true; // to process only the first pixel encountered of this color (from this neighbor plane)
                    intersection inter (&planes[idx_plane-1], &planes[idx_plane-1], delta_phi, delta_theta);
                    inter.isObject = true;
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
        std::cout<<"points number : "<<intersections[k].points.size()<<std::endl<<std::endl;
        std::cout<<"plane_ref : "<<intersections[k].plane_ref->index<<"  plane_neigh : "<<intersections[k].plane_neigh->index<<std::endl<<std::endl;
        std::cout<<"Number of pixels of this intersection : "<<intersections[k].pixels.size()<<std::endl;
        std::cout<<"Number of other_pixels of this intersection : "<<intersections[k].other_pixels.size()<<std::endl<<std::endl;

        //-------------------------------------------------------------------------------------------------------------------------------------------------
        // if opening
        if(intersections[k].isOpening)
        {
            std::cout<<"It is an opening"<<std::endl<<std::endl;

            if(intersections[k].other_pixels.size()>=min_number_pixels)
            {
                intersections[k].computeTheoriticalLineFeaturesOpening(); //RANSAC to detect the main line of the intersection

                if(intersections[k].pixels.size()>=min_number_pixels) // if the line found contains sufficient pixels
                {
                    std::cout<<"length : "<<intersections[k].length<<std::endl;
                    std::cout<<"number of pixels : "<<intersections[k].pixels.size()<<std::endl;
                    std::cout<<"number of other_pixels : "<<intersections[k].other_pixels.size()<<std::endl;
                    std::cout<<"normal : "<<intersections[k].normal.transpose()<<"   tangente : "<<intersections[k].tangente.transpose()<<"  distance : "<<intersections[k].distance<<"  initial point : "<<intersections[k].pt_mean.transpose()<<std::endl<<std::endl;
                    intersections[k].isLine = true;

                    //create new intersections to look for more lines in remaining points
                    if(intersections[k].other_pixels.size()>=min_number_pixels)
                    {
                        intersection inter(intersections[k].plane_ref, intersections[k].plane_neigh, delta_phi, delta_theta);
                        inter.other_points = intersections[k].other_points;
                        inter.other_pixels = intersections[k].other_pixels;
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
                    for(int idx_remaining_points = 0; idx_remaining_points < intersections[k].remaining_points.size(); ++idx_remaining_points)
                    {
                        inter_remaining.pixels.push_back(intersections[k].remaining_pixels[idx_remaining_points]);
                        inter_remaining.points.push_back(intersections[k].remaining_points[idx_remaining_points]);
                    }
                }
                else
                {
                    std::cout<<"intersection n°"<<k<<" could not be modelized by a line, the best line found only contains = "<<intersections[k].points.size()<<" points"<<std::endl<<std::endl;
                    if(intersections[k].pixels.size() + intersections[k].not_repeated.size()>=min_number_pixels) //if there remain enough points we consider this intersection is not modelisable by a line but it still exists
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
        //-----------------------------------------------------------------------------------------------------------------------------------
        else if(intersections[k].isObstruction) //if detected obstruction (obstruction baby)
        {
            std::cout<<"It is an obstruction baby between plane ref : "<<intersections[k].plane_ref->index<<" and plane neigh : "<<intersections[k].plane_neigh->index<<std::endl<<std::endl;

            if(intersections[k].other_pixels.size()>=min_number_pixels)
            {
                intersections[k].computeTheoriticalLineFeaturesObstruction(neighborPix2currentPix); //RANSAC to detect the main line of the intersection

                if(intersections[k].pixels.size()>=min_number_pixels) // if the line found contains sufficient pixels
                {
                    std::cout<<"length : "<<intersections[k].length<<std::endl;
                    std::cout<<"number of pixels : "<<intersections[k].pixels.size()<<std::endl;
                    std::cout<<"number of other_pixels : "<<intersections[k].other_pixels.size()<<std::endl;
                    std::cout<<"normal : "<<intersections[k].normal.transpose()<<"   tangente : "<<intersections[k].tangente.transpose()<<"  distance : "<<intersections[k].distance<<"  initial point : "<<intersections[k].pt_mean.transpose()<<std::endl<<std::endl;
                    processed_planes(intersections[k].plane_neigh->index, intersections[k].plane_ref->index) = true; // for following plane not to treat it
                    intersections[k].isLine = true;

                    //create new intersections to look for more lines in remaining points
                    if(intersections[k].other_pixels.size()>=min_number_pixels)
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
                            inter_remaining.pixels.push_back(intersections[k].other_pixels[*it_not_repeated]);
                            inter_remaining.points.push_back(intersections[k].other_points[*it_not_repeated]);
                        }
                    }

                    for(int idx_remaining_points = 0; idx_remaining_points < intersections[k].remaining_points.size(); ++idx_remaining_points)
                    {
                        inter_remaining.pixels.push_back(intersections[k].remaining_pixels[idx_remaining_points]);
                        inter_remaining.points.push_back(intersections[k].remaining_points[idx_remaining_points]);
                    }
                }
                else
                {
                    std::cout<<"intersection n°"<<k<<" does not contain enough points onto the plane studied"<<std::endl;
                    if(intersections[k].pixels.size() + intersections[k].not_repeated.size()>=min_number_pixels) //if there remain enough points we consider this intersection is not modelisable by a line but it still exists
                    {
                        std::cout<<"intersection n°"<<k<<" can not be modelize by a line because the best line found only contains = "<<intersections[k].points.size()<<" points"<<std::endl;
                        //make reappear points not used in all_boundaries_image for other walls to use it
                        //put remaining points to inter_remaining of this plane
                        for(auto it_not_repeated = intersections[k].not_repeated.begin(); it_not_repeated != intersections[k].not_repeated.end(); ++it_not_repeated)
                        {
                            std::pair<int,int> pix = intersections[k].other_pixels[*it_not_repeated];
                            all_boundaries_image(pix.first, pix.second) = 0;
                            intersections[k].pixels.push_back(intersections[k].other_pixels[*it_not_repeated]);
                            intersections[k].points.push_back(intersections[k].other_points[*it_not_repeated]);
                        }
                        intersections[k].isLine = false;
                    }
                    else // if there does not remain enough points we consider they are noise around the other intersection and we remove them
                        remove_intersection(&k);
                }
            }
        }
        //-----------------------------------------------------------------------------------------------------------------------------------
        else if(intersections[k].isObject) //if object
        {
            std::cout<<"It is a connection or an obstruction by object from plane ref : "<<intersections[k].plane_ref->index<<std::endl<<std::endl;

            if(intersections[k].other_pixels.size()>=min_number_pixels) // if not enough points at the begining points are forgetted for ever
            {
                intersections[k].computeTheoriticalLineFeaturesObject(neighborPix2currentPix); //RANSAC to detect the main line of the intersectio
                if(intersections[k].pixels.size()>=min_number_pixels) // if the line found contains sufficient pixels
                {
                    std::cout<<"length : "<<intersections[k].length<<std::endl;
                    std::cout<<"number of pixels : "<<intersections[k].pixels.size()<<std::endl;
                    std::cout<<"number of other_pixels : "<<intersections[k].other_pixels.size()<<std::endl;
                    std::cout<<"normal : "<<intersections[k].normal.transpose()<<"   tangente : "<<intersections[k].tangente.transpose()<<"  distance : "<<intersections[k].distance<<"  initial point : "<<intersections[k].pt_mean.transpose()<<std::endl<<std::endl;
                    intersections[k].isLine = true;

                    //create new intersections to look for more lines in remaining points
                    if(intersections[k].other_pixels.size()>=min_number_pixels)
                    {
                        intersection inter(intersections[k].plane_ref, intersections[k].plane_neigh, delta_phi, delta_theta);
                        inter.other_points = intersections[k].other_points;
                        inter.other_pixels = intersections[k].other_pixels;
                        inter.isObject = true;
                        inter.repeated = intersections[k].repeated;

                        inter.indices_self = intersections[k].indices_self;
                        inter.indices_sister = intersections[k].indices_sister;
                        inter.isOpening  = false;
                        inter.isLine = false;
                        intersections.push_back(inter);
                        std::cout<<"The baby is : "<<intersections.size()-1<<std::endl<<std::endl;
                        std::cout<<"It contains "<<intersections[k].repeated.size()<<" repeated points"<<std::endl<<std::endl;
                    }
                    else
                    {
                        if(intersections[k].other_pixels.size()<min_number_pixels)
                            std::cout<<"no baby added because not enough points"<<std::endl<<std::endl;
                        //make reappear points not used in all_boundaries_image for other walls to use it
                        //put remaining points to inter_remaining of this plane
                        for(auto it_not_repeated = intersections[k].not_repeated.begin(); it_not_repeated != intersections[k].not_repeated.end(); ++it_not_repeated)
                        {
                            std::pair<int,int> pix = intersections[k].other_pixels[*it_not_repeated];
                            all_boundaries_image(pix.first, pix.second) = 0;
                            inter_remaining.pixels.push_back(intersections[k].other_pixels[*it_not_repeated]);
                            inter_remaining.points.push_back(intersections[k].other_points[*it_not_repeated]);
                        }
                    }

                    for(int idx_remaining_points = 0; idx_remaining_points < intersections[k].remaining_points.size(); ++idx_remaining_points)
                    {
                        inter_remaining.pixels.push_back(intersections[k].remaining_pixels[idx_remaining_points]);
                        inter_remaining.points.push_back(intersections[k].remaining_points[idx_remaining_points]);
                    }
                }
                else
                {
                    std::cout<<"intersection n°"<<k<<" does not contain enough points onto the plane studied"<<std::endl;
                    if(intersections[k].pixels.size() + intersections[k].not_repeated.size()>=min_number_pixels) //if there remain enough points we consider this intersection is not modelisable by a line but it still exists
                    {
                        std::cout<<"intersection n°"<<k<<" can not be modelize by a line because the best line found only contains = "<<intersections[k].points.size()<<" points"<<std::endl;
                        //make reappear points not used in all_boundaries_image for other walls to use it
                        //put remaining points to inter_remaining of this plane
                        for(auto it_not_repeated = intersections[k].not_repeated.begin(); it_not_repeated != intersections[k].not_repeated.end(); ++it_not_repeated)
                        {
                            std::pair<int,int> pix = intersections[k].other_pixels[*it_not_repeated];
                            all_boundaries_image(pix.first, pix.second) = 0;
                            intersections[k].pixels.push_back(intersections[k].other_pixels[*it_not_repeated]);
                            intersections[k].points.push_back(intersections[k].other_points[*it_not_repeated]);
                        }
                        intersections[k].isLine = false;
                    }
                    else // if there does not remain enough points we consider they are noise around the other intersection and we remove them
                        remove_intersection(&k);
                }
            }
        }
        //-----------------------------------------------------------------------------------------------------------------------------------
        else if(!intersections[k].isConnection && !intersections[k].isObstruction && !intersections[k].isOpening && !intersections[k].isObject)           //begining : can be obstruction or connection
        {
            std::cout<<"this intersection is a connection or an obstruction which has not been detected yet between planes n°"<<intersections[k].plane_ref->index<<" and plane n°"<<intersections[k].plane_neigh->index<<std::endl;
            if(intersections[k].pixels.size()>=min_number_pixels)
            {
                intersections[k].computeTheoriticalLineFeaturesConnection(); // does not touch points/pixels other_points/other_pixels pt_mean -> just modify normal/tangente
                processed_planes(intersections[k].plane_neigh->index, intersections[k].plane_ref->index) = true;

                intersections[k].definePlaneConnection();               // define if the intersection is a real connection between planes or an obstruction (obstructed or obstructing) and if it is a line
                                                                        //move all points which do not belong to a connection from points to other_points
                //create new intersections with other_points to search obstructions

                if(intersections[k].other_pixels.size()>=min_number_pixels)
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
                    for(int idx_pixels = 0; idx_pixels < intersections[k].other_pixels.size(); ++idx_pixels)
                    {
                        all_boundaries_image(intersections[k].other_pixels[idx_pixels].first, intersections[k].other_pixels[idx_pixels].second) = 0;
                        inter_remaining.pixels.push_back(intersections[k].other_pixels[idx_pixels]);
                        inter_remaining.points.push_back(intersections[k].other_points[idx_pixels]);
                    }
                }

                if(intersections[k].isConnection)
                {
                    std::cout<<"Conclusion : it is a connection with plane n°"<<intersections[k].plane_neigh->index<<std::endl<<std::endl;
                    intersections[k].isLine = true;
                    std::vector<intersection> vec_sisters = intersections[k].export_sisters();

                    std::cout<<"Exporting connection pieces"<<std::endl<<std::endl;
                    std::cout<<"Number of pieces : "<<  vec_sisters.size() + 1<<std::endl<<std::endl;
                    for(int idx_vec_sisters = 0; idx_vec_sisters < vec_sisters.size(); ++idx_vec_sisters)
                    {
                        if(vec_sisters[idx_vec_sisters].points.size() >= min_number_points_on_line)
                        {
                            std::cout<<"baby connection is n°"<<intersections.size()<<std::endl<<std::endl;
                            intersections.push_back(vec_sisters[idx_vec_sisters]);
                        }
                        else
                        {
                            std::cout<<"baby connection not added :  only contains : "<<vec_sisters[idx_vec_sisters].points.size()<<" points"<<std::endl<<std::endl;
                            for(int idx_points_vec_sisters = 0; idx_points_vec_sisters<vec_sisters[idx_vec_sisters].points.size(); ++idx_points_vec_sisters)
                                inter_remaining.points.push_back(vec_sisters[idx_vec_sisters].points[idx_points_vec_sisters]);
                        }
                    }
                    if(intersections[k].points.size()< min_number_points_on_line)
                        remove_intersection(&k);
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
            auto it_before_end_edges = edges[plane_idx].end();
            --it_before_end_edges;
            for(auto it_edge = edges[plane_idx].begin(); it_edge != it_before_end_edges; ++it_edge)
            {
                auto it_edge_other_temp = it_edge;
                ++it_edge_other_temp;
                for(auto it_edge_other = it_edge_other_temp; it_edge_other!= edges[plane_idx].end(); ++it_edge_other)
                {
                    if(!all_edges[*it_edge].isConnection || !all_edges[*it_edge_other].isConnection)
                    {
                        if(abs(all_edges[*it_edge].tangente.dot(all_edges[*it_edge_other].tangente)) < not_parallel_dot_threshold)
                        {                                
                            corner c;
                            c.setLines(&all_edges[*it_edge], &all_edges[*it_edge_other]);
                            c.computePointFromLines(); // compute intersection point between lines on studied plane
                            c.isLineIntersection = true;

                            bool isCorner = replaceLim2(all_edges[*it_edge], all_edges[*it_edge_other], c.pt);

                            if(isCorner)
                                possible_corners.push_back(c);
                        }
                    }
                }
            }
        }
    }
}

void manager::computeParallelLinesIntersections()
{
    for(int plane_idx = 0; plane_idx<planes.size(); ++plane_idx)
    {
        if(edges[plane_idx].size()>=2)
        {
            auto it_before_end_edges = edges[plane_idx].end();
            --it_before_end_edges;
            for(auto it_edge = edges[plane_idx].begin(); it_edge != it_before_end_edges; ++it_edge)
            {
                auto it_edge_other_temp = it_edge;
                ++it_edge_other_temp;
                for(auto it_edge_other = it_edge_other_temp; it_edge_other!= edges[plane_idx].end(); ++it_edge_other)
                {
                    if(!all_edges[*it_edge].isConnection || !all_edges[*it_edge_other].isConnection)
                    {
                        double dist;
                        Eigen::Vector3d pt;
                        bool to_change;

                        if(all_edges[*it_edge].tangente.dot(all_edges[*it_edge_other].tangente)<0) // two lines oriented in opposite directions
                        {

                            double dist1 = (all_edges[*it_edge].end_pt - all_edges[*it_edge_other].end_pt).norm();
                            double dist2 = (all_edges[*it_edge].start_pt - all_edges[*it_edge_other].start_pt).norm();

                            if(dist1<dist2)
                            {
                                dist = dist1;
                                pt = (all_edges[*it_edge].end_pt + all_edges[*it_edge_other].end_pt) / 2;
                                to_change = !all_edges[*it_edge].end_changed && !all_edges[*it_edge_other].end_changed;
                            }
                            else
                            {
                                dist = dist2;
                                pt = (all_edges[*it_edge].start_pt + all_edges[*it_edge_other].start_pt) / 2;
                                to_change = !all_edges[*it_edge].start_changed && !all_edges[*it_edge_other].start_changed;
                            }
                        }
                        else // two lines oriented in the same direction
                        {
                            double dist1 = (all_edges[*it_edge].start_pt - all_edges[*it_edge_other].end_pt).norm();
                            double dist2 = (all_edges[*it_edge].end_pt - all_edges[*it_edge_other].start_pt).norm();

                            if(dist1<dist2)
                            {
                                dist = dist1;
                                pt = (all_edges[*it_edge].start_pt + all_edges[*it_edge_other].end_pt) / 2;
                                to_change = !all_edges[*it_edge].start_changed && !all_edges[*it_edge_other].end_changed;
                            }
                            else
                            {
                                dist = dist2;
                                pt = (all_edges[*it_edge].end_pt + all_edges[*it_edge_other].start_pt) / 2;
                                to_change = !all_edges[*it_edge].end_changed && !all_edges[*it_edge_other].start_changed;
                            }
                        }

                        bool lim_distance = dist < max_spat_dist_for_line_connection;

                        if(lim_distance && to_change)
                        {
                            corner c;
                            c.setLines(&all_edges[*it_edge], &all_edges[*it_edge_other]);
                            c.pt = pt;
                            c.isLineContinuity = true;

                            bool isCorner = replaceLim2(all_edges[*it_edge], all_edges[*it_edge_other], c.pt);
                            if(isCorner)
                            {
                                possible_corners.push_back(c);
                                std::cout<<"corner of continuity : "<<pt.transpose()<<" "<<"between lines : "<<*it_edge<<" and "<<*it_edge_other<<std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
}

bool manager::replaceLim2(intersection& inter1, intersection& inter2, Eigen::Vector3d potential_corner_pt)
{
    std::pair<int,int> potential_corner_pixel = pt2XY(potential_corner_pt, delta_phi, delta_theta, rot_axis, axis_init_phi);
    std::pair<int,int> start_pixel1 = pt2XY(inter1.start_pt, delta_phi, delta_theta, rot_axis, axis_init_phi);
    std::pair<int,int> end_pixel1 = pt2XY(inter1.end_pt, delta_phi, delta_theta, rot_axis, axis_init_phi);
    std::pair<int,int> start_pixel2 = pt2XY(inter2.start_pt, delta_phi, delta_theta, rot_axis, axis_init_phi);
    std::pair<int,int> end_pixel2 = pt2XY(inter2.end_pt, delta_phi, delta_theta, rot_axis, axis_init_phi);

    float start_diff_pix1 = sqrt( pow( start_pixel1.first - potential_corner_pixel.first,2 ) + pow( start_pixel1.second - potential_corner_pixel.second,2 ));
    float end_diff_pix1 = sqrt( pow( end_pixel1.first - potential_corner_pixel.first,2 ) + pow( end_pixel1.second - potential_corner_pixel.second,2 ));
    float start_diff_pix2 = sqrt( pow( start_pixel2.first - potential_corner_pixel.first,2 ) + pow( start_pixel2.second - potential_corner_pixel.second,2 ));
    float end_diff_pix2 = sqrt( pow( end_pixel2.first - potential_corner_pixel.first,2 ) + pow( end_pixel2.second - potential_corner_pixel.second,2 ));

    float start_diff_spat1 = (inter1.start_pt - potential_corner_pt).norm();
    float end_diff_spat1   = (inter1.end_pt - potential_corner_pt).norm();
    float start_diff_spat2 = (inter2.start_pt - potential_corner_pt).norm();
    float end_diff_spat2   = (inter2.end_pt - potential_corner_pt).norm();

    float dist_line1 = ((potential_corner_pt - inter1.pt_mean) - (potential_corner_pt - inter1.pt_mean).dot(inter1.tangente) * inter1.tangente).norm();
    float dist_line2 = ((potential_corner_pt - inter2.pt_mean) - (potential_corner_pt - inter2.pt_mean).dot(inter2.tangente) * inter2.tangente).norm();

    double eps = 0.00001;

    bool start1_replaced = ( start_diff_pix1 < end_diff_pix1 && ( start_diff_pix1-inter1.max_pixel_diff_start<eps || start_diff_spat1-inter1.max_spatial_diff_start<eps ) );
    bool end1_replaced = ( end_diff_pix1 < start_diff_pix1  && ( end_diff_pix1-inter1.max_pixel_diff_end<eps || end_diff_spat1-inter1.max_spatial_diff_end<eps ) );
    bool start2_replaced = ( start_diff_pix2 < end_diff_pix2 && ( start_diff_pix2-inter2.max_pixel_diff_start<eps || start_diff_spat2-inter2.max_spatial_diff_start<eps ) );
    bool end2_replaced = ( end_diff_pix2 < start_diff_pix2  && ( end_diff_pix2-inter2.max_pixel_diff_end<eps || end_diff_spat2-inter2.max_spatial_diff_end<eps ) );

    bool start1_can_be_replaced = ( dist_line1<max_line_distance && ( start_diff_pix1 < max_pix_dist_for_line_connection || start_diff_spat1 < max_spat_dist_for_line_connection) );
    bool end1_can_be_replaced = ( dist_line1<max_line_distance && ( end_diff_pix1 < max_pix_dist_for_line_connection || end_diff_spat1 < max_spat_dist_for_line_connection) );
    bool start2_can_be_replaced = ( dist_line2<max_line_distance && (start_diff_pix2 < max_pix_dist_for_line_connection || start_diff_spat2 < max_spat_dist_for_line_connection) );
    bool end2_can_be_replaced = ( dist_line2<max_line_distance && (end_diff_pix2 < max_pix_dist_for_line_connection || end_diff_spat2 < max_spat_dist_for_line_connection) );

    if( (start1_can_be_replaced || end1_can_be_replaced) && (start2_can_be_replaced || end2_can_be_replaced) ) // if one line extremity can be the other line extremity
    {
        if(start1_replaced) // if the current potential corner is better than previous ones (in pixels or in space)
        {
            inter1.new_start_pt = potential_corner_pt;
            inter1.max_pixel_diff_start = start_diff_pix1;
            inter1.max_spatial_diff_start  = start_diff_spat1;
            inter1.start_changed = true;
        }
        if(end1_replaced)
        {
            inter1.new_end_pt = potential_corner_pt;
            inter1.max_pixel_diff_end = end_diff_pix1;
            inter1.max_spatial_diff_end = end_diff_spat1;
            inter1.end_changed = true;
        }

        //--------------------------------------------------------------------------

        if(start2_replaced)
        {
            inter2.new_start_pt = potential_corner_pt;
            inter2.max_pixel_diff_start = start_diff_pix2;
            inter2.max_spatial_diff_start  = start_diff_spat2;
            inter2.start_changed = true;
        }
        if(end2_replaced)
        {
            inter2.new_end_pt = potential_corner_pt;
            inter2.max_pixel_diff_end = end_diff_pix2;
            inter2.max_spatial_diff_end  = end_diff_spat2;
            inter2.end_changed = true;
        }

        if( (start1_replaced || end1_replaced) || (start2_replaced || end2_replaced) )
            return true;
    }

    return false;
}


void manager::computeTheoriticalPlanesIntersections()
{
    std::vector<std::vector<std::vector<bool>>> already_treated(planes.size(), std::vector<std::vector<bool>>(planes.size(), std::vector<bool>(planes.size(), false)));
    for(int plane_idx = 0; plane_idx<planes.size(); ++plane_idx)
    {
        if(edges[plane_idx].size()>=2)
        {
            auto it_edge_end = edges[plane_idx].end();
            --it_edge_end;
            for(auto it_edge_first = edges[plane_idx].begin(); it_edge_first != it_edge_end; ++it_edge_first)
            {
                int i = *it_edge_first;
                if(all_edges[i].isConnection)
                {
                    auto it_edge_second_temp = it_edge_first;
                    ++it_edge_second_temp;
                    for(auto it_edge_second = it_edge_second_temp; it_edge_second != edges[plane_idx].end(); ++it_edge_second)
                    {
                        int j = *it_edge_second;
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
                            auto it_edge_third = edges[third_plane->index].begin();
                            for(it_edge_third = edges[third_plane->index].begin(); it_edge_third != edges[third_plane->index].end(); ++it_edge_third)
                            {
                                k = *it_edge_third;
                                if(all_edges[k].isConnection && (all_edges[k].plane_ref->index == second_plane->index || all_edges[k].plane_neigh->index == second_plane->index))
                                {
                                    if(!already_treated[first_plane->index][second_plane->index][third_plane->index] && it_edge_third != edges[third_plane->index].end())
                                    {
                                        corner c;
                                        c.setPlanes(first_plane, second_plane, third_plane);
                                        c.isPlaneIntersection = true;
                                        c.computePointFromPlanes();

                                        bool isCorner = replaceLim3(all_edges[i], all_edges[j], all_edges[k], c.pt);

                                        if(isCorner)
                                        {
                                            possible_corners.push_back(c);
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
            }
        }
    }
}

bool manager::replaceLim3(intersection& inter1, intersection& inter2, intersection& inter3, Eigen::Vector3d potential_corner_pt)
{
    std::pair<int,int> potential_corner_pixel = pt2XY(potential_corner_pt, delta_phi, delta_theta, rot_axis, axis_init_phi);
    std::pair<int,int> start_pixel1 = pt2XY(inter1.start_pt, delta_phi, delta_theta, rot_axis, axis_init_phi);
    std::pair<int,int> end_pixel1 = pt2XY(inter1.end_pt, delta_phi, delta_theta, rot_axis, axis_init_phi);
    std::pair<int,int> start_pixel2 = pt2XY(inter2.start_pt, delta_phi, delta_theta, rot_axis, axis_init_phi);
    std::pair<int,int> end_pixel2 = pt2XY(inter2.end_pt, delta_phi, delta_theta, rot_axis, axis_init_phi);
    std::pair<int,int> start_pixel3 = pt2XY(inter3.start_pt, delta_phi, delta_theta, rot_axis, axis_init_phi);
    std::pair<int,int> end_pixel3 = pt2XY(inter3.end_pt, delta_phi, delta_theta, rot_axis, axis_init_phi);

    float start_diff_pix1 = sqrt( pow( start_pixel1.first - potential_corner_pixel.first,2 ) + pow( start_pixel1.second - potential_corner_pixel.second,2 ));
    float end_diff_pix1 = sqrt( pow( end_pixel1.first - potential_corner_pixel.first,2 ) + pow( end_pixel1.second - potential_corner_pixel.second,2 ));
    float start_diff_pix2 = sqrt( pow( start_pixel2.first - potential_corner_pixel.first,2 ) + pow( start_pixel2.second - potential_corner_pixel.second,2 ));
    float end_diff_pix2 = sqrt( pow( end_pixel2.first - potential_corner_pixel.first,2 ) + pow( end_pixel2.second - potential_corner_pixel.second,2 ));
    float start_diff_pix3 = sqrt( pow( start_pixel3.first - potential_corner_pixel.first,2 ) + pow( start_pixel3.second - potential_corner_pixel.second,2 ));
    float end_diff_pix3 = sqrt( pow( end_pixel3.first - potential_corner_pixel.first,2 ) + pow( end_pixel3.second - potential_corner_pixel.second,2 ));

    float start_diff_spat1 = (inter1.start_pt - potential_corner_pt).norm();
    float end_diff_spat1   = (inter1.end_pt - potential_corner_pt).norm();
    float start_diff_spat2 = (inter2.start_pt - potential_corner_pt).norm();
    float end_diff_spat2   = (inter2.end_pt - potential_corner_pt).norm();
    float start_diff_spat3 = (inter3.start_pt - potential_corner_pt).norm();
    float end_diff_spat3   = (inter3.end_pt - potential_corner_pt).norm();

    float dist_line1 = ((potential_corner_pt - inter1.pt_mean) - (potential_corner_pt - inter1.pt_mean).dot(inter1.tangente) * inter1.tangente).norm();
    float dist_line2 = ((potential_corner_pt - inter2.pt_mean) - (potential_corner_pt - inter2.pt_mean).dot(inter2.tangente) * inter2.tangente).norm();
    float dist_line3 = ((potential_corner_pt - inter3.pt_mean) - (potential_corner_pt - inter3.pt_mean).dot(inter3.tangente) * inter3.tangente).norm();

    double eps = 0.00001;

    bool start1_replaced = ( start_diff_pix1 < end_diff_pix1 && ( start_diff_pix1-inter1.max_pixel_diff_start<eps || start_diff_spat1-inter1.max_spatial_diff_start<eps ) );
    bool end1_replaced = ( end_diff_pix1 < start_diff_pix1  && ( end_diff_pix1-inter1.max_pixel_diff_end<eps || end_diff_spat1-inter1.max_spatial_diff_end<eps ) );
    bool start2_replaced = ( start_diff_pix2 < end_diff_pix2 && ( start_diff_pix2-inter2.max_pixel_diff_start<eps || start_diff_spat2-inter2.max_spatial_diff_start<eps ) );
    bool end2_replaced = ( end_diff_pix2 < start_diff_pix2  && ( end_diff_pix2-inter2.max_pixel_diff_end<eps || end_diff_spat2-inter2.max_spatial_diff_end<eps ) );
    bool start3_replaced = ( start_diff_pix3 < end_diff_pix3 && ( start_diff_pix3-inter3.max_pixel_diff_start<eps || start_diff_spat3-inter3.max_spatial_diff_start<eps ) );
    bool end3_replaced = ( end_diff_pix3 < start_diff_pix3  && ( end_diff_pix3-inter3.max_pixel_diff_end<eps || end_diff_spat3-inter3.max_spatial_diff_end<eps ) );

    bool start1_can_be_replaced = ( dist_line1< max_line_distance && ( start_diff_pix1 < max_pix_dist_for_line_connection || start_diff_spat1 < max_spat_dist_for_line_connection) );
    bool end1_can_be_replaced = ( dist_line1< max_line_distance && (end_diff_pix1 < max_pix_dist_for_line_connection || end_diff_spat1 < max_spat_dist_for_line_connection) );
    bool start2_can_be_replaced = ( dist_line2< max_line_distance && (start_diff_pix2 < max_pix_dist_for_line_connection || start_diff_spat2 < max_spat_dist_for_line_connection) );
    bool end2_can_be_replaced = ( dist_line2< max_line_distance && (end_diff_pix2 < max_pix_dist_for_line_connection || end_diff_spat2 < max_spat_dist_for_line_connection) );
    bool start3_can_be_replaced = ( dist_line3< max_line_distance && (start_diff_pix3 < max_pix_dist_for_line_connection || start_diff_spat3 < max_spat_dist_for_line_connection) );
    bool end3_can_be_replaced = ( dist_line3< max_line_distance && ( end_diff_pix3 < max_pix_dist_for_line_connection || end_diff_spat3 < max_spat_dist_for_line_connection) );

    if( (start1_can_be_replaced || end1_can_be_replaced) && (start2_can_be_replaced || end2_can_be_replaced) && (start3_can_be_replaced || end3_can_be_replaced) ) // if one line extremity can be the other line extremity
    {
        if(start1_replaced) // if the current potential corner is better than previous ones (in pixels or in space)
        {
            inter1.new_start_pt = potential_corner_pt;
            inter1.max_pixel_diff_start = start_diff_pix1;
            inter1.max_spatial_diff_start  = start_diff_spat1;
            inter1.start_changed = true;
        }
        if(end1_replaced)
        {
            inter1.new_end_pt = potential_corner_pt;
            inter1.max_pixel_diff_end = end_diff_pix1;
            inter1.max_spatial_diff_end = end_diff_spat1;
            inter1.end_changed = true;
        }

        //--------------------------------------------------------------------------

        if(start2_replaced)
        {
            inter2.new_start_pt = potential_corner_pt;
            inter2.max_pixel_diff_start = start_diff_pix2;
            inter2.max_spatial_diff_start  = start_diff_spat2;
            inter2.start_changed = true;
        }
        if(end2_replaced)
        {
            inter2.new_end_pt = potential_corner_pt;
            inter2.max_pixel_diff_end = end_diff_pix2;
            inter2.max_spatial_diff_end  = end_diff_spat2;
            inter2.end_changed = true;
        }

        //--------------------------------------------------------------------------

        if(start3_replaced)
        {
            inter3.new_start_pt = potential_corner_pt;
            inter3.max_pixel_diff_start = start_diff_pix3;
            inter3.max_spatial_diff_start  = start_diff_spat3;
            inter3.start_changed = true;
        }
        if(end3_replaced)
        {
            inter3.new_end_pt = potential_corner_pt;
            inter3.max_pixel_diff_end = end_diff_pix3;
            inter3.max_spatial_diff_end  = end_diff_spat3;
            inter3.end_changed = true;
        }

        if( (start1_replaced || end1_replaced) || (start2_replaced || end2_replaced) || (start3_replaced || end3_replaced) )
            return true;
    }

    return false;
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

void manager::order_polygone_corners()
{
    for(int idx_plane; idx_plane < planes.size(); ++idx_plane)
    {
        std::cout<<"plane n°"<<idx_plane<<std::endl;
        std::vector<Eigen::Vector3d> ordered_corners;
        std::vector<int> processed_line_indices(all_edges.size(), 0);

        Eigen::Vector3d current_point;
        Eigen::Vector3d init_point;
        std::cout<<"selecting seed "<<std::endl;

        if(planes[idx_plane].intersections_indices.size()>1)
        {
            //select seed edge--------------------------------------------------------
            for(auto it_idx_lines = planes[idx_plane].intersections_indices.begin(); it_idx_lines != planes[idx_plane].intersections_indices.end(); ++it_idx_lines)
            {
                if( all_edges[*it_idx_lines].isLine)
                {
                    if(isLineConnectedInPlane(*it_idx_lines, idx_plane, 0))
                    {
                        current_point = all_edges[*it_idx_lines].new_start_pt;
                        init_point = all_edges[*it_idx_lines].new_end_pt;
                    }
                    else if(isLineConnectedInPlane(*it_idx_lines, idx_plane, 1))
                    {
                        current_point = all_edges[*it_idx_lines].new_end_pt;
                        init_point = all_edges[*it_idx_lines].new_start_pt;
                    }

                    std::cout<<"got seed intersection n°"<<*it_idx_lines<<std::endl;

                    ordered_corners.push_back(init_point);
                    ordered_corners.push_back(current_point);
                    processed_line_indices[*it_idx_lines] = 1;
                    break;
                }
            }

            std::cout<<"init point : "<<init_point.transpose()<<std::endl<<std::endl;
            std::cout<<"creating ordered "<<std::endl;
            double min_dist;
            int N_fact = 3;
            do
            {
                min_dist = N_fact*(current_point - init_point).norm();
                Eigen::Vector3d hypothetic_current_point;
                Eigen::Vector3d temp_point;
                int idx;
                for ( auto it_idx_lines_other = planes[idx_plane].intersections_indices.begin(); it_idx_lines_other != planes[idx_plane].intersections_indices.end(); ++it_idx_lines_other)
                {
                    std::cout<<"processed? "<<processed_line_indices[*it_idx_lines_other]<<std::endl;
                    if(processed_line_indices[*it_idx_lines_other] == 0 && all_edges[*it_idx_lines_other].isLine)
                    {
                        std::vector<Eigen::Vector3d> jonction(2);
                        jonction[0] = current_point;
                        jonction[1] = all_edges[*it_idx_lines_other].new_start_pt;

                        //check if the new link is minimal (less than any other links and less than going back to init) + check if it does not cross other previous lines
                        if( (all_edges[*it_idx_lines_other].new_start_pt - current_point).norm() < min_dist && !cross(jonction, ordered_corners, planes[idx_plane].normal) )
                        {
                            if( (all_edges[*it_idx_lines_other].new_start_pt - current_point).norm() < epsilon) // either it's the same point and lines are connected, or it is not the same point and the lines are not connected to anything else.
                            {
                                min_dist = (all_edges[*it_idx_lines_other].new_start_pt - current_point).norm();
                                hypothetic_current_point = all_edges[*it_idx_lines_other].new_end_pt;
                                temp_point = all_edges[*it_idx_lines_other].new_start_pt;
                                idx = *it_idx_lines_other;
                            }
                            else if (!isLineConnectedInPlane(*it_idx_lines_other, idx_plane, 0))
                            {
                                min_dist = (all_edges[*it_idx_lines_other].new_start_pt - current_point).norm();
                                hypothetic_current_point = all_edges[*it_idx_lines_other].new_end_pt;
                                temp_point = all_edges[*it_idx_lines_other].new_start_pt;
                                idx = *it_idx_lines_other;
                            }
                        }

                        jonction[0] = current_point;
                        jonction[1] = all_edges[*it_idx_lines_other].new_end_pt;
                        //check if the new link is minimal (less than any other links and less than going back to init) + check if it does not cross other previous lines
                        if( (all_edges[*it_idx_lines_other].new_end_pt - current_point).norm() < min_dist && !cross(jonction, ordered_corners, planes[idx_plane].normal))
                        {
                            if( (all_edges[*it_idx_lines_other].new_end_pt - current_point).norm() < epsilon) // either it's the same point and lines are connected with each other, or it is not the same point and the lines are not connected to anything else.
                            {
                                min_dist = (all_edges[*it_idx_lines_other].new_end_pt - current_point).norm();
                                hypothetic_current_point = all_edges[*it_idx_lines_other].new_start_pt;
                                temp_point = all_edges[*it_idx_lines_other].new_end_pt;
                                idx = *it_idx_lines_other;
                            }
                            else if (!isLineConnectedInPlane(*it_idx_lines_other, idx_plane, 1))
                            {
                                min_dist = (all_edges[*it_idx_lines_other].new_end_pt - current_point).norm();
                                hypothetic_current_point = all_edges[*it_idx_lines_other].new_start_pt;
                                temp_point = all_edges[*it_idx_lines_other].new_end_pt;
                                idx = *it_idx_lines_other;
                            }
                        }
                    }
                }

                if(abs(min_dist - N_fact*(current_point - init_point).norm()) > epsilon && (temp_point-current_point).norm() > 0.001) //check if one line was found and if it is not connected to the anterior
                    ordered_corners.push_back(temp_point);

                processed_line_indices[idx] = 1;
                current_point = hypothetic_current_point;

                if((current_point - init_point).norm()>epsilon)
                  ordered_corners.push_back(current_point);

                std::cout<<"current point : "<<current_point.transpose()<<std::endl;
            }
            while((current_point - init_point).norm()>epsilon && abs(min_dist - N_fact*(current_point - init_point).norm())>epsilon); // while a cycle of points is not found

            planes[idx_plane].ordered_corners = ordered_corners;
        }
    }
}

bool manager::isLineConnectedInPlane(int idx_line, int idx_plane, int end) // start = 0, end = 1
{
    for(auto it_idx_lines = planes[idx_plane].intersections_indices.begin(); it_idx_lines != planes[idx_plane].intersections_indices.end(); ++it_idx_lines)
    {
        if(idx_line != *it_idx_lines)
        {
            if(!end)
            {
                if(all_edges[idx_line].indices_connected_lines_start.find(*it_idx_lines) != all_edges[idx_line].indices_connected_lines_start.end())
                    return true;
            }
            else
            {
                if(all_edges[idx_line].indices_connected_lines_end.find(*it_idx_lines) != all_edges[idx_line].indices_connected_lines_start.end())
                    return true;
            }
        }

    }
    return false;
}

bool manager::cross(std::vector<Eigen::Vector3d> jonction, std::vector<Eigen::Vector3d> poly, Eigen::Vector3d normal)
{
    if(poly.size()>2 && (jonction[1]-jonction[0]).norm()>0.001)
    {
        //define features of lines in 2D
        Eigen::Affine3d rot = Eigen::Affine3d::Identity();

        float angle = acos(normal(2));   //angle between normal and z axis [0, pi]
        Eigen::Vector3d axis;
        if(angle>epsilon)
        {
            axis = normal.cross(Eigen::Vector3d(0,0,1)); // rotation axis to align normal onto z axis
            axis /= axis.norm();
            rot.rotate( Eigen::AngleAxisd(angle, axis) );
        }

        Eigen::Vector3d vec = jonction[1]-jonction[0];
        vec /= vec.norm();
        Eigen::Vector3d temp = rot.linear() * vec;
        Eigen::Vector2d tangente2D = {temp(0), temp(1)};
        temp = rot.linear() * jonction[0];
        Eigen::Vector2d pinit = {temp(0), temp(1)};
        temp = rot.linear() * jonction[1];
        Eigen::Vector2d pfin = {temp(0), temp(1)};
        Eigen::Vector2d normal2D = pfin - pfin.dot(tangente2D) * tangente2D;
        normal2D /= normal2D.norm();
        double distance2D = pfin.dot(normal2D);

        for(int i = 0; i<poly.size()-1; ++i)
        {
            Eigen::Vector3d vec_temp = poly[i+1]-poly[i];
            vec_temp /= vec_temp.norm();
            if(abs(vec.dot(vec_temp)) < 0.999) // vectors not parallel
            {
                Eigen::Vector2d tangente2D_temp = {(rot.linear() * vec_temp)(0), (rot.linear() * vec_temp)(1)};
                Eigen::Vector2d pinit_temp = {(rot.linear() * poly[i])(0), (rot.linear() * poly[i])(1)};
                Eigen::Vector2d pfin_temp = {(rot.linear() * poly[i+1])(0), (rot.linear() * poly[i+1])(1)};
                Eigen::Vector2d normal2D_temp = pfin_temp - pfin_temp.dot(tangente2D_temp) * tangente2D_temp;
                normal2D_temp /= normal2D_temp.norm();
                double distance2D_temp = pfin_temp.dot(normal2D_temp);

                Eigen::Matrix2d N ;
                N.row(0) = normal2D_temp;
                N.row(1) = normal2D;
                Eigen::Vector2d d = {distance2D_temp, distance2D};
                Eigen::Vector2d pt_intersect2D = N.colPivHouseholderQr().solve(d);
                if((pt_intersect2D-pinit).dot(tangente2D)>0 && (pt_intersect2D-pfin).dot(tangente2D)<0 && (pt_intersect2D-pinit_temp).dot(tangente2D_temp)>0 && (pt_intersect2D-pfin_temp).dot(tangente2D_temp)<0)
                    return true;
            }
        }
    }
    return false;
}

void manager::export_mesh()
{
    std::cout<<"Start exporting"<<std::endl<<std::endl;

    std::cout<<"Start ordering polygons"<<std::endl<<std::endl;
    order_polygone_corners();

    std::cout<<"Start writing obj"<<std::endl<<std::endl;
    std::ofstream obj_file ("mesh.obj");
    int idx_line_of_file = 1;
    std::vector<std::vector<int>> idx_of_line_of_file_per_plane(planes.size());
    int n_vertices = 0;
    int n_faces = 0;

    for(int idx_plane = 0; idx_plane<planes.size(); ++idx_plane)
    {
        if(planes[idx_plane].ordered_corners.size()>1)
        {
            n_vertices += planes[idx_plane].ordered_corners.size();
            ++n_faces;
        }
    }

    obj_file << n_vertices<<"\n";
    obj_file << n_faces<<"\n";
    for(int idx_plane = 0; idx_plane<planes.size(); ++idx_plane)
    {
        if(planes[idx_plane].ordered_corners.size()>1)
        {
            for(int idx_corner = 0; idx_corner<planes[idx_plane].ordered_corners.size(); ++idx_corner)
            {
                obj_file << "v "<<planes[idx_plane].ordered_corners[idx_corner].transpose()<<"\n";
                idx_of_line_of_file_per_plane[idx_plane].push_back(idx_line_of_file);
                ++idx_line_of_file;
            }
        }
    }
    for(int idx_plane = 0; idx_plane<planes.size(); ++idx_plane)
    {
        if(planes[idx_plane].ordered_corners.size()>1)
        {
            obj_file << "f";
            for(int idx_corner = 0; idx_corner<idx_of_line_of_file_per_plane[idx_plane].size(); ++idx_corner)
                obj_file <<" "<< idx_of_line_of_file_per_plane[idx_plane][idx_corner];
            obj_file <<"\n";
        }
    }
    obj_file.close();
}


void manager::fill_edges()
{
    edges.clear();
    edges.resize(planes.size()+1);
    for(int k = 0; k<all_edges.size(); ++k)
    {
        all_edges[k].index = k;
        if(all_edges[k].isLine)
        {
            edges[all_edges[k].plane_ref->index].insert(k);
            if(all_edges[k].isConnection)
                edges[all_edges[k].plane_neigh->index].insert(k);
        }
    }
}

void manager::actualizeChanged()
{
    for(int c = 0; c<possible_corners.size(); ++c)
    {
        if(possible_corners[c].isLineIntersection)
        {
            bool is_l0_start= (possible_corners[c].pt-possible_corners[c].lines[0]->new_start_pt).norm() < epsilon;
            bool is_l0_end = (possible_corners[c].pt-possible_corners[c].lines[0]->new_end_pt).norm() < epsilon;
            bool is_l1_start= (possible_corners[c].pt-possible_corners[c].lines[1]->new_start_pt).norm() < epsilon;
            bool is_l1_end = (possible_corners[c].pt-possible_corners[c].lines[1]->new_end_pt).norm() < epsilon;
            bool exist_in_line0= is_l0_start || is_l0_end;
            bool exist_in_line1= is_l1_start || is_l1_end;

            if( (!exist_in_line0 || !exist_in_line1) && (exist_in_line0 || exist_in_line1))
            {
                if(exist_in_line0)
                {
                    if(is_l0_start)
                        possible_corners[c].lines[0]->start_changed = false;
                    else
                        possible_corners[c].lines[0]->end_changed = false;
                }
                else
                {
                    if(is_l1_start)
                        possible_corners[c].lines[1]->start_changed = false;
                    else
                        possible_corners[c].lines[1]->end_changed = false;
                }
            }
        }
    }
}
