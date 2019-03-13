manager::manager(typename pcl::PointCloud<pcl::PointNormal>::Ptr c, int Nr, int Nc, double pa, double ta, int mc)
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
}

void manager::quantifyPhiTheta()
{
    int n = 0;
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
    std::cout<<"points non inserted in XY2PtN because same X and same Y (due to 3d noise affecting small theta): "<<n<<std::endl<<std::endl;
}

void manager::initSearchCluster(double rad, Eigen::Vector3d rot_axis)
{
    man3D.getNormals(rad);
    man3D.setRotAxis(rot_axis);
    man3D.createMap();
    PhiTheta2PtN = man3D.getMap();
    quantifyPhiTheta();
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
                    image_clusterized_indices(i,j) = planes.size()+1; // put indice of another cluster
            }
        }
    }
}

bool manager::isSeed(std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator it, int radius, double thresh_neigh_for_seed)
{
    std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator> neigh = getNeighbors(it, radius);
    Eigen::Vector3d mean_pt;
    for(int i = 0; i< neigh.size(); ++i)
        mean_pt += neigh[i]->second.first;

    mean_pt /= neigh.size();
    mean_pt = it->second.first;

    for(int i = 0; i< neigh.size(); ++i)
    {
        double dist_neigh = abs(neigh[i]->second.first.dot(neigh[i]->second.second));
        double dist = abs(mean_pt.dot(it->second.second));
        if( abs(dist- dist_neigh) > thresh_neigh_for_seed || abs(neigh[i]->second.second.dot(it->second.second))<normals_similarity_threshold_to_select_seed)
            return false;
    }
    it->second.first = mean_pt;
    return true;
//    return man3D.isSeed(it->second.first, it->second.second, 0.1, thresh_neigh_for_seed);
}

void manager::searchClusters(double thresh_plane_belonging, int radius, double thresh_neigh_for_seed) // create the object regions in which there is a vector of points with their associated mode
{
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> processed_seed = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> processed_growing = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();
    pcl::PointCloud<pcl::PointXYZ> pc_clus;

    for(auto it = XY2PtN.begin(); it != XY2PtN.end(); ++it)
    {
        if(!processed_seed(it->first.first, it->first.second)) // in order not to process one point which has already been selected in a cluster
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
                    bool all_neigh_with_non_coherent_normal = true;
                    for(int k = 0; k<neighbors.size(); ++k)
                    {
                        if(abs(neighbors[k]->second.second.dot(it->second.second))>normals_similarity_to_add_neighbors)
                        {
                            all_neigh_with_non_coherent_normal = false;
                            break;
                        }
                    }
                    if(!all_neigh_with_non_coherent_normal)
                    {
                        for(int k = 0; k<neighbors.size(); ++k)
                        {
                            if(!processed_growing(neighbors[k]->first.first, neighbors[k]->first.second)) //in order not to process twice one neighbor
                            {
                                if(abs((neighbors[k]->second.first-it->second.first).dot(it->second.second)) < thresh_plane_belonging)// && acos(neighbors[k]->second.dot(it->second)/(neighbors[k]->second.norm()*it->second.norm())) < thresh_angle )
                                {
                                    processed_growing(neighbors[k]->first.first, neighbors[k]->first.second) = true;
                                    region.push_back(neighbors[k]);

                                    pcl::PointXYZ pt_clus;
                                    pt_clus.x = neighbors[k]->second.first(0);
                                    pt_clus.y = neighbors[k]->second.first(1);
                                    pt_clus.z = neighbors[k]->second.first(2);
                                    pc_clus.points.push_back(pt_clus);
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
                    seeds.push_back(it);
//                    std::vector<Eigen::MatrixXi, Eigen::aligned_allocator<Eigen::MatrixXi>> test(3);
//                    test[0] = coeffWiseMulti(image_init[0], processed_seed.cast<int>());
//                    test[1] = coeffWiseMulti(image_init[1], processed_seed.cast<int>());
//                    test[2] = coeffWiseMulti(image_init[2], processed_seed.cast<int>());
//                    save_image_ppm("processed", "", test,max_col);
//                    save_image_pgm("processed", "", processed_growing.cast<int>(), 1);

//                    pc_clus.width = pc_clus.points.size();
//                    pc_clus.height = 1;
//                    pc_clus.is_dense=false;
//                    pcl::io::savePCDFileASCII ("cloud_clus.pcd", pc_clus);

//                    system("bash pcd2csv.sh cloud_clus.pcd");

//                    getchar();
//                    getchar();
                }
            }
        }
    }
}

void manager::cleanClusters() // from regions, keep all points in XY2Plane_clusterized and then clean points having various modes
{
    XY2Plane_clusterized.clear();
    XY2Plane_clusterized_clean.clear();
    doubt.resize(regions.size(), regions.size());
    for(auto it = regions.begin(); it!=regions.end(); ++it) // for each cluster
    {
        for(int i = 0; i<it->first.size(); ++i)
            XY2Plane_clusterized.insert(std::make_pair(it->first[i]->first, it->second));
    }

    std::cout<<"Number of points in XY2Plane_clusterized "<<XY2Plane_clusterized.size()<<std::endl<<std::endl;

    std::cout<<"Start cleaning points having various clusters "<<std::endl<<std::endl;
    int n =0;
    planes.resize(regions.size());
    for(auto it = XY2Plane_clusterized.begin(); it!=XY2Plane_clusterized.end(); ++it) // for each cluster
    {
        auto it_XY2PtN = XY2PtN.find(it->first);              // search in imageRGBf the initial normal of point-> { (X,Y), [(px,py,pz)(nx,ny,nz)]}
        int cnt = XY2Plane_clusterized.count(it->first);
        if(cnt>1) // points which belong to various clusters
        {
            //Count points which have multiple modes
            ++n;

            //Select the best cluster (mode most similar)
            auto it_multiple_limits = XY2Plane_clusterized.equal_range(it->first);
            std::map<double, std::multimap<std::pair<int,int>, Eigen::Vector3d>::iterator> dotty;
            for(auto it_multiple = it_multiple_limits.first; it_multiple != it_multiple_limits.second; ++it_multiple)
            {
                double dot = (it_multiple->second/it_multiple->second.norm()).dot(it_XY2PtN->second.second); // compute scalar product between real normal and region normal
                dotty.insert(std::make_pair(dot, it_multiple));                     // keep result in map to order scalar product and select the region normal which is more similar
            }
            auto it_map = dotty.end();
            --it_map;                                                               // it_map is the ptr of [ scalar_product, idx] with idx corresponding to the more similar region normal

            XY2Plane_clusterized_clean.insert(std::make_pair(it->first, it_map->second->second));

            //Actualize iterator
            for(int i =1; i!=cnt; ++i)
                ++it;
        }
        else // points which belong to one unique cluster
        {
            if( abs( ((it->second/it->second.norm()).dot(it_XY2PtN->second.second)) ) > normals_similarity_threshold_for_cleaning_when_one_cluster )
                XY2Plane_clusterized_clean.insert(std::make_pair(it->first, it->second));
        }

        // To recompute planes features with cleaned points :

        auto it_XY2Plane_clusterized_clean_found = XY2Plane_clusterized_clean.find(it->first);
        for(int i = 0; i<regions.size(); ++i)
        {
            if(it_XY2Plane_clusterized_clean_found != XY2Plane_clusterized_clean.end() && (it_XY2Plane_clusterized_clean_found->second-regions[i].second).norm()<epsilon)
            {
                planes[i].seed = seeds[i];
                planes[i].appendPoint(it_XY2PtN->second.first); // add points belonging to the plane i
                planes[i].appendPixel(it_XY2PtN->first);
                break;
            }
        }
    }

    std::cout<<"Stop cleaning points having various clusters "<<std::endl<<std::endl;

    std::cout<<"Start recomputing normal and distance "<<std::endl<<std::endl;
    //Recompute features

    for(int i = 0; i<planes.size(); ++i)
    {
        if(planes[i].points.rows()>min_number_of_pixels)    // keep only clusters with more than 100 points
        {
            planes[i].computeNormal();                      //recompute normal with points detected in plane
            planes[i].index = i;
        }
        else
        {
            planes.erase(planes.begin() + i);
            regions.erase(regions.begin() + i);
            --i;
        }
    }

    std::cout<<"Stop recomputing normal and distance "<<std::endl<<std::endl;

    // Use planes to correct XY2Plane_clusterized_clean planes
    for(auto it = XY2Plane_clusterized_clean.begin(); it!=XY2Plane_clusterized_clean.end();) // for each cluster
    {
        int i;
        for(i = 0; i<planes.size(); ++i)
        {
            if((it->second-regions[i].second).norm()<epsilon)
            {
                it->second = planes[i].normal * planes[i].distance;
                break;
            }
            if(i == planes.size()-1)
            {
                it = XY2Plane_clusterized_clean.erase(it);
                --it;
            }
        }
        ++it;
    }

    processed_planes = Eigen::MatrixXi::Zero(planes.size(), planes.size()).cast<bool>();
    std::cout<<"Number of points which have multiple modes : "<<n<<std::endl<<std::endl;
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

    int inity;
    int finity;
    if(it->first.second-rad>0)
        inity = it->first.second-rad;
    else
        inity = 0;

    if(it->first.second+rad<Ncol)
        finity = it->first.second+rad;
    else
        finity = Ncol-1;

    for(int i = initx; i != finitx; ++i)
    {
        if(i>=Nrow)
            i = -1;
        else
        {
            for(int j = inity; j != finity; ++j)
            {
                auto idx = XY2PtN.find(std::make_pair(i,j));
                if(idx != XY2PtN.end())
                    neighbors.push_back(idx);
            }
        }
    }
    return neighbors;
}

void manager::extractBoundImage()
{
    man2D.setImageClusterized(image_clusterized_indices); // image_clusterized_indices is a matrix which contains the indices of planes (+1 to let 0 be the non cluster)
    man2D.setMaxCol(max_col);
    pcl::PointCloud<pcl::PointXYZ> all_boundaries;
    std::vector<std::vector<intersection>> vec_openings(planes.size());
    edges.resize(planes.size()+1);
    system("rm results/binarized*.pgm");
    system("rm results/morpho*.pgm");
    system("rm results/boundary*.pgm");
    system("rm results/cloud_boundary_*.pcd");
    system("rm results/all_boundaries.pcd");
    system("rm results/*.csv");
    system("rm theoritical_line*.csv");
    system("rm theoritical_corner*.csv");
    system("rm cloud_boundary_*.pcd");
    system("rm cloud_boundary_*.csv");

    std::stringstream stm;

    for(int k = 0; k<planes.size(); ++k)
    {
        std::cout<<std::endl<<"PLANE NUMBER : "<<k<<std::endl;
        //binarize
        man2D.binarize(k);
        Eigen::Matrix<bool,Eigen::Dynamic, Eigen::Dynamic> binarized = man2D.getImBinarized();

        //apply closure
        Eigen::Matrix<bool, 5, 5> SE_morpho;
        SE_morpho<<true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true;
        man2D.closure(SE_morpho);
        Eigen::Matrix<bool,Eigen::Dynamic, Eigen::Dynamic> morpho = man2D.getImBinarizedMorpho();

        //compute boundary on image
        Eigen::Matrix<bool, 3, 3> SE_grad;
        SE_grad<<true, true, true, true, true, true, true, true, true;
        int margin = (int)(theta_margin/delta_theta); // where we start to take the thetas into account
        man2D.computeBoundaries(SE_grad, margin);
        boundary = man2D.getBoundaries();

        pcl::PointCloud<pcl::PointXYZ> cloud_boundary = extractBoundCloud();
        all_boundaries += cloud_boundary;

        std::cout<<"Start computing theoritical connexions"<<std::endl;
        computeConnections(k); // compute theoritical line corresponding to each pixel of the boundary + fill pixels corresponding to this intersection
        for(int i = 0; i<intersections.size(); ++i)
        {
            if(intersections[i].isLine)
                edges[intersections[i].plane_ref->index].push_back(intersections[i]);
            else
                edges[planes.size()].push_back(intersections[i]);
            if((intersections[i].isOpening || intersections[i].isObstruction) && intersections[i].isLine)
                vec_openings[intersections[i].plane_ref->index].push_back(intersections[i]);
        }
        std::cout<<"Number of openings edges in this plane : "<<vec_openings[k].size()<<std::endl<<std::endl;
        std::cout<<"Number of edges : "<<intersections.size()<<std::endl<<std::endl;
        std::cout<<"Stop computing connexions"<<std::endl<<std::endl;

        //visualization------------------------------------------------------------------------------------------------------------------
        stm.str("");
        stm<<"cloud_boundary_"<<k<<".pcd";
        pcl::io::savePCDFileASCII (stm.str(), cloud_boundary);

        std::stringstream stm2("");
        stm2<<"bash pcd2csv.sh "<<stm.str();
        std::string command = stm2.str();
        system(command.c_str());

        save_image_pgm("binarized", std::to_string(k), binarized.cast<int>()*max_col, max_col);
        save_image_pgm("morpho", std::to_string(k), morpho.cast<int>()*max_col, max_col);
        Eigen::MatrixXi boundary_image = man2D.getBoundariesImage();
        save_image_pgm("boundary", std::to_string(k), boundary_image, max_col);
        //--------------------------------------------------------------------------------------------------------------------------------
    }

    std::cout<<"Start computing theoritical corners"<<std::endl;
    computeCorners();
    computeTheoriticalLinesIntersections(vec_openings);
    std::cout<<"Stop computing theoritical corners"<<std::endl<<std::endl;

    std::vector< std::vector< Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > > LineIntersections;
    LineIntersections.resize(planes.size());

    float delta = 0.01;
    for(int k = 0; k<edges.size(); ++k)
    {
        std::cout<<"Plane number : "<<k<<std::endl;
        std::cout<<"-------------------------------------------------------------"<<std::endl<<std::endl;

        std::set<int> indices_possible_inter3D;
        for(int l = 0; l<edges[k].size(); ++l)
        {
            if(edges[k][l].isLine)
            {
                std::cout<<"Intersection number : "<<l<<std::endl<<std::endl;
                edges[k][l].computeLim(corners, possible_inter3D[k]);
                for(int m = 0; m<edges[k][l].theoriticalLineIntersection.size(); ++m)
                    indices_possible_inter3D.insert(edges[k][l].theoriticalLineIntersection[m]);

                stm.str("");
                stm<<"theoritical_line_"<<k<<"_"<<l<<".csv";
                std::ofstream file(stm.str());
                int n = 0;
                while (n*delta+edges[k][l].start<edges[k][l].end)
                {
                    file << ((n*delta+edges[k][l].start) * edges[k][l].tangente + edges[k][l].distance*edges[k][l].normal).transpose()<<"\n";
                    ++n;
                }
                file.close();
            }
        }

        for(auto it_indices_possible_inter3D = indices_possible_inter3D.begin(); it_indices_possible_inter3D != indices_possible_inter3D.end(); ++it_indices_possible_inter3D)
            LineIntersections[k].push_back(possible_inter3D[k][*it_indices_possible_inter3D]);
    }

    if(corners.size() != 0)
    {
        stm.str("");
        stm<<"theoritical_corners.csv";
        std::ofstream file(stm.str());

        for(int l = 0; l<corners.size(); ++l)
            file << corners[l].pt.transpose()<<"\n";
        file.close();
    }

    pcl::PointCloud<pcl::PointXYZ> LineIntersections_pc;

    for( int k = 0; k<LineIntersections.size(); ++k)
    {
        for( int l = 0; l<LineIntersections[k].size(); ++l)
        {
            pcl::PointXYZ pt;
            pt.x = LineIntersections[k][l](0);
            pt.y = LineIntersections[k][l](1);
            pt.z = LineIntersections[k][l](2);
            LineIntersections_pc.push_back(pt);
        }
    }

    LineIntersections_pc.width = LineIntersections_pc.points.size();
    LineIntersections_pc.height = 1;
    LineIntersections_pc.is_dense = false;
    pcl::io::savePCDFileASCII ("results/LineIntersections.pcd", LineIntersections_pc);
    pcl::io::savePCDFileASCII ("results/all_boundaries.pcd", all_boundaries);

    system("bash pcd2csv.sh results/LineIntersections.pcd");
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
                    std::cout<<"can't find this boundary point in 3D pointcloud"<<std::endl<<std::endl;
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

void manager::computeConnections(int idx_plane)
{
    ++idx_plane;
    intersections.resize(0);
    int intersections_number = 0;
    int idx_plane_neighbor;
    std::map<int, int> idx_plane2idx_intersection;

    std::cout<<"size boundary before :"<<boundary.size()<<std::endl<<std::endl;
    //1_ process neighboring planes
    for(auto it_boundary = boundary.begin(); it_boundary != boundary.end();)//++it_boundary)
    {
        idx_plane_neighbor = it_boundary->second;

        if(idx_plane_neighbor-1 != planes.size() && idx_plane_neighbor != idx_plane) // if(idx_plane == planes.size()), it means that there is an obstructing object
        {
            if(!processed_planes(idx_plane_neighbor-1, idx_plane-1)) // to not process intersection if the neighbor plane has already been processed
            {
                //Compute theoritical intersection with first pixel encountered
                if(!processed_planes(idx_plane-1, idx_plane_neighbor-1)) // to not process the pixels if the intersection has already been processed with other pixel
                {
                    processed_planes(idx_plane-1, idx_plane_neighbor-1) = true; // to process only the first pixel of this color
                    intersection inter (&planes[idx_plane-1], &planes[idx_plane_neighbor-1], delta_phi, delta_theta);
                    intersections.push_back(inter);
                    idx_plane2idx_intersection.insert(std::make_pair(idx_plane_neighbor-1, intersections_number)); //save which neighboring plane correspond to which intersection number
                    ++intersections_number;
                }

                //fill the points closest to the intersection
                // put current boundary point into the attribute "pixels" of the specific intersection
                intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].pixels.push_back(it_boundary->first);
                // put current boundary point into the attribute "points" of the specific intersection
                auto idx_found_XY2PtN = XY2PtN.find(it_boundary->first);
                if(idx_found_XY2PtN!=XY2PtN.end())
                    intersections[idx_plane2idx_intersection[idx_plane_neighbor-1]].points.push_back(idx_found_XY2PtN->second.first);
                else
                    std::cout<<"did not find this boundary point in XY2PtN"<<std::endl<<std::endl;
            }

            //erase points of the same studied plane to be sure not to detect false openings
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

    //2_ process remaining = openings
    bool first_time = true;
    for(auto it_boundary = boundary.begin(); it_boundary != boundary.end();)
    {
        idx_plane_neighbor = it_boundary->second;

        if(idx_plane_neighbor == idx_plane) // idx_plane == planes.size(), means that there is an obstructing object
        {
            if(!processed_planes(idx_plane_neighbor-1, idx_plane-1)) // to not process intersection if the neighbor plane has already been processed
            {
                //Compute theoritical intersection with first pixel encountered
                if(first_time) // to not process the pixels if the intersection has already been processed with other pixel
                {
                    intersection inter (&planes[idx_plane-1], &planes[idx_plane_neighbor-1], delta_phi, delta_theta);
                    intersections.push_back(inter);
                    idx_plane2idx_intersection.insert(std::make_pair(idx_plane_neighbor-1, intersections_number)); //save which neighboring plane correspond to which intersection number
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
                    std::cout<<"did not find this boundary point in XY2PtN"<<std::endl<<std::endl;
            }
            it_boundary = boundary.erase(it_boundary);
        }
        else
            ++it_boundary;
    }

    std::cout<<"size boundary after opening pixels :"<<boundary.size()<<std::endl<<std::endl;

    //Define the type of connection and compute limits of lines
    for (int k = 0; k<intersections.size(); ++k)
    {
        std::cout<<"Intersection number : "<<k<<std::endl<<std::endl;
        if(intersections[k].plane_ref->index == intersections[k].plane_neigh->index) // if opening or detected obstruction
        {
            if(!intersections[k].isObstruction)
                intersections[k].isOpening = true;

            if(!intersections[k].isObstruction)
                std::cout<<"It is an opening"<<std::endl<<std::endl;
            else
                std::cout<<"It is an obstruction"<<std::endl<<std::endl;

            if(intersections[k].other_points.size()>5)
            {
                intersections[k].computeTheoriticalLineFeatures(); //Hough or RANSAC
                if(intersections[k].points.size()>10) // if it is a line with sufficient points
                {
                    //Continue process of intersection
                    std::cout<<"normal : "<<intersections[k].normal.transpose()<<"   tangente : "<<intersections[k].tangente.transpose()<<"  distance : "<<intersections[k].distance<<"  initial point : "<<intersections[k].pt.transpose()<<std::endl<<std::endl;
                    processed_planes(intersections[k].plane_neigh->index, intersections[k].plane_ref->index) = true;
                    intersections[k].isLine = true;
                    //create new intersections to look for more lines in remaining points
                    intersection inter(intersections[k].plane_neigh, intersections[k].plane_ref, delta_phi, delta_theta);
                    inter.other_points = intersections[k].other_points;
                    inter.other_pixels = intersections[k].other_pixels;
                    inter.isObstruction = intersections[k].isObstruction;
                    inter.isOpening  = intersections[k].isOpening;
                    intersections.push_back(inter);
                    std::cout<<"The baby is : "<<intersections.size()-1<<std::endl<<std::endl;
                }
                else
                    remove_intersection(k);
            }
            else
               remove_intersection(k);
        }
        else             //if connection or not detected obstruction
        {
            if(intersections[k].points.size()>5)
            {
                intersections[k].isOpening = false;
                intersections[k].computeTheoriticalLineFeatures();
                std::cout<<"normal : "<<intersections[k].normal.transpose()<<"   tangente : "<<intersections[k].tangente.transpose()<<"  distance : "<<intersections[k].distance<<"  initial point : "<<intersections[k].pt.transpose()<<std::endl<<std::endl;
                processed_planes(intersections[k].plane_neigh->index, intersections[k].plane_ref->index) = true;
                intersections[k].definePlaneConnection();               // define if the intersection is a real connection between planes or an obstruction (obstructed or obstructing) and if it is a line
                if(intersections[k].isConnection)
                {
                    std::cout<<"It is a connection with plane n°"<<intersections[k].plane_ref->index<<std::endl<<std::endl;
                    intersections[k].isLine = true;
                }
                else
                {
                    std::cout<<"It is an obstruction"<<std::endl<<std::endl;
                    //create new intersections with same
                    intersection inter1 (intersections[k].plane_neigh, intersections[k].plane_neigh, delta_phi, delta_theta);
                    inter1.pixels = intersections[k].pixels;
                    inter1.other_points = intersections[k].points;
                    inter1.other_pixels = intersections[k].pixels;
                    inter1.isObstruction = true;
                    intersections.push_back(inter1);

                    intersection inter2 (intersections[k].plane_ref, intersections[k].plane_ref, delta_phi, delta_theta);
                    inter2.pixels = intersections[k].pixels;
                    inter2.other_points = intersections[k].points;
                    inter2.other_pixels = intersections[k].pixels;
                    inter2.isObstruction = true;
                    intersections.push_back(inter2);
                    intersections.erase(intersections.begin()+k);
                    --k;
                    std::cout<<"The baby obstructions are : "<<intersections.size()-2<<" and "<<intersections.size()-1<<std::endl<<std::endl;
                }
            }
            else
                remove_intersection(k);
        }
    }

    inter_remaining.pixels.clear();
    inter_remaining.points.clear();
    inter_remaining.setPlanes(&planes[idx_plane-1], &planes[idx_plane-1]);
    //3_ process remaining points (details + obstructing objects) and put it in a last intersection
    for(auto it_boundary = boundary.begin(); it_boundary != boundary.end(); ++it_boundary)
    {
        inter_remaining.pixels.push_back(it_boundary->first);
        // put current boundary point into the attribute "points" of the specific intersection
        auto idx_found_XY2PtN = XY2PtN.find(it_boundary->first);
        if(idx_found_XY2PtN!=XY2PtN.end())
            inter_remaining.points.push_back(idx_found_XY2PtN->second.first);
    }
    intersections.push_back(inter_remaining);
}

void manager::remove_intersection(int k)
{
    for(auto it_pixels = intersections[k].pixels.begin(); it_pixels != intersections[k].pixels.end(); ++it_pixels)
        inter_remaining.pixels.push_back(*it_pixels);
    for(int i = 0; i<intersections[k].points.size(); ++i)
        inter_remaining.points.push_back(intersections[k].points[i]);

    intersections.erase(intersections.begin() + k);
    --k;
}

void manager::computeTheoriticalLinesIntersections(std::vector<std::vector<intersection>> vec_openings)
{
    possible_inter3D.resize(planes.size());
    for(int plane_idx = 0; plane_idx<vec_openings.size(); ++plane_idx)
    {
        if(vec_openings[plane_idx].size()>=2)
        {
            for(int edge_idx = 0; edge_idx<vec_openings[plane_idx].size()-1; ++edge_idx)
            {
                for(int edge_idx_other = edge_idx+1; edge_idx_other<vec_openings[plane_idx].size(); ++edge_idx_other)
                {
                    if( abs(vec_openings[plane_idx][edge_idx].tangente2D.dot(vec_openings[plane_idx][edge_idx_other].tangente2D)) < 0.999)
                    {
                        Eigen::Vector2d Dn = vec_openings[plane_idx][edge_idx].distance2D * vec_openings[plane_idx][edge_idx].normal2D - vec_openings[plane_idx][edge_idx_other].distance2D * vec_openings[plane_idx][edge_idx_other].normal2D;
                        float X;

                        float ax = pow(vec_openings[plane_idx][edge_idx_other].tangente2D(0),2) - pow(vec_openings[plane_idx][edge_idx].tangente2D(0),2);
                        float bx = -2*vec_openings[plane_idx][edge_idx].tangente2D(0)*Dn(0);
                        float cx = (pow(vec_openings[plane_idx][edge_idx].distance2D,2)-pow(vec_openings[plane_idx][edge_idx_other].distance2D,2))*pow(vec_openings[plane_idx][edge_idx_other].tangente2D(0),2)-pow(Dn(0),2);

                        float ay = pow(vec_openings[plane_idx][edge_idx_other].tangente2D(1),2) - pow(vec_openings[plane_idx][edge_idx].tangente2D(1),2);
                        float by = -2*vec_openings[plane_idx][edge_idx].tangente2D(1)*Dn(1);
                        float cy = (pow(vec_openings[plane_idx][edge_idx].distance2D,2)-pow(vec_openings[plane_idx][edge_idx_other].distance2D,2))*pow(vec_openings[plane_idx][edge_idx_other].tangente2D(1),2)-pow(Dn(1),2);

                        X = (cx*(ay/ax)-cy)/(by-bx*(ay/ax));

                        Eigen::Vector2d pt_intersect2D = vec_openings[plane_idx][edge_idx].distance2D * vec_openings[plane_idx][edge_idx].normal2D + X * vec_openings[plane_idx][edge_idx].tangente2D;
                        Eigen::Vector3d pt_intersect3D = {pt_intersect2D(0), pt_intersect2D(1), 0};
                        pt_intersect3D = vec_openings[plane_idx][edge_idx].rot_inv.linear() * vec_openings[plane_idx][edge_idx].trans_z_inv*pt_intersect3D;

                        possible_inter3D[plane_idx].push_back(pt_intersect3D);
                    }
                }
            }
        }
        std::cout<<"for plane n°"<<plane_idx<<"    number of intersections of lines points : "<<possible_inter3D[plane_idx].size()<<std::endl<<std::endl;
    }
}

void manager::computeCorners()
{
    corners.resize(0);
    for(int i = 0; i<planes.size()-1; ++i)
    {
        for(int j = i+1; j<planes.size(); ++j)
        {
            if(planes[i].connected.find(planes[j].index) != planes[i].connected.end())
            {
                for(int k = j+1; k<planes.size(); ++k)
                {
                    if(planes[k].connected.find(planes[i].index) != planes[k].connected.end() && planes[k].connected.find(planes[j].index) != planes[k].connected.end())
                    {
                        corner c;
                        c.setPlanes(&planes[i], &planes[j], &planes[k]);
                        c.computePoint();
                        corners.push_back(c);
                    }
                }
            }
        }
    }
}

void manager::searchClusters2(double thresh_plane_belonging) // create the object regions in which there is a vector of points with their associated mode
{
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> processed_seed = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> processed_growing = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();
    pcl::PointCloud<pcl::PointXYZ> pc_clus;

    regions.resize(0);

    for(int planes_it = 0; planes_it<planes.size(); ++planes_it)
    {
        auto it = planes[planes_it].seed;
        pc_clus.points.resize(0);
        processed_growing = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();
        std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator> region; // vector of iterator
        region.push_back(it);

        for(int i = 0; i<region.size(); ++i)
        {
            std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3d, Eigen::Vector3d>>::iterator> neighbors = getNeighbors(region[i], pixel_radius_for_growing);
            bool all_neigh_with_non_coherent_normal = true;
            for(int k = 0; k<neighbors.size(); ++k)
            {
                if(abs(neighbors[k]->second.second.dot(planes[planes_it].normal))>normals_similarity_to_add_neighbors)
                {
                    all_neigh_with_non_coherent_normal = false;
                    break;
                }
            }
            if(!all_neigh_with_non_coherent_normal)
            {
                for(int k = 0; k<neighbors.size(); ++k)
                {
                    if(!processed_growing(neighbors[k]->first.first, neighbors[k]->first.second)) //in order not to process twice one neighbor
                    {
                        if(abs((neighbors[k]->second.first-it->second.first).dot(planes[planes_it].normal)) < thresh_plane_belonging)// && acos(neighbors[k]->second.dot(it->second)/(neighbors[k]->second.norm()*it->second.norm())) < thresh_angle )
                        {
                            processed_growing(neighbors[k]->first.first, neighbors[k]->first.second) = true;
                            region.push_back(neighbors[k]);

                            pcl::PointXYZ pt_clus;
                            pt_clus.x = neighbors[k]->second.first(0);
                            pt_clus.y = neighbors[k]->second.first(1);
                            pt_clus.z = neighbors[k]->second.first(2);
                            pc_clus.points.push_back(pt_clus);
                        }
                    }
                }
            }
        }

        if(region.size()>min_number_of_pixels)
        {
            processed_seed = processed_seed || processed_growing;
            Eigen::Vector3d cluster_normal = planes[planes_it].normal;
            double cluster_distance = abs(it->second.first.dot(planes[planes_it].normal));
            cluster_normal /= cluster_normal.norm();
            cluster_normal *= cluster_distance;
            regions.push_back(std::make_pair(region, cluster_normal));
            std::cout<<"New region number: "<< regions.size()-1 << std::endl<< "seed : "<<it->first.first<<" "<<it->first.second<<"-->"<<it->second.first.transpose()<<std::endl<<"normal : "<<it->second.second.transpose()<<std::endl<<"distance : "<<cluster_distance<<std::endl<<"number of points : "<<region.size()<<std::endl<<std::endl;
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
        }
    }
}
