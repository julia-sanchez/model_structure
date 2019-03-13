image_manager::image_manager(std::map<std::pair<float, float>, std::pair<Eigen::Vector3f, Eigen::Vector3f>> p, int Nr, int Nc, float pa, float ta, int mc)
{
    planes = p;
    Nrow = Nr;
    Ncol = Nc;
    phi_app = pa;
    theta_app = ta;
    max_col = mc;
}

void image_manager::makeImage()
{
    float maxi_dist = -100000;

    for(auto it = planes.begin(); it!=planes.end(); ++it)
    {
        float dist = abs(it->second.second.dot(it->second.first));
        if(maxi_dist<dist)
            maxi_dist = dist;
    }

    for(auto it = planes.begin(); it!=planes.end(); ++it)
    {
        Eigen::Vector3i rgb;
        Eigen::Vector3f temp;
        temp = (it->second.second + Eigen::Vector3f::Ones())/2;
        float dist = abs(it->second.second.dot(it->second.first));
        rgb = ( max_col*(1 - dist/maxi_dist)*temp ) .cast<int>();

        rgb = (rgb.cast<float>() * (1-dist/maxi_dist)).cast<int>();
        int X = (int)(it->first.first*Nrow/(phi_app+0.0001));
        int Y = (int)(it->first.second*Ncol/(theta_app+0.0001));
        imageRGB.insert(std::make_pair(std::make_pair(X,Y), rgb));
        imageRGBf.insert(std::make_pair(std::make_pair(X,Y), it->second));
    }

    save_image_ppm("ioi/testy", "", imageRGB, Nrow, Ncol, max_col);
}

void image_manager::CreateImageClusterized()
{
    float maxi_dist = -100000;
    for(auto it = regions.begin(); it!=regions.end(); ++it)
    {
        if(maxi_dist<it->second.norm())
            maxi_dist = it->second.norm();
    }

    for(auto it = imageRGB_clusterizedf_clean.begin(); it!=imageRGB_clusterizedf_clean.end(); ++it) // for each cluster
    {
        Eigen::Vector3i rgb;
        Eigen::Vector3f temp;
        temp = (it->second/it->second.norm() + Eigen::Vector3f::Ones())/2;
        rgb = ( max_col*(1 - it->second.norm()/maxi_dist)*temp ) .cast<int>();
        imageRGB_clusterized.insert(std::make_pair(it->first, rgb));
    }

    save_image_ppm("ioi/clusterized", "", imageRGB_clusterized, Nrow, Ncol, max_col);
}

void image_manager::CleanClusters()
{
    for(auto it = regions.begin(); it!=regions.end(); ++it) // for each cluster
    {
        for(int i = 0; i<it->first.size(); ++i)
            imageRGB_clusterizedf.insert(std::make_pair(it->first[i]->first, it->second));
    }
    int n =0;
    for(auto it = imageRGB_clusterizedf.begin(); it!=imageRGB_clusterizedf.end(); ++it) // for each cluster
    {
        if(imageRGB_clusterizedf.count(it->first)>1)
        {
            ++n;
            auto it_multiple_limits = imageRGB_clusterizedf.equal_range(it->first);
            std::map<float, std::multimap<std::pair<int,int>, Eigen::Vector3f>::iterator> dotty;
            for(auto it_multiple = it_multiple_limits.first; it_multiple != it_multiple_limits.second; ++it_multiple)
            {
                auto it_in_image = imageRGBf.find(it_multiple->first);              // search in imageRGBf the initial normal of point-> { (X,Y), [(px,py,pz)(nx,ny,nz)]}
                float dot = (it_multiple->second/it_multiple->second.norm()).dot(it_in_image->second.second); // compute scalar product between real normal and region normal
                dotty.insert(std::make_pair(dot, it_multiple));                     // keep result in map to order scalar product and select the region normal which is more similar
            }
            auto it_map = dotty.end();
            --it_map;                                                               // it_map is the ptr of [ scalar_product, idx] with idx corresponding to the more similar region normal
            imageRGB_clusterizedf_clean.insert(std::make_pair(it->first, it_map->second->second));
        }
        else
            imageRGB_clusterizedf_clean.insert(std::make_pair(it->first, it->second));
    }

    std::cout<<"Number of points which have multiple modes : "<<n<<std::endl<<std::endl;
}

////Eigen::MatrixXi image_manager::getImage()
////{
////    imagei = Eigen::MatrixXi::Zero(Nrow, Ncol);

////    int k = 0;
////    std::multimap<std::pair<float, float>, std::pair<Eigen::Vector3f, float>>::iterator it_best;
////    for(auto it = clusters.begin(); it != clusters.end(); ++it)
////    {
////        //if a point belongs to two clusters (edge) -> choose the closest one
////        int modes_nbr = clusters.count(it->first);
////        if(modes_nbr>1)
////        {
////            auto modes_it = clusters.equal_range(it->first);
////            k = 0;
////            float temp = 10000;
////            for( auto it_clus = modes_it.first; it_clus != modes_it.second; ++it_clus )
////            {
////                float dist = abs(planes[it->first].second-it_clus->second.second);
////                if(temp>dist)
////                {
////                    temp = dist;
////                    clusters.erase(it_best);
////                    it_best = it_clus;
////                }
////                else
////                {
////                    clusters.erase(it_clus);
////                }
////                ++k;
////            }
////        }

////        int col;
////        k=0;

////        for(auto it_clus = clusters_sizes.begin(); it_clus!=clusters_sizes.end(); ++it_clus)
////        {
////            ++k;
////            if(it_clus->second == it->second)
////            {
////                col = (int)(k*max_col/clusters_sizes.size());
////                break;
////            }
////        }

////        int X = (int)(it->first.first/(phi_app/Nrow));
////        int Y = (int)(it->first.second/(theta_app/Ncol));
////        imagei(Y,X) = col;

////    }
////    return imagei;
////}

////void image_manager::binarize(std::pair<Eigen::Vector3f, float> plane)
////{
////    image_bool.resize(Nrow, Ncol);
////    image_other_bool.resize(Nrow, Ncol);

////    for(auto it = clusters.begin(); it != clusters.end(); ++it)
////    {
////        //image_bool
////        if( it->second==plane )
////            image_bool(it->first.first,it->first.second) = true;
////        else
////            image_bool(it->first.first,it->first.second) = false;

////        //image_other_bool
////        if( it->second!=plane )
////            image_other_bool(it->first.first,it->first.second) = true;
////        else
////            image_other_bool(it->first.first,it->first.second) = false;
////    }
////}

////Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_manager::mynot(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_in)
////{
////    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_out = image_in;
////    for (int i = 0; i < Nrow; ++i)
////    {
////        for (int j = 0; j< Ncol; ++j)
////        {
////            if(image_in(i,j))
////                image_out(i,j) = false;
////            else
////                image_out(i,j) = true;
////        }
////    }

////    return image_out;
////}

////Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_manager::morpho(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_in, Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> SE)
////{
////    //Structure element is convolved with the image with an OR process

////    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_out = image_in;
////    int rad = (SE.rows()-1)/2;
////    for (int i = rad; i< Nrow-rad; ++i)
////    {
////        for (int j = rad; j< Ncol-rad; ++j)
////        {
////            if(image_in(i,j))
////                image_out.block<5,5>(i-rad,j-rad) = image_out.block<5,5>(i-rad,j-rad) || SE;
////        }
////    }
////    return image_out;
////}

////void image_manager::closure(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> SE)
////{
////    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_dilated = morpho(image_bool, SE); //dilate
////    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_closed = morpho(mynot(image_dilated), SE); //not(erode)
////    image_morpho = mynot(image_closed); // not
////    image_morpho = image_morpho - (image_other_bool && image_morpho); // remove other walls
////}

////Eigen::MatrixXi image_manager::getBinarized(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_in)
////{
////    Eigen::MatrixXi image_out = Eigen::MatrixXi::Zero(Nrow, Ncol);
////    for (int i = 0; i<Nrow; ++i)
////    {
////        for (int j = 0; j< Ncol; ++j)
////        {
////            if(image_in(i,j))
////                image_out(i,j) = max_col;
////        }
////    }

////    return image_out;
////}

////void image_manager::computeBoundaries(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> SE)
////{
////    std::cout<<"Start computing boundaries of all planes of one direction"<<std::endl<<std::endl;

////    image_grad.resize(Nrow,Ncol);
////    int rad = (SE.rows()-1)/2;

////    for (int i = rad; i<Nrow-rad; ++i)
////    {
////        for (int j = rad; j<Ncol-rad; ++j)
////        {
////            if(image_bool(i,j))
////            {
////                int sum = (image_morpho.block<3,3>(i-rad,j-rad) && SE).cast<int>().sum();
////                image_grad(i,j) = (sum!=SE.rows()*SE.rows());
////            }
////            else
////                image_grad(i,j) = false;
////        }
////    }

////    std::cout<<"Stop computing boundaries of all planes of one direction"<<std::endl<<std::endl;
////}

////void image_manager::searchClusters(float thresh_plane_belonging, float thresh_angle)
////{
////    float mean_dist;
////    Eigen::Vector3f mean_norm;
////    buildTree();
////    int n_neigh = 30;
////    float radius = 0.001;

////    for (auto it = planes.begin(); it!=planes.end(); ++it)
////    {
////        if(clusters.find(it->first)==clusters.end()) // if the new seed does not belong to any existing cluster
////        {
////            mean_dist = it->second.second;
////            mean_norm = it->second.first;
////            std::map<std::pair<float, float>, int> region;
////            std::vector<std::map<std::pair<float, float>, std::pair<Eigen::Vector3f, float>>::iterator> list; // element that are in region
////            std::vector<float> dis (n_neigh);
////            std::vector<int> neigh (n_neigh);
////            list.push_back(it);

////            auto t1 = std::chrono::high_resolution_clock::now();
////            for(int n = 0; n<list.size(); ++n)
////            {
////                SearchFLANNTree(tree, list[n]->first, neigh, dis, radius);
////                int i =1;
////                int k=0;
////                while (neigh[i] !=-1)
////                {
////                    std::map<std::pair<float, float>, std::pair<Eigen::Vector3f, float>>::iterator it_neigh = planes.begin();
////                    for (int j = 0; j<neigh[i]; ++j)
////                        ++it_neigh;

////                    if(region.find(it_neigh->first)==region.end()) // add it only if the new neighbor does not belong to the list already
////                    {
////                        float dot = it_neigh->second.first.dot(it->second.first);
////                        if(dot>1)
////                            dot = 1;

////                        if(abs( it_neigh->second.second - it->second.second) < thresh_plane_belonging && acos(dot) < thresh_angle)
////                        {
////                            list.push_back(it_neigh);
////                            mean_dist += it_neigh->second.second;
////                            mean_norm += it_neigh->second.first;
////                            ++k;
////                            region.insert(std::make_pair(it_neigh->first , k));
////                        }
////                    }
////                    ++i;
////                }
////            }
////            auto t2 = std::chrono::high_resolution_clock::now();
////            mean_dist /= list.size();
////            mean_norm /= list.size();
////            std::cout<<"NEW REGION"<<std::endl<<"mode:"<<mean_norm<<std::endl<<mean_dist<<std::endl<<"number of points :"<<list.size()<<std::endl<<std::endl;
////            std::cout<<"Time to find the region :" <<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<std::endl<<std::endl;
////            clusters_sizes.insert(std::make_pair(list.size(), it->second));
////            for(auto it_reg = list.begin(); it_reg!=list.end(); ++it_reg)
////                clusters.insert(std::make_pair((*it_reg)->first, std::make_pair(mean_norm, mean_dist)));
////        }
////    }
////}



bool image_manager::selectSeed(std::map<std::pair<int,int>, std::pair<Eigen::Vector3f, Eigen::Vector3f>>::iterator it_imageRGBf)
{
    bool ok_pt = false;
    
    std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3f, Eigen::Vector3f>>::iterator> neighbors = getNeighbors(it_imageRGBf, );
    return ok_pt;
}






void image_manager::searchClustersImage(float thresh_plane_belonging)
{
    Eigen::Vector3f mean_normal;
    float mean_distance;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> processed_seed = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> processed_growing = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();


    for(auto it = imageRGBf.begin(); it != imageRGBf.end(); ++it)
    {
        if(!processed_seed(it->first.first, it->first.second)) // in order not to process one point which has already been selected in a cluster
        {
            if(selectSeed(it))
            {
                processed_growing = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();
                mean_normal= it->second.second;
                mean_distance= it->second.first.dot(it->second.second);
                std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3f, Eigen::Vector3f>>::iterator> region; // vector of iterator
                region.push_back(it);
    
                for(int i = 0; i<region.size(); ++i)
                {
                    std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3f, Eigen::Vector3f>>::iterator> neighbors = getNeighbors(region[i], 30);
                    for(int k = 0; k<neighbors.size(); ++k)
                    {
                        if(!processed_growing(neighbors[k]->first.first, neighbors[k]->first.second)) //in order not to process twice one neighbor
                        {
    //                        std::cout<<neighbors[k]->second<<std::endl<<std::endl;
                            if(abs((neighbors[k]->second.first-it->second.first).dot(it->second.second)) < thresh_plane_belonging)// && acos(neighbors[k]->second.dot(it->second)/(neighbors[k]->second.norm()*it->second.norm())) < thresh_angle )
                            {
                                processed_seed(neighbors[k]->first.first, neighbors[k]->first.second) = true;
                                processed_growing(neighbors[k]->first.first, neighbors[k]->first.second) = true;
                                region.push_back(neighbors[k]);
                                mean_normal += neighbors[k]->second.second;
                                mean_distance += neighbors[k]->second.first.dot(neighbors[k]->second.second);
                            }
                        }
                    }
                }
            }

            regions.push_back(std::make_pair(region, mean_normal));
            mean_normal /= regions.size();
            mean_distance /= regions.size();
            mean_normal /= mean_normal.norm();
            mean_normal *= mean_distance;
            std::cout<<"New region : "<<std::endl<<"normal : "<<(mean_normal/mean_distance).transpose()<<std::endl<<"distance : "<<mean_distance<<std::endl<<"number of points : "<<region.size()<<std::endl<<std::endl;
        }
    }
}

std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3f, Eigen::Vector3f>>::iterator> image_manager::getNeighbors(std::map<std::pair<int,int>, std::pair<Eigen::Vector3f, Eigen::Vector3f>>::iterator it, int rad)
{
    std::vector<std::map<std::pair<int,int>, std::pair<Eigen::Vector3f, Eigen::Vector3f>>::iterator> neighbors;
    int initx;
    int finitx;
    if(it->first.first-rad>0)
        initx = it->first.first-rad;
    else
        initx = 0;

    if(it->first.first+rad<Nrow)
        finitx = it->first.first+rad;
    else
        finitx = Nrow;

    int inity;
    int finity;
    if(it->first.second-rad>0)
        inity = it->first.second-rad;
    else
        inity = 0;

    if(it->first.second+rad<Ncol)
        finity = it->first.second+rad;
    else
        finity = Ncol;

    for(int i = initx; i < finitx; ++i)
    {
        for(int j = inity; j< finity; ++j)
        {
            auto idx = imageRGBf.find(std::make_pair(i,j));
            if(idx != imageRGBf.end())
                neighbors.push_back(idx);
        }
    }
    return neighbors;
}

////void image_manager::SearchFLANNTree(flann::Index<flann::L2<float>>* index,
////                            std::pair<float,float> input,
////                            std::vector<int>& indices,
////                            std::vector<float>& dists,
////                            float radius)
////{
////    std::vector<float> query(2);
////    query[0] = input.first;
////    query[1] = input.second;
////    flann::Matrix<float> query_mat(&query[0], 1, 2);

////    flann::Matrix<int> indices_mat(&indices[0], 1, indices.size());
////    flann::Matrix<float> dists_mat(&dists[0], 1, indices.size());

////    index->radiusSearch(query_mat, indices_mat, dists_mat, radius, flann::SearchParams(128));
////}

////void image_manager::buildTree()
////{
////    std::vector<float> pc(planes.size() * 2);
////    flann::Matrix<float> flann_mat(&pc[0], planes.size(), 2);

////    int n = 0;

////    for (auto it =planes.begin(); it != planes.end(); ++it)
////    {
////        pc[n] = it->first.first;
////        ++n;
////        pc[n] = it->first.second;
////        ++n;
////    }

////    tree = new flann::Index<flann::L2<float>>(flann_mat, flann::KDTreeSingleIndexParams(15));
////    tree->buildIndex();
////}

