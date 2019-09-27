manager2D::manager2D()
{
}

void manager2D::setImageClusterized(Eigen::MatrixXi p, std::pair<int,int> lt)
{
    image_clusterized = p;
    Nrow = image_clusterized.rows();
    Ncol = image_clusterized.cols();
    image_grad = Eigen::MatrixXi::Zero(Nrow, Ncol);
    lim_theta = lt;
    image_bool = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>::Zero(Nrow, Ncol);
    image_morpho = Eigen::MatrixXi::Zero(Nrow, Ncol);
    image_other_bool = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>::Zero(Nrow, Ncol);
}

void manager2D::binarize(int p) // p is the index of the plane
{
    image_bool = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();
    image_other_bool = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();

    image_clusterized_idx = p+1;

    for(int i = 0; i< Nrow; ++i)
    {
        for(int j = 0; j< Ncol; ++j)
        {
            //image_bool
            if( image_clusterized(i,j)==(p+1))
            {
                image_bool(i,j) = true;
                image_other_bool(i,j) = false;
            }
            else if(image_clusterized(i,j)>0)
            {
                image_other_bool(i,j) = true;
                image_bool(i,j) = false;
            }
            else
            {
                image_bool(i,j) = false;
                image_other_bool(i,j) = false;
            }
        }
    }
}

Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> manager2D::mynot(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_in)
{
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_out = image_in;
    for (int i = 0; i < Nrow; ++i)
    {
        for (int j = 0; j< Ncol; ++j)
        {
            if(image_in(i,j))
                image_out(i,j) = false;
            else
                image_out(i,j) = true;
        }
    }

    return image_out;
}

void manager2D::morpho(int rad)
{
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_dilated = dilate(image_bool, rad); // dilate
    image_morpho = (image_dilated - (image_other_bool && image_dilated)).cast<int>();
//    segmentBiggest();
}

Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> manager2D::dilate(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_in, int rad)
{
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_out = image_in;
    reference.clear();
    //upper part
    for (int j = lim_theta.first+rad; j< lim_theta.second-rad; ++j)
    {
        for(int i = 0; i < rad; ++i)
        {
            if(image_in(i,j))
            {
                int min_ki = std::max(i-rad, 0);
                int max_ki = std::min(i+rad, Nrow-1);
                int min_kj = std::max(j-rad, lim_theta.first);
                int max_kj = std::min(j+rad, lim_theta.second);

                for(int ki = min_ki; ki<=max_ki; ++ki)
                {
                    for(int kj = min_kj; kj<=max_kj; ++kj)
                    {
                        if(!image_in(ki,kj))
                        {
                            image_out(ki,kj) = true;
                            reference.insert( std::make_pair(std::make_pair(ki,kj) , std::make_pair(i,j)) );
                        }
                    }
                }

                min_ki = Nrow - rad + i;
                max_ki = Nrow-1;
                min_kj = std::max(j-rad, lim_theta.first);
                max_kj = std::min(j+rad, lim_theta.second);

                for(int ki = min_ki; ki<=max_ki; ++ki)
                {
                    for(int kj = min_kj; kj<=max_kj; ++kj)
                    {
                        if(!image_in(ki,kj))
                        {
                            image_out(ki,kj) = true;
                            reference.insert( std::make_pair(std::make_pair(ki,kj) , std::make_pair(i,j)) );
                        }
                    }
                }
            }
        }
    }

    //bottom part

    for (int j = lim_theta.first+rad; j<= lim_theta.second-rad; ++j)
    {
        for(int i = Nrow-rad; i < Nrow; ++i)
        {
            if(image_in(i,j))
            {
                int min_ki = std::max(i-rad, 0);
                int max_ki = std::min(i+rad, Nrow-1);
                int min_kj = std::max(j-rad, lim_theta.first);
                int max_kj = std::min(j+rad, lim_theta.second);

                for(int ki = min_ki; ki<=max_ki; ++ki)
                {
                    for(int kj = min_kj; kj<=max_kj; ++kj)
                    {
                        if(!image_in(ki,kj))
                        {
                            image_out(ki,kj) = true;
                            reference.insert( std::make_pair(std::make_pair(ki,kj) , std::make_pair(i,j)) );
                        }
                    }
                }

                min_ki = 0;
                max_ki = rad + (Nrow - 1 - i);
                min_kj = std::max(j-rad, lim_theta.first);
                max_kj = std::min(j+rad, lim_theta.second);

                for(int ki = min_ki; ki<=max_ki; ++ki)
                {
                    for(int kj = min_kj; kj<=max_kj; ++kj)
                    {
                        if(!image_in(ki,kj))
                        {
                            image_out(ki,kj) = true;
                            reference.insert( std::make_pair(std::make_pair(ki,kj) , std::make_pair(i,j)) );
                        }
                    }
                }
            }
        }
    }

    //remaining part

    for (int i = rad; i< Nrow-rad; ++i)
    {
        for (int j = lim_theta.first+rad; j <= lim_theta.second-rad; ++j)
        {
            if(image_in(i,j))
            {
                int min_ki = std::max(i-rad, 0);
                int max_ki = std::min(i+rad, Nrow-1);
                int min_kj = std::max(j-rad, lim_theta.first);
                int max_kj = std::min(j+rad, lim_theta.second);

                for(int ki = min_ki; ki<=max_ki; ++ki)
                {
                    for(int kj = min_kj; kj<=max_kj; ++kj)
                    {
                        if(!image_in(ki,kj))
                        {
                            image_out(ki,kj) = true;
                            reference.insert( std::make_pair(std::make_pair(ki,kj) , std::make_pair(i,j)) );
                        }
                    }
                }
            }
        }
    }
    return image_out;
}

void manager2D::segmentBiggest()
{
   Eigen::MatrixXi region  = Eigen::MatrixXi::Zero(Nrow, Ncol);
   Eigen::MatrixXi region_temp = Eigen::MatrixXi::Zero(Nrow, Ncol);
   Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> processed = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>::Zero(Nrow, Ncol);
   std::set<std::pair<int,int>> list_pixels;
    //croissance de région à partir du premier pixel blanc
    //si la région est plus petite elle dégage
    for(int i = 0; i < Nrow; ++i)
    {
        for(int j = lim_theta.first; j < lim_theta.second; ++j)
        {
            if(image_morpho(i,j)>0 && !processed(i,j)) // seed
            {
                list_pixels.insert(std::make_pair(i,j));
                while (list_pixels.size() > 0)
                {
                    region_temp(list_pixels.begin()->first,list_pixels.begin()->second) = 1;
                    processed(list_pixels.begin()->first,list_pixels.begin()->second) = true;
                    std::vector<std::pair<int,int>> neighbors = getNeighbors(image_morpho, std::make_pair(list_pixels.begin()->first,list_pixels.begin()->second), 1, processed);
                    list_pixels.erase(list_pixels.begin());
                    for(int k = 0; k < neighbors.size(); ++k)
                        list_pixels.insert(neighbors[k]);
                }

                if(region_temp.sum()>region.sum())
                    region = region_temp;

                region_temp = Eigen::MatrixXi::Zero(Nrow, Ncol);
                list_pixels.clear();
            }
        }
    }

    image_morpho = region;
}

void manager2D::removeLittleObjects(int size_min)
{
   Eigen::MatrixXi region  = Eigen::MatrixXi::Zero(Nrow, Ncol);
   Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> processed = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>::Zero(Nrow, Ncol);
   std::set<std::pair<int,int>> list_pixels;
    //croissance de région à partir du premier pixel blanc
    //si la région est trop petite, elle dégage

   int max_value =  image_clusterized.maxCoeff();

    for(int i = 0; i < Nrow; ++i)
    {
        for(int j = lim_theta.first; j < lim_theta.second; ++j)
        {
            if(image_clusterized(i,j) == max_value && !processed(i,j)) // seed
            {
                list_pixels.insert(std::make_pair(i,j));
                while (list_pixels.size() > 0)
                {
                    region(list_pixels.begin()->first, list_pixels.begin()->second) = 1;
                    processed(list_pixels.begin()->first, list_pixels.begin()->second) = true;
                    std::vector<std::pair<int,int>> neighbors = getNeighbors(image_clusterized, std::make_pair(list_pixels.begin()->first,list_pixels.begin()->second), 1, processed);
                    list_pixels.erase(list_pixels.begin());
                    for(int k = 0; k < neighbors.size(); ++k)
                        list_pixels.insert(neighbors[k]);
                }

                if(region.sum()<=size_min)
                {
                    for(int region_i = 0; region_i < Nrow; ++region_i)
                    {
                        for(int region_j = 0; region_j < Ncol; ++region_j)
                        {
                            if(region(region_i, region_j))
                            {
                                int rad = 1;
                                bool found = false;
                                while (!found || rad>=10)
                                {
                                    int min_ki = std::max(region_i-rad, 0);
                                    int max_ki = std::min(region_i+rad, Nrow-1);
                                    int min_kj = std::max(region_j-rad, lim_theta.first);
                                    int max_kj = std::min(region_j+rad, lim_theta.second);

                                    for(int ri = min_ki; ri <= max_ki; ++ri)
                                    {
                                        for(int rj = min_kj; rj <= max_kj; ++rj)
                                        {
                                            if(!region(ri, rj))
                                            {
                                                image_clusterized(region_i, region_j) = image_clusterized(ri,rj);
                                                found = true;
                                                break;
                                            }
                                        }
                                        if(found)
                                            break;
                                    }
                                    ++rad;
                                }
                            }
                        }
                    }
                }

                region = Eigen::MatrixXi::Zero(Nrow, Ncol);
                list_pixels.clear();
            }
        }
    }
}

std::vector<std::pair<int,int>> manager2D::getNeighbors(Eigen::MatrixXi& image, std::pair<int,int> pixel, int rad, Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& processed)
{
    std::vector<std::pair<int,int>> neighbors;
    int i = pixel.first;
    int j = pixel.second;
    int min_ki = std::max(i-rad, 0);
    int max_ki = std::min(i+rad, Nrow-1);
    int min_kj = std::max(j-rad, lim_theta.first);
    int max_kj = std::min(j+rad, lim_theta.second);

    int value = image(pixel.first, pixel.second);

    for(int ki = min_ki; ki<=max_ki; ++ki)
    {
        for(int kj = min_kj; kj<=max_kj; ++kj)
        {
            if(!(i == ki && j == kj))
            {
                if(image(ki,kj) == value && !processed(ki,kj))
                    neighbors.push_back(std::make_pair(ki ,kj));
            }
        }
    }

    if(i<rad)
    {
        min_ki = Nrow - rad + i;
        max_ki = Nrow-1;
        min_kj = std::max(j-rad, lim_theta.first);
        max_kj = std::min(j+rad, lim_theta.second);

        for(int ki = min_ki; ki<=max_ki; ++ki)
        {
            for(int kj = min_kj; kj<=max_kj; ++kj)
            {
                if(!(i == ki && j == kj))
                {
                    if(image(ki,kj) == value && !processed(ki,kj))
                        neighbors.push_back(std::make_pair(ki ,kj));
                }
            }
        }
    }

    if(i+rad>Nrow)
    {
        min_ki = 0;
        max_ki = rad + (Nrow - 1 - i);
        min_kj = std::max(j-rad, lim_theta.first);
        max_kj = std::min(j+rad, lim_theta.second);

        for(int ki = min_ki; ki<=max_ki; ++ki)
        {
            for(int kj = min_kj; kj<=max_kj; ++kj)
            {
                if(!(i == ki && j == kj))
                {
                    if(image(ki,kj)>0  && !processed(ki,kj))
                        neighbors.push_back(std::make_pair(ki ,kj));
                }
            }
        }
    }

    return neighbors;
}

void manager2D::computeBoundaries(Eigen::MatrixXi& all_boundaries)
{
    int rad = 1;
    neighborPix2currentPix.clear();
    XY2planeIdx_boundary.clear();
    int max_value = image_clusterized.maxCoeff();
    for (int i = 0; i<Nrow; ++i)
    {
        for (int j = lim_theta.first+rad; j<=lim_theta.second-rad; ++j)
        {
            if(image_morpho(i,j)>0 && all_boundaries(i,j) == 0 )
            {
                int sum;
                if(i==0)
                {
                    sum = 0;
                    sum += image_morpho(Nrow-1,j-1);
                    sum += image_morpho(Nrow-1,j);
                    sum += image_morpho(Nrow-1,j+1);
                    sum += image_morpho(0,j-1);
                    sum += image_morpho(0,j);
                    sum += image_morpho(0,j+1);
                    sum += image_morpho(1,j-1);
                    sum += image_morpho(1,j);
                    sum += image_morpho(1,j+1);
                }
                else if (i==Nrow-1)
                {
                    sum = 0;
                    sum += image_morpho(Nrow-2,j-1);
                    sum += image_morpho(Nrow-2,j);
                    sum += image_morpho(Nrow-2,j+1);
                    sum += image_morpho(Nrow-1,j-1);
                    sum += image_morpho(Nrow-1,j);
                    sum += image_morpho(Nrow-1,j+1);
                    sum += image_morpho(0,j-1);
                    sum += image_morpho(0,j);
                    sum += image_morpho(0,j+1);
                }
                else
                {
                    sum = 0;
                    sum += image_morpho(i-1,j-1);
                    sum += image_morpho(i-1,j);
                    sum += image_morpho(i-1,j+1);
                    sum += image_morpho(i,j-1);
                    sum += image_morpho(i,j);
                    sum += image_morpho(i,j+1);
                    sum += image_morpho(i+1,j-1);
                    sum += image_morpho(i+1,j);
                    sum += image_morpho(i+1,j+1);
                }

                if(sum!=(2*rad+1)*(2*rad+1)) // at least one neighbor is different from (i,j)
                {
                    bool isOpenning = true;
                    //search who is the neighbor which leads to this boundary

                    int min_ki = std::max(i-rad, 0);
                    int max_ki = std::min(i+rad, Nrow-1);
                    int min_kj = std::max(j-rad, lim_theta.first+2);
                    int max_kj = std::min(j+rad, lim_theta.second-2);

                    int selected = image_clusterized_idx;

                    for(int ki = min_ki; ki<=max_ki; ++ki)
                    {
                        for(int kj = min_kj; kj<=max_kj; ++kj)
                        {
                            bool is_different = image_clusterized(ki, kj) != image_clusterized_idx;
                            bool is_not_black = image_clusterized(ki, kj)>0;
                            bool is_not_object = image_clusterized(ki, kj) != max_value;
                            if(is_different && is_not_black && is_not_object) //if pixels belong to different clusters and neighbor not black and neighbor not object
                            {
                                if(all_boundaries(ki, kj) == 0 )
                                {
                                    if(image_clusterized(ki, kj) == selected || selected==image_clusterized_idx)
                                    {
                                        selected = image_clusterized(ki, kj);
                                        if(image_clusterized(i, j) != 0)
                                        {
                                             XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(i, j), selected));
                                            neighborPix2currentPix.insert(std::make_pair(std::make_pair(ki, kj), std::make_pair(i, j)));
                                            all_boundaries(i,j) = 1;
                                        }
                                        else
                                        {
                                            auto pair_it_reference = reference.equal_range(std::make_pair(i,j));
                                            for(auto it_reference = pair_it_reference.first; it_reference!=pair_it_reference.second; ++it_reference)
                                            {
                                                if(all_boundaries(it_reference->second.first, it_reference->second.second) == 0 )
                                                {
                                                    XY2planeIdx_boundary.insert(std::make_pair(it_reference->second, selected));
                                                    neighborPix2currentPix.insert(std::make_pair(std::make_pair(ki, kj), it_reference->second));
                                                    all_boundaries(it_reference->second.first, it_reference->second.second) = 1;
                                                }
                                            }
                                        }
                                        XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(ki, kj), selected));
                                    }
                                }
                                isOpenning = false;
                            }
                        }
                    }

                    //-------------------------------------------------------------------------------------------------------------------------------------------------------
                    //search object neighbors

                    if(isOpenning)
                    {
                        for(int ki = min_ki; ki<=max_ki; ++ki)
                        {
                            for(int kj = min_kj; kj<=max_kj; ++kj)
                            {
                                bool is_different = image_clusterized(ki, kj) != image_clusterized_idx;
                                bool is_not_black = image_clusterized(ki, kj)>0;
                                if(is_different && is_not_black) //if pixels belong to different clusters and neighbor not black and neighbor not object
                                {
                                    if(all_boundaries(ki, kj) == 0 )
                                    {
                                        if(image_clusterized(ki, kj) == selected || selected==image_clusterized_idx)
                                        {
                                            selected = image_clusterized(ki, kj);
                                            if(image_clusterized(i, j) != 0)
                                            {
                                                 XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(i, j), selected));
                                                neighborPix2currentPix.insert(std::make_pair(std::make_pair(ki, kj), std::make_pair(i, j)));
                                                all_boundaries(i,j) = 1;
                                            }
                                            else
                                            {
                                                auto pair_it_reference = reference.equal_range(std::make_pair(i,j));
                                                for(auto it_reference = pair_it_reference.first; it_reference!=pair_it_reference.second; ++it_reference)
                                                {
                                                    if(all_boundaries(it_reference->second.first, it_reference->second.second) == 0 )
                                                    {
                                                        XY2planeIdx_boundary.insert(std::make_pair(it_reference->second, selected));
                                                        neighborPix2currentPix.insert(std::make_pair(std::make_pair(ki, kj), it_reference->second));
                                                        all_boundaries(it_reference->second.first, it_reference->second.second) = 1;
                                                    }
                                                }
                                            }
                                            XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(ki, kj), selected));
                                        }
                                    }
                                    isOpenning = false;
                                }
                            }
                        }
                    }

                    //---------------------------------------------------------------------------------------------------------------------------------------------------------
                    //check if it is not just a hole
                    if(isOpenning)
                    {
                        int rad1 = rad + 1;
                        min_ki = std::max(i-rad1, 0);
                        max_ki = std::min(i+rad1, Nrow-1);
                        min_kj = std::max(j-rad1, 0);
                        max_kj = std::min(j+rad1, Ncol-1);

                        int selected = image_clusterized_idx;

                        for(int ki = min_ki; ki<=max_ki; ++ki)
                        {
                            for(int kj = min_kj; kj<=max_kj; ++kj)
                            {
                                bool is_different = image_clusterized(ki, kj) != image_clusterized_idx;
                                bool is_not_black = image_clusterized(ki, kj)>0;
                                if(is_different && is_not_black) //if pixels belong to different clusters and neighbor not black
                                {
                                    if( all_boundaries(ki, kj) == 0 )
                                    {
                                        if(image_clusterized(ki, kj) == selected || selected==image_clusterized_idx)
                                        {
                                            selected = image_clusterized(ki, kj);
                                            if(image_clusterized(i, j)!=0)
                                            {
                                                XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(i, j), selected));
                                                neighborPix2currentPix.insert(std::make_pair(std::make_pair(ki, kj), std::make_pair(i, j)));
                                                all_boundaries(i,j) = 1;
                                            }
                                            else
                                            {
                                                auto pair_it_reference = reference.equal_range(std::make_pair(i,j));
                                                for(auto it_reference = pair_it_reference.first; it_reference!=pair_it_reference.second; ++it_reference)
                                                {
                                                    if(all_boundaries(it_reference->second.first, it_reference->second.second) == 0 )
                                                    {
                                                        XY2planeIdx_boundary.insert(std::make_pair(it_reference->second, selected));
                                                        neighborPix2currentPix.insert(std::make_pair(std::make_pair(ki, kj), it_reference->second));
                                                        all_boundaries(it_reference->second.first, it_reference->second.second) = 1;
                                                    }
                                                }
                                            }
                                            XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(ki, kj), selected));
                                        }
                                    }
                                    isOpenning = false;
                                }
                            }
                        }
                    }
                    //-----------------------------------------------------
                    //if the different neighbors remain black  = >  it is an opening =next to opening keep color of the current studied plane
                    if(isOpenning)
                    {
                        min_ki = std::max(i-rad, 0);
                        max_ki = std::min(i+rad, Nrow-1);
                        min_kj = std::max(j-rad, 0);
                        max_kj = std::min(j+rad, Ncol-1);

                        for(int ki = min_ki; ki<=max_ki; ++ki)
                        {
                            for(int kj = min_kj; kj<=max_kj; ++kj)
                            {
                                if( all_boundaries(ki, kj) == 0 )
                                {
                                    bool is_black_in_clusterized = image_clusterized(ki, kj)==0;
                                    bool is_black_in_morpho = image_morpho(ki, kj)==0;
                                    if(is_black_in_clusterized && is_black_in_morpho) //if pixels belong to different clusters and neighbor not black
                                    {
                                        if(image_clusterized(i, j) != 0)
                                        {
                                            XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(i, j), image_clusterized_idx));
                                            all_boundaries(i,j) = 1;
                                        }
                                        else
                                        {
                                            auto pair_it_reference = reference.equal_range(std::make_pair(i,j));
                                            for(auto it_reference = pair_it_reference.first; it_reference!=pair_it_reference.second; ++it_reference)
                                            {
                                                if(all_boundaries(it_reference->second.first, it_reference->second.second) == 0 )
                                                {
                                                    XY2planeIdx_boundary.insert(std::make_pair(it_reference->second, image_clusterized_idx));
                                                    all_boundaries(it_reference->second.first, it_reference->second.second) = 1;
                                                }
                                            }
                                        }
                                        isOpenning = true;
                                    }
                                }
                            }
                        }
                    }
                    //-----------------------------------------------------
                }
            }
        }
    }
}

Eigen::MatrixXi manager2D::getBoundariesImage()
{
    image_grad = Eigen::MatrixXi::Zero(Nrow, Ncol);
    if(XY2planeIdx_boundary.size()>0)
    {
        for (auto it_XY2planeIdx_boundary = XY2planeIdx_boundary.begin(); it_XY2planeIdx_boundary != XY2planeIdx_boundary.end(); ++it_XY2planeIdx_boundary )
            image_grad(it_XY2planeIdx_boundary->first.first, it_XY2planeIdx_boundary->first.second) = it_XY2planeIdx_boundary->second;
    }
    return image_grad;
}
