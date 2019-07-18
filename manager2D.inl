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

Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> manager2D::morpho(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_in, Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> SE)
{
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_out = image_in;
    int rad = (SE.rows()-1)/2;

    for (int j = lim_theta.first+rad; j< lim_theta.second-rad; ++j)
    {
        for(int i = 0; i < rad; ++i)
        {
            if(image_in(i,j))
            {
                image_out.block(0,j-rad,rad+i,5) = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>::Ones(rad+i,5);
                image_out.block(Nrow-rad+i,j-rad,rad-i,5) = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>::Ones(rad-i,5);
            }
        }
    }

    for (int j = lim_theta.first+rad; j<= lim_theta.second-rad; ++j)
    {
        for(int i = Nrow-rad; i < Nrow; ++i)
        {
            if(image_in(i,j))
            {
                image_out.block(i,j-rad, Nrow-i,5) = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>::Ones(Nrow-i,5);
                image_out.block(0,j-rad,5-Nrow-i,5) = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>::Ones(5-Nrow-i,5);
            }
        }
    }

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
//                image_out.block(i-rad,j-rad,5,5) = image_out.block(i-rad,j-rad,5,5) || SE;
        }
    }


    return image_out;
}

void manager2D::closure(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> SE)
{
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_dilated = morpho(image_bool, SE); //dilate
//    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_closed = morpho(mynot(image_dilated), SE); //dilate(not) = erode
//    image_morpho = mynot(image_closed); // not
//    image_morpho = image_morpho - (image_other_bool && image_morpho); // remove other walls
    image_morpho = image_dilated - (image_other_bool && image_dilated); // remove other walls
}

void manager2D::computeBoundaries(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> SE, Eigen::MatrixXi all_boundaries)
{
    int rad = (SE.rows()-1)/2;
    neighborPix2currentPix.clear();

    XY2planeIdx_boundary.clear();
    for (int i = 0; i<Nrow; ++i)
    {
        for (int j = lim_theta.first+rad; j<=lim_theta.second-rad; ++j)
        {
            if(image_morpho(i,j) && all_boundaries(i,j) == 0 )
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
                    sum = (image_morpho.block(i-rad, j-rad, SE.rows(), SE.cols()) && SE).cast<int>().sum();

                if(sum!=SE.rows()*SE.cols()) // at least one neighbor is different from (i,j)
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
                            if(all_boundaries(ki, kj) == 0 )
                            {
                                bool is_different = image_clusterized(ki, kj) != image_clusterized_idx;
                                bool is_not_black = image_clusterized(ki, kj)>0;
                                if(is_different && is_not_black) //if pixels belong to different clusters and neighbor not black
                                {
                                    if(image_clusterized(ki, kj) == selected || selected==image_clusterized_idx)
                                    {
                                        selected = image_clusterized(ki, kj);
                                        if(image_clusterized(i, j)!=0)
                                        {
                                            XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(i, j), selected));
                                            neighborPix2currentPix.insert(std::make_pair(std::make_pair(ki, kj), std::make_pair(i, j)));
                                        }
                                        else
                                        {
                                            auto pair_it_reference = reference.equal_range(std::make_pair(i,j));
                                            for(auto it_reference = pair_it_reference.first; it_reference!=pair_it_reference.second; ++it_reference)
                                            {
                                                XY2planeIdx_boundary.insert(std::make_pair(it_reference->second, selected));
                                                neighborPix2currentPix.insert(std::make_pair(std::make_pair(ki, kj), it_reference->second));
                                            }

                                        }
                                        XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(ki, kj), selected));
                                        isOpenning = false;
                                    }
                                }
                            }
                            else
                                isOpenning = false;
                        }
                    }

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
                                if(all_boundaries(ki, kj) == 0 )
                                {
                                    bool is_different = image_clusterized(ki, kj) != image_clusterized_idx;
                                    bool is_not_black = image_clusterized(ki, kj)>0;
                                    if(is_different && is_not_black) //if pixels belong to different clusters and neighbor not black
                                    {
                                        if(image_clusterized(ki, kj) == selected || selected==image_clusterized_idx)
                                        {
                                            selected = image_clusterized(ki, kj);
                                            if(image_clusterized(i, j)!=0)
                                            {
                                                XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(i, j), selected));
                                                neighborPix2currentPix.insert(std::make_pair(std::make_pair(ki, kj), std::make_pair(i, j)));
                                            }
                                            else
                                            {
                                                auto pair_it_reference = reference.equal_range(std::make_pair(i,j));
                                                for(auto it_reference = pair_it_reference.first; it_reference!=pair_it_reference.second; ++it_reference)
                                                {
                                                    XY2planeIdx_boundary.insert(std::make_pair(it_reference->second, selected));
                                                    neighborPix2currentPix.insert(std::make_pair(std::make_pair(ki, kj), it_reference->second));
                                                }
                                            }
                                            XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(ki, kj), selected));
                                            isOpenning = false;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //if the different neighbors remain black  = >  it is an opening =next to opening keep color of the current studied plane
                    if(isOpenning)
                    {
                        if(image_clusterized(i, j)!=0)
                            XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(i, j), image_clusterized_idx));
                        else
                        {
                            auto pair_it_reference = reference.equal_range(std::make_pair(i,j));
                            for(auto it_reference = pair_it_reference.first; it_reference!=pair_it_reference.second; ++it_reference)
                                XY2planeIdx_boundary.insert(std::make_pair(it_reference->second, image_clusterized_idx));
                        }
                    }
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
        image_grad = max_col * image_grad/image_grad.maxCoeff();
    }
    return image_grad;
}
