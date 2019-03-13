manager2D::manager2D()
{
}

void manager2D::setImageClusterized(Eigen::MatrixXi p)
{
    image_clusterized = p;
    Nrow = image_clusterized.rows();
    Ncol = image_clusterized.cols();
    image_grad = Eigen::MatrixXi::Zero(Nrow, Ncol);
}

void manager2D::binarize(int p) // p is the index of the plane
{
    image_bool = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();;
    image_other_bool = Eigen::MatrixXi::Zero(Nrow, Ncol).cast<bool>();;

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
    //Structure element is convolved with the image with an OR process

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_out = image_in;
    int rad = (SE.rows()-1)/2;
    for (int i = rad; i< Nrow-rad; ++i)
    {
        for (int j = rad; j< Ncol-rad; ++j)
        {
            if(image_in(i,j))
                image_out.block<5,5>(i-rad,j-rad) = image_out.block<5,5>(i-rad,j-rad) || SE;
        }
    }
    return image_out;
}

void manager2D::closure(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> SE)
{
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_dilated = morpho(image_bool, SE); //dilate
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> image_closed = morpho(mynot(image_dilated), SE); //not(erode)
    image_morpho = mynot(image_closed); // not
    image_morpho = image_morpho - (image_other_bool && image_morpho); // remove other walls
}

void manager2D::computeBoundaries(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> SE, int margin)
{
    int rad = (SE.rows()-1)/2;

    XY2planeIdx_boundary.clear();
    for (int i = rad; i<Nrow-rad; ++i)
    {
        for (int j = margin+2; j<Ncol-(margin+2); ++j)
        {
            if(image_bool(i,j))
            {
                int sum = (image_morpho.block<3,3>(i-rad,j-rad) && SE).cast<int>().sum();
                if(sum!=SE.rows()*SE.cols()) // at least one neighbor is different from (i,j)
                {
                    bool isOpenning = true;
                    //search who is the neighbor which leads to this boundary
                    for(int ki = 0; ki<SE.rows(); ++ki)
                    {
                        for(int kj = 0; kj<SE.cols(); ++kj)
                        {
                            bool is_different = image_clusterized(i-rad+ki, j-rad +kj) != image_clusterized(i,j);
                            bool is_not_black = image_clusterized(i-rad+ki, j-rad +kj)>0;
                            if(is_different && is_not_black) //if pixels belong to different clusters and neighbor not black
                            {
                                XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(i, j), image_clusterized(i-rad+ki, j-rad +kj)));
                                XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(i-rad+ki, j-rad +kj), image_clusterized(i-rad+ki, j-rad +kj)));
                                isOpenning = false;
                            }
                        }
                    }

                    //if the different neighbors are black  = >  next to opening keep color of the current studied plane
                    if(isOpenning)
                       XY2planeIdx_boundary.insert(std::make_pair(std::make_pair(i,j), image_clusterized(i, j)));
                }
            }
        }
    }
}

Eigen::MatrixXi manager2D::getBoundariesImage()
{
    image_grad = Eigen::MatrixXi::Zero(Nrow, Ncol);

    for (auto it_XY2planeIdx_boundary = XY2planeIdx_boundary.begin(); it_XY2planeIdx_boundary != XY2planeIdx_boundary.end(); ++it_XY2planeIdx_boundary )
        image_grad(it_XY2planeIdx_boundary->first.first, it_XY2planeIdx_boundary->first.second) = it_XY2planeIdx_boundary->second;
    image_grad = max_col * image_grad/image_grad.maxCoeff();
    return image_grad;
}
