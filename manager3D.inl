manager3D::manager3D()
{
    rot_axis = {0, 0, 1}; // rotation secondaire du scanner
    axis_init_phi = {1, 0, 0};
}

void manager3D::getNormals(double radius)
{
    std::cout<< "start computing normals"<<std::endl;

    if(cloud->points[10].normal_x > epsilon && cloud->points[10].normal_y == 0 && cloud->points[10].normal_z == 0) // if no normal is found in the file -> compute with PCA
    {
        pcl::NormalEstimationOMP<pcl::PointNormal, pcl::PointNormal> normal_estimation;
        normal_estimation.setSearchMethod(boost::make_shared<pcl::search::KdTree<pcl::PointNormal>>(tree));
        normal_estimation.setRadiusSearch(radius);
        normal_estimation.setViewPoint (0, 0, 0);
        normal_estimation.setInputCloud(cloud);
        normal_estimation.compute(*cloud);
    }
    else
    {
        orient();
    }

    std::vector<int> indices;
    pcl::removeNaNNormalsFromPointCloud(*cloud, *cloud, indices);

    std::cout<< "stop computing normals"<<std::endl;
    std::cout<< "source: points number after preprocessing : "<<cloud->size()<<std::endl;
}

void manager3D::orient()
{
    for (int i = 0; i<cloud->points.size(); ++i)
    {
        if( cloud->points[i].normal_x* cloud->points[i].x +  cloud->points[i].normal_y* cloud->points[i].y +  cloud->points[i].normal_z* cloud->points[i].z > 0 )
        {
            cloud->points[i].normal_x *= -1;
            cloud->points[i].normal_y *= -1;
            cloud->points[i].normal_z *= -1;
        }
    }
}

bool manager3D::isSeed(Eigen::Vector3d& pt, Eigen::Vector3d& normal,  double radius, double threshold)
{
    std::vector<int> pointIdxNKNSearch;
    std::vector<float> pointNKNSquaredDistance;
    Eigen::Vector3d mean = {0, 0, 0};
    pcl::PointNormal searchPoint;
    searchPoint.x = pt(0);
    searchPoint.y = pt(1);
    searchPoint.z = pt(2);

    if ( tree.radiusSearch (searchPoint, radius, pointIdxNKNSearch, pointNKNSquaredDistance) > 0 )
    {
        if(pointIdxNKNSearch.size ()<10)
            return false;
        for (size_t i = 0; i < pointIdxNKNSearch.size (); ++i)
        {
            Eigen::Vector3d neigh = {cloud->points[pointIdxNKNSearch[i]].x, cloud->points[pointIdxNKNSearch[i]].y, cloud->points[pointIdxNKNSearch[i]].z};
            mean += neigh;
        }
        mean /= pointIdxNKNSearch.size ();
        for (size_t i = 0; i < pointIdxNKNSearch.size (); ++i)
        {
            Eigen::Vector3d neigh = {cloud->points[pointIdxNKNSearch[i]].x, cloud->points[pointIdxNKNSearch[i]].y, cloud->points[pointIdxNKNSearch[i]].z};
            if(abs((neigh-pt).dot(normal))>threshold)
                return false;
        }
    }
    else
        return false;

    pt = mean;
    return true;
}

void manager3D::createMap()
{
    pcl::PointCloud<pcl::PointNormal> pc_init;
    /// depth_map : each color is a cluster on the gaussian image (main normal)
    int n = 0;

    for (int i =0; i<cloud->size(); ++i)
    {
        //compute point location phi/theta
        Eigen::Vector3d pt = {cloud->points[i].x, cloud->points[i].y, cloud->points[i].z};
        Eigen::Vector3d normal = {cloud->points[i].normal_x, cloud->points[i].normal_y, cloud->points[i].normal_z};
        double theta;
        double phi;
        if(pt.norm()!=0)
        {
            theta = acos(pt.dot(rot_axis)/pt.norm());
            phi = atan2((rot_axis.cross(axis_init_phi)).dot(pt), axis_init_phi.dot(pt));
            if (phi<0)
                phi += 2*M_PI;
        }

        if(theta>theta_margin && theta <M_PI-theta_margin) // take theta from 5 degrees because accuracy is insufficient before
        {
            pc_init.points.push_back(cloud->points[i]);
            auto idx_is_inserted = planes.insert(std::make_pair(std::make_pair(phi, theta), std::make_pair(pt,normal)));
            if(idx_is_inserted.second == false)
                ++n;
        }
    }

    pc_init.width = pc_init.points.size();
    pc_init.height = 1;
    pc_init.is_dense=false;
    pcl::io::savePCDFileASCII ("pc_init.pcd", pc_init);
    system("bash pcd2csv.sh pc_init.pcd");

    std::cout<<"points not added to PhiTheta2PtN because already exists : "<<n<<std::endl<<std::endl;
}
