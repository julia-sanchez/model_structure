void plane::computeNormal()
{
    Eigen::MatrixXd points = Eigen::MatrixXd::Zero(pts.size(), 3);
    for(int k = 0; k<pts.size(); ++k)
        points.row(k) = pts[k];
    Eigen::MatrixXd centered_points (points.rows(), 3);
    mean_point_ = points.colwise().mean();
    centered_points = points - Eigen::VectorXd::Ones(points.rows())*(mean_point_.transpose());
    Eigen::Matrix3d covariance = centered_points.transpose()*centered_points/points.rows();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(covariance);
    normal = es.eigenvectors().col(0);
    if(normal.dot(mean_point_)>0)
        normal *= -1;
    distance = abs(normal.dot(mean_point_));
}

void plane::computeMeanPoint()
{
    mean_point_ = Eigen::Vector3d::Zero();
    for(int i = 0; i<pts.size(); ++i)
        mean_point_ += pts[i];
    mean_point_ /= pts.size();
}


void plane::appendPoint(Eigen::Vector3d pt)
{
    pts.push_back(pt);
}

void plane::appendPixel(std::pair<int,int> p)
{
    pixels.push_back(p);
}

void plane::mean()
{
    //filter points to compute the mean independently from the sampling
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_to_sample (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_sampled (new pcl::PointCloud<pcl::PointXYZ>);
    cloud_to_sample->width = pts.size();
    cloud_to_sample->height = 1;
    cloud_to_sample->points.resize(pts.size());
    for(int i = 0; i<pts.size(); ++i)
    {
        cloud_to_sample->points[i].x = pts[i](0);
        cloud_to_sample->points[i].y = pts[i](1);
        cloud_to_sample->points[i].z = pts[i](2);
    }
    pcl::UniformSampling< pcl::PointXYZ > filter;
    filter.setInputCloud(cloud_to_sample);
    filter.setRadiusSearch(0.1);
    filter.filter(*cloud_sampled);

    std::vector<Eigen::Vector3d> sampled(cloud_sampled->size());

    for(int i = 0; i<cloud_sampled->size(); ++i)
    {
        sampled[i](0) = cloud_sampled->points[i].x;
        sampled[i](1) = cloud_sampled->points[i].y;
        sampled[i](2) = cloud_sampled->points[i].z;
    }

    Eigen::Vector3d mean_point = Eigen::Vector3d::Zero();

    for(int i = 0; i<sampled.size(); ++i)
        mean_point += sampled[i];

    mean_point /= sampled.size();

    //take the point the closest from the mean as the new seed

    int index = 0;
    float temp = 100000000;

    for(int i = 0; i<pts.size(); ++i)
    {
        if((pts[i]-mean_point).norm()<temp)
        {
            temp = (pts[i]-mean_point).norm();
            index = i;
        }
    }

    seed.first = pixels[index];
    seed.second = pts[index];
}
