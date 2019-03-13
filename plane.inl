void plane::computeNormal()
{
    Eigen::MatrixXd centered_points (points.rows(), 3);
    Eigen::Vector3d mean_point = points.colwise().mean();
    centered_points = points - Eigen::VectorXd::Ones(points.rows())*(mean_point.transpose());
    Eigen::Matrix3d covariance = centered_points.transpose()*centered_points/points.rows();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(covariance);
    normal = es.eigenvectors().col(0);
    if(normal.dot(mean_point)>0)
        normal *= -1;
    distance = abs(normal.dot(mean_point));
}


void plane::appendPoint(Eigen::Vector3d pt)
{
    points.conservativeResize(points.rows()+1, points.cols());
    points.row(points.rows()-1) = pt;
}

void plane::appendPixel(std::pair<int,int> p)
{
    pixels.insert(p);
}
