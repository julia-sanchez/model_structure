void corner::computePoint()
{
    Eigen::Matrix3d N ;
    N.row(0) = -planes[0]->normal;
    N.row(1) = -planes[1]->normal;
    N.row(2) = -planes[2]->normal;
    Eigen::Vector3d d = {planes[0]->distance, planes[1]->distance, planes[2]->distance};
    pt = N.colPivHouseholderQr().solve(d);
}

void corner::setPlanes(plane* p1, plane* p2, plane* p3)
{
    planes.resize(3);
    planes[0] = p1;
    planes[1] = p2;
    planes[2] = p3;
}
