void corner::computePointFromPlanes()
{
    Eigen::Matrix3d N ;
    N.row(0) = -planes[0]->normal;
    N.row(1) = -planes[1]->normal;
    N.row(2) = -planes[2]->normal;
    Eigen::Vector3d d = {planes[0]->distance, planes[1]->distance, planes[2]->distance};
    pt = N.colPivHouseholderQr().solve(d);
}

void corner::computePointFromLines()
{
    plane* p;
    if(!lines[0]->isConnection)
        p = lines[0]->plane_ref;
    else
        p = lines[1]->plane_ref;

    //define features of lines in 2D
    Eigen::Affine3d rot = Eigen::Affine3d::Identity();
    Eigen::Affine3d rot_inv = Eigen::Affine3d::Identity();

    float angle = acos(-p->normal(2));   //angle between normal and z axis [0, pi]
    Eigen::Vector3d axis;
    if(angle>eps)
    {
        axis = (-p->normal).cross(Eigen::Vector3d(0,0,1)); // rotation axis to align normal onto z axis
        axis /= axis.norm();
        rot.rotate( Eigen::AngleAxisd(angle, axis) );
        rot_inv.rotate( Eigen::AngleAxisd(-angle, axis) );
    }

    lines[0]->rot = rot;
    lines[1]->rot = rot;
    lines[0]->ProjectLineFeaturesOnPlane();
    lines[1]->ProjectLineFeaturesOnPlane();

    Eigen::Matrix2d N ;
    N.row(0) = lines[0]->normal2D;
    N.row(1) = lines[1]->normal2D;
    Eigen::Vector2d d = {lines[0]->distance2D, lines[1]->distance2D};
    Eigen::Vector2d pt_intersect2D = N.colPivHouseholderQr().solve(d);

    Eigen::Vector3d pt_intersect3D = {pt_intersect2D(0), pt_intersect2D(1), p->distance};
    pt = rot_inv.linear() *  pt_intersect3D;
}

void corner::setPlanes(plane* p1, plane* p2, plane* p3)
{
    planes_indices.clear();
    planes.resize(3);
    planes[0] = p1;
    planes_indices.insert(p1->index);
    planes[1] = p2;
    planes_indices.insert(p2->index);
    planes[2] = p3;
    planes_indices.insert(p3->index);
}

void corner::setLines(intersection* l1, intersection* l2)
{
    planes_indices.clear();
    lines.resize(2);
    lines[0] = l1;
    lines[1] = l2;

    planes.resize(1);
    planes[0] = l1->plane_ref;

    //add correspondant planes to corner
    auto ret = planes_indices.insert(l1->plane_neigh->index);
    if (ret.second)
        planes.push_back(l1->plane_neigh);
    //-------------------------------------
    ret = planes_indices.insert(l2->plane_ref->index);
    if (ret.second)
        planes.push_back(l2->plane_ref);
    //-------------------------------------
    ret = planes_indices.insert(l2->plane_neigh->index);
    if (ret.second)
        planes.push_back(l2->plane_neigh);
    //-------------------------------------
}
