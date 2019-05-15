void corner::computePointFromPlanes()
{
    Eigen::Matrix3d N ;
    N.row(0) = -planes[0]->normal;
    N.row(1) = -planes[1]->normal;
    N.row(2) = -planes[2]->normal;
    Eigen::Vector3d d = {planes[0]->distance, planes[1]->distance, planes[2]->distance};
    pt = N.colPivHouseholderQr().solve(d);
}

void corner::computePointFromLines(int plane_idx)
{
    if(lines[0]->isConnection)
        lines[0]->ProjectLineFeaturesOnPlane(plane_idx);
    if(lines[1]->isConnection)
        lines[1]->ProjectLineFeaturesOnPlane(plane_idx);

    Eigen::Vector2d Dn = lines[0]->distance2D * lines[0]->normal2D - lines[1]->distance2D * lines[1]->normal2D;
    float X;

    float ax = pow(lines[1]->tangente2D(0),2) - pow(lines[0]->tangente2D(0),2);
    float bx = -2*lines[0]->tangente2D(0)*Dn(0);
    float cx = (pow(lines[0]->distance2D,2)-pow(lines[1]->distance2D,2))*pow(lines[1]->tangente2D(0),2)-pow(Dn(0),2);

    float ay = pow(lines[1]->tangente2D(1),2) - pow(lines[0]->tangente2D(1),2);
    float by = -2*lines[0]->tangente2D(1)*Dn(1);
    float cy = (pow(lines[0]->distance2D,2)-pow(lines[1]->distance2D,2))*pow(lines[1]->tangente2D(1),2)-pow(Dn(1),2);

    X = (cx*(ay/ax)-cy)/(by-bx*(ay/ax));

    std::cout<<"tangente 1  "<<lines[0]->tangente2D.transpose()<<std::endl;
    std::cout<<"tangente 2 "<<lines[1]->tangente2D.transpose()<<std::endl;

    Eigen::Vector2d pt_intersect2D = lines[0]->distance2D * lines[0]->normal2D + X * lines[0]->tangente2D;
    Eigen::Vector3d pt_intersect3D = {pt_intersect2D(0), pt_intersect2D(1), 0};
    pt = lines[0]->rot_inv.linear() * lines[0]->trans_z_inv*pt_intersect3D;
}

void corner::setPlanes(plane* p1, plane* p2, plane* p3)
{
    planes.resize(3);
    planes[0] = p1;
    planes[1] = p2;
    planes[2] = p3;
}

void corner::setLines(intersection* l1, intersection* l2)
{
    lines.resize(2);
    lines[0] = l1;
    lines[1] = l2;

    planes.resize(1);
    planes[0] = l1->plane_ref;

    //add correspondant planes to corner

    bool already_exists;
    //-------------------------------------
    already_exists = false;
    for(int i = 0; i<planes.size(); ++i)
    {
        if(planes[i]->index == l1->plane_neigh->index)
            already_exists = true;
    }
    if(!already_exists)
        planes.push_back(l1->plane_neigh);
    //-------------------------------------
    already_exists = false;
    for(int i = 0; i<planes.size(); ++i)
    {
        if(planes[i]->index == l2->plane_ref->index)
            already_exists = true;
    }
    if(!already_exists)
        planes.push_back(l2->plane_ref);
    //-------------------------------------
    already_exists = false;
    for(int i = 0; i<planes.size(); ++i)
    {
        if(planes[i]->index == l2->plane_neigh->index)
            already_exists = true;
    }
    if(!already_exists)
        planes.push_back(l2->plane_neigh);
    //-------------------------------------
}
