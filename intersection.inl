void intersection::setDeltas(float dp, float dt)
{
    delta_phi = dp;
    delta_theta = dt;
    Nrow = (int)((2*M_PI+2*eps)/delta_phi);
    Ncol = (int)((M_PI+2*eps)/delta_theta);
}


void intersection::computeTheoriticalLineFeatures()
{
    if(!isOpening && !isObstruction)
    {
        std::cout<<"this intersection is an obstruction or a connection"<<std::endl;
        computeTheoriticalPhiTheta();
        Eigen::Vector3d normal1 =  -plane_neigh->normal;
        Eigen::Vector3d normal2 =  -plane_ref->normal;
        float distance1 = plane_neigh->distance;
        float distance2 = plane_ref->distance;
        tangente = normal1.cross(normal2);
        tangente /= tangente.norm();
        Eigen::Vector3d axis_x = {1,0,0};
        Eigen::Vector3d axis_y = {0,1,0};

        if(abs(normal1.dot(axis_x))<0.6 && abs(normal2.dot(axis_x))<0.6)
        {
            //as we will divide by normal1(2), we have to check it is not null to keep accuracy
            if(abs(normal1(2))<1e-3)
            {
                normal1 = -plane_ref->normal;
                normal2 = -plane_neigh->normal;
                distance1 = plane_ref->distance;
                distance2 = plane_neigh->distance;
            }
            pt(0) = 0; // we fix one element of the point
            float denom = normal2(1)-(normal2(2) * normal1(1)/normal1(2));
            float num = distance2-pt(0)*normal2(0) + (pt(0)*normal1(0)-distance1)*(normal2(2)/normal1(2));
            pt(1) = num/denom;
            pt(2) = (distance1-normal1(0)*pt(0)-pt(1)*normal1(1))/normal1(2);
        }
        else if (abs(normal1.dot(axis_y))<0.6 && abs(normal2.dot(axis_y))<0.6)
        {
            //as we will divide by normal1(2), we have to check it is not null to keep accuracy
            if(abs(normal1(2))<1e-3)
            {
                normal1 = -plane_ref->normal;
                normal2 = -plane_neigh->normal;
                distance1 = plane_ref->distance;
                distance2 = plane_neigh->distance;
            }
            pt(1) = 0; // we fix one element of the point
            float denom = normal2(0)-(normal2(2)*normal1(0)/normal1(2));
            float num = distance2-pt(1)*normal2(1) + (pt(1)*normal1(1)-distance1)*(normal2(2)/normal1(2));
            pt(0) = num/denom;
            pt(2) = (distance1-normal1(1)*pt(1)-pt(0)*normal1(0))/normal1(2);
        }
        else
        {
            //as we will divide by normal1(1), we have to check it is not null to keep accuracy
            if(abs(normal1(1))<1e-3)
            {
                normal1 = -plane_ref->normal;
                normal2 = -plane_neigh->normal;
                distance1 = plane_ref->distance;
                distance2 = plane_neigh->distance;
            }
            pt(2) = 0; // we fix one element of the point
            float denom = normal2(0)-(normal2(1) * normal1(0)/normal1(1));
            float num = distance2-pt(2)*normal2(2) + (pt(2)*normal1(2)-distance1)*(normal2(1)/normal1(1));
            pt(0) = num/denom;
            pt(1) = (distance1-normal1(2)*pt(2)-pt(0)*normal1(0))/normal1(1);
        }

        normal = pt-tangente.dot(pt)*tangente;
        normal /= normal.norm();
        if(pt.dot(normal)<0)
            normal *= -1;
        distance = pt.dot(normal);
    }
    else
    {
        std::cout<<"this intersection is an opening or an obstruction"<<std::endl;
        std::cout<<"It contains "<<other_points.size()<<" points"<<std::endl<<std::endl;
        ProjectBoundaryOnPlane();
//        Hough(2*M_PI/180, 0.03);
        RANSAC(0.03, 10000);
    }
}

void intersection::computeTheoriticalPhiTheta()
{
    // Fix Phi and compute theta
    for(int i = 0; i<Nrow; ++i)
    {
        Eigen::Vector3d normal1 = -plane_neigh->normal;
        Eigen::Vector3d normal2 = -plane_ref->normal;
        float phi = i*delta_phi + delta_phi/2;
        float A1 = normal1(0) * cos(phi);
        float B1 = normal1(1) * sin(phi);
        float C1 = normal1(2);
        float A2 = normal2(0) * cos(phi);
        float B2 = normal2(1) * sin(phi);
        float C2 = normal2(2);
        float D = plane_ref->distance/plane_neigh->distance;
        float theta = atan( (D*C1-C2)/(A2+B2 - D*(A1+B1)) );
        if(theta<0)
            theta += M_PI;
        theoritical_phi_theta.push_back(std::make_pair(phi,theta));
        int X = i;
        int Y = (int)(theta/delta_theta);
        if(theoritical_pixels.find(std::make_pair(X,Y)) == theoritical_pixels.end())
            theoritical_pixels.insert(std::make_pair(X,Y));
    }

    for(int j = 0; j<Ncol; ++j)
    {
        Eigen::Vector3d normal1 = -plane_neigh->normal;
        Eigen::Vector3d normal2 = -plane_ref->normal;
        float theta = j*delta_theta + delta_theta/2;
        float A1 = normal1(0) * sin(theta);
        float B1 = normal1(1) * sin(theta);
        float C1 = normal1(2) * cos(theta);
        float A2 = normal2(0) * sin(theta);
        float B2 = normal2(1) * sin(theta);
        float C2 = normal2(2) * cos(theta);
        float D = plane_ref->distance/plane_neigh->distance;
        float a = A2 - D*A1;
        float b = B2 - D*B1;
        float c = D*C1 - C2;
        float R = sqrt(a*a + b*b);
        if(abs(c/R)<1)
        {
            float alpha = atan2(b, a);
            float phi = acos(c/R)+alpha;
            Eigen::Vector3d pt_norm = {cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)};
            if(normal1.dot(pt_norm)<0 || normal2.dot(pt_norm)<0)
                phi = -acos(c/R)+alpha;

            if(phi<0)
                phi += 2*M_PI;
            theoritical_phi_theta.push_back(std::make_pair(phi,theta));
            int X = (int)(phi/delta_phi);
            int Y = j;
            if(theoritical_pixels.find(std::make_pair(X,Y)) == theoritical_pixels.end())
                theoritical_pixels.insert(std::make_pair(X,Y));
        }
    }
}

void intersection::computeLim(std::vector<corner>& corners, std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> possible_inter3D)
{
    lim.clear();

    //for connexions : the limits of the intersections are the corners---------------------------------------------------------------
    if(plane_ref->index != plane_neigh->index)
    {
        for(int i = 0; i<corners.size(); ++ i)
        {
            bool has_plane0 = corners[i].planes[0]->index == plane_ref->index || corners[i].planes[0]->index == plane_neigh->index;
            bool has_plane1 = corners[i].planes[1]->index == plane_ref->index || corners[i].planes[1]->index == plane_neigh->index;
            bool has_plane2 = corners[i].planes[2]->index == plane_ref->index || corners[i].planes[2]->index == plane_neigh->index;
            if( (has_plane0 && has_plane1) || (has_plane0 && has_plane2) || (has_plane1 && has_plane2))
                lim.insert(corners[i].pt.dot(tangente));
        }
    }

    //when planes are not connected to anything------------------------------------------------------------------------------------------------------------------------------

    std::vector<float> dot;

    for (auto it_points=points.begin(); it_points != points.end();)
    {
        if(abs(abs(it_points->dot(plane_ref->normal))-plane_ref->distance)<0.05)
        {
            dot.push_back(it_points->dot(tangente));
            ++it_points;
        }
        else
            it_points = points.erase(it_points);
    }

    for(auto it_lim = lim.begin(); it_lim != lim.end(); ++it_lim)
        dot.push_back(*it_lim);

    Eigen::Vector3d start_pt = points[std::distance(dot.begin(), std::min_element(dot.begin(), dot.end()) )];
    Eigen::Vector3d end_pt =   points[std::distance(dot.begin(), std::max_element(dot.begin(), dot.end()) )];
    start = dot[std::distance(dot.begin(), std::min_element(dot.begin(), dot.end()) )];
    end =   dot[std::distance(dot.begin(), std::max_element(dot.begin(), dot.end()) )];

    std::cout<<"visible start from data points = "<<start_pt.transpose()<<std::endl;
    std::cout<<"visible end from data points = "<<end_pt.transpose()<<std::endl<<std::endl;

    std::vector<float> diff_start;
    std::vector<float> diff_end;

    //for openings or obstructions (ie for edges which are not intersections of planes) we check if we can find intersections between 2d lines---------------------------------------------------------------
    if(plane_ref->index == plane_neigh->index)
    {
        for(int i = 0; i<possible_inter3D.size(); ++i)
        {
            diff_start.push_back( (possible_inter3D[i]-start_pt).norm() );
            diff_end.push_back( (possible_inter3D[i]-end_pt).norm() );
        }

        auto it_min_start = std::min_element(diff_start.begin(), diff_start.end());
        auto it_min_end = std::min_element(diff_end.begin(), diff_end.end());
        int idx_start = std::distance(diff_start.begin(), it_min_start);
        int idx_end = std::distance(diff_end.begin(), it_min_end);
        if(diff_start[idx_start] < 0.3)
        {
            theoriticalLineIntersection.push_back(idx_start);
            if(possible_inter3D[idx_start].dot(tangente)<start )
                start = possible_inter3D[idx_start].dot(tangente);
        }
        if(diff_end[idx_end] < 0.3)
        {
            theoriticalLineIntersection.push_back(idx_end);
            if(possible_inter3D[idx_end].dot(tangente)>end)
                end = possible_inter3D[idx_end].dot(tangente);
        }
    }
}

void intersection::definePlaneConnection()
{
//if sufficient quantity of current intersection pixels are really the theoritical connection, the intersection is a connection. if not it is an obstruction

    float connect = 0.0;
    for (auto it_pix = pixels.begin(); it_pix != pixels.end(); ++it_pix)
    {
        if(  theoritical_pixels.find(std::make_pair(it_pix->first, it_pix->second)) != theoritical_pixels.end() )
            ++connect;
    }

    std::cout<<"Number of pixels which contain the real intersection : "<<pixels.size()<<std::endl;

    connect = connect/(float)(pixels.size());

    std::cout<<"percentage of pixels which contain the theoritical intersection : "<<connect<<std::endl<<std::endl;

    isConnection = false;
    isObstruction = false;
    isOpening = false;
    if( connect > 0.25) // if more than half the quantity of pixels of the intersection correspond to a theoritical connection
    {
        isConnection = true;
        isLine = true;
        plane_neigh->connected.insert(plane_ref->index);
        plane_ref->connected.insert(plane_neigh->index);
    }
    else
    {
        isObstruction = true;
        float phi = pixels.begin()->first * delta_phi + delta_phi/2;
        float theta = pixels.begin()->second * delta_theta + delta_theta/2;
        Eigen::Vector3d pt_norm = {cos(phi)*sin(theta), sin(phi)* sin(theta), cos(theta)};
        float distance_pt1 = abs(plane_neigh->distance/(pt_norm.dot(plane_neigh->normal)));
        float distance_pt2 = abs(plane_ref->distance/(pt_norm.dot(plane_ref->normal)));
        if( distance_pt1 - distance_pt2 > 0 )
        {
            plane_ref->obstructed.insert(plane_neigh->index);
            plane_neigh->obstructing.insert(plane_ref->index);
        }
        else
        {
            plane_neigh->obstructing.insert(plane_ref->index);
            plane_ref->obstructed.insert(plane_neigh->index);
        }
    }
}





void intersection::ProjectBoundaryOnPlane()
{
    rot = Eigen::Affine3d::Identity();
    rot_inv = Eigen::Affine3d::Identity();

    float angle = acos(-plane_ref->normal(2));
    Eigen::Vector3d axis;
    if(angle>eps)
    {
        axis = (-plane_ref->normal).cross(Eigen::Vector3d(0,0,1));
        axis /= axis.norm();
        rot.rotate( Eigen::AngleAxisd(angle, axis) );
        rot_inv.rotate( Eigen::AngleAxisd(angle, -axis) );
    }

    points2D.resize(other_points.size());
    for(int i = 0; i<other_points.size(); ++i)
    {
        Eigen::Vector3d turned = rot.linear()*other_points[i];
        float dist_to_plane = turned(2)-plane_ref->distance;
        if(abs(dist_to_plane)<0.05) // points on the same plane with noise
        {
            points2D[i](0) = turned(0);
            points2D[i](1) = turned(1);
        }
        else if (dist_to_plane<0) // point on obstructing plane
        {
            turned *= ( 1 + ((plane_ref->distance/turned(2))-1) );
            points2D[i](0) = turned(0);
            points2D[i](1) = turned(1);
        }
        else
        {
            turned *= ( 1 - (1-(plane_ref->distance/turned(2))));
            points2D[i](0) = turned(0);
            points2D[i](1) = turned(1);
        }
    }

//    pt_mean = Eigen::Vector2d::Zero();
//    for(int i = 0; i<points2D.size(); ++i)
//        pt_mean += points2D[i];

//    pt_mean /= points2D.size();

//    for(int i = 0; i<points2D.size(); ++i)
//        points2D[i] -= pt_mean;

    trans_z = Eigen::Affine3d::Identity();
    trans_z_inv = Eigen::Affine3d::Identity();
    trans_z.translate(Eigen::Vector3d(0.0, 0.0, -plane_ref->distance));
    trans_z_inv.translate(Eigen::Vector3d(0.0, 0.0, plane_ref->distance));

    trans_xy = Eigen::Affine3d::Identity();
//    trans_xy_inv = Eigen::Affine3d::Identity();
//    trans_xy.translate(-Eigen::Vector3d(pt_mean(0), pt_mean(1), 0));
//    trans_xy_inv.translate(Eigen::Vector3d(pt_mean(0), pt_mean(1), 0));




    pcl::PointCloud<pcl::PointXYZ> points_in_line_cloud0;
    points_in_line_cloud0.width = other_points.size();
    points_in_line_cloud0.height = 1;
    points_in_line_cloud0.is_dense = false;
    points_in_line_cloud0.points.resize(other_points.size());

    for(int k = 0; k<other_points.size(); ++k)
    {
        points_in_line_cloud0.points[k].x = points2D[k](0);
        points_in_line_cloud0.points[k].y = points2D[k](1);
        points_in_line_cloud0.points[k].z = 0;
    }

    pcl::io::savePCDFileASCII ("initial_projected.pcd", points_in_line_cloud0);
    system("bash pcd2csv.sh initial_projected.pcd");

    pcl::PointCloud<pcl::PointXYZ> points_in_line_cloud1;
    points_in_line_cloud1.width = other_points.size();
    points_in_line_cloud1.height = 1;
    points_in_line_cloud1.is_dense = false;
    points_in_line_cloud1.points.resize(other_points.size());

    for(int k = 0; k<other_points.size(); ++k)
    {
        points_in_line_cloud1.points[k].x = (rot.linear()*other_points[k])(0);
        points_in_line_cloud1.points[k].y = (rot.linear()*other_points[k])(1);
        points_in_line_cloud1.points[k].z = (rot.linear()*other_points[k])(2);
    }

    pcl::io::savePCDFileASCII ("initial.pcd", points_in_line_cloud1);
    system("bash pcd2csv.sh initial.pcd");

    std::cout<<"plane index : "<<plane_ref->index<<std::endl<<std::endl;
}

void intersection::RANSAC(double error, int tests)
{
    std::cout<<"number of points at the begining: "<<other_points.size()<<std::endl<<std::endl;
    int idx1, idx2;
    int max_points = 0;
    std::vector<int> indices_on_line;

    std::srand (std::time(NULL));

    Eigen::Vector2d Vector_tested;

    for(int i = 0; i< tests; ++i)
    {
        idx1 = rand() % other_points.size();// [0; other_points.size()[
        idx2 = rand() % other_points.size();

        while(idx2 == idx1)
            idx2 = rand() % other_points.size();

        Eigen::Vector2d tangente2D_temp = (points2D[idx1]-points2D[idx2]);
        tangente2D_temp /= tangente2D_temp.norm();
        Eigen::Vector2d normal2D_temp = points2D[idx1]-points2D[idx1].dot(tangente2D_temp)*tangente2D_temp;
        normal2D_temp /= normal2D_temp.norm();
        float distance2D_temp = points2D[idx1].dot(normal2D_temp);

        std::vector<int> indices_on_line_temp;

        for(int j = 0; j< points2D.size(); ++j)
        {
            if(abs(points2D[j].dot(normal2D_temp)-distance2D_temp) < error)
                indices_on_line_temp.push_back(j);
        }

        if(indices_on_line_temp.size()>max_points)
        {
            Vector_tested = (points2D[idx1]-points2D[idx2]);
            normal2D = normal2D_temp;
            tangente2D = tangente2D_temp;
            distance2D = distance2D_temp;
            indices_on_line = indices_on_line_temp;
            max_points = indices_on_line_temp.size();
        }
    }

    std::cout<<"normal2D : "<<normal2D.transpose()<<std::endl<<"tangente2D : "<<tangente2D.transpose()<<std::endl<<"distance2D : "<<distance2D<<std::endl<<"number of points : "<<max_points<<std::endl<<std::endl;

    //Erase , points of line and its pixels neighbors from other_points + add pixels neighbors to points
//-----------------------------------------------------------------------------------------------------------

    std::set<int> to_erase;
    int rad = 1;
    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> pt2D;
    points.clear();
    pixels.clear();
    for(int k = 0; k<indices_on_line.size(); ++k)
    {
        int X = other_pixels[indices_on_line[k]].first;
        int Y = other_pixels[indices_on_line[k]].second;
        int imin = std::max(0, X-rad);
        int imax = std::min(Nrow-1, X+rad);
        int jmin = std::max(0, Y-rad);
        int jmax = std::min(Ncol-1, Y+rad);

        for(int i = imin; i<=imax; ++i)
        {
            for(int j = jmin; j<=jmax; ++j)
            {
                auto it_found_pixels = std::find(other_pixels.begin(), other_pixels.end(), std::make_pair(i, j));
                // Get index of element from iterator

                if(it_found_pixels != other_pixels.end())
                {
                    int index = std::distance(other_pixels.begin(), it_found_pixels);
                    to_erase.insert(index);
                    points.push_back(other_points[index]);
                    pixels.push_back(other_pixels[index]);
                    pt2D.push_back(points2D[index]);
                }
            }
        }
    }

    points2D = pt2D;

    auto it_to_erase_start = to_erase.end();
    --it_to_erase_start;
    auto it_to_erase_end = to_erase.begin();
    --it_to_erase_end;
    for(auto it_to_erase = it_to_erase_start; it_to_erase!=it_to_erase_end; --it_to_erase)
    {
        other_points.erase(other_points.begin() + *it_to_erase);
        other_pixels.erase(other_pixels.begin() + *it_to_erase);
    }

//        least square for accuracy
//---------------------------------------------------------------------------------------------
    float theta = acos(-normal2D(0));

   Eigen::MatrixXd JTJ(2,2);
    Eigen::MatrixXd JTr(2, 1);
    Eigen::MatrixXd J(2,1);
    double r;
    int nbr_iterations = 10;
    for(int k =0; k<nbr_iterations; ++k)
    {
        JTJ.setZero();
        JTr.setZero();
        for(int i = 0; i< points2D.size(); ++i)
        {
            J(0,0) = points2D[i].dot(Eigen::Vector2d(sin(theta), cos(theta)));
            J(1,0) = -1;
            r = points2D[i].dot(normal2D)-distance2D;
            JTJ += J * J.transpose();
            JTr += J * r;
        }

        Eigen::MatrixXd result(2, 1);
        result = -JTJ.llt().solve(JTr);
        theta += result(0);
        normal2D = {-cos(theta), sin(theta)};
        distance2D += result(1);
    }
    tangente2D = {sin(theta), cos(theta)};
    if(distance2D<0)
    {
        distance2D *= -1;
        normal2D *= -1;
    }

    //convert to 3D
//---------------------------------------------------------------------------------------------

    tangente(0) = tangente2D(0);
    tangente(1) = tangente2D(1);
    tangente(2) = 0;
    tangente = rot_inv.linear()*tangente;
    normal(0) = normal2D(0);
    normal(1) = normal2D(1);
    normal(2) = 0;
    normal = rot_inv.linear()*normal;
    pt = plane_ref->distance*(-plane_ref->normal) + normal * distance2D;
    normal = pt - pt.dot(tangente)*tangente;
    normal /= normal.norm();
    distance = pt.dot(normal);

    std::cout<<"number of remaining points at the end: "<<other_points.size()<<std::endl<<std::endl;
}





void intersection::Hough(double error_angle, double error_distance)
{
    std::cout<<"Hough : Number of initial points in which we seek a line: "<<other_points.size()<<std::endl<<std::endl;
    double max_theta = 2*M_PI;
    double min_theta = 0;
    double delta_angle = error_angle*2;
    int Nhist_theta = std::ceil( (max_theta-min_theta + eps)/delta_angle );
    double delta_distance = error_distance * 2;
    std::vector<double> thetas(Nhist_theta);
    std::vector<double> distances(points2D.size() * Nhist_theta);

    normal2D.Zero();

    for (int n = 0; n<Nhist_theta; ++n)
        thetas[n] = min_theta + n*delta_angle + delta_angle/2;

    //define all possible distances for one point and one theta
    for (int k = 0; k<points2D.size(); ++k)
    {
        for (int n = 0; n<Nhist_theta; ++n)
        {
            normal2D = {-sin(thetas[n]), cos(thetas[n])};
            distances[k*Nhist_theta +n] = points2D[k].dot(normal2D);
        }
    }

    //fill histogram with theta/distance
    double max_distance = *(std::max_element(distances.begin(), distances.end()));
    double min_distance = *(std::min_element(distances.begin(), distances.end()));
    int Nhist_distance = std::ceil((max_distance - min_distance + eps)/delta_distance);

    Eigen::MatrixXi hist = Eigen::MatrixXi::Zero(Nhist_theta, Nhist_distance);
    std::multimap<std::pair<int, int>, int> angle_distance2boundary_points;
    for (int k = 0; k<points2D.size(); ++k)
    {
        for (int n = 0; n<Nhist_theta; ++n)
        {
            if(distances[k*Nhist_theta + n]>0)
            {
                ++hist((int)((thetas[n]-min_theta)/delta_angle), (int)((distances[k*Nhist_theta + n]-min_distance)/delta_distance));
                angle_distance2boundary_points.insert(std::make_pair(std::make_pair((int)((thetas[n]-min_theta)/delta_angle), (int)((distances[k*Nhist_theta + n]-min_distance)/delta_distance)), k));
            }
        }
    }

//    save_image_pgm("Boundary_Hough", "", hist, 100);

    //search interesting bins in histogram

    int i;
    int j;
    std::vector<int> indices_on_line;

    //put points corresponding to line in points_in_line
    auto points_found_idx = angle_distance2boundary_points.equal_range(std::make_pair(i,j));

    for(auto it = points_found_idx.first; it !=points_found_idx.second; ++it)
        indices_on_line.push_back(it->second);

    //---------------------------------------------------------------------------------------------------------------------------------

    // compute line features (theta distance) with bin features (coord i and coord j)
    double theta = i * delta_angle + min_theta + delta_angle/2;
    distance2D = j * delta_distance + min_distance + delta_distance/2;
    tangente2D = {cos(theta), sin(theta)};
    normal2D = {-sin(theta), cos(theta)};
//    if(normal2D.dot(points2D[indices_on_line[0]])<0)
//        normal2D *= -1;

    points.resize(indices_on_line.size());
    pixels.resize(indices_on_line.size());
    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> pt2D(indices_on_line.size());
    for(int i = 0; i < indices_on_line.size(); ++i)
    {
        points[i] = other_points[indices_on_line[i]];
        pixels[i] = other_pixels[indices_on_line[i]];
        pt2D[i] = points2D[indices_on_line[i]];
    }

    std::cout<<"Number of points found in line : "<<indices_on_line.size()<<std::endl<<std::endl;

    //Erase , points of line and its pixels neighbors from other_points + add pixels neighbors to points
    std::set<int> to_erase;
    int rad = 1;
    for(int k = 0; k<indices_on_line.size(); ++k)
    {
        int X = other_pixels[indices_on_line[k]].first;
        int Y = other_pixels[indices_on_line[k]].second;
        int imin = std::max(0, X-rad);
        int imax = std::min(Nrow-1, X+rad);
        int jmin = std::max(0, Y-rad);
        int jmax = std::min(Ncol-1, Y+rad);

        for(int i = imin; i<=imax; ++i)
        {
            for(int j = jmin; j<=jmax; ++j)
            {
                auto it_found_pixels = std::find(other_pixels.begin(), other_pixels.end(), std::make_pair(i, j));
                // Get index of element from iterator
                if(it_found_pixels != other_pixels.end())
                {
                    int index = std::distance(other_pixels.begin(), it_found_pixels);
                    to_erase.insert(index);
                    points.push_back(other_points[index]);
                    pt2D.push_back(points2D[index]);
                }
            }
        }
    }

    points2D = pt2D;

    auto it_to_erase_start = to_erase.end();
    --it_to_erase_start;
    auto it_to_erase_end = to_erase.begin();
    --it_to_erase_end;
    for(auto it_to_erase = it_to_erase_start; it_to_erase!=it_to_erase_end; --it_to_erase)
    {
        other_points.erase(other_points.begin() + *it_to_erase);
        other_pixels.erase(other_pixels.begin() + *it_to_erase);
    }

    //---------------------------------------------------------------------------------------------
   Eigen::MatrixXd JTJ(2,2);
    Eigen::MatrixXd JTr(2, 1);
    Eigen::MatrixXd J(2,1);
    double r;
    int nbr_iterations = 10;
    for(int k =0; k<nbr_iterations; ++k)
    {
        JTJ.setZero();
        JTr.setZero();
        for(int i = 0; i< points2D.size(); ++i)
        {
            J(0,0) = points2D[i].dot(Eigen::Vector2d(-cos(theta), -sin(theta)));
            J(1,0) = -1;
            r = points2D[i].dot(normal2D)-distance2D;
            JTJ += J * J.transpose();
            JTr += J * r;
        }

        Eigen::MatrixXd result(2, 1);
        result = -JTJ.llt().solve(JTr);
        theta += result(0);
        normal2D = {-sin(theta), cos(theta)};
        distance2D += result(1);
    }
    tangente2D = {cos(theta), sin(theta)};
    if(distance2D<0)
    {
        distance2D *= -1;
        normal2D *= -1;
    }

    //---------------------------------------------------------------------------------------------

    //convert to 3D

    //retransform in xy:
    distance2D = (distance2D*normal2D+pt_mean).dot(normal2D);
    tangente(0) = tangente2D(0);
    tangente(1) = tangente2D(1);
    tangente(2) = 0;
    tangente = rot_inv.linear()*tangente;
    normal(0) = normal2D(0);
    normal(1) = normal2D(1);
    normal(2) = 0;
    normal = rot_inv.linear()*normal;
    pt = plane_ref->distance*(-plane_ref->normal) + normal * distance2D;
    normal = pt - pt.dot(tangente)*tangente;
    normal /= normal.norm();
    distance = pt.dot(normal);

    std::cout<<"number of remaining points at the end: "<<other_points.size()<<std::endl<<std::endl;
}
