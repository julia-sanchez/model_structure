void intersection::setDeltas(float dp, float dt)
{
    delta_phi = dp;
    delta_theta = dt;
    Nrow = (int)((2*M_PI+2*eps)/delta_phi);
    Ncol = (int)((M_PI+2*eps)/delta_theta);
}

std::pair<int,int> pt2XY(Eigen::Vector3d pt, float delta_phi, float delta_theta, Eigen::Vector3d rot_axis, Eigen::Vector3d axis_init_phi)
{
    float theta = acos(pt.dot(rot_axis)/pt.norm());
    float phi = atan2((rot_axis.cross(axis_init_phi)).dot(pt), axis_init_phi.dot(pt));
    if (phi<0)
        phi += 2*M_PI;

    return std::make_pair((int)(phi/delta_phi), (int)(theta/delta_theta));
}


void intersection::computeTheoriticalLineFeatures()
{
    if(!isOpening && !isObstruction)
    {
        computeTheoriticalPhiTheta();
        Eigen::Vector3d normal1 =  -plane_neigh->normal;
        Eigen::Vector3d normal2 =  -plane_ref->normal;
        float distance1 = plane_neigh->distance;
        float distance2 = plane_ref->distance;
        tangente = normal1.cross(normal2);
        tangente /= tangente.norm();
        Eigen::Vector3d axis_x = {1,0,0};
        Eigen::Vector3d axis_y = {0,1,0};

        double thresh_to_define_axis = 0.6;

        if(abs(normal1.dot(axis_x))<thresh_to_define_axis && abs(normal2.dot(axis_x))<thresh_to_define_axis)
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
        else if (abs(normal1.dot(axis_y))<thresh_to_define_axis && abs(normal2.dot(axis_y))<thresh_to_define_axis)
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
        ProjectBoundaryOnPlane(other_points, 100000);
//        Hough(2*M_PI/180, 0.03);
        RANSAC(max_line_distance, ransac_iterations);
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

void intersection::selectCorners(int rad, float delta_phi, float delta_theta, Eigen::Vector3d rot_axis, Eigen::Vector3d axis_init_phi)
{
    //compare with theoretical corners

    std::map<float, int> diff_start;
    std::map<float, int> diff_end;

    if(theoreticalLim.size() != 0)
    {
        for(int i = 0; i<theoreticalLim.size(); ++i)
        {
            diff_start.insert( std::make_pair((theoreticalLim[i]-start_pt).norm(), i) );
            diff_end.insert(   std::make_pair((theoreticalLim[i]-end_pt).norm(),   i) );
        }
        //_______________________________________________________________________________________________________________________
        Eigen::Vector3d best_theo_corner_start = theoreticalLim[diff_start.begin()->second];
        std::pair<int, int> theo_pixel_start = pt2XY(best_theo_corner_start, delta_phi, delta_theta, rot_axis, axis_init_phi);

        //_____________________________________________________________________________________________________________________
        Eigen::Vector3d best_theo_corner_end = theoreticalLim[diff_end.begin()->second];
        std::pair<int, int> theo_pixel_end = pt2XY(best_theo_corner_end, delta_phi, delta_theta, rot_axis, axis_init_phi);

        //_____________________________________________________________________________________________________________________
        for(int ki = -rad; ki<=rad; ++ki)
        {
            for(int kj = -rad; kj<=rad; ++kj)
            {
                std::pair<int,int> test_pixel = std::make_pair(theo_pixel_start.first+ki, theo_pixel_start.second+kj);
                if(test_pixel == start_pixel)
                {
                    start_pt = best_theo_corner_start;
                    start = best_theo_corner_start.dot(tangente);
                }
            }
        }
        //_____________________________________________________________________________________________________________________
        for(int ki = -rad; ki<=rad; ++ki)
        {
            for(int kj = -rad; kj<=rad; ++kj)
            {
                std::pair<int,int> test_pixel = std::make_pair(theo_pixel_end.first+ki, theo_pixel_end.second+kj);
                if(test_pixel == end_pixel)
                {
                    end_pt = best_theo_corner_end;
                    end = best_theo_corner_end.dot(tangente);
                }
            }
        }
    }
}


void intersection::computeLim()
{
    std::vector<float> dot;

    //compute real current limits
    //take the farthest points as current limits--------------------------------------------------------------------

    for (auto it_points=points.begin(); it_points != points.end();)
    {
        if(abs(abs(it_points->dot(plane_ref->normal))-plane_ref->distance)<max_plane_distance) // for obstruction cases : we only keep points onto the current plane
        {
            dot.push_back(it_points->dot(tangente)); // project points of intersection onto tangent line
            ++it_points;
        }
        else
        {
            it_points = points.erase(it_points);
            int idx = std::distance(points.begin(), it_points );
            pixels.erase(pixels.begin()+idx);
        }
    }

    int min_idx = std::distance(dot.begin(), std::min_element(dot.begin(), dot.end()) );
    int max_idx = std::distance(dot.begin(), std::max_element(dot.begin(), dot.end()) );

    start = dot[min_idx];
    end =   dot[max_idx];

    start_pt = points[min_idx];
    end_pt =   points[max_idx];

    start_pixel = pixels[min_idx];
    end_pixel =   pixels[max_idx];

//    start_pt = distance * normal + start * tangente;
//    end_pt =   distance * normal + end * tangente;

//    start_pixel = pt2XY(start_pt, delta_phi, delta_theta, rot_axis, axis_init_phi);
//    end_pixel = pt2XY(end_pt, delta_phi, delta_theta, rot_axis, axis_init_phi);

    std::cout<<"visible start from data points = "<<start_pt.transpose()<<std::endl;
    std::cout<<"visible end from data points = "<<end_pt.transpose()<<std::endl<<std::endl;
    std::cout<<"visible start pixel from data points = "<<start_pixel.first<<" "<<start_pixel.second<<std::endl;
    std::cout<<"visible end pixel from data points = "<<end_pixel.first<<" "<<end_pixel.second<<std::endl<<std::endl;
}

void intersection::definePlaneConnection()
{
//if sufficient quantity of current intersection pixels are really the theoritical connection, the intersection is a connection. if not it is an obstruction

    float connect = 0.0;
    int rad = 2;
    int nbr_pts = pixels.size();
    std::vector<int> to_erase;

    for (int i = 0; i < pixels.size(); ++i)
    {
        bool isFound = false;
        int min_ki = std::max(pixels[i].first-rad, 0);
        int max_ki = std::min(pixels[i].first+rad, Nrow-1);
        int min_kj = std::max(pixels[i].second-rad, 0);
        int max_kj = std::min(pixels[i].second+rad, Ncol-1);


        for(int ki = min_ki; ki<=max_ki; ++ki)
        {
            for(int kj = min_kj; kj<=max_kj; ++kj)
            {
                if( theoritical_pixels.find(std::make_pair(ki, kj)) != theoritical_pixels.end() )
                {
                    ++connect;
                    isFound = true;
                    break;
                }
            }
            if(isFound)
                break;
        }

        if(!isFound)
        {
            other_pixels.push_back(pixels[i]);
            other_points.push_back(points[i]);
            to_erase.push_back(i);
        }
    }

    std::cout<<"Number of pixels which contain the real intersection : "<<pixels.size()<<std::endl;

    connect = connect/(float)(nbr_pts);

    std::cout<<"percentage of pixels which contain the theoritical intersection : "<<connect<<std::endl<<std::endl;

    isConnection = false;
    isObstruction = false;
    isOpening = false;
    if( connect > perc_pixels_belonging_to_theoretical) // if more than 0.1 the quantity of pixels of the intersection correspond to a theoretical connection
    {
        //erase points not found in real connection (may be obstructions from the same neighbor plane (example : door -> one connection + one obstruction))
        for (int k = to_erase.size()-1; k >= 0; --k)
        {
            pixels.erase(pixels.begin() + to_erase[k]);
            points.erase(points.begin() + to_erase[k]);
        }
        isConnection = true;
        isLine = true;
        plane_neigh->connected.insert(plane_ref->index);
        plane_ref->connected.insert(plane_neigh->index);
    }
    else
    {
        other_pixels = pixels;
        other_points = points;
        isObstruction = true;    
        //to know which one is obstructing the other
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

void intersection::ProjectLineFeaturesOnPlane(int plane_idx)
{
    //define features of lines in 2D
    ProjectBoundaryOnPlane(points, plane_idx);
    normal2D(0) = (rot * normal)(0);
    normal2D(1) = (rot * normal)(1);
    tangente2D(0) = (rot * tangente)(0);
    tangente2D(1) = (rot * tangente)(1);
    Eigen::Vector2d pt2D = {(rot*pt)(0), (rot*pt)(1)};
    if(pt2D.dot(normal2D)<0)
        normal2D *=-1;
    distance2D = pt2D.dot(normal2D);
}


void intersection::ProjectBoundaryOnPlane(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> pts_to_projet, int plane_idx)
{
    rot = Eigen::Affine3d::Identity();
    rot_inv = Eigen::Affine3d::Identity();

    plane* plane_to_project;

    if(plane_ref->index == plane_idx || plane_idx == 100000)
        plane_to_project = plane_ref;
    else
        plane_to_project = plane_neigh;

    float angle = acos(-plane_to_project->normal(2));
    Eigen::Vector3d axis;
    if(angle>eps)
    {
        axis = (-plane_to_project->normal).cross(Eigen::Vector3d(0,0,1));
        axis /= axis.norm();
        rot.rotate( Eigen::AngleAxisd(angle, axis) );
        rot_inv.rotate( Eigen::AngleAxisd(angle, -axis) );
    }

    points2D.resize(pts_to_projet.size());
    for(int i = 0; i<pts_to_projet.size(); ++i)
    {
        Eigen::Vector3d turned = rot.linear()*pts_to_projet[i];
        turned = turned * (plane_to_project->distance/turned.dot(rot.linear()*(-plane_to_project->normal))); // project points on plane
        points2D[i](0) = turned(0);
        points2D[i](1) = turned(1);
    }

    trans_z = Eigen::Affine3d::Identity();
    trans_z_inv = Eigen::Affine3d::Identity();
    trans_z.translate(Eigen::Vector3d(0.0, 0.0, -plane_to_project->distance));
    trans_z_inv.translate(Eigen::Vector3d(0.0, 0.0, plane_to_project->distance));

    //display stuff -v

//    pcl::PointCloud<pcl::PointXYZ> points_in_line_cloud0;
//    points_in_line_cloud0.width = pts_to_projet.size();
//    points_in_line_cloud0.height = 1;
//    points_in_line_cloud0.is_dense = false;
//    points_in_line_cloud0.points.resize(pts_to_projet.size());

//    for(int k = 0; k<pts_to_projet.size(); ++k)
//    {
//        points_in_line_cloud0.points[k].x = pts_to_projet[k](0);
//        points_in_line_cloud0.points[k].y = pts_to_projet[k](1);
//        points_in_line_cloud0.points[k].z = pts_to_projet[k](2);
//    }

//    pcl::io::savePCDFileASCII ("initial_projected.pcd", points_in_line_cloud0);
//    system("bash pcd2csv.sh initial_projected.pcd");

//    pcl::PointCloud<pcl::PointXYZ> points_in_line_cloud1;
//    points_in_line_cloud1.width = pts_to_projet.size();
//    points_in_line_cloud1.height = 1;
//    points_in_line_cloud1.is_dense = false;
//    points_in_line_cloud1.points.resize(pts_to_projet.size());

//    for(int k = 0; k<pts_to_projet.size(); ++k)
//    {
//        points_in_line_cloud1.points[k].x = points2D[k](0);
//        points_in_line_cloud1.points[k].y = points2D[k](1);
//        points_in_line_cloud1.points[k].z = 0;
//    }

//    pcl::io::savePCDFileASCII ("initial.pcd", points_in_line_cloud1);
//    system("bash pcd2csv.sh initial.pcd");

//    std::cout<<"plane index : "<<plane_to_project->index<<std::endl<<std::endl;
//    getchar();
//    getchar();
//    getchar();
}

void intersection::RANSAC(double error, int tests)
{
    std::cout<<"number of points at the begining: "<<other_points.size()<<std::endl<<std::endl;
    int idx1, idx2;
    int max_points = 0;
    std::map<double, int> proj;

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

        std::map<double, int> proj_temp; // map to keep the link point projection on line / index of point and to order depending on projection

        for(int j = 0; j< points2D.size(); ++j)
        {
            if(abs(points2D[j].dot(normal2D_temp)-distance2D_temp) < error)
                proj_temp.insert(std::make_pair(points2D[j].dot(tangente2D_temp), j));
        }

        //filter indices on line to remove points far from the others----------------------------------------------------------------------------------

        if(proj_temp.size()>3)
        {
            auto last_it = proj_temp.end();
            --last_it;

            for(auto it_proj = proj_temp.begin(); it_proj != last_it; ++it_proj)
            {
                auto it_proj_after = it_proj;
                ++it_proj_after;
                //dist_after is the distance acceptable in pixels
                float dist_after = sqrt(pow(other_pixels[it_proj_after->second].first - other_pixels[it_proj->second].first, 2) + pow(other_pixels[it_proj_after->second].second - other_pixels[it_proj->second].second, 2));
                if(dist_after>20)
                {
                    proj_temp.erase(it_proj_after, proj_temp.end());
                    break;
                }
            }
        }

        //check if the line is the best--------------------------------------------------------------------------------------------------------------------------------------------------------------------

        if(proj_temp.size()>max_points)
        {
            Vector_tested = (points2D[idx1]-points2D[idx2]);
            normal2D = normal2D_temp;
            tangente2D = tangente2D_temp;
            distance2D = distance2D_temp;
            max_points = proj_temp.size();
            proj = proj_temp;
        }
    }
//    std::cout<<"normal2D : "<<normal2D.transpose()<<std::endl<<"tangente2D : "<<tangente2D.transpose()<<std::endl<<"distance2D : "<<distance2D<<std::endl<<"number of points : "<<max_points<<std::endl<<std::endl;

    //Erase , points of line and its pixels neighbors from other_points + add pixels neighbors to points
//-----------------------------------------------------------------------------------------------------------

    auto it_proj_end = proj.end();
    --it_proj_end;
    if(proj.size()>min_number_points_on_line && (it_proj_end->first - proj.begin()->first)>min_line_length)
    {
        struct compare {
            bool operator() (const int a, const int b) const {
                return a > b;
            }
        };
        std::set<int, compare> to_erase;
        std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> pt2D;
        points.clear();
        pixels.clear();
        int index;
        auto it_proj = proj.begin();
        // for first indices of line (limit of line) -> just push the point to the line features
        int rad = 2;
        while(it_proj != proj.end())
        {
            index = it_proj->second;
            int X = other_pixels[index].first;
            int Y = other_pixels[index].second;
            int imin = std::max(0, X-rad);
            int imax = std::min(Nrow-1, X+rad);
            int jmin = std::max(0, Y-rad);
            int jmax = std::min(Ncol-1, Y+rad);

            for(int i = imin; i<=imax; ++i)
            {
                for(int j = jmin; j<=jmax; ++j)
                {
                    auto it_found_pixels = std::find(pixels.begin(), pixels.end(), std::make_pair(i, j));
                    if(it_found_pixels == pixels.end())
                    {
                        auto it_found_other_pixels = std::find(other_pixels.begin(), other_pixels.end(), std::make_pair(i, j));
                        if(it_found_other_pixels != other_pixels.end())
                        {
                            int index = std::distance(other_pixels.begin(), it_found_other_pixels);
                            if(abs(points2D[index].dot(tangente2D)-proj.begin()->first)>line_margin && abs(points2D[index].dot(tangente2D)-it_proj_end->first)>line_margin)
                                to_erase.insert(index);
                            points.push_back(other_points[index]);
                            pixels.push_back(other_pixels[index]);
                            pt2D.push_back(points2D[index]);
                        }
                    }
                }
            }

            ++it_proj;
        }

        std::ofstream f("test_line.csv");
        for(int i = 0; i<other_points.size(); ++i)
                f << other_points[i](0)<<","<<other_points[i](1)<<","<<other_points[i](2)<<"\n";
        f.close();

        f.open("test_to_erase.csv");
        for(auto it_to_erase = to_erase.begin(); it_to_erase!=to_erase.end(); ++it_to_erase)
                f << other_points[*it_to_erase](0)<<","<<other_points[*it_to_erase](1)<<","<<other_points[*it_to_erase](2)<<"\n";
        f.close();

        f.open("points2D.csv");
        for(int i = 0; i<points2D.size(); ++i)
                f << points2D[i](0)<<","<<points2D[i](1)<<","<<0<<"\n";
        f.close();

        f.open("points2D_to_erase.csv");
        for(auto it_to_erase = to_erase.begin(); it_to_erase!=to_erase.end(); ++it_to_erase)
                f << points2D[*it_to_erase](0)<<","<<points2D[*it_to_erase](1)<<","<<0<<"\n";
        f.close();

        f.open("points2D_of_line.csv");
        for(int i = 0; i<pt2D.size(); ++i)
                f << pt2D[i](0)<<","<<pt2D[i](1)<<","<<0<<"\n";
        f.close();


        tangente = {tangente2D(0), tangente2D(1), 0};
        normal = {normal2D(0), normal2D(1), 0};
        distance = distance2D;

        std::ofstream file("test_tangente.csv");
        int n = 0;
        float delta = 0.001;
        while ( (n*delta + proj.begin()->first)<it_proj_end->first)
        {
            file << ((n*delta + proj.begin()->first) * tangente + distance*normal).transpose()<<"\n";
            ++n;
        }
        file.close();


        points2D = pt2D;
        for(auto it_to_erase = to_erase.begin(); it_to_erase!=to_erase.end(); ++it_to_erase)
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

//        //ultimate filter to remove points too far from corrected line
//        //---------------------------------------------------------------------------------------------
//            to_erase.clear();
//            for(auto pt_ptr = points2D.begin(); pt_ptr != points2D.end(); ++pt_ptr)
//            {
//                if(abs(abs((*pt_ptr).dot(normal2D))-distance2D)>error)
//                {
//                    int idx_pt = std::distance(points2D.begin(), pt_ptr);
//                    other_pixels.push_back(pixels[idx_pt]);
//                    other_points.push_back(points[idx_pt]);
//                    to_erase.insert(idx_pt);
//                }
//            }

//            if(to_erase.size()>0)
//            {
//                for(auto it_to_erase = to_erase.begin(); it_to_erase!=to_erase.end(); ++it_to_erase)
//                {
//                    points2D.erase(points2D.begin()  + *it_to_erase);
//                    pixels.erase(pixels.begin() + *it_to_erase);
//                    points.erase(points.begin() + *it_to_erase);
//                }
//            }

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
    }

    std::cout<<"number of remaining points at the end: "<<other_points.size()<<std::endl<<std::endl;
}


//------------------------------------------------Hough-----------------------------------------------------------------------------------------------------

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
