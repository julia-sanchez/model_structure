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

void intersection::SeparatePoints(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& points_of_plane_ref, std::vector<std::pair<int,int>>& pixels_of_plane_ref, std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& points_of_plane_neigh, std::vector<std::pair<int,int>>& pixels_of_plane_neigh)
{
    std::cout<<"Separate points : "<<std::endl<<"indices on ref plane : "<<indices_self.size()<<std::endl;
    std::cout<<"indices on neigh plane : "<<indices_sister.size()<<std::endl<<std::endl;

    int n = 0;
    points_of_plane_ref.resize(indices_self.size());
    pixels_of_plane_ref.resize(indices_self.size());
    for (auto it_indices_self = indices_self.begin(); it_indices_self != indices_self.end(); ++it_indices_self)
    {
        points_of_plane_ref[n] = other_points[*it_indices_self];
        pixels_of_plane_ref[n] = other_pixels[*it_indices_self];
        ++n;
    }

    n = 0;
    points_of_plane_neigh.resize(indices_sister.size());
    pixels_of_plane_neigh.resize(indices_sister.size());
    for (auto it_indices_sister = indices_sister.begin(); it_indices_sister != indices_sister.end(); ++it_indices_sister)
    {
        points_of_plane_neigh[n] = other_points[*it_indices_sister];
        pixels_of_plane_neigh[n] = other_pixels[*it_indices_sister];
        ++n;
    }
}

void intersection::computeTheoriticalLineFeaturesObstruction(std::multimap<std::pair<int,int>,std::pair<int,int>> neighborPix2currentPix)
{
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> points_of_plane_ref;
    std::vector<std::pair<int,int>> pixels_of_plane_ref;
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> points_of_plane_neigh;
    std::vector<std::pair<int,int>> pixels_of_plane_neigh;

    SeparatePoints(points_of_plane_ref, pixels_of_plane_ref, points_of_plane_neigh, pixels_of_plane_neigh);

    std::cout<<"number of points of ref : "<<points_of_plane_ref.size()<<std::endl;
    std::cout<<"number of points of neigh : "<<points_of_plane_neigh.size()<<std::endl;

    if(points_of_plane_ref.size()<min_number_points_on_line || points_of_plane_neigh.size()<min_number_points_on_line)
        return;

    points.clear();
    pixels.clear();
    other_points.clear();
    other_pixels.clear();
    indices_sister.clear();
    indices_self.clear();

    //-----------------------------------------------------------------------------------------------------------------------
    ProjectBoundaryOnPlane(points_of_plane_neigh, plane_neigh);
    std::set<int> indices_line_neigh, repeated_neigh;
    RANSAC(points_of_plane_neigh, pixels_of_plane_neigh, max_line_distance, ransac_iterations, indices_line_neigh, repeated_neigh);

    if(indices_line_neigh.size()<min_number_points_on_line)
    {
        points.clear();
        indices_sister.clear();
        indices_self.clear();
        pixels.clear();
        other_points.clear();
        other_pixels.clear();
        std::cout<<"exit computing obstruction"<<std::endl;
        return;
    }

    auto it_indices_line_neigh = indices_line_neigh.begin();
    auto it_repeated_neigh = repeated_neigh.begin();
    int n = 0;
    indices_sister.clear();
    for(int k = 0; k < pixels_of_plane_neigh.size(); ++k)
    {
        if(k==*it_indices_line_neigh)
        {
            points.push_back(points_of_plane_neigh[k]);
            pixels.push_back(pixels_of_plane_neigh[k]);
            ++it_indices_line_neigh;
            if(it_indices_line_neigh == indices_line_neigh.end())
                --it_indices_line_neigh;
        }
        else
        {
            not_repeated.insert(other_points.size());
            other_points.push_back(points_of_plane_neigh[k]);
            other_pixels.push_back(pixels_of_plane_neigh[k]);
            indices_sister.insert(n);
            ++n;
        }

        if (k==*it_repeated_neigh)
        {
            repeated.insert(other_points.size());
            other_points.push_back(points_of_plane_neigh[k]);
            other_pixels.push_back(pixels_of_plane_neigh[k]);
            ++it_repeated_neigh;
            if(it_repeated_neigh == repeated_neigh.end())
                --it_repeated_neigh;
            indices_sister.insert(n);
            ++n;
        }
    }

    normal_sister = normal;
    tangente_sister = tangente;
    pt_mean_sister = pt_mean;
    distance_sister = distance;

    std::cout<<"neighbor plane of obstruction : points number of line : "<<points.size()<<std::endl;
    std::cout<<"neighbor plane of obstruction : points number remaining : "<<other_points.size()<<std::endl;

    //-----------------------------------------------------------------------------------------------------------------------

    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> points_of_plane_ref_temp;
    std::vector<std::pair<int,int>> pixels_of_plane_ref_temp;
    std::set<int> indices_pixels_of_plane_ref;
    //for all pixels of neigh, add its relative pixels of ref
    for(auto it_neigh_pixel = pixels.begin(); it_neigh_pixel != pixels.end(); ++it_neigh_pixel)
    {
        auto pair_of_it_neighborPix2currentPix = neighborPix2currentPix.equal_range(*it_neigh_pixel);
        for(auto it_neighborPix2currentPix = pair_of_it_neighborPix2currentPix.first; it_neighborPix2currentPix != pair_of_it_neighborPix2currentPix.second; ++it_neighborPix2currentPix )
        {
            auto it_pixels_of_plane_ref = std::find(pixels_of_plane_ref.begin(), pixels_of_plane_ref.end(), it_neighborPix2currentPix->second);
            if(it_pixels_of_plane_ref != pixels_of_plane_ref.end())
            {
                int idx = std::distance(pixels_of_plane_ref.begin(), it_pixels_of_plane_ref);
                indices_pixels_of_plane_ref.insert(idx); // just put it in set because it can appear various times and we just want one index of each
            }
            else
                std::cout<<"correspondance not found in plane_ref points"<<std::endl;
        }
    }

    if(indices_pixels_of_plane_ref.size() == 0)
    {
        points.clear();
        pixels.clear();
        other_points.clear();
        other_pixels.clear();
        indices_sister.clear();
        indices_self.clear();
        return;
    }
    indices_self.clear();

    std::cout<<"number of corresponding  points in reference : "<<indices_pixels_of_plane_ref.size()<<std::endl;
    auto it_indices_pixels_of_plane_ref = indices_pixels_of_plane_ref.begin();
    for(int k = 0; k < pixels_of_plane_ref.size(); ++k)
    {
        if(k == *it_indices_pixels_of_plane_ref)
        {
            pixels_of_plane_ref_temp.push_back(pixels_of_plane_ref[k]);
            points_of_plane_ref_temp.push_back(points_of_plane_ref[k]);
            ++it_indices_pixels_of_plane_ref;
            if(it_indices_pixels_of_plane_ref == indices_pixels_of_plane_ref.end())
                --it_indices_pixels_of_plane_ref;
        }
        else
        {
            other_points.push_back(points_of_plane_ref[k]);
            other_pixels.push_back(pixels_of_plane_ref[k]);
            indices_self.insert(n);
            ++n;
        }
    }

    pixels_of_plane_ref = pixels_of_plane_ref_temp;
    points_of_plane_ref = points_of_plane_ref_temp;

    if(points_of_plane_ref.size()<min_number_points_on_line)
    {
        points.clear();
        pixels.clear();
        other_points.clear();
        other_pixels.clear();
        indices_sister.clear();
        indices_self.clear();
        return;
    }

    std::cout<<"number of points of ref close to neigh : "<<points_of_plane_ref.size()<<std::endl<<std::endl;

    ProjectBoundaryOnPlane(points_of_plane_ref, plane_ref);
    std::set<int> indices_line_ref, repeated_ref;
    RANSAC(points_of_plane_ref, pixels_of_plane_ref, max_line_distance, ransac_iterations, indices_line_ref, repeated_ref);

    if(indices_line_ref.size()<min_number_points_on_line)
    {
        points.clear();
        pixels.clear();
        other_points.clear();
        other_pixels.clear();
        indices_sister.clear();
        indices_self.clear();
        std::cout<<"exit computing obstruction"<<std::endl;
        return;
    }

    auto it_indices_line_ref = indices_line_ref.begin();
    auto it_repeated_ref = repeated_ref.begin();
    for(int k = 0; k < pixels_of_plane_ref.size(); ++k)
    {
        if(k==*it_indices_line_ref)
        {
            points.push_back(points_of_plane_ref[k]);
            pixels.push_back(pixels_of_plane_ref[k]);
            ++it_indices_line_ref;
            if(it_indices_line_ref == indices_line_ref.end())
                --it_indices_line_ref;

        }
        else
        {
            not_repeated.insert(other_points.size());
            other_points.push_back(points_of_plane_ref[k]);
            other_pixels.push_back(pixels_of_plane_ref[k]);
            indices_self.insert(n);
            ++n;
        }

        if (k==*it_repeated_ref)
        {
            repeated.insert(other_points.size());
            other_points.push_back(points_of_plane_ref[k]);
            other_pixels.push_back(pixels_of_plane_ref[k]);
            indices_self.insert(n);
            ++n;
            ++it_repeated_ref;
            if(it_repeated_ref == repeated_ref.end())
                --it_repeated_ref;
        }
    }

    std::cout<<"obstruction tangente : "<<tangente.transpose()<<std::endl;
    std::cout<<"obstruction normal : "<<normal.transpose()<<std::endl;
    std::cout<<"obstruction distance : "<<distance<<std::endl<<std::endl;

    std::cout<<"obstruction sister tangente : "<<tangente_sister.transpose()<<std::endl;
    std::cout<<"obstruction sister normal : "<<normal_sister.transpose()<<std::endl;
    std::cout<<"obstruction sister distance : "<<distance_sister<<std::endl<<std::endl;
}

intersection intersection::export_sister()
{
    intersection inter(plane_neigh, plane_ref, delta_phi, delta_theta);
    inter.isObstruction = true;
    inter.isOpening = false;
    inter.isLine = true;
    inter.other_points = other_points;
    inter.other_pixels = other_pixels;
    inter.points = points;
    inter.pixels = pixels;
    inter.normal = normal_sister;
    inter.tangente = tangente_sister;
    inter.pt_mean = pt_mean_sister;
    inter.distance = distance_sister;

    return inter;
}

void intersection::computeTheoriticalLineFeaturesOpening()
{
    ProjectBoundaryOnPlane(other_points, plane_ref);
    std::set<int> indices_line, repeated;
    RANSAC(other_points, other_pixels, max_line_distance, ransac_iterations, indices_line, repeated);

    if(indices_line.size()==0)
    {
        points.clear();
        pixels.clear();
        other_points.clear();
        other_pixels.clear();
        return;
    }

    auto it_indices_line = indices_line.begin();
    auto it_repeated = repeated.begin();
    points.clear();
    pixels.clear();
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> other_points_temp;
    std::vector<std::pair<int,int>> other_pixels_temp;
    std::set<int> repeated_temp;

    for(int k = 0; k < other_pixels.size(); ++k)
    {
        if(k == *it_indices_line)
        {
            points.push_back(other_points[k]);
            pixels.push_back(other_pixels[k]);
            ++it_indices_line;
            if(it_indices_line == indices_line.end())
                --it_indices_line;
        }
        else
        {
            not_repeated.insert(other_points_temp.size());
            other_points_temp.push_back(other_points[k]);
            other_pixels_temp.push_back(other_pixels[k]);
        }

        if (k==*it_repeated)
        {
            repeated_temp.insert(other_points_temp.size());
            other_points_temp.push_back(other_points[k]);
            other_pixels_temp.push_back(other_pixels[k]);
            ++it_repeated;
            if(it_repeated == repeated.end())
                --it_repeated;
        }
    }
    other_points = other_points_temp;
    other_pixels = other_pixels_temp;

    std::cout<<"opening tangente : "<<tangente.transpose()<<std::endl;
    std::cout<<"opening normal : "<<normal.transpose()<<std::endl;
    std::cout<<"opening distance : "<<distance<<std::endl<<std::endl;
}

void intersection::computeTheoriticalLineFeaturesConnection()
{
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
    distance = abs(pt.dot(normal));

    std::cout<<"connection tangente : "<<tangente.transpose()<<std::endl<<std::endl;
}

void intersection::definePlaneConnection()
{
//if sufficient quantity of current intersection pixels are really the theoritical connection, the intersection is a connection. if not it is an obstruction
    computeTheoriticalPhiTheta();
    isOpening = false;
    float connect = 0.0;
    int rad = 2;
    int nbr_pts = pixels.size();
    std::set<int> to_erase;
    other_points.clear();
    other_pixels.clear();

    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> points_of_line;
    std::vector<std::pair<int,int>> pixels_of_line;
    std::set<int> indices_of_line;

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
                    isFound = true;
                    break;
                }
            }
            if(isFound)
                break;
        }

        if(isFound)
        {
            ++connect;
            indices_of_line.insert(i);
        }
    }

    std::cout<<"Number of pixels which contain the real intersection : "<<indices_of_line.size()<<std::endl;

    connect = connect/(float)(nbr_pts);

    std::cout<<"percentage of pixels which contain the theoritical intersection : "<<connect*100<<"%"<<std::endl<<std::endl;

    isConnection = false;
    isObstruction = false;
    isOpening = false;
    if( connect > perc_pixels_belonging_to_theoretical) // if more than 0.1 the quantity of pixels of the intersection correspond to a theoretical connection
    {
        std::cout<<"Erase pixels not belonging to theoretical lines from pixels"<<std::endl<<std::endl;
        //erase points not found in real connection (may be obstructions from the same neighbor plane (example : door -> one connection + various obstructions))

        points_of_line.resize(indices_of_line.size());
        pixels_of_line.resize(indices_of_line.size());
        int n = 0;
        int n_other = 0;
        auto it_indices_of_line = indices_of_line.begin();
        auto it_indices_self = indices_self.begin();
        auto it_indices_sister = indices_sister.begin();
        std::set<int> indices_self_temp;
        std::set<int> indices_sister_temp;

        for(int k = 0; k< pixels.size(); ++k)
        {
            if(k == *it_indices_of_line)
            {
                points_of_line[n] = points[k];
                pixels_of_line[n] = pixels[k];
                ++n;
                ++it_indices_of_line;
                if(it_indices_of_line == indices_of_line.end())
                    --it_indices_of_line;
            }
            else
            {
                other_points.push_back(points[k]);
                other_pixels.push_back(pixels[k]);

                if(k == *it_indices_sister)
                    indices_sister_temp.insert(n_other);
                else if(k == *it_indices_self)
                    indices_self_temp.insert(n_other);
                else
                    std::cout<<"no indice for neigh or ref"<<std::endl;
                ++n_other;
            }

            if(k == *it_indices_sister)
            {
                ++it_indices_sister;
                if(it_indices_sister == indices_sister.end())
                    --it_indices_sister;
            }
            else if(k == *it_indices_self)
            {
                ++it_indices_self;
                if(it_indices_self == indices_self.end())
                    --it_indices_self;
            }
        }

        points = points_of_line;
        pixels = pixels_of_line;
        indices_self = indices_self_temp;
        indices_sister = indices_sister_temp;

        std::cout<<"points remaining of neghbor plane : "<<indices_sister.size()<<std::endl;
        std::cout<<"points remaining of current ref plane : "<<indices_self.size()<<std::endl;

        pt_mean = Eigen::Vector3d::Zero();
        for (int i = 0; i < points.size(); ++i)
            pt_mean += points[i];

        pt_mean /= points.size();

        isConnection = true;
        isLine = true;
        plane_neigh->connected.insert(plane_ref->index);
        plane_ref->connected.insert(plane_neigh->index);
        std::vector<float> proj;

        for (int i = 0; i < points.size(); ++i)
            proj.push_back(points[i].dot(tangente));

        float min_proj = *std::min(proj.begin(), proj.end());
        float max_proj = *std::max(proj.begin(), proj.end());

        //remove closest points from points of line ---> to change : je veux faire le for sur other_pixels (il doit y en avoir moins) si ya un voisin dans pixel je veux remettre le other_pixel dans pixel
        for (int i = 0; i < pixels.size(); ++i)
        {
            int min_ki = std::max(pixels[i].first-rad, 0);
            int max_ki = std::min(pixels[i].first+rad, Nrow-1);
            int min_kj = std::max(pixels[i].second-rad, 0);
            int max_kj = std::min(pixels[i].second+rad, Ncol-1);

            for(int ki = min_ki; ki<=max_ki; ++ki)
            {
                for(int kj = min_kj; kj<=max_kj; ++kj)
                {
                    auto it_other_pixel_found = std::find(other_pixels.begin(), other_pixels.end(), std::make_pair(ki,kj));
                    if(it_other_pixel_found != other_pixels.end())
                    {
                        int idx = std::distance(other_pixels.begin(), it_other_pixel_found);
                        if(abs(proj[idx]-min_proj)>line_margin && abs(proj[idx]-max_proj)>line_margin)
                            to_erase.insert(idx);
                    }
                }
            }
        }

        std::cout<<"Add neighbors pixels of current selected pixels to pixels vector"<<std::endl<<std::endl;

        if(to_erase.size() > 0 && indices_self.size()>0 && indices_sister.size()>0)
        {
            n_other = 0;
            it_indices_self = indices_self.begin();
            it_indices_sister = indices_sister.begin();
            indices_self_temp.clear();
            indices_sister_temp.clear();
            auto it_to_erase = to_erase.begin();
            std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> other_points_temp;
            std::vector<std::pair<int,int>> other_pixels_temp;

            for(int k = 0; k < other_points.size(); ++k)
            {
                if(k == *it_to_erase)
                {
                     points.push_back(other_points[k]);
                     pixels.push_back(other_pixels[k]);
                     ++it_to_erase;
                     if(it_to_erase == to_erase.end())
                         --it_to_erase;
                }
                else
                {
                    other_points_temp.push_back(other_points[k]);
                    other_pixels_temp.push_back(other_pixels[k]);
                    if(k == *it_indices_sister)
                        indices_sister_temp.insert(n_other);
                    if(k == *it_indices_self)
                        indices_self_temp.insert(n_other);
                    ++n_other;
                }

                if(k == *it_indices_sister)
                {
                    ++it_indices_sister;
                    if(it_indices_sister == indices_sister.end())
                        --it_indices_sister;
                }
                if(k == *it_indices_self)
                {
                    ++it_indices_self;
                    if(it_indices_self == indices_self.end())
                        --it_indices_self;
                }
            }

            other_points = other_points_temp;
            other_pixels = other_pixels_temp;
            indices_self = indices_self_temp;
            indices_sister = indices_sister_temp;

            if(indices_sister_temp.size() < min_number_points_on_line || indices_self_temp.size() < min_number_points_on_line)
            {
                other_points.clear();
                other_pixels.clear();
                indices_self.clear();
                indices_sister.clear();
            }
        }
    }
    else
    {
        other_pixels = pixels;
        other_points = points;
        pixels.clear();
        points.clear();
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

void intersection::computeTheoriticalPhiTheta()
{
    Eigen::Vector3d normal1 = -plane_ref->normal;
    Eigen::Vector3d normal2 = -plane_neigh->normal;
    std::cout<<"normal ref : "<<normal1.transpose()<< " distance ref : "<<plane_ref->distance<<std::endl;
    std::cout<<"normal neigh : "<<normal2.transpose()<< " distance neigh : "<<plane_neigh->distance<<std::endl<<std::endl;
    float D = plane_ref->distance/plane_neigh->distance;
    // Fix phi and compute theta
    for(int i = 0; i<Nrow; ++i)
    {
        float phi = i*delta_phi + delta_phi/2;
        float A1 = normal1(0) * cos(phi);
        float B1 = normal1(1) * sin(phi);
        float C1 = normal1(2);
        float A2 = normal2(0) * cos(phi);
        float B2 = normal2(1) * sin(phi);
        float C2 = normal2(2);
        float theta = atan( (D*C2-C1)/ (A1+B1 - D*(A2+B2)) );

        while(theta<0)
            theta += 2*M_PI;

        Eigen::Vector3d pt_norm = {cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)};
        float distance1 = plane_ref->distance/normal1.dot(pt_norm);
        float distance2 = plane_neigh->distance/normal2.dot(pt_norm);

        if(abs(distance1 - distance2) < 0.001)
        {
            if(theta<M_PI)
            {
                theoritical_phi_theta.push_back(std::make_pair(phi,theta));
                int X = i;
                int Y = (int)(theta/delta_theta);
                theoritical_pixels.insert(std::make_pair(X,Y));
            }
        }

        theta -= M_PI;
        while(theta<0)
            theta += 2*M_PI;

        pt_norm = {cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)};
        distance1 = plane_ref->distance/normal1.dot(pt_norm);
        distance2 = plane_neigh->distance/normal2.dot(pt_norm);

        if(abs(distance1 - distance2) < 0.001)
        {
            if(theta<M_PI)
            {
                theoritical_phi_theta.push_back(std::make_pair(phi,theta));
                int X = i;
                int Y = (int)(theta/delta_theta);
                theoritical_pixels.insert(std::make_pair(X,Y));
            }
        }
    }

    // Fix theta and compute phi
    float alpha = atan2(normal1(1)-D*normal2(1),normal1(0)-D*normal2(0));
    for(int j = 0; j<Ncol; ++j)
    {
        float theta = j*delta_theta + delta_theta/2;
        float A1 = normal1(0) * sin(theta);
        float B1 = normal1(1) * sin(theta);
        float C1 = normal1(2) * cos(theta);
        float A2 = normal2(0) * sin(theta);
        float B2 = normal2(1) * sin(theta);
        float C2 = normal2(2) * cos(theta);
        float a = A1 - D*A2;
        float b = B1 - D*B2;
        float c = D*C2 - C1;
        float R = sqrt(a*a + b*b);
        if(abs(c/R)<1)
        {
            float phi = acos(c/R)+alpha;
            Eigen::Vector3d pt_norm = {cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)};
            float distance1 = plane_ref->distance/normal1.dot(pt_norm);
            float distance2 = plane_neigh->distance/normal2.dot(pt_norm);

            if(abs(distance1 - distance2) < 0.001)
            {
                while(phi<0)
                    phi += 2*M_PI;
                theoritical_phi_theta.push_back(std::make_pair(phi,theta));
                int X = (int)(phi/delta_phi);
                int Y = j;
                theoritical_pixels.insert(std::make_pair(X,Y));
            }

            phi = -acos(c/R)+alpha;

            pt_norm = {cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)};
            distance1 = plane_ref->distance/normal1.dot(pt_norm);
            distance2 = plane_neigh->distance/normal2.dot(pt_norm);

            if(abs(distance1 - distance2) < 0.001)
            {
                while(phi<0)
                    phi += 2*M_PI;
                theoritical_phi_theta.push_back(std::make_pair(phi,theta));
                int X = (int)(phi/delta_phi);
                int Y = j;
                theoritical_pixels.insert(std::make_pair(X,Y));
            }
        }
    }

//    if(plane_ref->index == 2 && plane_neigh->index == 8)
//    {
//        Eigen::MatrixXi image_test = Eigen::MatrixXi::Zero(Nrow, Ncol);
//        for(auto it = theoritical_pixels.begin(); it != theoritical_pixels.end(); ++it)
//            image_test(it->first, it->second) = 1;
//        save_image_pgm("theoretical_line", "", image_test, 1);
//        getchar();
//        getchar();
//        getchar();
//    }
}

void intersection::selectCorners(int rad, float delta_phi, float delta_theta, Eigen::Vector3d rot_axis, Eigen::Vector3d axis_init_phi)
{
    //compare with theoretical corners

    std::map<float, int> diff_start;
    std::map<float, int> diff_end;

    replaced_start = false;
    replaced_end = false;

    std::cout<<"intersection between planes :"<<plane_ref->index <<" and "<<plane_neigh->index<<" initial limits are : "<<start_pt.transpose()<<"      "<<end_pt.transpose()<<std::endl;

    if(theoreticalLim.size() != 0)
    {
        for(int i = 0; i<theoreticalLim.size(); ++i)
        {
            diff_start.insert( std::make_pair((theoreticalLim[i]-start_pt).norm(), i) );   //difference between potential lims and start_pt
            diff_end.insert(   std::make_pair((theoreticalLim[i]-end_pt).norm(),   i) );
        }

        //we extract the closest potential limit

        //_______________________________________________________________________________________________________________________
        int corner_start_idx = possible_corner_indices[diff_start.begin()->second];
        Eigen::Vector3d best_theo_corner_start = theoreticalLim[diff_start.begin()->second];
        std::pair<int, int> theo_pixel_start = pt2XY(best_theo_corner_start, delta_phi, delta_theta, rot_axis, axis_init_phi);

        //_____________________________________________________________________________________________________________________
        int corner_end_idx = possible_corner_indices[diff_end.begin()->second];
        Eigen::Vector3d best_theo_corner_end = theoreticalLim[diff_end.begin()->second];
        std::pair<int, int> theo_pixel_end = pt2XY(best_theo_corner_end, delta_phi, delta_theta, rot_axis, axis_init_phi);

        //_____________________________________________________________________________________________________________________
        bool start_has_close_theoretical = true;
        bool end_has_close_theoretical = true;
        if(corner_start_idx == corner_end_idx)  // if both bounds have the same theoretical correspondance
        {
            if(diff_start.begin()->first<diff_end.begin()->first)
            {
                start_has_close_theoretical = true;
                end_has_close_theoretical = false;
            }
            else
            {
                start_has_close_theoretical = false;
                end_has_close_theoretical = true;
            }
        }

        if((theoreticalLim[diff_start.begin()->second]-end_pt).norm()<diff_start.begin()->first) // if end_pt is closer to the supposed correspondance of start_pt than start_pt
            start_has_close_theoretical = false;

        if((theoreticalLim[diff_end.begin()->second]-start_pt).norm()<diff_end.begin()->first) // if end_pt is closer to the supposed correspondance of start_pt than start_pt
            end_has_close_theoretical = false;

        //take the decision if the closest potential lim is the real lim : if the potential lim is a close pixel

        if(start_has_close_theoretical)
        {
            for(int ki = -rad; ki<=rad; ++ki)
            {
                int kj;
                for(kj = -rad; kj<=rad; ++kj)
                {
                    std::pair<int,int> test_pixel = std::make_pair(theo_pixel_start.first+ki, theo_pixel_start.second+kj);
                    if(std::find(pixels.begin(), pixels.end(), test_pixel) != pixels.end() || std::find(other_pixels.begin(), other_pixels.end(), test_pixel) != other_pixels.end())
                    {
                        start_pt = best_theo_corner_start;
                        start = best_theo_corner_start.dot(tangente);
                        replaced_start = true;
                        corner_indices.push_back(corner_start_idx);
                        break;
                    }
                }
                if(replaced_start)
                    break;
            }
        }
        //_____________________________________________________________________________________________________________________
        if(end_has_close_theoretical)
        {
            for(int ki = -rad; ki<=rad; ++ki)
            {
                int kj;
                for(kj = -rad; kj<=rad; ++kj)
                {
                    std::pair<int,int> test_pixel = std::make_pair(theo_pixel_end.first+ki, theo_pixel_end.second+kj);
                    if(std::find(pixels.begin(), pixels.end(), test_pixel) != pixels.end() || std::find(other_pixels.begin(), other_pixels.end(), test_pixel) != other_pixels.end())
                    {
                        end_pt = best_theo_corner_end;
                        end = best_theo_corner_end.dot(tangente);
                        replaced_end = true;
                        corner_indices.push_back(corner_end_idx);
                        break;
                    }
                }
                if(replaced_end)
                    break;
            }
        }

        if(start>end)
            std::cout<<"intersection between planes :"<<plane_ref->index <<" and "<<plane_neigh->index<<" problem is : start > end"<<std::endl;
        if(start==end)
            std::cout<<"intersection between planes :"<<plane_ref->index <<" and "<<plane_neigh->index<<" problem is : start = end"<<std::endl;
        if(corner_start_idx==corner_end_idx)
            std::cout<<"intersection between planes :"<<plane_ref->index <<" and "<<plane_neigh->index<<" the two lims of line have the same correspondance with theoretical corner"<<std::endl;

    }

    std::cout<<"intersection between planes :"<<plane_ref->index <<" and "<<plane_neigh->index<<" final limits are : "<<start_pt.transpose()<<"      "<<end_pt.transpose()<<std::endl<<std::endl;
}


void intersection::computeLim()
{
    std::set<float> dot;
    //compute real current limits
    //take the farthest points as current limits--------------------------------------------------------------------

    for (auto it_points=points.begin(); it_points != points.end(); ++it_points)
    {
        if(abs(abs(it_points->dot(plane_ref->normal))-plane_ref->distance)<max_plane_distance) // for obstruction cases : we only keep points onto the current plane
            dot.insert(it_points->dot(tangente)); // project points of intersection onto tangent line
    }

    if(dot.size() == 0)
    {
        for(auto it_pix = pixels.begin(); it_pix != pixels.end(); ++it_pix)
            std::cout<<it_pix->first<<" "<<it_pix->second<<std::endl;
        std::cout<<"no points of intersection on plane"<<std::endl<<std::endl;
    }

    auto dot_end = dot.end();
    --dot_end;

    start = *dot.begin();
    end =   *dot_end;

    start_pt = distance * normal + start * tangente;
    end_pt =   distance * normal + end * tangente;

    std::cout<<"visible start from data points = "<<start_pt.transpose()<<std::endl;
    std::cout<<"visible end from data points = "<<end_pt.transpose()<<std::endl<<std::endl;
}


void intersection::ProjectLineFeaturesOnPlane(plane* p)
{
    //define features of lines in 2D
    ProjectBoundaryOnPlane(points, p);
    normal2D(0) = (rot * normal)(0);
    normal2D(1) = (rot * normal)(1);
    tangente2D(0) = (rot * tangente)(0);
    tangente2D(1) = (rot * tangente)(1);
    Eigen::Vector2d pt2D = {(rot*pt)(0), (rot*pt)(1)};
    if(pt2D.dot(normal2D)<0)
        normal2D *=-1;
    distance2D = pt2D.dot(normal2D);
}


void intersection::ProjectBoundaryOnPlane(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> pts_to_projet, plane* plane_to_project)
{
    rot = Eigen::Affine3d::Identity();
    rot_inv = Eigen::Affine3d::Identity();

    float angle = acos(-plane_to_project->normal(2));   //angle between normal and z axis
    Eigen::Vector3d axis;
    if(angle>eps)
    {
        axis = (-plane_to_project->normal).cross(Eigen::Vector3d(0,0,1)); // rotation axis to align normal onto z axis
        axis /= axis.norm();
        rot.rotate( Eigen::AngleAxisd(angle, axis) );
        rot_inv.rotate( Eigen::AngleAxisd(angle, -axis) );
    }

//    pcl::PointCloud<pcl::PointXYZ> points_in_line_cloud0;
//    points_in_line_cloud0.width = pts_to_projet.size();
//    points_in_line_cloud0.height = 1;
//    points_in_line_cloud0.is_dense = false;
//    points_in_line_cloud0.points.resize(pts_to_projet.size());

    points2D.resize(pts_to_projet.size());
    for(int i = 0; i<pts_to_projet.size(); ++i)
    {
        Eigen::Vector3d turned = pts_to_projet[i] * (plane_to_project->distance/pts_to_projet[i].dot(-plane_to_project->normal)); // project points on plane
        turned = rot.linear()*turned; //turn it
        points2D[i](0) = turned(0);
        points2D[i](1) = turned(1);

//        points_in_line_cloud0.points[i].x = turned(0);
//        points_in_line_cloud0.points[i].y = turned(1);
//        points_in_line_cloud0.points[i].z = turned(2);
    }

    trans_z = Eigen::Affine3d::Identity();
    trans_z_inv = Eigen::Affine3d::Identity();
    trans_z.translate(Eigen::Vector3d(0.0, 0.0, -plane_to_project->distance));
    trans_z_inv.translate(Eigen::Vector3d(0.0, 0.0, plane_to_project->distance));

    //display stuff -v

//    if(plane_to_project->index == 4)
//    {
////        pcl::PointCloud<pcl::PointXYZ> points_in_line_cloud0;
////        points_in_line_cloud0.width = pts_to_projet.size();
////        points_in_line_cloud0.height = 1;
////        points_in_line_cloud0.is_dense = false;
////        points_in_line_cloud0.points.resize(pts_to_projet.size());

////        for(int k = 0; k<pts_to_projet.size(); ++k)
////        {
////            points_in_line_cloud0.points[k].x = pts_to_projet[k](0);
////            points_in_line_cloud0.points[k].y = pts_to_projet[k](1);
////            points_in_line_cloud0.points[k].z = pts_to_projet[k](2);
////        }

//        pcl::io::savePCDFileASCII ("intersection_for_RANSAC.pcd", points_in_line_cloud0);
//        system("bash pcd2csv.sh intersection_for_RANSAC.pcd");

//        pcl::PointCloud<pcl::PointXYZ> points_in_line_cloud1;
//        points_in_line_cloud1.width = pts_to_projet.size();
//        points_in_line_cloud1.height = 1;
//        points_in_line_cloud1.is_dense = false;
//        points_in_line_cloud1.points.resize(pts_to_projet.size());

//        for(int k = 0; k<pts_to_projet.size(); ++k)
//        {
//            points_in_line_cloud1.points[k].x = points2D[k](0);
//            points_in_line_cloud1.points[k].y = points2D[k](1);
//            points_in_line_cloud1.points[k].z = 0;
//        }

//        pcl::io::savePCDFileASCII ("projection.pcd", points_in_line_cloud1);
//        system("bash pcd2csv.sh projection.pcd");

//        std::cout<<"plane index : "<<plane_to_project->index<<std::endl<<std::endl;

//        getchar();
//        getchar();
//        getchar();
//    }
}

void intersection::RANSAC(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& points_init, std::vector<std::pair<int,int>>& pixels_init, double error, int tests, std::set<int>& indices_line, std::set<int>& repeated )
{
    std::cout<<"number of points at the begining: "<<points_init.size()<<std::endl<<std::endl;
    int idx1, idx2;
    int max_points = 0;
    std::map<double, int> proj;

    std::srand (std::time(NULL));

    Eigen::Vector2d Vector_tested;

    //to limit number of tests if few points-----
    int combinaisons = (int)((float)(points_init.size() * points_init.size()-1)/2.0);
    tests = std::min(5*combinaisons, tests);
    //-----//-----//-----//-----//-----//-----//-----

    for(int i = 0; i< tests; ++i)
    {
        idx1 = rand() % points_init.size();// [0; other_points.size()[
        idx2 = rand() % points_init.size();

        while(idx2 == idx1)
            idx2 = rand() % points_init.size();

        if((points2D[idx1]-points2D[idx2]).norm()<1e-3)
            continue;
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
                //dist_after is the distance acceptable in pixels->
                float dist_after = sqrt(pow(pixels_init[it_proj_after->second].first - pixels_init[it_proj->second].first, 2) + pow(pixels_init[it_proj_after->second].second - pixels_init[it_proj->second].second, 2));
                if(dist_after>max_dist_between_pixels_in_line)
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
    auto start_proj = proj.begin();
    auto end_proj = proj.end();
    --end_proj;
    length = abs(end_proj->first - start_proj->first);

    if(proj.size()>min_number_points_on_line && length>min_line_length)
    {
        std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> pt2D;
        int proj_idx;
        auto it_proj = proj.begin();
        int rad = 2;

        //on remplit indices_line et repeated
        while(it_proj != proj.end())
        {
            proj_idx = it_proj->second;
            int X = pixels_init[proj_idx].first;
            int Y = pixels_init[proj_idx].second;
            int imin = std::max(0, X-rad);
            int imax = std::min(Nrow-1, X+rad);
            int jmin = std::max(0, Y-rad);
            int jmax = std::min(Ncol-1, Y+rad);

            for(int i = imin; i<=imax; ++i)
            {
                for(int j = jmin; j<=jmax; ++j)
                {
                    auto it_found_pixels_init = std::find(pixels_init.begin(), pixels_init.end(), std::make_pair(i, j));
                    if(it_found_pixels_init != pixels_init.end())
                    {
                        int pixels_init_idx = std::distance(pixels_init.begin(), it_found_pixels_init);
                        indices_line.insert(pixels_init_idx);
                        if(abs(points2D[pixels_init_idx].dot(tangente2D)-start_proj->first)<line_margin || abs(points2D[pixels_init_idx].dot(tangente2D)-end_proj->first)<line_margin)
                            repeated.insert(pixels_init_idx);
                    }
                }
            }
            ++it_proj;
        }


        for(auto it_indices_line = indices_line.begin(); it_indices_line!=indices_line.end(); ++it_indices_line)
            pt2D.push_back(points2D[*it_indices_line]); //not in previous loop because indices would have been repeated when looking for neighbors...

    //        least square for accuracy
    //---------------------------------------------------------------------------------------------
        float theta = acos(-normal2D(0));

        Eigen::MatrixXd JTJ(2,2);
        Eigen::MatrixXd JTr(2, 1);
        Eigen::MatrixXd J(2,1);
        double r;
        int nbr_iterations = 10;
        bool too_thin = true;
        Eigen::Vector2d pt2D_mean;
        for(int k =0; k<nbr_iterations; ++k)
        {
            JTJ.setZero();
            JTr.setZero();
            too_thin = true;
            int n = 0;
            pt2D_mean = Eigen::Vector2d::Zero();
            for(int i = 0; i< pt2D.size(); ++i)
            {
                if(abs(pt2D[i].dot(tangente2D)-proj.begin()->first)>line_margin && abs(pt2D[i].dot(tangente2D)-end_proj->first)>line_margin)
                {
                    J(0,0) = pt2D[i].dot(Eigen::Vector2d(sin(theta), cos(theta)));
                    J(1,0) = -1;
                    r = pt2D[i].dot(normal2D)-distance2D;
                    JTJ += J * J.transpose();
                    JTr += J * r;
                    too_thin = false;
                    pt2D_mean += pt2D[i];
                    ++n;
                }
            }
            pt2D_mean /= n;

            if(too_thin == false)
            {
                Eigen::MatrixXd result(2, 1);
                result = -JTJ.llt().solve(JTr);
                theta += result(0);
                normal2D = {-cos(theta), sin(theta)};
                distance2D += result(1);
                tangente2D = {sin(theta), cos(theta)};
                if(distance2D<0)
                {
                    distance2D *= -1;
                    normal2D *= -1;
                }
            }
            else
            {
                indices_line.clear();
                std::cout<<"points only in the 2cm extremities"<<std::endl<<std::endl;
                return;
            }
        }


        //convert to 3D
        //---------------------------------------------------------------------------------------------

        tangente(0) = tangente2D(0);
        tangente(1) = tangente2D(1);
        tangente(2) = 0;
        tangente = rot_inv.linear()*tangente;
        pt_mean = {pt2D_mean(0), pt2D_mean(1), 0};
        pt_mean = rot_inv.linear() * (pt_mean + trans_z_inv.translation());
        normal = pt_mean - pt_mean.dot(tangente)*tangente;
        normal /= normal.norm();
        distance = abs(pt_mean.dot(normal));
        points2D = pt2D;
    }
    else
    {
        std::cout<< "line too thin or not enough points"<<std::endl;
        return;
    }

    std::cout<<"number of remaining points at the end: "<<repeated.size() + (points_init.size() - indices_line.size())<<std::endl<<std::endl;
    std::cout<<"number of remaining points at the end without repeated: "<<points_init.size() - indices_line.size()<<std::endl<<std::endl;
}


//------------------------------------------------Hough-----------------------------------------------------------------------------------------------------

//void intersection::Hough(double error_angle, double error_distance)
//{
//    std::cout<<"Hough : Number of initial points in which we seek a line: "<<other_points.size()<<std::endl<<std::endl;
//    double max_theta = 2*M_PI;
//    double min_theta = 0;
//    double delta_angle = error_angle*2;
//    int Nhist_theta = std::ceil( (max_theta-min_theta + eps)/delta_angle );
//    double delta_distance = error_distance * 2;
//    std::vector<double> thetas(Nhist_theta);
//    std::vector<double> distances(points2D.size() * Nhist_theta);

//    normal2D.Zero();

//    for (int n = 0; n<Nhist_theta; ++n)
//        thetas[n] = min_theta + n*delta_angle + delta_angle/2;

//    //define all possible distances for one point and one theta
//    for (int k = 0; k<points2D.size(); ++k)
//    {
//        for (int n = 0; n<Nhist_theta; ++n)
//        {
//            normal2D = {-sin(thetas[n]), cos(thetas[n])};
//            distances[k*Nhist_theta +n] = points2D[k].dot(normal2D);
//        }
//    }

//    //fill histogram with theta/distance
//    double max_distance = *(std::max_element(distances.begin(), distances.end()));
//    double min_distance = *(std::min_element(distances.begin(), distances.end()));
//    int Nhist_distance = std::ceil((max_distance - min_distance + eps)/delta_distance);

//    Eigen::MatrixXi hist = Eigen::MatrixXi::Zero(Nhist_theta, Nhist_distance);
//    std::multimap<std::pair<int, int>, int> angle_distance2boundary_points;
//    for (int k = 0; k<points2D.size(); ++k)
//    {
//        for (int n = 0; n<Nhist_theta; ++n)
//        {
//            if(distances[k*Nhist_theta + n]>0)
//            {
//                ++hist((int)((thetas[n]-min_theta)/delta_angle), (int)((distances[k*Nhist_theta + n]-min_distance)/delta_distance));
//                angle_distance2boundary_points.insert(std::make_pair(std::make_pair((int)((thetas[n]-min_theta)/delta_angle), (int)((distances[k*Nhist_theta + n]-min_distance)/delta_distance)), k));
//            }
//        }
//    }

////    save_image_pgm("Boundary_Hough", "", hist, 100);

//    //search interesting bins in histogram

//    int i;
//    int j;
//    std::vector<int> indices_on_line;

//    //put points corresponding to line in points_in_line
//    auto points_found_idx = angle_distance2boundary_points.equal_range(std::make_pair(i,j));

//    for(auto it = points_found_idx.first; it !=points_found_idx.second; ++it)
//        indices_on_line.push_back(it->second);

//    //---------------------------------------------------------------------------------------------------------------------------------

//    // compute line features (theta distance) with bin features (coord i and coord j)
//    double theta = i * delta_angle + min_theta + delta_angle/2;
//    distance2D = j * delta_distance + min_distance + delta_distance/2;
//    tangente2D = {cos(theta), sin(theta)};
//    normal2D = {-sin(theta), cos(theta)};
////    if(normal2D.dot(points2D[indices_on_line[0]])<0)
////        normal2D *= -1;

//    points.resize(indices_on_line.size());
//    pixels.resize(indices_on_line.size());
//    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> pt2D(indices_on_line.size());
//    for(int i = 0; i < indices_on_line.size(); ++i)
//    {
//        points[i] = other_points[indices_on_line[i]];
//        pixels[i] = other_pixels[indices_on_line[i]];
//        pt2D[i] = points2D[indices_on_line[i]];
//    }

//    std::cout<<"Number of points found in line : "<<indices_on_line.size()<<std::endl<<std::endl;

//    //Erase , points of line and its pixels neighbors from other_points + add pixels neighbors to points
//    std::set<int> to_erase;
//    int rad = 1;
//    for(int k = 0; k<indices_on_line.size(); ++k)
//    {
//        int X = other_pixels[indices_on_line[k]].first;
//        int Y = other_pixels[indices_on_line[k]].second;
//        int imin = std::max(0, X-rad);
//        int imax = std::min(Nrow-1, X+rad);
//        int jmin = std::max(0, Y-rad);
//        int jmax = std::min(Ncol-1, Y+rad);

//        for(int i = imin; i<=imax; ++i)
//        {
//            for(int j = jmin; j<=jmax; ++j)
//            {
//                auto it_found_pixels = std::find(other_pixels.begin(), other_pixels.end(), std::make_pair(i, j));
//                // Get index of element from iterator
//                if(it_found_pixels != other_pixels.end())
//                {
//                    int other_pixels_idx = std::distance(other_pixels.begin(), it_found_pixels);
//                    to_erase.insert(other_pixels_idx);
//                    points.push_back(other_points[other_pixels_idx]);
//                    pixels.push_back(other_pixels[other_pixels_idx]);
//                    pt2D.push_back(points2D[other_pixels_idx]);
//                }
//            }
//        }
//    }

//    points2D = pt2D;

//    auto it_to_erase_start = to_erase.end();
//    --it_to_erase_start;
//    auto it_to_erase_end = to_erase.begin();
//    --it_to_erase_end;
//    for(auto it_to_erase = it_to_erase_start; it_to_erase!=it_to_erase_end; --it_to_erase)
//    {
//        other_points.erase(other_points.begin() + *it_to_erase);
//        other_pixels.erase(other_pixels.begin() + *it_to_erase);
//    }

//    //---------------------------------------------------------------------------------------------
//   Eigen::MatrixXd JTJ(2,2);
//    Eigen::MatrixXd JTr(2, 1);
//    Eigen::MatrixXd J(2,1);
//    double r;
//    int nbr_iterations = 10;
//    for(int k =0; k<nbr_iterations; ++k)
//    {
//        JTJ.setZero();
//        JTr.setZero();
//        for(int i = 0; i< points2D.size(); ++i)
//        {
//            J(0,0) = points2D[i].dot(Eigen::Vector2d(-cos(theta), -sin(theta)));
//            J(1,0) = -1;
//            r = points2D[i].dot(normal2D)-distance2D;
//            JTJ += J * J.transpose();
//            JTr += J * r;
//        }

//        Eigen::MatrixXd result(2, 1);
//        result = -JTJ.llt().solve(JTr);
//        theta += result(0);
//        normal2D = {-sin(theta), cos(theta)};
//        distance2D += result(1);
//    }
//    tangente2D = {cos(theta), sin(theta)};
//    if(distance2D<0)
//    {
//        distance2D *= -1;
//        normal2D *= -1;
//    }

//    //---------------------------------------------------------------------------------------------

//    //convert to 3D

//    //retransform in xy:
//    distance2D = (distance2D*normal2D+pt_mean).dot(normal2D);
//    tangente(0) = tangente2D(0);
//    tangente(1) = tangente2D(1);
//    tangente(2) = 0;
//    tangente = rot_inv.linear()*tangente;
//    normal(0) = normal2D(0);
//    normal(1) = normal2D(1);
//    normal(2) = 0;
//    normal = rot_inv.linear()*normal;
//    pt = plane_ref->distance*(-plane_ref->normal) + normal * distance2D;
//    normal = pt - pt.dot(tangente)*tangente;
//    normal /= normal.norm();
//    distance = pt.dot(normal);

//    std::cout<<"number of remaining points at the end: "<<other_points.size()<<std::endl<<std::endl;
//}

void intersection::correctObstructions(intersection& sister)
{
    std::cout<<"normal of self : "<<normal.transpose()<<std::endl;
    std::cout<<"normal of sister : "<<sister.normal.transpose()<<std::endl;

    Eigen::Vector3d normal_plane_ref = -plane_ref->normal;
    Eigen::Vector3d normal_plane_neigh = -plane_neigh->normal;

    //compute mean point projection on each plane
    Eigen::Vector3d pt_mean_sister_proj = ( plane_ref->distance/(sister.pt_mean.dot(normal_plane_ref)) ) * sister.pt_mean;
    Eigen::Vector3d pt_mean_self_proj = ( plane_neigh->distance/(     pt_mean.dot(normal_plane_neigh)) ) * pt_mean;

    //compute mean point + tangente projection on each plane
    Eigen::Vector3d pt_mean_and_tangente_sister_proj = ( plane_ref->distance/( (sister.pt_mean+sister.tangente).dot(normal_plane_ref)) ) * ( sister.pt_mean+sister.tangente );
    Eigen::Vector3d pt_mean_and_tangente_self_proj = ( plane_neigh->distance/( (pt_mean + tangente).dot(normal_plane_neigh)) ) * ( pt_mean + tangente );

    //deduce tangente projection on each plane
    Eigen::Vector3d tangente_sister_proj = pt_mean_and_tangente_sister_proj - pt_mean_sister_proj;
    tangente_sister_proj /= tangente_sister_proj.norm();
    Eigen::Vector3d tangente_self_proj = pt_mean_and_tangente_self_proj - pt_mean_self_proj;
    tangente_self_proj /= tangente_self_proj.norm();

    if(acos(tangente.dot(tangente_sister_proj))>M_PI/2)
        tangente_sister_proj *=-1;
    if(acos(sister.tangente.dot(tangente_self_proj))>M_PI/2)
        tangente_self_proj *=-1;

    //compute tangentes as mean of current and projected
    tangente = (tangente + tangente_sister_proj)/2;
    tangente /= tangente.norm();
    sister.tangente = (sister.tangente + tangente_self_proj)/2;
    sister.tangente /= sister.tangente.norm();

    //compute new mean points as mean of current and projected
    pt_mean = (pt_mean + pt_mean_sister_proj)/2;
    sister.pt_mean = (sister.pt_mean + pt_mean_self_proj) / 2;

    //deduce normal of lines
    normal = pt_mean - pt_mean.dot(tangente)*tangente;
    normal /= normal.norm();
    sister.normal = sister.pt_mean - sister.pt_mean.dot(sister.tangente)*sister.tangente;
    sister.normal /= sister.normal.norm();

    //deduce distance of lines
    distance = abs(pt_mean.dot(normal));
    sister.distance = abs(sister.pt_mean.dot(sister.normal));

    std::cout<<"normal of corrected self : "<<normal.transpose()<<std::endl<<std::endl;
}
