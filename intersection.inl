void intersection::setDeltas(float dp, float dt)
{
    delta_phi = dp;
    delta_theta = dt;
    Nrow = (int)((2*M_PI+2*eps)/delta_phi);
    Ncol = (int)((M_PI+2*eps)/delta_theta);
}

void intersection::SeparatePoints(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& points_of_plane_ref, std::vector<std::pair<int,int>>& pixels_of_plane_ref, std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& points_of_plane_neigh, std::vector<std::pair<int,int>>& pixels_of_plane_neigh, std::set<int>& repeated_ref, std::set<int>& repeated_neigh)
{
    std::cout<<"Separate points"<<std::endl<<std::endl;

    int n = 0;
    points_of_plane_ref.resize(indices_self.size());
    pixels_of_plane_ref.resize(indices_self.size());
    for (auto it_indices_self = indices_self.begin(); it_indices_self != indices_self.end(); ++it_indices_self)
    {
        points_of_plane_ref[n] = other_points[*it_indices_self];
        pixels_of_plane_ref[n] = other_pixels[*it_indices_self];
        if(repeated.find(*it_indices_self) != repeated.end())
            repeated_ref.insert(n);
        ++n;
    }

    n = 0;
    points_of_plane_neigh.resize(indices_sister.size());
    pixels_of_plane_neigh.resize(indices_sister.size());
    for (auto it_indices_sister = indices_sister.begin(); it_indices_sister != indices_sister.end(); ++it_indices_sister)
    {
        points_of_plane_neigh[n] = other_points[*it_indices_sister];
        pixels_of_plane_neigh[n] = other_pixels[*it_indices_sister];
        if(repeated.find(*it_indices_sister) != repeated.end())
            repeated_neigh.insert(n);
        ++n;
    }
}

void intersection::computeTheoriticalLineFeaturesObstruction(std::multimap<std::pair<int,int>,std::pair<int,int>> neighborPix2currentPix)
{
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> points_of_plane_ref;
    std::vector<std::pair<int,int>> pixels_of_plane_ref;
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> points_of_plane_neigh;
    std::vector<std::pair<int,int>> pixels_of_plane_neigh;
    std::set<int> repeated_neigh;
    std::set<int> repeated_ref;

    std::set<int> repeated_temp;
    std::set<int> not_repeated_temp;

    SeparatePoints(points_of_plane_ref, pixels_of_plane_ref, points_of_plane_neigh, pixels_of_plane_neigh, repeated_neigh, repeated_ref);

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

    std::set<int> repeated_ref_temp;
    //-----------------------------------------------------------------------------------------------------------------------
    //compute RANSAC and remove points if RANSAC line not good enough

    std::set<int> indices_line_neigh;
    searchLine( points_of_plane_neigh, pixels_of_plane_neigh, plane_neigh, indices_line_neigh, repeated_neigh );

    //--------------------------

    if(indices_line_neigh.size() < min_number_points_on_line)
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

    int n = 0;
    indices_sister.clear();

    if(repeated_neigh.size() != 0)
    {
        for(int k = 0; k < pixels_of_plane_neigh.size(); ++k)
        {
            auto it_repeated_neigh = repeated_neigh.begin();
            if (k!=*it_repeated_neigh && k!=*it_indices_line_neigh)
            {
                not_repeated_temp.insert(other_points.size());
                other_points.push_back(points_of_plane_neigh[k]);
                other_pixels.push_back(pixels_of_plane_neigh[k]);
                indices_sister.insert(n);
                ++n;
            }

            if(k==*it_indices_line_neigh)
            {
                points.push_back(points_of_plane_neigh[k]);
                pixels.push_back(pixels_of_plane_neigh[k]);
                ++it_indices_line_neigh;
                if(it_indices_line_neigh == indices_line_neigh.end())
                    --it_indices_line_neigh;
            }

            if(k==*it_repeated_neigh)
            {
                repeated_temp.insert(other_points.size());
                other_points.push_back(points_of_plane_neigh[k]);
                other_pixels.push_back(pixels_of_plane_neigh[k]);
                ++it_repeated_neigh;
                if(it_repeated_neigh == repeated_neigh.end())
                    --it_repeated_neigh;
                indices_sister.insert(n);
                ++n;
            }
        }
    }
    else
    {
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
                not_repeated_temp.insert(other_points.size());
                other_points.push_back(points_of_plane_neigh[k]);
                other_pixels.push_back(pixels_of_plane_neigh[k]);
                indices_sister.insert(n);
                ++n;
            }
        }
    }

    normal_sister = normal;
    tangente_sister = tangente;
    pt_mean_sister = pt_mean;
    distance_sister = distance;

    std::cout<<"neighbor plane of obstruction : points number of line : "<<points.size()<<std::endl;
    std::cout<<"neighbor plane of obstruction : points number remaining : "<<other_points.size()<<std::endl;

    //-----------------------------------------------------------------------------------------------------------------------
    //keep self points/pixels close to first line of neigh found

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
        }
    }

    if(indices_pixels_of_plane_ref.size() < min_number_points_on_line)
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
    repeated_ref_temp.clear();

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
            if(repeated_ref.find(k) != repeated_ref.end())
                repeated_ref_temp.insert(n);
            ++n;
        }
    }

    repeated_ref = repeated_ref_temp;

    pixels_of_plane_ref = pixels_of_plane_ref_temp;
    points_of_plane_ref = points_of_plane_ref_temp;

    std::cout<<"number of points of ref close to neigh : "<<points_of_plane_ref.size()<<std::endl<<std::endl;

    //-----------------------------------------------------------------------------------------------------------------------
    //find line in self
    std::set<int> indices_line_ref;
    searchLine( points_of_plane_ref, pixels_of_plane_ref, plane_ref, indices_line_ref, repeated_ref );

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
    if(repeated_ref.size() != 0)
    {
        auto it_repeated_ref = repeated_ref.begin();
        for(int k = 0; k < pixels_of_plane_ref.size(); ++k)
        {
            if(k!=*it_repeated_ref && k!=*it_indices_line_ref)
            {
                not_repeated_temp.insert(other_points.size());
                other_points.push_back(points_of_plane_ref[k]);
                other_pixels.push_back(pixels_of_plane_ref[k]);
                indices_self.insert(n);
                ++n;
            }

            if(k==*it_indices_line_ref)
            {
                points.push_back(points_of_plane_ref[k]);
                pixels.push_back(pixels_of_plane_ref[k]);
                ++it_indices_line_ref;
                if(it_indices_line_ref == indices_line_ref.end())
                    --it_indices_line_ref;

            }

            if (k==*it_repeated_ref)
            {
                repeated_temp.insert(other_points.size());
                other_points.push_back(points_of_plane_ref[k]);
                other_pixels.push_back(pixels_of_plane_ref[k]);
                indices_self.insert(n);
                ++n;
                ++it_repeated_ref;
                if(it_repeated_ref == repeated_ref.end())
                    --it_repeated_ref;
            }
        }
    }
    else
    {
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
                not_repeated_temp.insert(other_points.size());
                other_points.push_back(points_of_plane_ref[k]);
                other_pixels.push_back(pixels_of_plane_ref[k]);
                indices_self.insert(n);
                ++n;
            }
        }
    }

    repeated = repeated_temp;
    not_repeated = not_repeated_temp;

    std::cout<<"number of points : "<<points.size()<<std::endl;
    std::cout<<"number of other_points : "<<other_points.size()<<std::endl;
    std::cout<<"obstruction tangente : "<<tangente.transpose()<<std::endl;
    std::cout<<"obstruction normal : "<<normal.transpose()<<std::endl;
    std::cout<<"obstruction distance : "<<distance<<std::endl<<std::endl;

    std::cout<<"obstruction sister tangente : "<<tangente_sister.transpose()<<std::endl;
    std::cout<<"obstruction sister normal : "<<normal_sister.transpose()<<std::endl;
    std::cout<<"obstruction sister distance : "<<distance_sister<<std::endl<<std::endl;
}


void intersection::searchLine( std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& points_of_plane, std::vector<std::pair<int,int>>& pixels_of_plane, plane* p, std::set<int>& indices_line_plane, std::set<int>& repeated_plane) //repeated_plane enter full and is completed
{
    ProjectBoundaryOnPlane(points_of_plane, p);
    std::set<int> repeated_current_plane_temp;
    RANSAC(points_of_plane, pixels_of_plane, ransac_iterations, indices_line_plane, repeated_current_plane_temp);

    while((indices_line_plane.size() < min_number_points_on_line || !has_points_after_ls) && points_of_plane.size() >= 2)
    {
        if(!has_points_after_ls)
            std::cout<<"has no point after ls : try again"<<std::endl<<std::endl;

        if(indices_line_plane.size() < min_number_points_on_line)
            std::cout<<"not enough points for line : try again"<<std::endl<<std::endl;

        for(auto it_indices_line_plane = indices_line_plane.begin(); it_indices_line_plane != indices_line_plane.end(); ++it_indices_line_plane)
        {
            if(repeated_plane.find(*it_indices_line_plane) != repeated_plane.end())
            {
                std::set<int> repeated_plane_temp;
                for(auto it_repeated_plane = repeated_plane.begin(); it_repeated_plane != repeated_plane.end(); ++it_repeated_plane )
                {
                    if(*it_repeated_plane > *it_indices_line_plane)
                        repeated_plane_temp.insert(*it_repeated_plane-1);
                    else if (*it_repeated_plane < *it_indices_line_plane)
                        repeated_plane_temp.insert(*it_repeated_plane);
                }
                repeated_plane = repeated_plane_temp;
            }
            else
            {
                remaining_points.push_back(points_of_plane[*it_indices_line_plane]);
                remaining_pixels.push_back(pixels_of_plane[*it_indices_line_plane]);
            }
            points_of_plane.erase(points_of_plane.begin() + *it_indices_line_plane);
            pixels_of_plane.erase(pixels_of_plane.begin() + *it_indices_line_plane);
        }
        std::cout<<"number of points remaining for line search : "<<points_of_plane.size()<<std::endl<<std::endl;
        if(points_of_plane.size() >= 2)
        {
            ProjectBoundaryOnPlane(points_of_plane, p);
            indices_line_plane.clear();
            repeated_current_plane_temp.clear();

            RANSAC(points_of_plane, pixels_of_plane, ransac_iterations, indices_line_plane, repeated_current_plane_temp);
        }
        else
        {
            indices_line_plane.clear();
            repeated_current_plane_temp.clear();
            return;
        }
    }

    //Repeated points on other_points/other_pixels regarding to the new line found
    if(repeated_current_plane_temp.size()>0)
    {
        for(auto it_repeated_current_plane_temp = repeated_current_plane_temp.begin(); it_repeated_current_plane_temp != repeated_current_plane_temp.end(); ++it_repeated_current_plane_temp)
            repeated_plane.insert(*it_repeated_current_plane_temp);
    }

    //remove points from repeated if all repeated used in one line
    int counter = 0;
    for(auto it_indices_line_plane = indices_line_plane.begin(); it_indices_line_plane != indices_line_plane.end(); ++it_indices_line_plane)
    {
        if(repeated_plane.find(*it_indices_line_plane) != repeated_plane.end())
            ++counter;
    }

    if(counter == indices_line_plane.size())
    {
        std::cout<<"all points of line are repeated"<<std::endl;
        std::cout<<"I remove them from repeated to not compute it again"<<std::endl;
        for(auto it_indices_line_plane = indices_line_plane.begin(); it_indices_line_plane != indices_line_plane.end(); ++it_indices_line_plane)
        {
            auto it_repeated_found = repeated_plane.find(*it_indices_line_plane);
            if(it_repeated_found != repeated_plane.end())
                repeated_plane.erase(it_repeated_found);
        }
    }
    std::cout<<"number of repeated points : "<<repeated_plane.size()<<std::endl;
}

intersection intersection::export_sister()
{
    intersection inter(plane_neigh, plane_ref, delta_phi, delta_theta);
    inter.isObstruction = isObstruction;
    inter.isObject = isObject;
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




void intersection::computeTheoriticalLineFeaturesObject(std::multimap<std::pair<int,int>,std::pair<int,int>> neighborPix2currentPix)
{
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> points_of_plane_ref;
    std::vector<std::pair<int,int>> pixels_of_plane_ref;
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> points_of_plane_neigh;
    std::vector<std::pair<int,int>> pixels_of_plane_neigh;
    std::set<int> repeated_neigh;
    std::set<int> repeated_ref;

    std::set<int> not_repeated_temp;
    std::set<int> repeated_temp;

    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> other_points_temp = other_points;
    std::vector<std::pair<int,int>> other_pixels_temp = other_pixels;

    SeparatePoints(points_of_plane_ref, pixels_of_plane_ref, points_of_plane_neigh, pixels_of_plane_neigh, repeated_neigh, repeated_ref);

    std::cout<<"number of points of ref : "<<points_of_plane_ref.size()<<std::endl;
    std::cout<<"number of points of neigh : "<<points_of_plane_neigh.size()<<std::endl;

    if(points_of_plane_ref.size()<min_number_points_on_line || points_of_plane_neigh.size()<min_number_points_on_line)
    {
        points_of_plane_ref = other_points_temp;
        pixels_of_plane_ref = other_pixels_temp;
        points.clear();
        pixels.clear();
        other_points.clear();
        other_pixels.clear();
        indices_sister.clear();
        indices_self.clear();
        points_of_plane_neigh.clear();
        pixels_of_plane_neigh.clear();
        has_sister = false;
        std::cout<<"Neighbor line does not exist : not enough points, Try self line with all points : number of points_of_plane_ref : "<<points_of_plane_ref.size()<<std::endl<<std::endl;
    }

    int n = 0;

    //-----------------------------------------------------------------------------------------------------------------------

    if(points_of_plane_neigh.size() >= min_number_points_on_line)
    {
        std::cout<<"Computing neighbor line"<<std::endl;
        has_sister = true;
        std::set<int> indices_line_neigh;
        std::cout<<"number of init points : "<<pixels_of_plane_neigh.size()<<std::endl;
        searchLine( points_of_plane_neigh, pixels_of_plane_neigh, plane_neigh, indices_line_neigh, repeated_neigh); //repeated_temp is her just to adjust points_of_plane_neigh if not enough point in ransac line. It can be outputed empty if the first ransac line is accepted and no point is erased

        if(indices_line_neigh.size() < min_number_points_on_line)
        {
            points_of_plane_ref = other_points_temp;
            pixels_of_plane_ref = other_pixels_temp;
            points.clear();
            pixels.clear();
            other_points.clear();
            other_pixels.clear();
            indices_sister.clear();
            indices_self.clear();
            points_of_plane_neigh.clear();
            pixels_of_plane_neigh.clear();
            has_sister = false;
            std::cout<<"Neighbor line does not exist : not enough points, Try self line with all points : number of points_of_plane_ref : "<<points_of_plane_ref.size()<<std::endl<<std::endl;
        }
        else
        {
            points.clear();
            pixels.clear();
            other_points.clear();
            other_pixels.clear();
            indices_sister.clear();
            indices_self.clear();

            auto it_indices_line_neigh = indices_line_neigh.begin();
            n = 0;
            indices_sister.clear();

            if(repeated_neigh.size() != 0)
            {
                auto it_repeated_neigh = repeated_neigh.begin();
                for(int k = 0; k < pixels_of_plane_neigh.size(); ++k)
                {
                    if (k!=*it_repeated_neigh && k!=*it_indices_line_neigh)
                    {
                        not_repeated_temp.insert(other_points.size());
                        other_points.push_back(points_of_plane_neigh[k]);
                        other_pixels.push_back(pixels_of_plane_neigh[k]);
                        indices_sister.insert(n);
                        ++n;
                    }

                    if(k==*it_indices_line_neigh)
                    {
                        points.push_back(points_of_plane_neigh[k]);
                        pixels.push_back(pixels_of_plane_neigh[k]);
                        ++it_indices_line_neigh;
                        if(it_indices_line_neigh == indices_line_neigh.end())
                            --it_indices_line_neigh;
                    }

                    if (k==*it_repeated_neigh)
                    {
                        repeated_temp.insert(other_points.size());
                        other_points.push_back(points_of_plane_neigh[k]);
                        other_pixels.push_back(pixels_of_plane_neigh[k]);
                        ++it_repeated_neigh;
                        if(it_repeated_neigh == repeated_neigh.end())
                            --it_repeated_neigh;
                        indices_sister.insert(n);
                        ++n;
                    }
                }
            }
            else
            {
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
                        not_repeated_temp.insert(other_points.size());
                        other_points.push_back(points_of_plane_neigh[k]);
                        other_pixels.push_back(pixels_of_plane_neigh[k]);
                        indices_sister.insert(n);
                        ++n;
                    }
                }
            }

            normal_sister = normal;
            tangente_sister = tangente;
            pt_mean_sister = pt_mean;
            distance_sister = distance;

            std::cout<<"object sister tangente : "<<tangente_sister.transpose()<<std::endl;
            std::cout<<"object sister normal : "<<normal_sister.transpose()<<std::endl;
            std::cout<<"object sister distance : "<<distance_sister<<std::endl<<std::endl;

            std::cout<<"neighbor plane of object : points number of line : "<<points.size()<<std::endl;
            std::cout<<"neighbor plane of object : points number remaining : "<<other_points.size()<<std::endl;

            //-----------------------------------------------------------------------------------------------------------------------
            //searching closest points in ref of neighbors points found

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
                }
            }

            if(indices_pixels_of_plane_ref.size() < min_number_points_on_line)
            {
                points_of_plane_ref = other_points_temp;
                pixels_of_plane_ref = other_pixels_temp;
                points_of_plane_neigh.clear();
                pixels_of_plane_neigh.clear();
                points.clear();
                pixels.clear();
                other_points.clear();
                other_pixels.clear();
                indices_sister.clear();
                indices_self.clear();
                repeated_temp.clear();
                not_repeated_temp.clear();
                has_sister = false;
            }
            else
            {
                indices_self.clear();
                std::set<int> repeated_ref_temp;

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
                        if(repeated_ref.find(k) != repeated_ref.end())
                                repeated_ref_temp.insert(n);
                        ++n;
                    }
                }

                repeated_ref = repeated_ref_temp;

                pixels_of_plane_ref = pixels_of_plane_ref_temp;
                points_of_plane_ref = points_of_plane_ref_temp;

                std::cout<<"number of points of ref close to neigh : "<<points_of_plane_ref.size()<<std::endl<<std::endl;
            }
        }
    }
    //-----------------------------------------------------------------------------------------------------------------------

    std::cout<<"Computing self line"<<std::endl;
    std::set<int> indices_line_ref;
    std::cout<<"number of init points : "<<pixels_of_plane_ref.size()<<std::endl;
    searchLine( points_of_plane_ref, pixels_of_plane_ref, plane_ref, indices_line_ref, repeated_ref );

    if(indices_line_ref.size() < min_number_points_on_line)
    {
        points.clear();
        pixels.clear();
        other_points.clear();
        other_pixels.clear();
        indices_sister.clear();
        indices_self.clear();
        std::cout<<"exit computing obstruction : not enough points on line"<<std::endl;
        return;
    }

    auto it_indices_line_ref = indices_line_ref.begin();

    if(repeated_ref.size() != 0)
    {
        auto it_repeated_ref = repeated_ref.begin();
        for(int k = 0; k < pixels_of_plane_ref.size(); ++k)
        {
            if(k!=*it_repeated_ref && k!=*it_indices_line_ref)
            {
                not_repeated_temp.insert(other_points.size());
                other_points.push_back(points_of_plane_ref[k]);
                other_pixels.push_back(pixels_of_plane_ref[k]);
                indices_self.insert(n);
                ++n;
            }

            if(k==*it_indices_line_ref)
            {
                points.push_back(points_of_plane_ref[k]);
                pixels.push_back(pixels_of_plane_ref[k]);
                ++it_indices_line_ref;
                if(it_indices_line_ref == indices_line_ref.end())
                    --it_indices_line_ref;

            }

            if (k==*it_repeated_ref)
            {
                repeated_temp.insert(other_points.size());
                other_points.push_back(points_of_plane_ref[k]);
                other_pixels.push_back(pixels_of_plane_ref[k]);
                indices_self.insert(n);
                ++n;
                ++it_repeated_ref;
                if(it_repeated_ref == repeated_ref.end())
                    --it_repeated_ref;
            }
        }
    }
    else
    {
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
                not_repeated_temp.insert(other_points.size());
                other_points.push_back(points_of_plane_ref[k]);
                other_pixels.push_back(pixels_of_plane_ref[k]);
                indices_self.insert(n);
                ++n;
            }
        }
    }

    repeated = repeated_temp;
    not_repeated = not_repeated_temp;

    std::cout<<"number of points : "<<points.size()<<std::endl;
    std::cout<<"number of other_points : "<<other_points.size()<<std::endl;
    std::cout<<"obstruction tangente : "<<tangente.transpose()<<std::endl;
    std::cout<<"obstruction normal : "<<normal.transpose()<<std::endl;
    std::cout<<"obstruction distance : "<<distance<<std::endl<<std::endl;
}




void intersection::computeTheoriticalLineFeaturesOpening()
{
    std::set<int> indices_line;
    std::set<int> repeated_temp;
    std::set<int> not_repeated_temp;

    std::cout<<"number of other_points before searchLine : "<<other_points.size()<<std::endl<<std::endl;
    searchLine( other_points, other_pixels, plane_ref, indices_line, repeated );
    std::cout<<"number of other_points after searchLine : "<<other_points.size()<<std::endl<<std::endl;

    if(indices_line.size()==0)
    {
        points.clear();
        pixels.clear();
        other_points.clear();
        other_pixels.clear();
        return;
    }

    auto it_indices_line = indices_line.begin();

    points.clear();
    pixels.clear();
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> other_points_temp;
    std::vector<std::pair<int,int>> other_pixels_temp;

    if(repeated.size() != 0)
    {
        auto it_repeated = repeated.begin();
        for(int k = 0; k < other_pixels.size(); ++k)
        {
            if(k!=*it_repeated && k != *it_indices_line)
            {
                not_repeated_temp.insert(other_points_temp.size());
                other_points_temp.push_back(other_points[k]);
                other_pixels_temp.push_back(other_pixels[k]);
            }

            if(k == *it_indices_line)
            {
                points.push_back(other_points[k]);
                pixels.push_back(other_pixels[k]);
                ++it_indices_line;
                if(it_indices_line == indices_line.end())
                    --it_indices_line;
            }

            if(k==*it_repeated)
            {
                repeated_temp.insert(other_points_temp.size());
                other_points_temp.push_back(other_points[k]);
                other_pixels_temp.push_back(other_pixels[k]);
                ++it_repeated;
                if(it_repeated == repeated.end())
                    --it_repeated;
            }
        }
    }
    else
    {
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
                not_repeated_temp.insert(other_points_temp.size());
                other_points_temp.push_back(other_points[k]);
                other_pixels_temp.push_back(other_pixels[k]);
            }
        }
    }

    repeated = repeated_temp;
    not_repeated = not_repeated_temp;

    other_points = other_points_temp;
    other_pixels = other_pixels_temp;

    std::cout<<"number of points : "<<points.size()<<std::endl;
    std::cout<<"number of other_points : "<<other_points.size()<<std::endl;
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
    Eigen::Vector3d pt = Eigen::Vector3d::Zero();

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
    pt_mean = pt;

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
    isObject = false;
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

        //actualize pt_mean :
        pt_mean = Eigen::Vector3d::Zero();
        for(int idx_points = 0; idx_points < points.size(); ++idx_points)
            pt_mean += points[idx_points];
        pt_mean /= points.size();
        pt_mean = distance * normal + pt_mean.dot(tangente)*tangente;

        std::cout<<"points remaining of neghbor plane : "<<indices_sister.size()<<std::endl;
        std::cout<<"points remaining of current ref plane : "<<indices_self.size()<<std::endl;

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

std::vector<intersection> intersection::export_sisters()
{
    // divide line into various if  large space in the middle (example : door opening in wall)

    std::map<double, int> proj_temp; // map to keep the link point projection on line / index of point and to order depending on projection

    for(int j = 0; j< points.size(); ++j)
        proj_temp.insert(std::make_pair(points[j].dot(tangente), j));

    std::vector<std::map<double, int>::iterator> iterators_lim;

    auto last_it = proj_temp.end();
    --last_it;

    if(proj_temp.size()>3)
    {
        for(auto it_proj = proj_temp.begin(); it_proj != last_it; ++it_proj)
        {
            auto it_proj_after = it_proj;
            ++it_proj_after;
            //dist_after is the distance acceptable in pixels->
            float dist_after_pix = sqrt( pow( pixels[it_proj_after->second].first - pixels[it_proj->second].first , 2 ) + pow( pixels[it_proj_after->second].second - pixels[it_proj->second].second , 2 ) );
            float dist_after_spat = (points[it_proj_after->second] - points[it_proj->second]).norm();
            if(dist_after_spat > max_dist_between_points_in_line && dist_after_pix > max_dist_between_pixels_in_line)
                iterators_lim.push_back(it_proj);
        }
    }

    std::vector<intersection> vec_sisters;
    if(iterators_lim.size()>0)
    {
        iterators_lim.push_back(proj_temp.end());
        intersection inter (plane_ref, plane_neigh, delta_phi, delta_theta);
        inter.normal = normal;
        inter.tangente = tangente;
        inter.pt_mean = pt_mean;
        inter.distance = distance;
        inter.isConnection = true;
        inter.isObstruction = false;
        inter.isObject = false;
        inter.isOpening = false;
        inter.isLine = true;

        int n = 1;
        auto it_start = iterators_lim[0];
        ++it_start;
        for(auto it_proj_temp = it_start; it_proj_temp != proj_temp.end(); ++it_proj_temp)
        {
            inter.points.push_back(points[it_proj_temp->second]);
            inter.pixels.push_back(pixels[it_proj_temp->second]);
            if(it_proj_temp == iterators_lim[n])
            {
                vec_sisters.push_back(inter);
                ++n;
                inter.points.clear();
                inter.pixels.clear();
            }
        }

        vec_sisters.push_back(inter);

        //fill this intersection with first line piece found
        std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> points_temp;
        std::vector<std::pair<int,int>> pixels_temp;

        auto it_end = iterators_lim[0];
        ++it_end;
        for(auto it_proj_temp =  proj_temp.begin(); it_proj_temp != it_end; ++it_proj_temp)
        {
            points_temp.push_back(points[it_proj_temp->second]);
            pixels_temp.push_back(pixels[it_proj_temp->second]);
        }

        points = points_temp;
        pixels = pixels_temp;
    }

    return vec_sisters;
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
}

bool intersection::computeLim()
{
    std::set<float> dot;
    //compute real current limits
    //take the farthest points as current limits--------------------------------------------------------------------

    for (auto it_points=points.begin(); it_points != points.end(); ++it_points)
    {
        if(abs(abs(it_points->dot(plane_ref->normal))-plane_ref->distance)<max_plane_distance) // for obstruction cases : we only keep points onto the current plane
            dot.insert(it_points->dot(tangente)); // project points of intersection onto tangent line
    }

    if(dot.size() < 2)
    {
        for(auto it_pix = pixels.begin(); it_pix != pixels.end(); ++it_pix)
            std::cout<<it_pix->first<<" "<<it_pix->second<<std::endl;
        std::cout<<"no points of intersection on plane"<<std::endl<<std::endl;
        return false;
    }

    auto dot_end = dot.end();
    --dot_end;

    start = *dot.begin();
    end =   *dot_end;

    start_pt = distance * normal + start * tangente;
    end_pt =   distance * normal + end * tangente;

    new_start_pt = start_pt;
    new_end_pt = end_pt;

    length = (start_pt - end_pt).norm();

    std::cout<<"visible start from data points = "<<start_pt.transpose()<<std::endl;
    std::cout<<"visible end from data points = "<<end_pt.transpose()<<std::endl<<std::endl;
    return true;
}


void intersection::ProjectLineFeaturesOnPlane()
{
    tangente2D = {(rot.linear() * tangente)(0), (rot.linear() * tangente)(1)};
    Eigen::Vector2d pt2D = {(rot.linear() * pt_mean)(0), (rot.linear() * pt_mean)(1)};
    normal2D = pt2D - pt2D.dot(tangente2D) * tangente2D;
    normal2D /= normal2D.norm();
    distance2D = pt2D.dot(normal2D);
}


void intersection::ProjectBoundaryOnPlane(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> pts_to_projet, plane* plane_to_project)
{
    rot = Eigen::Affine3d::Identity();
    rot_inv = Eigen::Affine3d::Identity();

    float angle = acos(-plane_to_project->normal(2));   //angle between normal and z axis [0, pi]
    Eigen::Vector3d axis;
    if(angle>eps)
    {
        axis = (-plane_to_project->normal).cross(Eigen::Vector3d(0,0,1)); // rotation axis to align normal onto z axis
        axis /= axis.norm();
        rot.rotate( Eigen::AngleAxisd(angle, axis) );
        rot_inv.rotate( Eigen::AngleAxisd(angle, -axis) );
    }

    points2D.resize(pts_to_projet.size());
    for(int i = 0; i<pts_to_projet.size(); ++i)
    {
        Eigen::Vector3d turned = pts_to_projet[i] * (plane_to_project->distance/pts_to_projet[i].dot(-plane_to_project->normal)); // project points on plane
        turned = rot.linear()*turned; //turn it
        points2D[i](0) = turned(0);
        points2D[i](1) = turned(1);
    }

    trans_z = Eigen::Affine3d::Identity();
    trans_z_inv = Eigen::Affine3d::Identity();
    trans_z.translate(Eigen::Vector3d(0.0, 0.0, -plane_to_project->distance));
    trans_z_inv.translate(Eigen::Vector3d(0.0, 0.0, plane_to_project->distance));
}

bool intersection::isNearCorner(std::vector<Eigen::Vector3d> corners_pt)
{
    for(int k = 0; k < corners_pt.size(); ++k)
    {
        if( (pt_mean-corners_pt[k]).norm() < 0.1)
            return true;
    }
    return false;
}

bool intersection::isDoubled(std::vector<intersection> all_edges)
{
    for(int k = 0; k<all_edges.size(); ++k)
    {
        if(all_edges[k].isLine && all_edges[k].plane_ref->index == plane_ref->index || (all_edges[k].isConnection && all_edges[k].plane_neigh->index == plane_ref->index))
        {
            double eps = 0.00001;
            bool same = (all_edges[k].tangente-tangente).norm() < eps && (all_edges[k].normal-normal).norm() < eps && (all_edges[k].pt_mean-pt_mean).norm() < eps;
            if(!same)
            {
                bool parallel = acos(abs(all_edges[k].tangente.dot(tangente))) < 30*M_PI/180;
                bool lines_distance = ((all_edges[k].pt_mean - pt_mean) - (all_edges[k].pt_mean - pt_mean).dot(all_edges[k].tangente)*all_edges[k].tangente).norm()  < 0.05;

                Eigen::Vector3d tangente_temp;

                if(tangente.dot(all_edges[k].tangente)<0)
                    tangente_temp = -tangente;
                else
                    tangente_temp = tangente;

                tangente_temp = (all_edges[k].tangente + tangente_temp)/2;
                tangente_temp/=tangente_temp.norm();

                double start1 = start_pt.dot(tangente_temp);
                double end1 = end_pt.dot(tangente_temp);
                double start2 = all_edges[k].start_pt.dot(tangente_temp);
                double end2 = all_edges[k].end_pt.dot(tangente_temp);

                bool superimposed = !( (start1<start2 && end1<start2) || (start2<start1 && end2<start1) );

                bool is_the_smallest = points.size()<all_edges[k].points.size();

                if(parallel && lines_distance && superimposed && is_the_smallest)
                    return true;
            }
        }
    }
    return false;
}

void intersection::RANSAC(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& points_init, std::vector<std::pair<int,int>>& pixels_init, int tests, std::set<int>& indices_line, std::set<int>& repeated_indices )
{
    int idx1, idx2;
    int max_points = 0;
    std::map<double, int> proj;

    std::srand (std::time(NULL));

    Eigen::MatrixXi processed = Eigen::MatrixXi::Zero(points_init.size(), points_init.size());

    for(int i = 0; i< tests; ++i)
    {
        idx1 = rand() % points_init.size();// [0; points_init.size()[
        idx2 = rand() % points_init.size();

        int n = 0;
        while((idx2 == idx1 || processed(idx1, idx2)) && n < ransac_iterations)
        {
            ++n;
            idx1 = rand() % points_init.size();
            idx2 = rand() % points_init.size();
        }

        if(n == ransac_iterations)
        {
            std::cout<<"can not find a new pair of points, have all been tested : exiting RANSAC"<<std::endl<<std::endl;
            break;
        }

        processed(idx1, idx2) = 1;
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
            if(abs(points2D[j].dot(normal2D_temp)-distance2D_temp) < max_line_distance)
                proj_temp.insert(std::make_pair(points2D[j].dot(tangente2D_temp), j));
        }

        //filter indices on line to remove points far from the others----------------------------------------------------------------------------------

        std::vector<std::map<double, int>::iterator> iterators_lim;
        iterators_lim.push_back(proj_temp.begin());
        if(proj_temp.size()>3)
        {
            auto last_it = proj_temp.end();
            --last_it;

            for(auto it_proj = proj_temp.begin(); it_proj != last_it; ++it_proj)
            {
                auto it_proj_after = it_proj;
                ++it_proj_after;
                float dist_after_pix = sqrt( pow( pixels_init[it_proj_after->second].first - pixels_init[it_proj->second].first , 2 ) + pow( pixels_init[it_proj_after->second].second - pixels_init[it_proj->second].second , 2 ) );
                float dist_after_spat = (points_init[it_proj_after->second] - points_init[it_proj->second]).norm();

                if(dist_after_spat > max_dist_between_points_in_line && dist_after_pix > max_dist_between_pixels_in_line) // pixel distance is added to counter problems of sampling on a line far from sensor
                    iterators_lim.push_back(it_proj_after);
            }
        }

        if(iterators_lim.size()>1)
        {
            int idx =0;
            int num_points = 0;
            for(int k = 1; k < iterators_lim.size(); ++k)
            {
                if(std::distance(iterators_lim[k-1], iterators_lim[k]) > num_points)
                {
                    idx = k;
                    num_points = std::distance(iterators_lim[k-1], iterators_lim[k]);
                }
            }

            std::map<double, int> proj_temp_temp;

            for(auto it_proj_temp = iterators_lim[idx-1]; it_proj_temp != iterators_lim[idx]; ++it_proj_temp)
                proj_temp_temp.insert(*it_proj_temp);

            proj_temp = proj_temp_temp;
        }

        //check if the line is the best--------------------------------------------------------------------------------------------------------------------------------------------------------------------

        if(proj_temp.size()>= max_points)
        {
            auto start_proj_temp = proj_temp.begin();
            auto end_proj_temp = proj_temp.end();
            --end_proj_temp;
            double length_temp = abs(end_proj_temp->first - start_proj_temp->first);
            length = length_temp;
            normal2D = normal2D_temp;
            tangente2D = tangente2D_temp;
            distance2D = distance2D_temp;
            max_points = proj_temp.size();
            proj = proj_temp; 
        }
    }

    //Erase , points of line and its pixels neighbors from other_points + add pixels neighbors to points
//-----------------------------------------------------------------------------------------------------------

    //-----------------
//    Eigen::MatrixXi image_test = Eigen::MatrixXi::Zero(Nrow, Ncol);
//    if(plane_ref->index == 4 && plane_neigh->index == 4)
//    {
//        for (auto it_proj = proj.begin(); it_proj != proj.end(); ++it_proj)
//            image_test(pixels_init[it_proj->second].first, pixels_init[it_proj->second].second) = 1;
//        save_image_pgm("inter","",image_test,1);
//        getchar();
//        getchar();
//    }
    //--------------------
    //build indices_line from proj elements

    std::map<double, int>::iterator start_proj;
    std::map<double, int>::iterator end_proj;

    if(proj.size()>0)
    {
        std::cout<<"max points on line : "<<max_points<<std::endl;
        std::cout<<"tangente2D : "<<tangente2D.transpose()<<std::endl<<std::endl;

        start_proj = proj.begin();
        end_proj = proj.end();
        --end_proj;

        int proj_idx;
        int rad = 2;

        repeated_indices.clear();

        //on remplit indices_line et repeated
        for (auto it_proj = proj.begin(); it_proj != proj.end(); ++it_proj)
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
                            repeated_indices.insert(pixels_init_idx);
                    }
                }
            }
        }
    }

    //if enough points and long enough -> least square
    if(proj.size()>=min_number_points_on_line && length > min_line_length)
    {
        std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> pt2D;
        for(auto it_indices_line = indices_line.begin(); it_indices_line!=indices_line.end(); ++it_indices_line)
            pt2D.push_back(points2D[*it_indices_line]); //not in previous loop because indices would have been repeated when looking for neighbors...

    //        least square for accuracy
    //---------------------------------------------------------------------------------------------
        float theta = atan2(normal2D(1),normal2D(0));
        tangente2D = {sin(theta), -cos(theta)};

        Eigen::MatrixXd JTJ(2,2);
        Eigen::MatrixXd JTr(2, 1);
        Eigen::MatrixXd J(2,1);
        double r;
        int nbr_iterations = 10;
        Eigen::Vector2d pt2D_mean;

        int n;
        for(int k =0; k<nbr_iterations; ++k)
        {
            n = 0;
            for(int i = 0; i< pt2D.size(); ++i)
            {
                if(abs(pt2D[i].dot(tangente2D)-proj.begin()->first)>line_margin * length && abs(pt2D[i].dot(tangente2D)-end_proj->first)>line_margin * length)
                    ++n;
            }

            if(n>2)
            {
                std::cout<<"Starting least square : "<<std::endl<<std::endl;
                JTJ.setZero();
                JTr.setZero();
                n = 0;
                pt2D_mean = Eigen::Vector2d::Zero();

                for(int i = 0; i< pt2D.size(); ++i)
                {
                    if(abs(pt2D[i].dot(tangente2D)-proj.begin()->first)>line_margin * length && abs(pt2D[i].dot(tangente2D)-end_proj->first)>line_margin * length)
                    {
                        J(0,0) = pt2D[i].dot(Eigen::Vector2d(-sin(theta), cos(theta)));
                        J(1,0) = -1;
                        r = pt2D[i].dot(normal2D)-distance2D;
                        JTJ += J * J.transpose();
                        JTr += J * r;
                        pt2D_mean += pt2D[i];
                        ++n;
                    }
                }

                pt2D_mean /= n;
                pt2D_mean = distance2D * normal2D + pt2D_mean.dot(tangente2D)*tangente2D;

                Eigen::MatrixXd result(2, 1);
                result = -JTJ.llt().solve(JTr);
                theta += result(0);

                normal2D = {cos(theta), sin(theta)};
                distance2D += result(1);

                if(distance2D<0)
                {
                    distance2D *= -1;
                    normal2D *= -1;
                    theta += M_PI;
                }

                tangente2D = {sin(theta), -cos(theta)};
                proj.clear();
                for(int j = 0; j< pt2D.size(); ++j)
                    proj.insert(std::make_pair(pt2D[j].dot(tangente2D), j));

                end_proj = proj.end();
                --end_proj;

                has_points_after_ls = true;

                std::cout<< "number of points : "<<n<<std::endl;
                std::cout<< "normal : "<<normal2D.transpose()<<std::endl;
                std::cout<< "distance 2D : "<<distance2D<<std::endl;
                std::cout<<  "theta = "<<theta<<std::endl<<std::endl;
            }
            else
            {
                has_points_after_ls = false;
                std::cout<<"not enough points on line after ls process"<<std::endl; // can happen if just one point on line for example (if there were 4 and 2 are on the boundary) and then, normal is computed randomly
                return;
            }
        }

        //rechoose points
        //---------------------------------------------------------------------------------------------------------

//        proj.clear();
//        repeated_indices.clear();
//        indices_line.clear();
//        for (int j = 0; j < points2D.size(); ++j)
//        {
//            if(abs(points2D[j].dot(normal2D)-distance2D) < max_line_distance)
//                proj.insert(std::make_pair(points2D[j].dot(tangente2D), j));
//        }

//        auto last_it = proj.end();
//        --last_it;

//        std::vector<std::map<double, int>::iterator> iterators_lim;
//        iterators_lim.push_back(proj.begin());
//        for(auto it_proj = proj.begin(); it_proj != last_it; ++it_proj)
//        {
//            auto it_proj_after = it_proj;
//            ++it_proj_after;
//            float dist_after_pix = sqrt( pow( pixels_init[it_proj_after->second].first - pixels_init[it_proj->second].first , 2 ) + pow( pixels_init[it_proj_after->second].second - pixels_init[it_proj->second].second , 2 ) );
//            float dist_after_spat = (points_init[it_proj_after->second] - points_init[it_proj->second]).norm();
//            if(dist_after_spat > max_dist_between_points_in_line && dist_after_pix > max_dist_between_pixels_in_line) // pixel distance is added to counter problems of sampling on a line far from sensor
//                iterators_lim.push_back(it_proj_after);
//        }

//        if(iterators_lim.size()>1)
//        {
//            int idx =0;
//            int num_points = 0;
//            for(int k = 1; k < iterators_lim.size(); ++k)
//            {
//                if(std::distance(iterators_lim[k-1], iterators_lim[k]) > num_points)
//                {
//                    idx = k;
//                    num_points = std::distance(iterators_lim[k-1], iterators_lim[k]);
//                }
//            }

//            std::map<double, int> proj_temp;

//            for(auto it_proj_temp = iterators_lim[idx-1]; it_proj_temp != iterators_lim[idx]; ++it_proj_temp)
//                proj_temp.insert(*it_proj_temp);

//            proj = proj_temp;
//        }

//        start_proj = proj.begin();
//        end_proj = proj.end();
//        --end_proj;

//        for (auto it_proj = proj.begin(); it_proj != proj.end(); ++it_proj)
//        {
//            int X = pixels_init[it_proj->second].first;
//            int Y = pixels_init[it_proj->second].second;
//            int imin = std::max(0, X-rad);
//            int imax = std::min(Nrow-1, X+rad);
//            int jmin = std::max(0, Y-rad);
//            int jmax = std::min(Ncol-1, Y+rad);

//            for(int i = imin; i<=imax; ++i)
//            {
//                for(int j = jmin; j<=jmax; ++j)
//                {
//                    auto it_found_pixels_init = std::find(pixels_init.begin(), pixels_init.end(), std::make_pair(i, j));
//                    if(it_found_pixels_init != pixels_init.end())
//                    {
//                        int pixels_init_idx = std::distance(pixels_init.begin(), it_found_pixels_init);
//                        indices_line.insert(pixels_init_idx);
//                        if(abs(points2D[pixels_init_idx].dot(tangente2D)-start_proj->first)<line_margin || abs(points2D[pixels_init_idx].dot(tangente2D)-end_proj->first)<line_margin)
//                            repeated_indices.insert(pixels_init_idx);
//                    }
//                }
//            }
//        }


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
    else  //if not enough points or not long enough -> remove current points and try again
    {
        std::cout<< "line too thin or not enough points"<<std::endl;
        std::cout<<"number of points found by RANSAC + dilatation : "<<proj.size()<<std::endl;
        std::cout<<"length of line found by RANSAC : "<<length<<std::endl<<std::endl;
        return;
    }

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
