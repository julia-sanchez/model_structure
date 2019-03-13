void save_image_ppm(std::string file_name, std::string complement, std::vector<Eigen::MatrixXi, Eigen::aligned_allocator<Eigen::MatrixXi>> image, int max_col)
{
    std::stringstream sstm;
    std::string file_name1;
    size_t lastindex_point = file_name.find_last_of(".");
    size_t lastindex_slash = file_name.find_last_of("/");
    if (lastindex_slash==std::string::npos)
    {
       lastindex_slash = -1;
    }

    file_name1 = file_name.substr(lastindex_slash+1, lastindex_point-(lastindex_slash+1));
    sstm.str("");
    sstm<<"results/"<<file_name1<<complement<<".ppm";
    std::string file_name_tot = sstm.str();
    std::ofstream file (file_name_tot, std::ios::trunc);
    file<<"P3\n";
    file<<image[0].cols()<<" "<<image[0].rows()<<"\n";
    file<<max_col<<"\n";

    sstm.str("");

    for(int i = 0; i< image[0].rows(); ++i)
    {
        for(int j = 0; j< image[0].cols(); ++j)
        {
            sstm<<image[0](i,j)<<" "<<image[1](i,j)<<" "<<image[2](i,j)<<" ";
        }
        file<<sstm.str();
        sstm.str("");
    }

    file.close();

}
