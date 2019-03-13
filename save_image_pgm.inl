void save_image_pgm(std::string file_name, std::string complement, Eigen::MatrixXi image, int max_col)
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
    sstm<<"results/"<<file_name1<<complement<<".pgm";
    std::string file_name_tot = sstm.str();
    std::ofstream file (file_name_tot, std::ios::trunc);
    file<<"P2\n";
    file<<image.cols()<<" "<<image.rows()<<"\n";
    file<<max_col<<"\n";
    file<<image;
    file.close();

}
