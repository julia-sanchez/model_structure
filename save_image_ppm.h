#ifndef SAVE_IMAGE_PPM
#define SAVE_IMAGE_PPM

#include <iostream>
#include <string>
#include <fstream>

#include </usr/local/include/eigen3/Eigen/Dense>
#include </usr/local/include/eigen3/Eigen/Core>
#include </usr/local/include/eigen3/Eigen/Eigen>


void save_image_ppm(std::string file_name, std::string complement, std::vector<Eigen::MatrixXi, Eigen::aligned_allocator<Eigen::MatrixXi>> image, int max_col);
#include "save_image_ppm.inl"

#endif // SAVE_IMAGE_PPM
