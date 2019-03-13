#ifndef SAVE_IMAGE_PGM
#define SAVE_IMAGE_PGM

#include <iostream>
#include <string>
#include <fstream>

#include </usr/local/include/eigen3/Eigen/Dense>
#include </usr/local/include/eigen3/Eigen/Core>
#include </usr/local/include/eigen3/Eigen/Eigen>

void save_image_pgm(std::string file_name, std::string complement, Eigen::MatrixXi image, int max_col);

#include "save_image_pgm.inl"

#endif // SAVE_IMAGE_PGM
