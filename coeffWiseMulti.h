#ifndef COEFFWISEMULTI
#define COEFFWISEMULTI

#include <iostream>
#include <string>
#include <fstream>

#include </usr/local/include/eigen3/Eigen/Dense>
#include </usr/local/include/eigen3/Eigen/Core>
#include </usr/local/include/eigen3/Eigen/Eigen>

Eigen::MatrixXi coeffWiseMulti(Eigen::MatrixXi a, Eigen::MatrixXi b)
{
    Eigen::MatrixXi c = a;
    for (int i = 0; i<a.rows(); ++i)
    {
        for (int j = 0; j<a.cols(); ++j)
        {
            c(i,j) = a(i,j)*b(i,j);
        }
    }
    return c;
}

#endif // COEFFWISEMULTI
