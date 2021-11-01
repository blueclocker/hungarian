/*
 * @Author: your name
 * @Date: 2021-10-18 09:51:24
 * @LastEditTime: 2021-10-21 14:09:35
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /hungarian-algorithm-cpp/test.cpp
 */
#include "hungarian.h"
#include <chrono>

int main()
{
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    Eigen::MatrixXd matrixin(4, 4);
    matrixin << 82, 83, 69, 92 ,77, 37, 49, 92 , 11, 69, 5, 86 , 8, 9, 98, 23;
    //matrixin << 82, 83, 69, 77, 37, 49, 11, 69, 5, 8, 9, 98;
    hungarian mat(matrixin);
    double costprint;
    //Eigen::VectorXi result(mat.solve(costprint));
    std::cout << mat.solve(costprint) << std::endl;
    std::cout << "cost=" << costprint << std::endl;
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout<<"time cost= "<< (time_used.count() * 1000) <<" ms."<<std::endl;
    return 0;
}