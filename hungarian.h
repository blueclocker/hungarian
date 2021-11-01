/*
 * @Author: your name
 * @Date: 2021-10-18 09:12:17
 * @LastEditTime: 2021-10-21 13:53:00
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /hungarian-algorithm-cpp/hungarian.h
 */

#include <Eigen/Dense>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <stdlib.h>

class hungarian
{
private:
    int Raws;
    int Cols;
    int maxdim;
    int mindim;
    double cost;
    Eigen::MatrixXd dismatrix;//保留
    Eigen::MatrixXd dismatrixFlow;//动态计算
    Eigen::VectorXi assignment;//结果，匹配失败则为-1
    bool init(const Eigen::MatrixXd &initMatrixin);
    void step12();//对于矩阵每行，减去该行最小的一个数;对于矩阵每列，减去该列最小的一个数
    bool step3(bool* coveredraws, bool* coveredcols);//用最少的横线或竖线覆盖所有0元素
    void step4(bool* rawsflag, bool* colsflag);//对于每个未覆盖的元素都减去这个数，对于每个覆盖了两次的元素加上这个数
    double computecost();//计算损失
    void outputassignment();//输出结果
public:
    Eigen::VectorXi solve(double &costout);
    hungarian(const Eigen::MatrixXd &initMatrix);
    ~hungarian();
};
