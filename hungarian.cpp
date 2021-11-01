/*
 * @Author: your name
 * @Date: 2021-10-18 09:11:23
 * @LastEditTime: 2021-10-21 13:56:07
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /hungarian-algorithm-cpp/hungarian.cpp
 */

#include "hungarian.h"

hungarian::hungarian(const Eigen::MatrixXd &initMatrix)
{
    /*if(init(initMatrix))
    {
        std::cout << "init finished" << std::endl;
    }*/
    init(initMatrix);
}

hungarian::~hungarian()
{
}

bool hungarian::init(const Eigen::MatrixXd &initMatrixin)
{
    this -> Raws = initMatrixin.rows();
    this -> Cols = initMatrixin.cols();
    if(Raws >= Cols)
    {
        this -> maxdim = Raws;
        this -> mindim = Cols;
    }else{
        this -> maxdim = Cols;
        this -> mindim = Raws;
    }
    this -> dismatrix.resize(Raws, Cols);
    //this -> dismatrixFlow.resize(maxdim, maxdim);
    this -> dismatrixFlow = Eigen::MatrixXd::Zero(maxdim, maxdim);
    for(int i = 0; i < Raws; i++)
    {
        for(int j = 0; j < Cols; j++)
        {
            dismatrix(i, j) = initMatrixin(i, j);
            dismatrixFlow(i, j) = initMatrixin(i, j);
        }
    }
    this -> assignment.resize(Raws);
    for(int i = 0; i < Raws; i++)
    {
        assignment(i) = -1;
    }
    return true;
}

Eigen::VectorXi hungarian::solve(double &costout)
{
    step12();
    bool coveredraws[maxdim] = {0};
    bool coveredcols[maxdim] = {0};
    //std::cout << step3(&coveredraws[0], &coveredcols[0]) << std::endl;
    //bool a[4] = {0, 1, 0, 1};
    //bool b[4] = {0, 0, 1, 0};
    while (!step3(&coveredraws[0], &coveredcols[0]))
    {
        step4(&coveredraws[0], &coveredcols[0]);
    }
    outputassignment();
    computecost();
    //std::cout << dismatrixFlow << std::endl;
    //std::cout << "cost= " << computecost() << std::endl;
    costout = cost;
    return assignment;
}

void hungarian::step12()
{
    //对于矩阵每行，减去该行最小的一个数
    for(int i = 0; i < Raws; i++)
    {
        double min = dismatrixFlow.row(i).minCoeff();
        for(int j = 0; j < Cols; j++)
        {
            dismatrixFlow(i, j) -= min;
        }
    }
    //对于矩阵每列，减去该列最小的一个数
    for(int i = 0; i < Cols; i++)
    {
        double min = dismatrixFlow.col(i).minCoeff();
        for(int j = 0; j < Raws; j++)
        {
            dismatrixFlow(j, i) -= min;
        }
    }
    //std::cout << dismatrixFlow << std::endl;
}

bool hungarian::step3(bool* coveredraws, bool* coveredcols)
{
    //用最少的横线或竖线覆盖所有0元素
    int count = 0;
    int maxraw, maxcol;
    Eigen::VectorXi count0raws = Eigen::VectorXi::Zero(maxdim);
    Eigen::VectorXi count0cols = Eigen::VectorXi::Zero(maxdim);
    //每行，每列0个数
    for(int i = 0; i < maxdim; i++)
    {
        for(int j = 0; j < maxdim; j++)
        {
            if(std::fabs(dismatrixFlow(i, j)) < DBL_EPSILON) 
                count0raws(i)++;
            if(std::fabs(dismatrixFlow(j, i)) < DBL_EPSILON) 
                count0cols(i)++;
        }
    }
    //std::cout << count0raws << std::endl;
    //std::cout << count0cols << std::endl;
    //贪心算法
    while (count0raws.maxCoeff() > 0 || count0cols.maxCoeff() > 0)
    {
        if(count0raws.maxCoeff(&maxraw) >= count0cols.maxCoeff(&maxcol))
        {
            coveredraws[maxraw] = true;
            count0raws(maxraw) = -1;
            for(int i = 0; i < maxdim; i++)
            {
                if(std::fabs(dismatrixFlow(maxraw, i)) < DBL_EPSILON)
                    count0cols(i)--;
            }
        }else{
            coveredcols[maxcol] = true;
            count0cols(maxcol) = -1;
            for(int i = 0; i < maxdim; i++)
            {
                if(std::fabs(dismatrixFlow(i, maxcol)) < DBL_EPSILON)
                    count0raws(i)--;
            }
        }
        count++;
        //std::cout << "count= " << count << std::endl;
    }
    return count == maxdim;
}

void hungarian::step4(bool* rawsflag, bool* colsflag)
{
    //找未标记的最小值
    double min = DBL_MAX;
    for(int i = 0; i < Raws; i++)
    {
        if(!rawsflag[i])
        {
            for(int j = 0; j < Cols; j++)
            {
                if(!colsflag[j])
                {
                    if(dismatrixFlow(i, j) < min) min = dismatrixFlow(i, j);
                }
            }
        }
    }
    for(int i = 0; i < maxdim; i++)
    {
        //对未覆盖的行减去min
        if(!rawsflag[i])
        {
            for(int j = 0; j < maxdim; j++)
            {
                dismatrixFlow(i, j) -= min;
            }
        }
        //对覆盖的列加上min
        if(colsflag[i])
        {
            for(int k = 0; k < maxdim; k++)
            {
                dismatrixFlow(k, i) += min;
            }
        }
    }
    //std::cout << dismatrixFlow << std::endl;
}

double hungarian::computecost()
{
    cost = 0;
    for(int i = 0; i < Raws; i++)
    {
        if(assignment(i) != -1)
            cost += dismatrix(i, assignment(i));
    }
    return cost;
}

void hungarian::outputassignment()
{
    int minraw, mincol;
    Eigen::VectorXi count0raws_ = Eigen::VectorXi::Zero(maxdim);
    Eigen::VectorXi count0cols_ = Eigen::VectorXi::Zero(maxdim);
    for(int i = 0; i < maxdim; i++)
    {
        for(int j = 0; j < maxdim; j++)
        {
            if(dismatrixFlow(i, j) < DBL_EPSILON)
                count0raws_(i)++;
            if(dismatrixFlow(j, i) < DBL_EPSILON)
                count0cols_(i)++;
        }
    }
    bool coveredcols[maxdim] = {0};
    int countcols = 0;
    while(count0raws_.minCoeff(&minraw) < Raws)
    {
        if(countcols == Cols) break;
        for(int j = 0; j < Cols; j++)
        {
            if(dismatrixFlow(minraw, j) < DBL_EPSILON && !coveredcols[j])
            {
                assignment[minraw] = j;
                coveredcols[j] = true;
                countcols++;
                break;
            }
        }
        count0raws_(minraw) = maxdim;
    }
    
}