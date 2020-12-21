/*
    euler.h: Euler methods.

    Copyright (C) 2020 Xinran Wei <weixr0605@sina.com>
    Full LICENCE: ./LICENCE
*/

#ifndef _EULER_H
#define _EULER_H

#include <functional>
#include "mpmat.h"

class EulerMethod {
public:
    EulerMethod(std::function<MpMat(const MpMat&)> F);
    virtual ~EulerMethod();
public:
    MpMat run(const MpMat& mpm_Y_0);
    virtual void init(const MpMat& mpm_M, const MpMat& L, double x_0, double x_n, double boundary) = 0;
    virtual MpMat step(const MpMat& mpm_Y_n) = 0;
protected:
    std::function<MpMat(const MpMat&)> F_;
    int FLOAT_BITS_;
    int NUM_STEPS_;
    mp_num_t mp_h_;
};


class ForwardEuler : public EulerMethod
{
public:
    ForwardEuler(
        std::function<MpMat(const MpMat&)> F,
        const MpMat& mpm_M, const MpMat& L,
        double x_0, double x_n, double boundary
    );
    ~ForwardEuler();
    void init(const MpMat& mpm_M, const MpMat& L, double x_0, double x_n, double boundary);
    MpMat step(const MpMat& mpm_Y_n);
};

class TransformEuler : public EulerMethod
{
public:
    TransformEuler(
        std::function<MpMat(const MpMat&)> F,
        const MpMat& mpm_M, const MpMat& L,
        double x_0, double x_n, double boundary
    );
    ~TransformEuler();
    void init(const MpMat& mpm_M, const MpMat& L, double x_0, double x_n, double boundary);
    MpMat step(const MpMat& mpm_Y_n);
};

const int EQ_NUM = 4;
typedef Eigen::Matrix<double, EQ_NUM, EQ_NUM> eigen_mat;
typedef Eigen::Matrix<double, EQ_NUM, 1> eigen_vec;

#endif //_EULER_H