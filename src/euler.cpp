/*
    euler.cpp: Euler methods.

    Copyright (C) 2020 Xinran Wei <weixr0605@sina.com>
    Full LICENCE: ./LICENCE
*/

#include <Eigen/Core>
#include "euler.h"
#include "mpmat_utils.h"

/********* Euler method framework *********/

MpMat EulerMethod::run(const MpMat& mpm_Y_0) {

    int num = NUM_STEPS_;
    MpMat mpm_Y_n = mpm_Y_0.copy();
    //uint32_t Y_size = mpm_Y_0.row_num();
    //Eigen::MatrixXd res(Y_size, num);
    for (int i = 0; i < num; i++) {
        mpm_Y_n = step(mpm_Y_n);
        //res.block(0, i, Y_size, 1) = Y_n.to_matrix();
    }
    return mpm_Y_n;
}

EulerMethod::EulerMethod(std::function<MpMat(const MpMat&)> F) {
    F_ = F;
    mp_num_init(mp_h_);
}

EulerMethod::~EulerMethod() {
    mp_num_clear(mp_h_);
}

/*

    // calculation parameters
    mp_num_t x_0, x_n;
    mp_num_init(x_0);
    mp_num_init(x_n);
    mp_num_set_d(x_n, 0);
    mp_num_set_d(x_n, x_n);
*/


/********* Forward Euler method *********/

ForwardEuler::ForwardEuler(
    std::function<MpMat(const MpMat&)> F,
    const MpMat& mpm_M, const MpMat& mpm_L,
    double x_0, double x_n) :
    EulerMethod(F)
{
    init(mpm_M, mpm_L, x_0, x_n);
}

MpMat ForwardEuler::step(const MpMat& mpm_Y_n) {
    return mpm_Y_n + F_(mpm_Y_n) * mp_h_;
}

void ForwardEuler::init(const MpMat& mpm_M, const MpMat& mpm_L, double x_0, double x_n) {

    // get t
    double t = x_n - x_0;

    // prepare
    uint32_t Y_size = mpm_M.row_num();
    Eigen::MatrixXd I(Y_size, Y_size);
    I.setIdentity(Y_size, Y_size);
    Eigen::MatrixXd tmp(Y_size, Y_size);

    // calc n
    Eigen::MatrixXd M = mpm_M.to_matrix();
    Eigen::MatrixXd tM = M * t;
    Eigen::MatrixXd e_tM = tM.exp();
    Eigen::MatrixXd M_ = M.inverse();
    tmp = (e_tM - I) * M_;


    // calc h
    double h = t / NUM_STEPS_;
    mp_num_set_d(mp_h_, h);

    // set precision
    mpf_set_default_prec(FLOAT_BITS_);
    printf("Using float number of %d bits.\n", mpf_get_default_prec());

}

ForwardEuler::~ForwardEuler() {
}

/********* Transform Euler method *********/

TransformEuler::TransformEuler(
    std::function<MpMat(const MpMat&)> F,
    const MpMat& mpm_M, const MpMat& mpm_L,
    double x_0, double x_n) :
    EulerMethod(F)
{
    init(mpm_M, mpm_L, x_0, x_n);
}

void TransformEuler::init(const MpMat& mpm_M, const MpMat& mpm_L, double x_0, double x_n) {

}

MpMat TransformEuler::step(const MpMat& mpm_Y_n) {
    mp_num_t mp_half_h;
    mp_num_init(mp_half_h);
    mp_num_div_ui(mp_half_h, mp_h_, 2);

    MpMat mpm_Y_n_bar = mpm_Y_n + F_(mpm_Y_n) * mp_half_h;
    return mpm_Y_n + F_(mpm_Y_n_bar) * mp_h_;
}

TransformEuler::~TransformEuler() {
}

