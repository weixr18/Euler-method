/*
    euler.cpp: Euler methods.

    Copyright (C) 2020 Xinran Wei <weixr0605@sina.com>
    Full LICENCE: ./LICENCE
*/

#include <Eigen/Core>
#include <Eigen/Dense>
#include "euler.h"
#include "mpmat_utils.h"

/********* Euler method framework *********/

MpMat EulerMethod::run(const MpMat& mpm_Y_0) {

    int num = NUM_STEPS_;
    MpMat mpm_Y_n = mpm_Y_0.copy();
    for (int i = 0; i < num; i++) {
        mpm_Y_n = step(mpm_Y_n);
        if (i % 10000 == 0) {
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\bStep: %8d", i);
        }
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

/********* Forward Euler method *********/

ForwardEuler::ForwardEuler(
    std::function<MpMat(const MpMat&)> F,
    const MpMat& mpm_M, const MpMat& mpm_L,
    double x_0, double x_n, double boundary) :
    EulerMethod(F)
{
    init(mpm_M, mpm_L, x_0, x_n, boundary);
}

MpMat ForwardEuler::step(const MpMat& mpm_Y_n) {
    return mpm_Y_n + F_(mpm_Y_n) * mp_h_;
}

void ForwardEuler::init(const MpMat& mpm_M, const MpMat& mpm_L, double x_0, double x_n, double boundary) {

    // get t
    double t = x_n - x_0;

    // prepare
    uint32_t Y_size = mpm_M.row_num();
    eigen_mat I;
    I.setIdentity(Y_size, Y_size);

    // calc NUM_STEPS_
    eigen_mat M = mpm_M.to_matrix();
    eigen_mat tM = M * t;
    //eigen_mat e_tM = tM.exp();
    eigen_mat e_tM = I + tM * (I + tM * 3 / 8 * (I + tM / 6 * (I + tM / 4)));
    eigen_mat M_ = M.inverse();
    eigen_mat tmp = (e_tM - I) * M_;
    eigen_vec ns = tmp * (t / boundary) * mpm_L.to_matrix();
    NUM_STEPS_ = (int)ns.maxCoeff() + 1;
    printf("Total step numbers: %d\n", NUM_STEPS_);
    printf("Step size: %f\n", (double)(t / NUM_STEPS_));

    // calc precision
    eigen_vec ones = eigen_vec::Ones();
    eigen_vec tmp_2 = tmp * ones;
    eigen_vec tmp_3 = tmp_2 * NUM_STEPS_ / (double)t / (double)boundary;
    FLOAT_BITS_ = (int)log2(tmp_2.maxCoeff());
    printf("Least bits of float: %d\n", FLOAT_BITS_);

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
    double x_0, double x_n, double boundary) :
    EulerMethod(F)
{
    init(mpm_M, mpm_L, x_0, x_n, boundary);
}

void TransformEuler::init(const MpMat& mpm_M, const MpMat& mpm_L, double x_0, double x_n, double boundary) {

}

MpMat TransformEuler::step(const MpMat& mpm_Y_n) {
    mp_num_t mp_half_h;
    mp_num_init(mp_half_h);
    mp_num_div_ui(mp_half_h, mp_h_, 2);

    MpMat mpm_Y_n_bar = mpm_Y_n + F_(mpm_Y_n) * mp_half_h;
    mp_num_clear(mp_half_h);
    return mpm_Y_n + F_(mpm_Y_n_bar) * mp_h_;
}

TransformEuler::~TransformEuler() {
}

