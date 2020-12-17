/*
    euler.cpp: Euler methods.
    Copyright (C) 2020 Xinran Wei.

    Full LICENCE: ./LICENCE
*/

#include "euler.h"
#include "mpmat_utils.h"

/********* Euler method framework *********/

MpMat EulerMethod::run(const MpMat& Y_0, const mp_num_t& x_0, const mp_num_t& x_n, const mp_num_t& n) {

    // verify n
    int num = (int)mp_num_get_d(n);
    if (num <= 0) {
        char num_str[1024] = { 0 };
        mp_num_get_str(num_str, NULL, 10, N_DIGITS, n);
        throw RangeException(num_str, "should be positive integer.");
    }

    // calc h
    mp_num_t h, tmp;
    mp_num_init(h);
    mp_num_init(tmp);
    mp_num_sub(tmp, x_n, x_0);
    mp_num_div(h, tmp, n);

    // steps
    MpMat Y_n = Y_0.copy();
    for (int i = 0; i < num; i++) {
        Y_n = step(Y_n, h);
    }

    mp_num_clear(h);
    mp_num_clear(tmp);

    return Y_n;
}


EulerMethod::EulerMethod(std::function<MpMat(const MpMat&)> F) {
    F_ = F;
}

EulerMethod::~EulerMethod()
{
}


/********* Forward Euler method *********/

MpMat ForwardEuler::step(const MpMat& Y_n, const mp_num_t& h) {
    return Y_n + F_(Y_n) * h;
}

ForwardEuler::~ForwardEuler()
{
}

/********* Transform Euler method *********/

MpMat TransformEuler::step(const MpMat& Y_n, const mp_num_t& h) {

    mp_num_t half_h;
    mp_num_init(half_h);
    mp_num_div_ui(half_h, h, 2);

    MpMat Y_n_and_half = Y_n + F_(Y_n) * half_h;
    return Y_n + F_(Y_n_and_half) * h;
}

TransformEuler::~TransformEuler()
{
}

