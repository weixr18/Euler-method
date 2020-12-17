/*
    mpmat.h: A multiple precision 2-D matrix implementation.

    Copyright (C) 2020 Xinran Wei <weixr0605@sina.com>
    Full LICENCE: ./LICENCE
*/

#ifndef _MPMAT_H
#define _MPMAT_H

#include <iostream>
#include <exception>
#include <Eigen/Core>
#include <gmp.h>

//#define mp_num mpf_

#define FLOAT_BITS 256
#define N_DIGITS 30

#define mp_num_t mpf_t
#define mp_num_init(x) mpf_init2(x, FLOAT_BITS)
#define mp_num_inits mpf_inits
#define mp_num_set mpf_set
#define mp_num_set_d mpf_set_d
#define mp_num_set_str mpf_set_str
#define mp_num_clear mpf_clear
#define mp_num_get_d mpf_get_d
#define mp_num_get_str mpf_get_str
#define mp_num_add mpf_add
#define mp_num_sub mpf_sub
#define mp_num_mul mpf_mul
#define mp_num_div mpf_div
#define mp_num_div_ui mpf_div_ui

class MpMat {
public:

    MpMat(uint32_t row, uint32_t col);
    MpMat(const MpMat& mat);
    ~MpMat();

    void init();
    void init(char const*** M);
    MpMat copy() const;
    static MpMat zeros_like(const MpMat& mat);
    MpMat& operator= (const MpMat& mat);
    Eigen::MatrixXd to_matrix() const;

    MpMat transpose() const;
    MpMat T() const;
    MpMat operator* (const MpMat& mat) const;
    MpMat operator* (const mp_num_t& q) const;
    MpMat operator+ (const MpMat& mat) const;
    MpMat operator- (const MpMat& mat) const;

    static void print(const MpMat& mat);
    uint32_t row_num() const;
    uint32_t col_num() const;


private:
    mp_num_t** v_ = nullptr;
    uint32_t row_;
    uint32_t col_;
};


#endif //_MPMAT_H