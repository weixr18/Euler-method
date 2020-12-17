/*
    mpmat.cpp: A multiple precision 2-D matrix implementation.
    Copyright (C) 2020 Xinran Wei.

    Full LICENCE: ./LICENCE
*/

#include "mpmat.h"
#include "mpmat_utils.h"

/*
* A matrix consists of mpq numbers.
*/
MpMat::MpMat(uint32_t row, uint32_t col) {
    row_ = row;
    col_ = col;
    v_ = new mp_num_t * [row_];
    for (int i = 0; i < row_; i++) {
        v_[i] = new mp_num_t[col_];
    }
}

/*
* copy construct function
*/
MpMat::MpMat(const MpMat& mat) {
    MpMat(mat.row_, mat.col_);
    for (int i = 0; i < row_; i++) {
        for (int j = 0; j < col_; j++) {
            mp_num_init(v_[i][j]);
            mp_num_set(v_[i][j], mat.v_[i][j]);
        }
    }
}



/*
* Initialize a mpq matrix by numbers.
*/
void MpMat::init(char const*** string_matrix) {

    this->init();
    for (int i = 0; i < row_; i++) {
        const char** m_row = (const char**)string_matrix + col_ * i;
        for (int j = 0; j < col_; j++) {
            const char* m_col = *(m_row + j);
            mp_num_set_str(v_[i][j], m_col, 10);
        }
    }
}


/*
* Initialize a mpq matrix.
*/
void MpMat::init() {
    for (int i = 0; i < row_; i++) {
        for (int j = 0; j < col_; j++) {
            mp_num_init(v_[i][j]);
            mp_num_set_d(v_[i][j], 0);
        }
    }
}


/*
* set value operator
*/
MpMat& MpMat::operator= (const MpMat& mat) {
    if (this != &mat)
    {
        if (v_ != nullptr) {
            this->~MpMat();
        }

        this->row_ = mat.row_;
        this->col_ = mat.col_;
        v_ = new mp_num_t * [row_];
        for (int i = 0; i < row_; i++) {
            v_[i] = new mp_num_t[col_];
        }
        for (int i = 0; i < row_; i++) {
            for (int j = 0; j < col_; j++) {
                mp_num_init(v_[i][j]);
                mp_num_set(v_[i][j], mat.v_[i][j]);
            }
        }
    }
    return *this;
}


/*
* Get a copy of the matrix.
*/
MpMat MpMat::copy() const {
    MpMat res(row_, col_);
    res.init();
    for (int i = 0; i < row_; i++) {
        for (int j = 0; j < col_; j++) {
            mp_num_set(res.v_[i][j], v_[i][j]);
        }
    }
    return res;
}

/*
* Get a copy of the matrix.
*/
MpMat MpMat::zeros_like(const MpMat& mat) {
    MpMat res(mat.row_, mat.col_);
    res.init();
    return res;
}

/*
* Print a matrix to stdout.
*/
void MpMat::print(const MpMat& mat) {
    std::cout << "-------" << mat.row_ << " * " << mat.col_ << "------" << std::endl;
    for (int i = 0; i < mat.row_; i++) {
        for (int j = 0; j < mat.col_; j++) {
            mp_num_t& a = mat.v_[i][j];
            mpf_out_str(stdout, 10, N_DIGITS, mat.v_[i][j]);
            printf("\t");
        }
        printf("\n");
    }
    printf("\n");
}



/*
* Get transpose of a mpq matrix.
*/
MpMat MpMat::T() const {
    return transpose();
}


/*
* Get transpose of a mpq matrix.
*/
MpMat MpMat::transpose() const {

    MpMat res(col_, row_);
    res.init();
    for (int i = 0; i < row_; i++) {
        for (int j = 0; j < col_; j++) {
            mp_num_set(res.v_[j][i], v_[i][j]);
        }
    }
    return res;
}


/*
* Get the sum of two matrix.
*/
MpMat MpMat::operator+(const MpMat& mat) const {

    if (this->row_ != mat.row_ || this->col_ != mat.col_) {
        throw SizeException(
            MpMatOprType::MpMatAdd,
            this->row_, this->col_,
            mat.row_, mat.col_
        );
    }

    MpMat res(row_, col_);
    res.init();
    for (int i = 0; i < row_; i++) {
        for (int j = 0; j < col_; j++) {
            mp_num_add(res.v_[i][j], this->v_[i][j], mat.v_[i][j]);
        }
    }
    return res;
}


/*
* Get the difference of two matrix.
*/
MpMat MpMat::operator-(const MpMat& mat) const {
    if (this->row_ != mat.row_ || this->col_ != mat.col_) {
        throw SizeException(
            MpMatOprType::MpMatSub,
            this->row_, this->col_,
            mat.row_, mat.col_
        );
    }
    MpMat res(row_, col_);
    res.init();
    for (int i = 0; i < row_; i++) {
        for (int j = 0; j < col_; j++) {
            mp_num_sub(res.v_[i][j], this->v_[i][j], mat.v_[i][j]);
        }
    }
    return res;
}


/*
* Get the product of two matrix.
*/
MpMat MpMat::operator*(const MpMat& mat) const {
    if (this->col_ != mat.row_) {
        throw SizeException(
            MpMatOprType::MpMatMultiply,
            this->row_, this->col_,
            mat.row_, mat.col_
        );
    }
    MpMat res(this->row_, mat.col_);
    res.init();

    mp_num_t tmp;
    mp_num_init(tmp);
    for (int i = 0; i < this->row_; i++) {
        for (int j = 0; j < mat.col_; j++) {
            mp_num_set_d(res.v_[i][j], 0);
            for (int k = 0; k < this->col_; k++) {
                mp_num_mul(tmp, this->v_[i][k], mat.v_[k][j]);
                mp_num_add(res.v_[i][j], tmp, res.v_[i][j]);
            }
        }
    }
    mp_num_clear(tmp);
    return res;
}


/*
* Get the product of a matrix and a number.
*/
MpMat MpMat::operator* (const mp_num_t& q) const {
    MpMat res(row_, col_);
    res.init();
    for (int i = 0; i < row_; i++) {
        for (int j = 0; j < col_; j++) {
            mp_num_mul(res.v_[i][j], this->v_[i][j], q);
        }
    }
    return res;
}


MpMat::~MpMat() {
    for (int i = 0; i < row_; i++) {
        for (int j = 0; j < col_; j++) {
            mp_num_clear(v_[i][j]);
        }
    }
    for (int i = 0; i < row_; i++) {
        delete v_[i];
    }
    delete v_;
    v_ = nullptr;
}
