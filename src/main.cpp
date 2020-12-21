/*
    main.cpp: Entrance of program.

    Copyright (C) 2020 Xinran Wei <weixr0605@sina.com>

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
    USA
*/


#include <functional>
#include <iostream>
#include <string.h>

#include "mpmat.h"
#include "mpmat_utils.h"
#include "euler.h"

char const* A[EQ_NUM][EQ_NUM] = {
        "0", "0", "1", "0",
        "0", "0", "0", "0",
        "1", "0", "0", "0",
        "0", "0", "0", "0",
};
char const* B[EQ_NUM][1] = {
    "-1e-5",
    "1e-5",
    "0",
    "0",
};
char const* C[EQ_NUM][EQ_NUM] = {
    "0", "0", "0", "0",
    "0", "-0.4", "0", "0",
    "0", "0.4", "-0.5", "0",
    "0", "0", "0.5", "0",
};
char const* Y[EQ_NUM][1] = {
    "8000",
    "2000",
    "0",
    "0"
};

char const* M[EQ_NUM][EQ_NUM] = {
    "0.1", "0", "0.1", "1e-9",
    "0.1", "0.4", "0.1", "0",
    "0", "0.4", "0.5", "0",
    "0", "0", "0.5", "0",
};
char const* L[EQ_NUM][1] = {
    "6000",
    "3400",
    "6500",
    "2000",
};

void get_args(int argc, char** argv, int& boundary, int& x_n) {

    bool args_format_right = true;
    boundary = -1;
    x_n = -1;
    for (int i = 1; i + 1 < argc; i += 2) {
        if (strncmp(argv[i], "-b", 10) == 0) {
            boundary = atoi(argv[i + 1]);
        }
        else if (strncmp(argv[i], "-x", 10) == 0) {
            x_n = atoi(argv[i + 1]);
        }
        else {
            args_format_right = false;
        }
    }

    if (boundary <= 0 || x_n < 0 || !args_format_right) {
        printf("Usage:\n");
        printf("\t (1) ./main -b <boundary> -x <value>\n");
        printf("Arguments:\n");
        printf("\t -b <boundary> Integer. Error boundary.\n");
        printf("\t -x <value> Integer. Terminal value of independant variable x.\n");
        exit(0);
    }
}

/*
    Input:  x_0=0, x_n, b
    Set:    F(·), Y_0, M, L
    Calc:   h, NUM_STEPS, FLOAT_BITS
    Output: Y(x_n)
*/

int main(int argc, char** argv) {

    // get args 
    int x_0 = 0,    // x_0
        x_n,        // target x
        boundary;   // boudary
    get_args(argc, argv, boundary, x_n);

    // differencial equation
    MpMat mpm_A(EQ_NUM, EQ_NUM), mpm_B(EQ_NUM, 1),
        mpm_C(EQ_NUM, EQ_NUM), mpm_Y(EQ_NUM, 1),
        mpm_M(EQ_NUM, EQ_NUM), mpm_L(4, 1);
    mpm_A.init((char const***)A);
    mpm_B.init((char const***)B);
    mpm_C.init((char const***)C);
    mpm_Y.init((char const***)Y);
    mpm_M.init((char const***)M);
    mpm_L.init((char const***)L);
    auto quardic = [&](const MpMat& mpm_Y) {return mpm_B * mpm_Y.T() * mpm_A * mpm_Y + mpm_C * mpm_Y;};

    // run
    printf("Running forward Euler algorithm.\n");
    ForwardEuler forward_euler(quardic, mpm_M, mpm_L, x_0, x_n, boundary);
    MpMat res = forward_euler.run(mpm_Y);
    printf("Result:\n");
    MpMat::print(res);

    /*
    printf("Running transform Euler algorithm.\n");
    TransformEuler transform_euler(quardic, mpm_M, mpm_L, x_0, x_n, boundary);
    res = transform_euler.run(mpm_Y);
    printf("Result:\n");
    MpMat::print(res);
    */

    printf("Mission completed.\n");
    return 0;
}
