/*
    main.cpp: Entrance of program.

    Copyright (C) 2020 Xinran Wei.

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


int main(int argc, char** argv) {

    // get args
    bool args_format_right = false;
    int i_x_n;
    int i_steps;
    args_format_right = (argc == 5)
        && (strncmp(argv[1], "-s", 10) == 0)
        && (strncmp(argv[3], "-t", 10) == 0)
        && ((i_steps = atoi(argv[2])) > 0)
        && ((i_x_n = atoi(argv[4])) > 0);

    if (!args_format_right) {
        printf("Arguments:\n");
        printf("\t -s <step-number> integer\n");
        printf("\t -t <end-time> integer\n");
        return 0;
    }

    mpf_set_default_prec(FLOAT_BITS);
    printf("default precision: %d\n", mpf_get_default_prec());

    // differencial equation
    char const* A[4][4] = {
        "0", "0", "1", "0",
        "0", "0", "0", "0",
        "1", "0", "0", "0",
        "0", "0", "0", "0",
    };
    char const* B[4][1] = {
        "-5e-7",
        "5e-7",
        "0",
        "0",
    };
    char const* C[4][4] = {
        "0", "0", "0", "0",
        "0", "-0.4", "0", "0",
        "0", "0.4", "-0.5", "0",
        "0", "0", "0.5", "0",
    };
    char const* Y[4][1] = {
        "8000",
        "2000",
        "0",
        "0"
    };

    MpMat mA(4, 4), mB(4, 1), mC(4, 4), mY(4, 1);
    mA.init((char const***)A);
    mB.init((char const***)B);
    mC.init((char const***)C);
    mY.init((char const***)Y);

    auto quardic = [&](const MpMat& Y) {return mB * Y.T() * mA * Y + mC * Y;};

    // calculation parameter
    mp_num_t x_0, x_n, step_num;
    mp_num_init(x_0);
    mp_num_init(x_n);
    mp_num_init(step_num);
    mp_num_set_d(x_n, i_x_n);
    mp_num_set_d(step_num, i_steps);
    printf("precision of x_0, x_n, step_num: %d, %d, %d\n",
        mpf_get_prec(x_0),
        mpf_get_prec(x_n),
        mpf_get_prec(step_num)
    );

    // run
    ForwardEuler forward_euler(quardic);
    MpMat res = forward_euler.run(mY, x_0, x_n, step_num);
    printf("ForwardEuler Result:\n");
    MpMat::print(res);

    TransformEuler transform_euler(quardic);
    res = transform_euler.run(mY, x_0, x_n, step_num);
    printf("TransformEuler Result:\n");
    MpMat::print(res);

    printf("Done.\n");
    mp_num_clear(x_0);
    mp_num_clear(x_n);
    mp_num_clear(step_num);
    return 0;
}