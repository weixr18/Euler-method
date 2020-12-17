/*
    mpmat_utils.cpp: Utilities for MpMat class.

    Copyright (C) 2020 Xinran Wei <weixr0605@sina.com>
    Full LICENCE: ./LICENCE
*/

#include "mpmat_utils.h"
#include <stdio.h>
#include <stdlib.h>

SizeException::SizeException(MpMatOprType type, int row0, int col0, int row1, int col1) {
    switch (type)
    {
    case MpMatOprType::MpMatAdd:
        sprintf(reason_, \
            "Size error while add two Mpq matrixes, shape 0 is: (%d, %d) and shape 1 is: (%d, %d)", \
            row0, col0, row1, col1);
        break;
    case MpMatOprType::MpMatSub:
        sprintf(reason_, \
            "Size error while sub two Mpq matrixes, shape 0 is: (%d, %d) and shape 1 is: (%d, %d)", \
            row0, col0, row1, col1);
        break;

    case MpMatOprType::MpMatMultiply:
        sprintf(reason_, \
            "Size error while multiply two Mpq matrixes, shape 0 is: (%d, %d) and shape 1 is: (%d, %d)", \
            row0, col0, row1, col1);
        break;

    default:
        break;
    }
}

RangeException::RangeException(const char* value, const char* limination) {
    sprintf(reason_, \
        "Range Error: value %s should meet following requirement: %s", \
        value, limination);
}