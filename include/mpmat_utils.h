/*
    mpmat_utils.h: Utilities for MpMat class.

    Copyright (C) 2020 Xinran Wei <weixr0605@sina.com>
    Full LICENCE: ./LICENCE
*/

#ifndef _MPMAT_UTILS_H
#define _MPMAT_UTILS_H

#include <exception>
#include <gmp.h>

enum MpMatOprType {
    MpMatAdd,
    MpMatSub,
    MpMatMultiply,
};

struct MpMatException : public std::exception
{
};

struct SizeException : MpMatException
{
    char reason_[1024] = "Size Error.";
    SizeException(MpMatOprType type, int row0, int col0, int row1, int col1);
    const char* what() const throw ()
    {
        return reason_;
    }
};

struct RangeException : MpMatException
{
    char reason_[1024] = "Size Error.";
    RangeException(const char* value, const char* limination);
    const char* what() const throw ()
    {
        return reason_;
    }
};
#endif //_MPMAT_UTILS_H