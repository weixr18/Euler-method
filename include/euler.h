#ifndef _EULER_H
#define _EULER_H

#include "mpmat.h"
#include <functional>

class EulerMethod {
public:
    EulerMethod(std::function<MpMat(const MpMat&)> F);
    virtual ~EulerMethod();
public:
    std::function<MpMat(const MpMat&)> F_;
    MpMat run(const MpMat& Y_0, const mp_num_t& t0, const mp_num_t& t, const mp_num_t& n);
    virtual MpMat step(const MpMat& Y_n, const mp_num_t& h) = 0;
};


class ForwardEuler : public EulerMethod
{
public:
    using EulerMethod::EulerMethod;
    ~ForwardEuler();
    MpMat step(const MpMat& Y_n, const mp_num_t& h);
};

class TransformEuler : public EulerMethod
{
public:
    using EulerMethod::EulerMethod;
    ~TransformEuler();
    MpMat step(const MpMat& Y_n, const mp_num_t& h);
};


#endif //_EULER_H