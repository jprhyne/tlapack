/// @file cdiv.hpp
/// @author Johnathan Rhyne, CU Denver, USA
/// Adapted from @see https://netlib.org/eispack/hqr2.f
//
// Copyright (c) 2013-2022, University of Colorado Denver. All rights reserved.
//
// This file is part of <T>LAPACK.
// <T>LAPACK is free software: you can redistribute it and/or modify it under
// the terms of the BSD 3-Clause license. See the accompanying LICENSE file.

#ifndef TLAPACK_CDIV_HH
#define TLAPACK_CDIV_HH
namespace tlapack
{
    /**
     * @brief Performs complex division using real numbers storing the result in cr, ci.
     *
     * This function performs complex division on 2 ordered pairs of real numbers (ar, ai) and (br, bi)
     * and stores the result in the ordered pair (cr,ci)
     * IE, on output (cr, ci) = (ar, ai) / (br, bi)
     *
     * @tparam real_t Some kind of real number data type like float or double
     *
     * @param[in] ar 
     *      The real part of the numerator
     * @param[in] ai
     *      The imaginary part of the numerator
     * @param[in] br 
     *      The real part of the denominator
     * @param[in] bi
     *      The imaginary part of the denominator
     * @param[out] cr
     *      On output, contains the real part of (ar + ai * i) / (br + bi * i)
     * @param[out] ci
     *      On output, contains the imaginary part of (ar + ai * i) / (br + bi * i)
     */
    template <class real_t>
    void cdiv(real_t ar, real_t ai, real_t br, real_t bi, real_t &cr, real_t &ci) {
        real_t s,ars,ais,brs,bis;
        s = tlapack::abs(br) + tlapack::abs(bi);
        ars = ar/s;
        ais = ai/s;
        brs = br/s;
        bis = bi/s;
        s = brs * brs + bis * bis; 
        cr = (ars*brs + ais*bis) / s;
        ci = (ais*brs - ars*bis) / s;
    }
}
#endif // TLAPACK_HQR_HH
