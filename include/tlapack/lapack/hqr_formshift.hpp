/// @file hqr_formshift.hpp
/// @author Johnathan Rhyne, CU Denver, USA
/// Adapted from @see https://netlib.org/eispack/hqr2.f
//
// Copyright (c) 2013-2022, University of Colorado Denver. All rights reserved.
//
// This file is part of <T>LAPACK.
// <T>LAPACK is free software: you can redistribute it and/or modify it under
// the terms of the BSD 3-Clause license. See the accompanying LICENSE file.

#ifndef TLAPACK_HQR_FORMSHIFT_HH
#define TLAPACK_HQR_FORMSHIFT_HH

#include <functional>


namespace tlapack
{
    template <class matrix_t>
    int hqr_formShift(
        size_type<matrix_t> low,
        matrix_t &A,
        size_type<matrix_t> its,
        size_type<matrix_t> itn,
        size_type<matrix_t> en,
        size_type<matrix_t> l,
        real_type<type_t<matrix_t>> &s,
        real_type<type_t<matrix_t>> &t,
        real_type<type_t<matrix_t>> &x,
        real_type<type_t<matrix_t>> &y,
        real_type<type_t<matrix_t>> &w )
    {
        using TA = type_t<matrix_t>;
        using idx_t = size_type<matrix_t>;
        // Not really sure what this is, however seems like it is asking what the type of the real 
        // components of the elements of A, Not sure if we need this as this algorithm is only for 
        // real matrices
        using real_t = real_type<TA>; 
        using pair = std::pair<idx_t,idx_t>;

        // Grab the number of columns of A, we only work on square matrices
        const idx_t n = ncols(A);

        // Perform the checks for our arguments
        // Why is the convention to use a 'check false' as opposed to a 'check true'?
        tlapack_check(n == nrows(A));


        x = A(en, en);
        if (l == en)
            return 1;
        y = A(en - 1, en - 1);
        w = A(en, en - 1) * A(en - 1, en);
        if (l == en - 1)
            return 2;
        if (itn == 0)
            return 3;
        if ((its != 10) and (itn != 20))
            return 0;
        t += x;
        for (idx_t i = low; i <= en; i++)
            A(i,i) = A(i,i) -  x;
        s = tlapack::abs(A(en, en - 1)) + tlapack::abs(A(en - 1, en - 2));
        x = 0.75 * s;
        y = x;
        w = -0.4375 * s * s;
        return 0;

    }

} // lapack

#endif // TLAPACK_HQR_FORMSHIFT_HH
