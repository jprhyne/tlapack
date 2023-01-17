/// @file hqr_doubleSubDiagonalSearch.hpp
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
    int hqr_doubleSubDiagonalSearch(
        matrix_t &A,
        size_type<matrix_t> en,
        size_type<matrix_t> l,
        real_type<type_t<matrix_t>> *s,
        real_type<type_t<matrix_t>> x,
        real_type<type_t<matrix_t>> y,
        real_type<type_t<matrix_t>> w,
        real_type<type_t<matrix_t>> *p,
        real_type<type_t<matrix_t>> *q,
        real_type<type_t<matrix_t>> *r,
        real_type<type_t<matrix_t>> *zz )
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
        tlapack_check_false(n != nrows(A));

        idx_t m;
        for (m = en - 2; m >= 2; m--) {
            *zz = A(m,m);
            *r = x - *zz;
            *s = y - *zz;
            *p = (*r * *s - w) / A(m + 1, m) + A(m, m + 1);
            *q = A(m+1,m+1) - *zz - *r - *s;
            *r = A(m+2,m+1);
            *s = fabs(*p) + fabs(*q) + fabs(*r);
            *p = *p / *s;
            *q = *q / *s;
            *r = *r / *s;
            if (m == l) break;
            real_t tst1 = fabs(*p) * (fabs(b0(m-1, m-1)) + fabs(*zz) + fabs(b0(m+1,m+1)));
            real_t tst2 = tst1 + fabs(b0(m,m-1)) * (fabs(*q) + fabs(*r));
            if (tst1 == tst2) break;

        }
        for (idx_t i = m + 2; i <= en; i++) {
            A(i, i - 2) = 0.0;
            if (i == m + 2)
                continue;
            A(i, i - 3) = 0.0;
        }
        return m;
    }

} // lapack

#endif // TLAPACK_HQR_FORMSHIFT_HH
