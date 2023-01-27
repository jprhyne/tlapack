/// @file hqr_doubleSubDiagonalSearch.hpp
/// @author Johnathan Rhyne, CU Denver, USA
/// Adapted from @see https://netlib.org/eispack/hqr2.f
//
// Copyright (c) 2013-2022, University of Colorado Denver. All rights reserved.
//
// This file is part of <T>LAPACK.
// <T>LAPACK is free software: you can redistribute it and/or modify it under
// the terms of the BSD 3-Clause license. See the accompanying LICENSE file.

#ifndef TLAPACK_HQR_DOUBLESUBDIAGSEARCH_HH
#define TLAPACK_HQR_DOUBLESUBDIAGSEARCH_HH



namespace tlapack
{
    /** 
     * TODO: Determine what this function does
     */
    template <class matrix_t>
    int hqr_doubleSubDiagonalSearch(
        matrix_t &A,
        size_type<matrix_t> en,
        size_type<matrix_t> l,
        real_type<type_t<matrix_t>> &s,
        real_type<type_t<matrix_t>> x,
        real_type<type_t<matrix_t>> y,
        real_type<type_t<matrix_t>> w,
        real_type<type_t<matrix_t>> &p,
        real_type<type_t<matrix_t>> &q,
        real_type<type_t<matrix_t>> &r,
        real_type<type_t<matrix_t>> &zz )
    {
        using TA = type_t<matrix_t>;
        using idx_t = size_type<matrix_t>;
        using real_t = real_type<TA>; 

        // Grab the number of columns of A, we only work on square matrices
        const idx_t n = ncols(A);

        // Perform the checks for our arguments
        tlapack_check(n == nrows(A));

        idx_t m;
        for (m = en - 2; m <= en && m >= l; m--) {
            zz = A(m,m);
            r = x - zz;
            s = y - zz;
            p = (r * s - w) / A(m + 1, m) + A(m, m + 1);
            q = A(m+1,m+1) - zz - r - s;
            r = A(m+2,m+1);
            s = tlapack::abs(p) + tlapack::abs(q) + tlapack::abs(r);
            p = p / s;
            q = q / s;
            r = r / s;
            if (m == l) break;
            real_t tst1 = tlapack::abs(p) * (tlapack::abs(A(m-1, m-1)) + tlapack::abs(zz) + tlapack::abs(A(m+1,m+1)));
            real_t tst2 = tst1 + tlapack::abs(A(m,m-1)) * (tlapack::abs(q) + tlapack::abs(r));
            if (tst1 == tst2) break;

        }
        for (idx_t i = m + 2; i <= en; i++) {
            A(i, i - 2) = real_t(0.0);
            if (i == m + 2)
                continue;
            A(i, i - 3) = real_t(0.0);
        }
        return m;
    }

} // lapack

#endif // TLAPACK_HQR_DOUBLESUBDIAGSEARCH_HH
