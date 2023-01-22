/// @file hqr.hpp
/// @author Johnathan Rhyne, CU Denver, USA
/// Adapted from @see https://netlib.org/eispack/hqr2.f
//
// Copyright (c) 2013-2022, University of Colorado Denver. All rights reserved.
//
// This file is part of <T>LAPACK.
// <T>LAPACK is free software: you can redistribute it and/or modify it under
// the terms of the BSD 3-Clause license. See the accompanying LICENSE file.

#ifndef TLAPACK_HQR_SUBDIAGSEARCH_HH
#define TLAPACK_HQR_SUBDIAGSEARCH_HH

namespace tlapack
{
    /**
     * @brief search for the first (coming from the bottom) index where there is a small subdiagonal element of A
     */
    template <class matrix_t>
    int hqr_subDiagonalSearch(
        size_type<matrix_t> low,
        matrix_t &A,
        size_type<matrix_t> en,
        real_type<type_t<matrix_t>> norm,
        real_type<type_t<matrix_t>> &s )
    {
        using TA = type_t<matrix_t>;
        using idx_t = size_type<matrix_t>;
        using real_t = real_type<TA>; 

        // Grab the number of columns of A, we only work on square matrices
        const idx_t n = ncols(A);

        // Perform the checks for our arguments
        tlapack_check(n == nrows(A));


        for (idx_t l = en; l <= en && l > low; l--) {
            s = tlapack::abs(A(l - 1, l - 1)) + tlapack::abs(A(l, l));
            if (s == 0.0)
                s = norm;
            real_t tst1 = s;
            real_t tst2 = tst1 + tlapack::abs(A(l,l - 1));
            if (tst1 == tst2)
                return l;
        }
        return low;
    }
} // lapack

#endif // TLAPACK_HQR_SUBDIAGSEARCH_HH
