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

#include <functional>


namespace tlapack
{
    template <class matrix_t>

    int hqr_subDiagonalSearch(
        size_type<matrix_t> low,
        matrix_t &A,
        size_type<matrix_t> en,
        real_type<type_t<matrix_t>> norm,
        real_type<type_t<matrix_t>> *s )
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


        for (idx_t l = en; l > low; l--) {
            *s = fabs(A(l - 1, l - 1)) + fabs(A(l, l));
            if (*s == 0.0)
                *s = norm;
            real_t tst1 = *s;
            real_t tst2 = tst1 + fabs(A(l,l - 1));
            if (tst1 == tst2)
                return l;
        }
        return low;
    }
} // lapack

#endif // TLAPACK_HQR_SUBDIAGSEARCH_HH
