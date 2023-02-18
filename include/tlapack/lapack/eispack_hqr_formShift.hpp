/// @file eispack_hqr_formshift.hpp
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

namespace tlapack
{
    /**
     *  @brief perform a form shift on A
     *
     */
    template <class matrix_t>
    int eispack_hqr_formShift(
        size_type<matrix_t> low,
        matrix_t &A,
        size_type<matrix_t> en,
        real_type<type_t<matrix_t>> &s,
        real_type<type_t<matrix_t>> &t,
        real_type<type_t<matrix_t>> &x,
        real_type<type_t<matrix_t>> &y,
        real_type<type_t<matrix_t>> &w )
    {
        using TA = type_t<matrix_t>;
        using idx_t = size_type<matrix_t>;
        using real_t = real_type<TA>; 

        real_t const constant1 = real_t(0.75);
        real_t const constant2 = real_t(-0.4375);

        // Grab the number of columns of A, we only work on square matrices
        const idx_t n = ncols(A);

        // Perform the checks for our arguments
        tlapack_check(n == nrows(A));

        t += x;
        for (idx_t i = low; i <= en; i++)
            A(i,i) = A(i,i) -  x;
        s = tlapack::abs(A(en, en - 1)) + tlapack::abs(A(en - 1, en - 2));
        x = constant1 * s;
        y = x;
        w = constant2 * s * s;
        return 0;

    }

} // lapack

#endif // TLAPACK_HQR_FORMSHIFT_HH
