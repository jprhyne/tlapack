/// @file eispack_hqr_doubleSubDiagonalSearch.hpp
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
     * @brief Find our lower index for our active sub window for us to shift on and compute the
     *  first column of the desired shifted matrix.
     *
     *  We use the matrix:
     *      C = H^2 - t H + d I where t is the trace and d the determinant of the bottom 
     *      2x2 sub-matrix given by:
     *      [
     *          A(en-1,en-1) A(en-1, en)
     *          A(en,en-1)   A(en,en)
     *      ]
     *
     *      This is called a Francis shift, and is shifting based on the eigenvalues of this
     *      bottom block regardless of them being complex conjugates or not.
     *
     * This function searches for an index m that satisfies one of the following conditions
     * 1) m is an index such that A(m,m-1) is "small" with respect to the first column of
     *      
     * 2) m = l, which is the top of our active window. This only happens when 1 fails, which means
     *      that Hess(C) is an irreducible Hessenberg matrix so Implicit Q applies on the entire 
     *      window
     *
     *
     * The rationale behind doing this is as follows.
     * If we would reach a point that would be small with respect to our shift, 
     * we would deflate in our matrix C (given above). Since we would no longer have a reduced 
     * Hessenberg matrix we can no longer guarantee the new Hessenberg after qr iteration would 
     * be the same as if we did this explicitly, which would mean our implicitly computed Q may differ
     * which could make this entire algorithm useless.
     *
     * While it is not directly obvious why starting at column 
     *
     * A good writeup of the Implicit Q Theorem can be found at: 
     *  https://www.cs.utexas.edu/users/flame/laff/alaff-beta/chapter10-implicit-q-theorem.html
     *
     *  Note: Not all iterations on this theorem require the sub-diagonal elements to be positive,
     *  just non-zero
     *
     * Input:   A   -- Our matrix
     *          en  -- The index of the eigenvalue we are searching for. This is the bottom
     *              of our active sub window
     *          l   -- The top of our active sub window.
     *          x   -- A(en,en) along with any accumulated shifts or the value computed by exceptional shift
     *          y   -- A(en-1,en-1) or the value computed by exceptional shift
     *          w   -- A(en-1,en) * A(en,en-1) or the value computed by exceptional shift
     *
     * Output:  p   -- The first element of the first column of the shifted matrix given above
     *          q   -- The second element of the first column of the shifted matrix given above
     *          r   -- The third element of the first column of the shifted matrix given above
     */
    template <class matrix_t>
    int eispack_hqr_doubleSubDiagonalSearch(
        matrix_t &A,
        size_type<matrix_t> en,
        size_type<matrix_t> l,
        real_type<type_t<matrix_t>> x,
        real_type<type_t<matrix_t>> y,
        real_type<type_t<matrix_t>> w,
        real_type<type_t<matrix_t>> &p,
        real_type<type_t<matrix_t>> &q,
        real_type<type_t<matrix_t>> &r
        )
    {
        using TA = type_t<matrix_t>;
        using idx_t = size_type<matrix_t>;
        using real_t = real_type<TA>; 

        // Grab the number of columns of A, we only work on square matrices
        const idx_t n = ncols(A);

        // Perform the checks for our arguments
        tlapack_check(n == nrows(A));

        idx_t m; // Return index
        real_t s;// Local variable
        /*
         * [p,q,r] appears to be some kind of reflector/shift and we are trying to find 
         * an index where this would result in "0" prematurely. This initial scaled p,q,r
         * will then be used as the first step of our QR Iteration. 
         */
        for (m = en - 2; m <= en && m >= l; m--) {
            // This section is creating the vector
            // c from page 120 of Matrix Algorithms Vol
            // II: Eigensystems. (Algorithm 3.1)
            //
            // Note: We do not scale like the given algorithm does immediately nor does it 
            // generate the householder vector u, alpha.
            r = x - A(m,m); 
            s = y - A(m,m); 
            p = (r * s - w) / A(m + 1, m) + A(m, m + 1);
            q = A(m+1,m+1) - A(m,m) - r - s;
            r = A(m+2,m+1);
            s = tlapack::abs(p) + tlapack::abs(q) + tlapack::abs(r);
            p = p / s;
            q = q / s;
            r = r / s;
            if (m == l) break; // We can't go past l anyway, so we stop.
            real_t tst1 = tlapack::abs(p) * (tlapack::abs(A(m-1, m-1)) + tlapack::abs(A(m,m)) + tlapack::abs(A(m+1,m+1)));
            real_t tst2 = tst1 + tlapack::abs(A(m,m-1)) * (tlapack::abs(q) + tlapack::abs(r));
            if (tst1 == tst2) break; // Means the sub diagonal elements are "small", so m is our desired index

        }
        // Sets values below the first
        // to be zero past m + 1 up to en
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
