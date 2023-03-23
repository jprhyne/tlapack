/// @file eispack_hqr_qrIteration.hpp
/// @author Johnathan Rhyne, CU Denver, USA
/// Adapted from @see https://netlib.org/eispack/hqr2.f
//
// Copyright (c) 2013-2022, University of Colorado Denver. All rights reserved.
//
// This file is part of <T>LAPACK.
// <T>LAPACK is free software: you can redistribute it and/or modify it under
// the terms of the BSD 3-Clause license. See the accompanying LICENSE file.

#ifndef TLAPACK_HQR_QRITERATION_HH
#define TLAPACK_HQR_QRITERATION_HH

#include <stdbool.h>
#include "tlapack/lapack/lapy2.hpp"

namespace tlapack
{
    /**
     * @brief perform a double implicit QR iteration on A
     */
    template <class matrix_t>
    int eispack_hqr_qrIteration(
        matrix_t &A,
        size_type<matrix_t> en,
        size_type<matrix_t> l,
        real_type<type_t<matrix_t>> &s,
        real_type<type_t<matrix_t>> &x,
        real_type<type_t<matrix_t>> &y,
        real_type<type_t<matrix_t>> &p,
        real_type<type_t<matrix_t>> &q,
        real_type<type_t<matrix_t>> &r,
        real_type<type_t<matrix_t>> &zz,
        size_type<matrix_t> m,
        bool want_Q,
        size_type<matrix_t> low,
        size_type<matrix_t> igh,
        matrix_t &Q)
    {
        using TA = type_t<matrix_t>;
        using idx_t = size_type<matrix_t>;
        using real_t = real_type<TA>; 

        real_t zero = real_t(0);

        // Grab the number of columns of A, we only work on square matrices
        const idx_t n = ncols(A);

        // Perform the checks for our arguments
        tlapack_check(n == nrows(A));

        idx_t i,j,k;
        bool notLast;
        for (k = m; k <= en - 1; k++) {
            notLast = k != en - 1;
            if (k != m) {
                p = A(k,k-1);
                q = A(k+1,k-1);
                r = zero;
                if (notLast) {
                    r = A(k+2,k-1);
                }
                x = tlapack::abs(p) + tlapack::abs(q) + tlapack::abs(r);
                if (x == zero) {
                    continue;
                }
                p = p / x;
                q = q / x;
                r = r / x;
            }
            real_t insideSqrt = p * p + q * q + r * r;
            if (p >= zero) {
                s = sqrt(insideSqrt);
            } else {
                s = -sqrt(insideSqrt);
            }
            if (k != m) {
                A(k,k-1) = -s * x;
            } else if (l != m) {
                A(k,k-1) = -A(k,k-1);
            }
            p = p + s;
            x = p / s;
            y = q / s;
            zz = r / s;
            q = q / p;
            r = r / p;
            if (notLast) {
                idx_t upperBound = (want_Q) ? n - 1 : en;
    //c     .......... row modification ..........
                for (j = k; j <= upperBound; j++) {
                    p = A(k,j) + q * A(k+1,j) + r * A(k+2,j);
                    A(k,j) = A(k,j) - p * x;
                    A(k+1,j) = A(k+1,j) - p * y;
                    A(k+2,j) = A(k+2,j) - p * zz;
                }
                if (en <= k+3) {
                    j = en;
                } else {
                    j = k+3;
                }
    //c     .......... column modification ..........
                idx_t lowerBound = (want_Q) ? 0 : l;
                for (i = lowerBound; i <= j; i++) {
                    p = x * A(i,k) + y * A(i,k+1) + zz * A(i,k+2);
                    A(i,k) = A(i,k) - p;
                    A(i,k+1) = A(i,k+1) - p * q;
                    A(i,k+2) = A(i,k+2) - p * r;
                }
                if (want_Q) {
    //c     .......... accumulate transformations  ..........
                    for (i = low; i <= igh; i++) {
                        p = x * Q(i,k) + y * Q(i,k+1) + zz * Q(i,k+2);
                        Q(i,k) = Q(i,k) - p;
                        Q(i, k + 1) = Q(i,k + 1) - p * q;
                        Q(i, k + 2) = Q(i,k + 2) - p * r;
                    }
                }
            } else {
    //c     .......... row modification ..........
                idx_t upperBound = (want_Q) ? n - 1 : en;
                for (j = k; j <= upperBound; j++) {
                    p = A(k,j) + q * A(k+1,j);
                    A(k,j) = A(k,j) - p * x;
                    A(k+1,j) = A(k+1,j) - p * y; 
                }
                if (en <= k+3) {
                    j = en;
                } else {
                    j = k+3;
                }
    //c     .......... column modification ..........
                idx_t lowerBound = (want_Q) ? 0 : l;
                for (i = lowerBound; i <= j; i++) {
                    p = x * A(i,k) + y * A(i,k+1);
                    A(i,k) = A(i,k) - p;
                    A(i,k+1) = A(i,k+1) - p * q;
                }
                if (want_Q) {
    //c     ...........accumulate transformations  ..........
                    for (i = low; i <= igh; i++) {
                        p = x * Q(i,k) + y * Q(i,k+1);
                        Q(i,k) = Q(i,k) - p;
                        Q(i, k + 1) = Q(i,k+1) - p * q;
                    }
                }
            }
        }
        return 0;
    }

    template <class matrix_t, class vector_t>
    int eispack_comqr_qrIteration(
        matrix_t &A,
        size_type<matrix_t> en,
        size_type<matrix_t> l,
        complex_type<type_t<matrix_t>> &s,
        complex_type<type_t<matrix_t>> &x,
        complex_type<type_t<matrix_t>> &y,
        complex_type<type_t<matrix_t>> &zz,
        vector_t &eigs,
        bool want_Q,
        size_type<matrix_t> low,
        size_type<matrix_t> igh,
        matrix_t &Q)
    {
        using TA = type_t<matrix_t>;
        using idx_t = size_type<matrix_t>;
        using real_t = real_type<TA>; 
        using complex_t = complex_type<TA>; 

        real_t rZero = real_t(0);
        complex_t cZero = complex_t(0);

        // Grab the number of columns of A, we only work on square matrices
        const idx_t n = ncols(A);
        real_t norm;

        // Perform the checks for our arguments
        tlapack_check(n == nrows(A));

        // Reduce to triangle (rows)
        for (idx_t i = l + 1; i <= en; i++) {
            s = complex_t(A(i, i - 1).real(), s.imag());
            norm = lapy2(tlapack::abs(A(i - 1, i - 1)), s.real());
            x = A(i - 1, i - 1) / norm; // Never checks if A is the zero matrix...
            eigs[i-1] = x;
            A(i - 1, i - 1) = complex_t(norm,0);
            A(i, i - 1) = complex_t(0, s.real()/norm);
            idx_t upperBound = (want_Q) ? (n - 1) : (en);
            for (idx_t j = i; j <= upperBound; j++) {
                y = A(i - 1, j);
                zz = A(i,j);
                real_t hr = x.real() * y.real() + x.imag() * y.imag() + A(i,i-1).imag() * zz.real();
                real_t hi = x.real() * y.imag() - x.imag() * y.real() + A(i,i-1).imag() * zz.imag();
                A(i-1,j) = complex_t(hr,hi);
                hr = x.real() * zz.real() - x.imag() * zz.imag() - A(i,i-1).imag() * y.real();
                hi = x.real() * zz.imag() + x.imag() * zz.real() - A(i,i-1).imag() * y.imag();
                A(i,j) = complex_t(hr,hi);
                //A(i,j) = x * zz - complex_t(A(i, i - 1).imag() * y.real(), A(i, i - 1).imag() * y.imag());
            }
        }

        s = complex_t(s.real(), A(en,en).imag());
        if (s.imag() != rZero) {
            norm = tlapack::abs(A(en,en));
            s = A(en,en) / norm;
            A(en,en) = complex_t(norm, rZero);
            if (want_Q && en != n-1) {
                for (idx_t j = en + 1; j < n; j++) {
                    y = A(en,j);
                    A(en,j) = s * conj(y);
                }
            }
        }


        // Inverse operation (columns)
        for (idx_t j = l + 1; j <= en; j++) {
            x = eigs[j - 1];
            idx_t lowerBound = (want_Q) ? (0) : (l);
            for (idx_t i = lowerBound; i <= j; i++) {
                y = (i == j) ? (complex_t(A(i,j-1).real(), 0)) : (A(i,j-1));
                zz = A(i,j);
                real_t hr = x.real() * y.real() - x.imag() * y.imag() + A(j,j-1).imag() * zz.real();
                real_t hi = (i == j) ? (A(i,j-1).imag()) : (x.real() * y.imag() + x.imag() * y.real() + A(j,j-1).imag() * zz.imag());
                A(i,j-1) = complex_t(hr,hi);
                hr = x.real() * zz.real() + x.imag() * zz.imag() - A(j,j-1).imag() * y.real();
                hi = x.real() * zz.imag() - x.imag() * zz.real() - A(j,j-1).imag() * y.imag();
                A(i,j) = complex_t(hr,hi);
            }
            if (want_Q) {
                for (idx_t i = low; i <= igh; i++) {
                    y = Q(i,j-1);
                    zz = Q(i,j);
                    real_t zr = x.real() * y.real() - x.imag() * y.imag() + A(j,j-1).imag() * zz.real();
                    real_t zi = x.real() * y.imag() + x.imag() * y.real() + A(j,j-1).imag() * zz.imag();
                    Q(i,j-1) = complex_t(zr,zi);
                    zr = x.real() * zz.real() + x.imag() * zz.imag() - A(j,j-1).imag() * y.real();
                    zi = x.real() * zz.imag() - x.imag() * zz.real() - A(j,j-1).imag() * y.imag();
                    Q(i,j) = complex_t(zr,zi);
                }
            }
        }
        if (s.imag() != rZero) {
            idx_t lowerBound = (want_Q) ? (0) : (l);
            for (idx_t i = lowerBound; i <= en; i++) {
                y = A(i,en);
                A(i,en) = s * y;
            }
            if (want_Q) {
                for (idx_t i = low; i <= igh; i++) {
                    y = Q(i,en);
                    Q(i,en) = s * y;
                }
            }
        }
        return 0;
    }
} // lapack

#endif // TLAPACK_HQR_QRITER_HH
