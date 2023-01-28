/// @file hqr_schurToEigen.hpp
/// @author Johnathan Rhyne, CU Denver, USA
/// Adapted from @see https://netlib.org/eispack/hqr2.f
//
// Copyright (c) 2013-2022, University of Colorado Denver. All rights reserved.
//
// This file is part of <T>LAPACK.
// <T>LAPACK is free software: you can redistribute it and/or modify it under
// the terms of the BSD 3-Clause license. See the accompanying LICENSE file.

#ifndef TLAPACK_HQR_SCHURTOEIGEN_HH
#define TLAPACK_HQR_SCHURTOEIGEN_HH

#include "tlapack/lapack/cdiv.hpp"

namespace tlapack
{
    template <class matrix_t, class vector_t>
    int hqr_schurToEigen(
        matrix_t &U,
        size_type<matrix_t> low,
        size_type<matrix_t> igh,
        vector_t wr,
        vector_t wi,
        matrix_t &Z,
        real_type<type_t<matrix_t>> norm )
    {
        using TA = type_t<matrix_t>;
        using idx_t = size_type<matrix_t>;
        using real_t = real_type<TA>; 

        // Grab the number of columns of A, we only work on square matrices
        const idx_t n = ncols(U);

        // Perform the checks for our arguments
        tlapack_check(n == nrows(U));
        tlapack_check((idx_t)size(wr) == n);
        tlapack_check((idx_t)size(wi) == n);

        tlapack_check((n == ncols(Z)));

        // Consider adding some checks for trivial cases if needed

        // Now, we actually start porting the behavior

        // Initialize some variables at the beginning.
        idx_t en,m;
        real_t p,q,w,r,zz,s,t,tst1,tst2,x,y,ra,sa,vr,vi;
        real_t zero = real_t(0);
        real_t one = real_t(1);
        real_t two = real_t(2);
        real_t scalingFactor = real_t(0.01);

        if (norm < zero) {
            tlapack_error(1,"Norm cannot be negative. Check the value then call again");
            return 1;
        } 
        if (norm == zero) {
            tlapack_error(2,"The norm is 0. Ensure U is actually the 0 matrix. If so, nothing needs to be done");
            return 2;
        }
        for (en = n - 1; en >= 0 && en <= n - 1; en--) {
            p = wr[en];
            q = wi[en];
            // Note: We do nothing if the imaginary part is positive. This comes from our 
            // construction of the eigenvectors. See the above documentation for the structure.
            if (q == zero) {
                // This means we have a real eigenvalue
                m = en;
                U(en,en) = one;
                if (en != 0) {
                    for (idx_t i = en - 1; i >= 0 && i <= n; i--) {
                        w = U(i,i) - p;
                        r = zero;
                        for (idx_t j = m; j <= en; j++)
                            r += U(i,j) * U(j,en);
                        if (wi[i] < zero) {
                            zz = w;
                            s = r;
                            continue;
                        } 
                        m = i;
                        if (wi[i] == zero) {
                            t = w;
                            if (t == zero) {
                                tst1 = norm;
                                t = tst1;
                                do {
                                    t *= scalingFactor;
                                    tst2 = norm + t;
                                } while (tst2 >= tst1);
                            }
                            U(i,en) = -r / t;
                        } else {
                            x = U(i, i + 1);
                            y = U(i + 1, i);
                            q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
                            t = (x * s - zz * r) / q;
                            U(i, en) = t;
                            if (tlapack::abs(x) > tlapack::abs(zz))
                                U(i + 1, en) = (-r - w * t) / x;
                            else 
                                U(i + 1, en) = (-s - y * t) / zz;
                        }
                        // Overflow control
                        t = tlapack::abs(U(i, en));
                        if ( t != zero ) {
                            tst1 = t;
                            tst2 = tst1 + one / tst1;
                            if (tst2 <= tst1) {
                                for (idx_t j = i; j <= en; j++)
                                    U(j, en) /= t;
                            }
                        }
                    }
                }
            } else if (q < zero) {
                m = en - 1;
                //last vector component chosen imaginary so that the eigenvector
                // matrix is triangular
                if (tlapack::abs(U(en, en - 1)) <= tlapack::abs(U(en - 1, en))){
                    real_t a = zero;
                    real_t b = zero;
                    cdiv(zero,  -U(en - 1, en),  U(en - 1, en - 1) - p, q, a, b);
                    U(en - 1, en - 1) = a;
                    U(en - 1, en) = b;
                } else {
                    U(en - 1,en - 1) = q / U(en,en - 1);
                    U(en - 1,en) = -(U(en,en) - p) / U(en,en - 1);
                }
                U(en,en - 1) = zero;
                U(en,en) = one;
                if (en != 1) {
                    for (idx_t i = en - 2; i >= 0 && i <= n; i--) {
                        w = U(i,i) - p;
                        ra = zero;
                        sa = zero;
                        for (idx_t j = m; j <= en; j++) {
                            ra = ra + U(i,j) * U(j,en - 1);
                            sa = sa + U(i,j) * U(j,en);
                        }
                        if (wi[i] < zero) {
                            zz = w;
                            r = ra;
                            s = sa;
                            continue;
                        }
                        m = i;
                        if (wi[i] == zero) {
                            real_t a = zero;
                            real_t b = zero;
                            cdiv(-ra, -sa, w, q, a, b);
                            U(i,en - 1) = a;
                            U(i,en) = b;
                        } else {
                            // solve complex equations
                            x = U(i, i + 1);
                            y = U(i + 1, i);
                            vr = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i] - q * q;
                            vi = (wr[i] - p) * two * q;
                            if (vr == zero && vi == zero) {
                                tst1 = norm * (tlapack::abs(w) + tlapack::abs(q) + tlapack::abs(x) + tlapack::abs(y) + tlapack::abs(zz));
                                vr = tst1;
                                do {
                                    vr = scalingFactor * vr;
                                    tst2 = tst1 + vr;
                                } while (tst2 > tst1);
                            }
                            real_t a = zero;
                            real_t b = zero;
                            cdiv(x * r - zz * ra + q * sa, x * s - zz * sa - q * ra, vr, vi,  a, b);
                            U(i,en - 1) = a;
                            U(i,en) = b;
                            if (tlapack::abs(x) > tlapack::abs(zz) + tlapack::abs(q)) {
                                U(i + 1, en - 1) = (-ra - w * U(i,en - 1) + q * U(i,en)) / x;
                                U(i + 1, en) = (-sa - w * U(i,en) - q * U(i,en - 1)) / x;
                            } else {
                                a = zero;
                                b = zero;
                                cdiv(-r - y * U(i, en - 1), -s - y * U(i,en), zz, q, a, b);
                                U(i + 1, en - 1) = a;
                                U(i + 1, en) = b;
                            }
                        }
                        // overflow control
                        t = tlapack::abs(U(i, en - 1));
                        if (t < tlapack::abs(U(i, en)))
                            t = tlapack::abs(U(i, en));
                        if ( t != zero ) {
                            tst1 = t;
                            tst2 = tst1 + one / tst1;
                            if (tst2 <= tst1) {
                                for (idx_t j = i; j <= en; j++) {
                                    U(j,en - 1) = U(j,en - 1) / t;
                                    U(j,en) = U(j,en) / t;
                                }
                            }
                        }
                    }
                }
            }
        }

        // This section produces the vectors of isolated roots. 
        // In our testing we do nothing with this, however this is kept in the event
        // that balance or a similar function is ported and/or used
        for (idx_t i = 0; i < n && (i < low || i > igh); i++) {
            for (idx_t j = 0; j < n; j++)
                Z(i,j) = U(i,j);
        }
        // This may be able to be refactored into smarter matrix multiplication
        for (idx_t j = n - 1; j >= low && j <= n - 1; j--) {
            m = j;
            if (m > igh)
                m = igh;
            for (idx_t i = low; i <= igh; i++) {
                zz = zero;
                for (idx_t k = low; k <= m; k++) 
                    zz = zz + Z(i,k) * U(k,j);
                Z(i,j) = zz;
            }
        }

        return 0;
    }

} // lapack

#endif // TLAPACK_HQR_SCHURTOEIGEN_HH
