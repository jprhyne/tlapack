/// @file hqr.hpp
/// @author Johnathan Rhyne, CU Denver, USA
/// Adapted from @see https://netlib.org/eispack/hqr2.f
//
// Copyright (c) 2013-2022, University of Colorado Denver. All rights reserved.
//
// This file is part of <T>LAPACK.
// <T>LAPACK is free software: you can redistribute it and/or modify it under
// the terms of the BSD 3-Clause license. See the accompanying LICENSE file.

#ifndef TLAPACK_HQR_HH
#define TLAPACK_HQR_HH


#include <functional>

#include "tlapack/lapack/hqr_subDiagonalSearch.hpp"
#include "tlapack/lapack/hqr_formShift.hpp"
#include "tlapack/lapack/hqr_doubleSubDiagonalSearch.hpp"
#include "tlapack/lapack/hqr_qrIteration.hpp"


namespace tlapack
{
    template <class matrix_t, class vector_t>

        // The final parameter *norm is supposed to store the norm of A on 
        // exit. This is to make sure we can return an actual exit
        // code instead of the hack I initially did for hqr.c
    int hqr(
        bool want_q,
        size_type<matrix_t> low,
        size_type<matrix_t> igh,
        matrix_t &A,
        vector_t &wr,
        vector_t &wi,
        matrix_t &Q,
        real_type<type_t<matrix_t>> &norm )
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
        tlapack_check((idx_t)size(wr) == n);
        tlapack_check((idx_t)size(wi) == n);

        if (want_q) {
            // If we want the Schur Vectors, we need to make sure that Q is the right size
            tlapack_check((n == ncols(Q)));
        }

        // Consider adding some checks for trivial cases if needed

        // Some other functions create a workspace, see if this is desirable
        // to preserve A.

        // Now, we actually start porting the behavior

        // Initialize some variables at the beginning. We use the following
        // general conversion schema
        // double   -> real_t
        // int      -> idx_t
        // double*  -> matrix_t || vector_t (context dependent)
        // The * is because, I ```believe``` &norm is a ptr to a variable named norm
        norm = 0.0; 
        idx_t i,j,en,m,itn,its,l,retVal;
        real_t p,q,r,s,t,w,x,y,z,zz;
        // Construct the sum of the absolute elements of A and check for isolated
        // eigenvalues
        for (i = 0; i < n; i++) {
            for (j = (i == 0) ? (i):(i-1); j <= n; j++) {
                norm += tlapack::abs(A(i,j)); // See if fabs is the correct absolute value function. May need to change it
            }
            if (i >= low && i <= igh)
                continue;
            wr[i] = A(i, i);
            wi[i] = 0.0;
        }

        // Initialize some more variables
        // en denotes the index of wr/wi we are trying to find. We start by
        //      trying to find the eigenvalues from the bottom to the top
        // t is a shifting variable that I am not too sure the theory behind
        //      why it works
        // itn is the total number of QR iterations we allow for the entire matrix
        //      The heuristic chosen was 30 * n
        en = igh;
        t = 0.0;
        itn = 30 * n;

        // This is a hack that we can possibly refactor out. We
        // set a flag to determine if we did a QR step, and if so we
        // reset its to be 0 as we have found another eigenvalue
        bool didQRStep = false;
        while (en >= low) {
            if (!didQRStep)
                its = 0;
            didQRStep = false;

            // Perform a subdiagonal search to determine where the first
            // small subdiagonal element is
            l = hqr_subDiagonalSearch(low,A,en,norm,&s);
            // Perform a form shift as well as check if we have found an eigenvalue
            retVal = hqr_formShift(low, A, its, itn, en, l, &s, &t, &x, &y, &w);
            // In order to emulate the behavior of the fortran code, instead 
            // of jumping to the right code inside there, we instead set a
            // return value and check what it is on exit
            // if retVal did not change from 0, we went through the entire
            // form shift section
            // if retVal is 1, then we found a single root
            // if retVal is 2, then we found a double root
            // if retVal is 3, then we did not converge and terminate with error
            // Any other value means there was an error inside formShift, and is
            // not supported by this implementation
            switch (retVal) {
                case 0:
                    // Full termination, so we need to do a QR Step
                    its += 1;
                    itn -= 1;
                    m = tlapack::hqr_doubleSubDiagonalSearch(A, en, l, &s, x, y, w, &p, &q, &r, &zz);
                    // double qr step
                    hqr_qrIteration(A, en, l, &s, &x, &y, &p, &q, &r, &zz, m, want_q, low, igh, Q);
                    didQRStep = true;
                    break;
                case 1:
                    // Found a single root
                    if (want_q)
                        A(en, en) = x + t;
                    wr[en] = x + t;
                    wi[en] = 0.0;
                    en -= 1;
                    break;
                case 2:
                    // We have found a double root, so we need to now
                    // determine if it's a complex or real pair then modify A and Q
                    // if Q is desired. (Note we only compute the Schur 
                    // form of A if Q is desired. This may be modifiable if 
                    // desired)
                    p = ( y - x ) / 2.0;
                    q = p * p + w;
                    zz = sqrt(fabs(q));
                    if (want_q) {
                        A(en,en) = x + t;
                        A(en - 1, en - 1) = y + t;
                    }
                    x = x + t;
                    if ( q < 0.0 ) {
                        // This means we found a complex pair
                        wr[en - 1] = x + p;
                        wr[en] = x + p;
                        wi[en - 1] = zz;
                        wi[en] = -zz;
                    } else {
                        // Otherwise we have a real pair
                        if (p >= 0.0) 
                            zz = p + zz;
                        else
                            zz = p - zz;
                        wr[en - 1] = x + zz;
                        wr[en] = (zz == 0.0) ? (wr[en - 1]) : (x - w / zz);
                        wi[en - 1] = 0.0;
                        wi[en] = 0.0;
                        if (want_q) {
                            x = A(en, en - 1);
                            s = fabs(x) + fabs(zz);
                            p = x / s;
                            q = zz / s;
                            r = sqrt(p*p + q*q);
                            p /= r;
                            q /= r;
                            // Row modifications
                            for (j = en - 1; j < n; j++) {
                                zz = A(en - 1, j);
                                A(en - 1, j) = q * zz + p * A(en, j);
                                A(en, j) = q * A(en, j) - p * zz;
                            } 
                            // Column modifications
                            for ( i = 0; i <= en; i++) {
                                zz = A(i, en - 1);
                                A(i, en - 1) = q * zz + p * A(i, en);
                                A(i, en) = q * A(i, en) - p * zz;
                            }
                            // Accumulate transformations
                            for ( i = low; i <= igh; i++) {
                                zz = Q(i, en - 1);
                                Q(i, en - 1) = q * zz + p * Q(i, en);
                                Q(i, en) = q * Q(i, en) - p * zz;
                            }
                        }
                    }
                    en -= 2;
                    break;
                case 3:
                    // This means we failed to converge
                    return 1;
            }
        }
        return 0;
    }

} // lapack

#endif // TLAPACK_HQR_HH
