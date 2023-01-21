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
    /**
     *  @brief Computes the eigenvalues of the real upper Hessenberg matrix A and optionally the Schur Vectors
     *
     *  hqr takes in a real upper Hessenberg matrix A and computes the eigenvalues of it, and then if desired
     *  we also compute the quasi-Schur form of A ie, we find Q,T s.t. A = Q * T * Q^T.
     *  On termination, A will be destroyed so if further computations are needed based on the
     *  original matrix A, it must be stored beforehand.
     *
     *  @tparam matrix_t the type of matrix of A, should be real
     *  @tparam vector_t the type of vector that we are storing the eigenvalues in
     *
     *  @param[in,out] A
     *      On input, A is our upper Hessenberg data matrix
     *      On output, 
     *          if want_Q, then A is the quasi-Schur (block upper triangular) of the
     *              original data matrix.
     *          if !want_Q, then A is destroyed and has no meaning
     *  @param[in] low
     *      The lower bound of the indices we are considering. A must be upper triangular from 
     *          columns 0:low. If this is unknown then choose 0.
     *  @param[in] igh
     *      the upper bound of the indices we are considering. A must be upper triangular from
     *          columns igh:(ncols(A)-1). If this is unknown, choose ncols(A) - 1.
     *  @param[out] wr
     *      The real parts of the eigenvalues of A
     *  @param[out] wi
     *      The imaginary parts of the eigenvalues of A
     *  @param[in] want_Q
     *      boolean value determining if the Schur vectors (stored in Q) are desired
     *  @param[in,out] Q
     *      On input: Q must either be the identity matrix or a unitary similarity transformation on some
     *          matrix to get A. IE, starting with B, A and Q must satisfy B = Q*A*Q^T like the
     *          Hessenberg reduction
     *      On output: If Q was eye(ncols(A)), then Q contains the Schur vectors of A.
     *          If Q was a similarity transformation, then the new Q will contain the 
     *          Schur vectors of the original B (see on input for what B is)
     *  @param[out] norm
     *      The sum of the absolute elements of the upper hessenberg portion of A used for
     *      computing the eigenvectors of A after computing the Schur vectors
     */
    template <class matrix_t, class vector_t>
    int hqr(
        matrix_t &A,
        size_type<matrix_t> low,
        size_type<matrix_t> igh,
        vector_t &wr,
        vector_t &wi,
        bool want_Q,
        matrix_t &Q,
        real_type<type_t<matrix_t>> &norm )
    {
        using TA = type_t<matrix_t>;
        using idx_t = size_type<matrix_t>;
        using real_t = real_type<TA>; 

        // Grab the number of columns of A, we only work on square matrices
        const idx_t n = ncols(A);

        // Perform the checks for our arguments
        tlapack_check(n == nrows(A));
        tlapack_check((idx_t)size(wr) == n);
        tlapack_check((idx_t)size(wi) == n);

        if (want_Q) {
            // If we want the Schur Vectors, we need to make sure that Q is the right size
            tlapack_check((n == ncols(Q)));
        }

        // Consider adding some checks for trivial cases if needed

        // Now, we actually start porting the behavior

        // Initialize some variables at the beginning. We use the following
        // general conversion schema
        // double   -> real_t
        // int      -> idx_t
        norm = 0.0; 
        idx_t i,j,en,m,itn,its,l;
        real_t p,q,r,s,t,w,x,y,z,zz;
        // Construct the sum of the absolute elements of A and check for isolated
        // eigenvalues
        for (i = 0; i < n; i++) {
            for (j = (i == 0) ? (i):(i-1); j < n; j++) {
                real_t tmp = tlapack::abs(A(i,j));
                norm += tlapack::abs(A(i,j));
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
        while (en >= low && en <= igh) {
            if (!didQRStep)
                its = 0;
            didQRStep = false;

            // We want to check here if we have exhausted our iterations ie if itn == 0. If so, we immediately
            // terminate and return the index of failure
            if (itn == 0)
                return en;

            // Creating a flag to help determine if we want to do a QR Step. We only want to do
            // this if we did not find an eigenvalue.
            bool foundEigenValue = false;

            // Perform a subdiagonal search to determine where the first
            // small subdiagonal element is
            l = hqr_subDiagonalSearch(low,A,en,norm,s);
            // Check to see if we found eigenvalues based on l from the subDiagonalSearch
            // if not, perform our form shift
            x = A(en, en);
            if (l == en) {
                // This means we found a single root
                if (want_Q)
                    A(en, en) = x + t;
                wr[en] = x + t;
                wi[en] = 0.0;
                en -= 1;
                foundEigenValue = true;
            } else {
                // Since in the old behavior this
                // is only done if we have not
                // found a single root
                y = A(en - 1, en - 1);
                w = A(en, en - 1) * A(en - 1, en);
            }
            // The above if block and this one will never execute at the same time.
            // See: if l = en we decrement en so l = en + 1 after the block, which is not en - 1
            // Even if en - 1 underflows we would have l = 0 which is not the maximal value for any 
            // reasonable unsigned data type
            if (l  == en - 1) {
                // This means we found a double root
                // We have found a double root, so we need to now
                // determine if it's a complex or real pair then modify A and Q
                // if Q is desired. (Note we only compute the Schur 
                // form of A if Q is desired. This may be modifiable if 
                // desired)
                p = ( y - x ) / 2.0;
                q = p * p + w;
                zz = sqrt(tlapack::abs(q));
                if (want_Q) {
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
                    if (want_Q) {
                        x = A(en, en - 1);
                        s = tlapack::abs(x) + tlapack::abs(zz);
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
                foundEigenValue = true;
            } 
            if (!foundEigenValue) {
                // We want to potentially shift if we did not find an eigenvalue and our check
                // must be done before the QR Step
                // We can put this check inside of formShift if we want to support shifts on 
                // every iteration or on a different about of iterations
                if (its % 10 == 0)
                    hqr_formShift(low, A, en, s, t, x, y, w);
                // We only do a QR step if we did not find an eigenvalue
                its += 1;
                itn -= 1;
                m = hqr_doubleSubDiagonalSearch(A, en, l, s, x, y, w, p, q, r, zz);
                // double qr step
                hqr_qrIteration(A, en, l, s, x, y, p, q, r, zz, m, want_Q, low, igh, Q);
                didQRStep = true;    
            }

        }
        return 0;
    }

} // lapack

#endif // TLAPACK_HQR_HH
