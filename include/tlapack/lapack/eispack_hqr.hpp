/// @file eispack_hqr.hpp
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

#include "tlapack/lapack/eispack_subDiagonalSearch.hpp"
#include "tlapack/lapack/eispack_formShift.hpp"
#include "tlapack/lapack/eispack_doubleSubDiagonalSearch.hpp"
#include "tlapack/lapack/eispack_qrIteration.hpp"


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
    template <class matrix_t, 
             class vector_t,
             enable_if_t<
                 is_complex<type_t<vector_t>>::value, 
                int> = 0
                >
    int eispack_hqr(
        matrix_t &A,
        size_type<matrix_t> low,
        size_type<matrix_t> igh,
        vector_t &eigs,
        bool want_T,
        bool want_Q,
        matrix_t &Q,
        real_type<type_t<matrix_t>> &norm )
    {
        using TA = type_t<matrix_t>;
        using idx_t = size_type<matrix_t>;
        using real_t = real_type<TA>; 
        using complex_t = type_t<vector_t>;
        // Constants
        real_t zero = real_t(0);

        // Grab the number of columns of A, we only work on square matrices
        const idx_t n = ncols(A);

        // Perform the checks for our arguments
        tlapack_check(n == nrows(A));
        tlapack_check((idx_t)size(eigs) == n);

        if (want_Q) {
            // If we want the Schur Vectors, we need to make sure that Q is the right size
            tlapack_check((n == ncols(Q)));
        }

        // Consider adding some checks for trivial cases if needed

        // Now, we actually start porting the behavior

        // Initialize some variables at the beginning.
        norm = zero; 
        idx_t i,j,en,m,itn,its,l;
        real_t p,q,r,s,t,w,x,y,z,zz;
        // Construct the sum of the absolute elements of A and check for isolated
        // eigenvalues
        for (i = 0; i < n; i++) {
            for (j = (i == 0) ? (i):(i-1); j < n; j++) {
                norm += tlapack::abs(A(i,j));
            }
            if (i >= low && i <= igh)
                continue;
            eigs[i] = complex_t(A(i,i), zero);
        }

        // Initialize some more variables
        // en denotes the index of wr/wi we are trying to find. We start by
        //      trying to find the eigenvalues from the bottom to the top
        // t is a shifting variable that I am not too sure the theory behind
        //      why it works
        // itn is the total number of QR iterations we allow for the entire matrix
        //      The heuristic chosen was 30 * n
        en = igh;
        t = zero;
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
            l = eispack_hqr_subDiagonalSearch(low,A,en,norm,s);
            // Check to see if we found eigenvalues based on l from the subDiagonalSearch
            // if not, perform our form shift
            x = A(en, en);
            if (l == en) {
                // This means we found a single root
                if (want_T)
                    A(en, en) = x + t;
                eigs[en] = complex_t(x+t, zero);
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
                p = ( y - x ) / real_t(2.0);
                q = p * p + w;
                zz = sqrt(tlapack::abs(q));
                if (want_T) {
                    A(en,en) = x + t;
                    A(en - 1, en - 1) = y + t;
                }
                x = x + t;
                if ( q < zero ) {
                    // This means we found a complex pair
                    eigs[en - 1] = complex_t(x + p, zz);
                    eigs[en] = complex_t(x + p, -zz);
                } else {
                    // Otherwise we have a real pair
                    if (p >= zero) 
                        zz = p + tlapack::abs(zz);
                    else
                        zz = p - tlapack::abs(zz);
                    eigs[en - 1] = complex_t(x+zz, zero);
                    if (zz == zero)
                        eigs[en] = eigs[en - 1];
                    else
                        eigs[en] = complex_t(x - w/zz, zero);
                        
                    if (want_T || want_Q) {
                        x = A(en, en - 1);
                        s = tlapack::abs(x) + tlapack::abs(zz);
                        p = x / s;
                        q = zz / s;
                        r = sqrt(p*p + q*q);
                        p = p / r;
                        q = q / r;
                        if (want_T) {
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
                        }
                        if (want_Q) {
                            // Accumulate transformations
                            for ( i = low; i <= igh; i++) {
                                zz = Q(i, en - 1);
                                Q(i, en - 1) = q * zz + p * Q(i, en);
                                Q(i, en) = q * Q(i, en) - p * zz;
                            }
                        }
                    }
                }
                en -= 2;
                foundEigenValue = true;
            } 
            if (!foundEigenValue) {
                // This is performing an "exceptional"
                // shift, so we sometimes perform a single
                // shift if we are exceptionally unlucky. We do an extra single shift.
                if (its % 10 == 0 && its > 0)
                    eispack_hqr_formShift(low, A, en, s, t, x, y, w);
                // We only do a QR step if we did not find an eigenvalue
                its += 1;
                itn -= 1;
                m = eispack_hqr_doubleSubDiagonalSearch(A, en, l, s, x, y, w, p, q, r, zz);
                // double qr step
                eispack_hqr_qrIteration(A, en, l, s, x, y, p, q, r, zz, m, want_T, want_Q, low, igh, Q);
                didQRStep = true;    
            }

        }
        return 0;
    }

    /*
     * Note: We do not actually use the norm value here,
     * however we keep it in the interface to make it the same as above
     */
    template <class matrix_t, 
             class vector_t,
             enable_if_t<
                 is_complex<type_t<vector_t>>::value, 
                int> = 0
                >
    int eispack_comqr(
        matrix_t &A,
        size_type<matrix_t> low,
        size_type<matrix_t> igh,
        vector_t &eigs,
        bool want_T,
        bool want_Q,
        matrix_t &Q,
        real_type<type_t<matrix_t>> &norm )
    {
        using TA = type_t<matrix_t>;
        using idx_t = size_type<matrix_t>;
        using real_t = real_type<TA>; 
        using complex_t = type_t<vector_t>;
        // Grab the number of columns of A, we only work on square matrices
        const idx_t n = ncols(A);

        // Perform the checks for our arguments
        tlapack_check(n == nrows(A));
        tlapack_check((idx_t)size(eigs) == n);

        if (want_Q) {
            // If we want the Schur Vectors, we need to make sure that Q is the right size
            tlapack_check((n == ncols(Q)));
        }

        // Consider adding some checks for trivial cases if needed

        // Now, we actually start porting the behavior
        // Constructing real sub diagonal elements
        for (idx_t i = low + 1; i <= igh; i++) {
            idx_t l = (i + 1 < igh) ? (i + 1) : (igh);
            if (A(i,i-1).imag() == 0)
                continue;
            norm = tlapack::abs(A(i,i-1));
            complex_t y = A(i, i - 1)/norm;
            A(i, i - 1) = complex_t(norm, 0);
            idx_t upperBound = (want_Q) ? (n - 1) : (igh);
            for (idx_t j = i; j <= upperBound; j++)
                A(i,j) = A(i,j) * conj(y);
            idx_t lowerBound = (want_Q) ? (0) : (low);
            for (idx_t j = lowerBound; j <= l; j++)
                A(j,i) = A(j,i) * y;
            if (want_Q) 
                for (idx_t j = low; j <= igh; j++)
                    Q(j,i) = Q(j,i) * y;
        }
        // Store the isolated roots
        for (idx_t i = 0; i <= n - 1; i++)
            if (i < low || i > igh)
                eigs[i] = A(i,i);

        complex_t s,x,y,zz;
        idx_t l;
        idx_t en = igh;
        complex_t t = 0;
        idx_t itn = 30 * n;
        idx_t its = 0;
        bool didQRStep = false;
        while (en >= low && en <= igh) {
            if (!didQRStep)
                its = 0;
            didQRStep = false;

            // We want to check here if we have exhausted our iterations ie if itn == 0. If so, we immediately
            // terminate and return the index of failure
            if (itn == 0)
                return en;

            // Perform a subdiagonal search to determine where the first
            // small subdiagonal element is
            l = eispack_comqr_subDiagonalSearch(low,A,en);
            if (l == en) {
                // Means we have found a root
                eigs[en] = A(en,en) + t;
                if (want_T)
                    A(en,en) = eigs[en];
                en--;
                continue;
            }
            // Compute the shift and do it on the diagonal of A. In addition, accumulate that shift
            eispack_comqr_formShift(low, A, its, en, s, t, x, y, zz);
            // Update the number of iterations we have done so far
            its++;
            itn--;
            eispack_comqr_qrIteration(A,en,l,s,x,y,zz,eigs,want_T, want_Q,low,igh,Q);
            didQRStep = true;
        }
        return 0;
    }

} // lapack

#endif // TLAPACK_HQR_HH
