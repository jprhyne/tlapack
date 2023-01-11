/// @file test_hqr.cpp
/// @brief Test HQR. 
//
// Copyright (c) 2022, University of Colorado Denver. All rights reserved.
//
// This file is part of <T>LAPACK.
// <T>LAPACK is free software: you can redistribute it and/or modify it under
// the terms of the BSD 3-Clause license. See the accompanying LICENSE file.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "testutils.hpp"

#include <tlapack/lapack/hqr.hpp>

// Auxiliary routines

using namespace tlapack;

TEMPLATE_TEST_CASE("schur form is backwards stable", "[hqr][schur]", TLAPACK_TYPES_TO_TEST)
{
    srand(1);

    rand_generator gen;

    using matrix_t = TestType;
    using T = type_t<matrix_t>;
    using idx_t = size_type<matrix_t>;
    using range = std::pair<idx_t, idx_t>;
    typedef real_type<T> real_t;

    const real_t eps = uls<real_t>();
    const real_t tol = real_t(n * eps);

    // Function
    Create<matrix_t> new_matrix;

    const T zero(0);
    const T one(1);

    idx_t n;

    n = GENERATE(10, 15, 20, 30);

    // Create matrices
    std::vector<T> A_; auto A = new_matrix( A_, n, n);
    std::vector<T> U_; auto U = new_matrix( U_, n, n);
    std::Vector<T> Q_; auto Q = new_matrix( Q_, n, n);
    std::vector<T> wr;
    std::vector<T> wi;
    //Populate A and U with random numbers
    for(idx_t i = 0; i < n; i++) {
        idx_t start = 0;
        if (i - 1 > 0)
            start = i - 1;
 	    for(int j = start; j < n; j++) {
            T val = rand_helper<T>(gen)
	        A(i,j) = val; 
	        U(i,j) = val; 
        }
    }

    // Start Q as eye(n)
    for (idx_t i = 0; i < n; i++)
        Q(i,i) = 1;

    //Call hqr
    T norm = 0.0;
    // Do we not need something like <matrix_t, vector_t> to tell
    // the call what we are using? Or does it work automatically?
    // Some testing functions I have seen do not include this
    int retCode = tlapack::hqr(true, 0, n-1, T, wr, wi, Q, &norm)
    if (retCode != 0) {
        // This means that hqr did not converge to at some index,
        // so we print it out and terminate execution as our Schur
        // vectors will not be correct
        printf("Did not converge at index: %e\n",-norm);
        return 1;
    }

    // Getting here means that we have successfully ran all of 
    // hqr and got an answer, so now we check if our Schur vectors are correct
    //  check || Q' * Q - I ||_F
    T orthZ, tmp;
    orthZ = 0e+00;
    for (idx_t i = 0; i < n; i++) {
        for (idx_t j = 0; j < n; j++) {
            tmp = ( i == j ) ? 1.0e+00 : 0.00e+00;
            for (idx_t k = 0; k < n; k++) {
                tmp -= Q(k,i)*Q(k,j);
            }
            orthZ += tmp * tmp;
        }
    }
    orthZ = sqrt( orthZ );
    CHECK(orthZ <= tol);
    // Zero out below quasi diagonal elements of T
    // First, zero out everything below the 1st subdiagonal
    for (idx_t i = 0; i < n; i++) 
        for (idx_t j = 0; j < i - 1; j++) 
            U(i,j) = 0;
    // if eigValsImag[k]  = 0 then the sub diagonal elements need to be 0
    // If eigValsImag[k] != 0 then we have a schur block  
    idx_t k;
    for (k = 0; k < n-1; k++) {
        if (wi[k] == 0) {
            U(k+1,k) = 0;
        } else if (k < n-2){
            // This means we are in a schur block, so the next sub diagonal
            // element must be 0
            U(k+2,k+1) = 0;
            k++;
        }
    }
    //  check || A - Z * T * Z^T ||_F / || A ||_F
    std::Vector<T> Qt_; auto Qt = new_matrix( Qt_, n, n);
    T normR, normA;
    normR = 0.0e+00;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Qt(i,j) = Q(j,i);
    // Not sure how matrix mult is done in this codebase. Wait for meeting to 
    // properly implement. The difference will be stored in a matrix called 
    // ans = A - QUQ'
    for (idx_t i = 0; i < n; i++)
        for (idx_t j = 0; j < n; j++)
            normR += ans(i,j) * ans(i,j);
    normR = sqrt( normR );
    normA = 0.0e+00;
    for (idx_t i = 0; i < n; i++) {
        for (idx_t j = 0; j < n; j++) {
            normA += A(i,j) * A(i,j) ;
        }
    }
    normA = sqrt( normA );

    CHECK(normR <= tol * normA);
}
