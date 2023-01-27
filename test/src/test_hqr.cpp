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

#include <tlapack/lapack/lacpy.hpp>
#include <tlapack/lapack/lange.hpp>

// Auxiliary routines

using namespace tlapack;

TEMPLATE_TEST_CASE("schur form is backwards stable", "[hqr][schur]", TLAPACK_REAL_TYPES_TO_TEST)
{
    srand(1);

    rand_generator gen;
    

    using matrix_t = TestType;
    using T = type_t<matrix_t>;
    using idx_t = size_type<matrix_t>;
    using range = std::pair<idx_t, idx_t>;
    typedef real_type<T> real_t;

    idx_t n;

    n = GENERATE(5, 10, 30, 50, 100, 125, 150, 250,  300, 400, 500);
    const real_t eps = uroundoff<real_t>(); 
    const T tol = T( 10 * n *  eps);

    T zero = T(0);
    T one = T(1);

    // Function
    Create<matrix_t> new_matrix;


    // Create matrices
    std::vector<T> A_; auto A = new_matrix( A_, n, n);
    std::vector<T> U_; auto U = new_matrix( U_, n, n);
    std::vector<T> Q_; auto Q = new_matrix( Q_, n, n);
    std::vector<T> wr(n);
    std::vector<T> wi(n);
    for (idx_t i = 0; i < n; i++) {
        for(idx_t j = 0; j < n; j++) {
            A(i,j) = zero;
            U(i,j) = zero;
        }
    }
    //Populate A and U with random numbers
    for(idx_t i = 0; i < n; i++) {
        idx_t start = (i >= 1) ? (i-1): 0;
        for(int j = start; j < n; j++) {
            T val = rand_helper<T>(gen);
            A(i,j) = val; 
            U(i,j) = val; 
        }
    }

    // Start Q as eye(n)
    for (idx_t i = 0; i < n; i++)
        Q(i,i) = one;

    //Call hqr
    real_t norm = real_t(zero);
    int retCode = tlapack::hqr(U, 0, n - 1, wr, wi, true, Q, norm);
    CHECK(retCode == 0);

    // Getting here means that we have successfully ran all of 
    // hqr and got an answer, so now we check if our Schur vectors are correct
    //  check || Q' * Q - I ||_F
    std::vector<T> res_; auto res = new_matrix( res_, n, n );
    auto orthZ = check_orthogonality(Q,res);
    
    CHECK(orthZ <= tol);
    // Zero out below quasi diagonal elements of T
    // First, zero out everything below the 1st subdiagonal
    for (idx_t i = 0; i < n; i++) 
        if (i != 0)
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
    //  check || A - Q * U * Q^T ||_F / || A ||_F
    auto normA = lange(tlapack::frob_norm, A);
    auto sim_res_norm = check_similarity_transform(A,Q,U);
    CHECK(sim_res_norm <= tol * normA);
}
