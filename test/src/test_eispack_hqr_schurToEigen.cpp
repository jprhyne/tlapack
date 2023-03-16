/// @file test_eispack_hqr_schurToEigen.cpp
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
#include <tlapack/lapack/lacpy.hpp>
#include <tlapack/lapack/lange.hpp>
#include <tlapack/lapack/gehrd.hpp>
#include <tlapack/lapack/unghr.hpp>
#include <tlapack/lapack/getri.hpp>
#include <tlapack/lapack/getrf.hpp>
#include <tlapack/blas/gemm.hpp>

#include <tlapack/lapack/eispack_hqr.hpp>
#include <tlapack/lapack/eispack_hqr_schurToEigen.hpp>
#include <tlapack/lapack/multishift_qr.hpp>

// Auxiliary routines

using namespace tlapack;

TEMPLATE_TEST_CASE("schur form is backwards stable", "[hqr][schur]", TLAPACK_REAL_TYPES_TO_TEST)
{
    srand(1);
    rand_generator gen;
    

    using matrix_t = TestType;
    using T = type_t<matrix_t>;
    using idx_t = size_type<matrix_t>;
    typedef real_type<T> real_t;
    // Complex numbers are constructed as complex_t(real,imaginary)
    typedef complex_type<T> complex_t; 

    idx_t n;
    // Need to review how to do this
    //auto matrix_type = GENERATE(as<std::string>{}, "Full Matrix", "Inner Window");
    auto matrix_type = GENERATE(as<std::string>{}, "Full Matrix");
    auto schur_type = GENERATE(as<std::string>{}, "Thijs");
    INFO("schur_type=" << schur_type);
    idx_t seed = GENERATE(123,623,134,5); // Numbers generated by my shell's random
    n = GENERATE(1,2,3,4,5, 10, 30, 32, 33, 41);
    gen.seed(seed);
    INFO("n=" << n);
    const real_t eps = uroundoff<real_t>(); 
    const T tol = T( 10 * n *  eps);

    const T one = T(1);
    const T zero = T(0);

    // Function
    Create<matrix_t> new_matrix; // For Real matrices
    // There may be an easier way of writing this, however it follows the
    // structure of a matrix typing.
    Create<legacyMatrix<std::complex<real_t>,std::size_t,Layout::ColMajor>> new_matrixC; // For Complex matrices


    // Create matrices
    std::vector<T> A_; auto A = new_matrix( A_, n, n);
    std::vector<T> H_; auto H = new_matrix( H_, n, n);
    std::vector<T> Z_; auto Z = new_matrix( Z_, n, n);
    std::vector<complex_t> s(n);
    for (idx_t i = 0; i < n; i++) {
        s[i] = complex_t(rand_helper<real_t>(gen), rand_helper<real_t>(gen));
    }

    // Generate our matrix A as a full matrix and then reduce it to 
    // hessenberg form
    for (idx_t i = 0; i < n; i++) {
        for (idx_t j = 0; j < n; j++) { 
            T val = rand_helper<T>(gen);
            A(i,j) = val; 
            H(i,j) = val;
        }
    }
    // Perform Hessenberg reduction on A.
    std::vector<T> tau(n);
    gehrd(0,n, H, tau);
    lacpy(Uplo::General, H, Z);
    unghr(0, n, Z, tau);

    // zero out the parts of H that represent Z
    // IE the 'reflectors'
    // After running tests, this is necessary for our function to run.
    for (idx_t i = 1; i < n; i++)
        for (idx_t j = 0; j < i - 1; j++) 
            H(i,j) = zero;
    real_t norm = 0;
    idx_t retCode;
    if (schur_type == "Thijs") {
        //Call Thijs' function to ensure we
        //work for other schur form computations
        idx_t ns = 4;
        idx_t nw = 4;
        francis_opts_t<idx_t> opts;
        opts.nshift_recommender = [ns](idx_t n, idx_t nh) -> idx_t 
        {
            return ns;
        };
        opts.deflation_window_recommender = [nw](idx_t n, idx_t nh) -> idx_t
        {
            return nw;
        };
        opts.nmin = 15;

        // Compute the norm of H as this is needed for schurToEigen
        for (idx_t i = 0; i < n; i++) 
            for (idx_t j = 0; j < n; j++) 
                norm += tlapack::abs(H(i,j));
        retCode = multishift_qr(true, true, 0, n, H, s, Z, opts);
        // Ensure we actually finished hqr.
        CHECK(retCode == 0);
    } else {
        //Call hqr
        retCode = eispack_hqr(H, 0, n - 1, s, true, Z, norm);
        //CHECK(retCode == 0);
    }
    // TODO: Add a way to determine if we are testing real or complex numbers then 
    // call hqr and comqr respectively
    // Getting here means that we have successfully ran all of hqr
    retCode = eispack_hqr_schurToEigen(H, 0, n - 1, s, Z, norm);
    CHECK(retCode == 0);
    // TODO: Add a way to determine if we are testing real or complex numbers then 
    // if complex always do #1 and construct our matrices more smartly
    // We have two kinds of testing schemes.
    // 1) Complex numbers. This is more of a legacy check and is kept in case our real tests break.
    //      Can remove at some point if we need to
    // 2) Real numbers. This way treats our matrices as their real meanings. This is because
    //      even though we get complex eigenvalues, the real and imaginary parts represent a 
    //      kind of invariant subspace.
    auto testingScheme = "Complex Multiplication";
    // If wanting to test with complex numbers
    if (testingScheme == "Complex Multiplication") {
        // Now, currently we are only testing matrices that are supposed to be diagonalizable
        // We will test representativity by constructing a matrix Zc such that Zc contains 
        // the eigenvectors stored in Z but as complex numbers
        // Similarly, do so for Dc containing the eigenvalues stored as complex numbers
        std::vector<complex_t> Zc_; auto Zc = new_matrixC( Zc_, n, n);
        for (idx_t j = 0; j < n; j++) { // For each column
            for (idx_t i = 0; i < n; i++) { // Grab the ith row
                // For more information on how the eigenvectors are constructed
                // see the documentation for hqr_schurToEigen
                // If wi[j] is 0, then we have a real eigenvector, so we only need to copy the current column
                if (s[j].imag() == zero)
                    Zc(i,j) = complex_t(Z(i,j), 0);
                // If wi[j] is positive, then we have an eigenvector of the form Z[:,j] + Z[:,j+1]*i
                else if (s[j].imag() > zero)
                    Zc(i,j) = complex_t(Z(i,j), Z(i,j + 1));
                // Otherwise, we found the conjugate pair so we need Z[:,j-1] - Z[:,j]*i
                else
                    Zc(i,j) = complex_t(Z(i,j - 1), -Z(i,j));
            }
        }
        std::vector<complex_t> Zi_; auto Zi = new_matrixC( Zi_, n, n);
        lacpy(Uplo::General, Zc, Zi);
        std::vector<T> Piv(n);
        // Perform LU Decomp of Zi
        retCode = getrf(Zi,Piv);
        CHECK(retCode == 0); // Ensure we properly computed the LU
        // Now compute the inverse of Zi
        retCode = getri(Zi,Piv);
        CHECK(retCode == 0); // Ensure we properly computed the Inverse
        // Now test VDV^{-1} - A
        std::vector<complex_t> lhs_; auto lhs = new_matrixC( lhs_, n, n);
        std::vector<complex_t> Dc_; auto Dc = new_matrixC( Dc_, n, n);
        for (idx_t i = 0; i < n; i++)
            Dc(i,i) = s[i];
        // We need to also construct Ac which is just a complex equivalent of A
        // with 0 for all imaginary parts
        std::vector<complex_t> Ac_; auto Ac = new_matrixC( Ac_, n, n);
        for (idx_t i = 0; i < n; i++)
            for (idx_t j = 0; j < n; j++)
                Ac(i,j) = complex_t(A(i,j),0);
    
        gemm(Op::NoTrans, Op::NoTrans, real_t(1), Zc, Dc, lhs);
        // Note: This overwrites A, but we already have the norm of A saved from hqr
        gemm(Op::NoTrans, Op::NoTrans, real_t(1), lhs, Zi, real_t(-1), Ac);
        // Compute the frobenius norm of the residual
        real_t normR = lange(tlapack::frob_norm, Ac);
        real_t normZ = lange(tlapack::frob_norm, Zc);
        real_t normZi = lange(tlapack::frob_norm, Zi);
        real_t normD = lange(tlapack::frob_norm, Dc);
        CHECK(normR <= tol * normZ * normZi * normD);
    } else {
        // This means our testing Scheme is going to be using only real arithmetic
        // So we will have to construct our diagonal matrix much more carefully
        std::vector<T> D_; auto D = new_matrix(D_, n, n);
        for (idx_t j = 0; j < n; j++) {
            // From here we need to determine if our eigenvalue is real or complex
            // if it is complex, then we must create a 2x2 block of the form
            // ( wr[j] wi[j]
            //  -wi[j] wr[j] )
            // Then skip over the next eigenvalue as wr and wi contain the conjugate pairs
            if (s[j].imag() != zero) {
                D(j,j) = s[j].real();
                D(j + 1,j + 1) = s[j].real();
                D(j, j + 1) = s[j].imag();
                D(j + 1, j) = -s[j].imag();
                j++;
            } else {
                // Otherwise we just put the real eigenvalue on the diagonal
                D(j,j) = s[j].real();
            }
        }
        // Now we can do basically what is above to test if we have the correct form
        // Now, construct Zi 
        std::vector<T> Zi_; auto Zi = new_matrix( Zi_, n, n);
        lacpy(Uplo::General, Z, Zi);
        std::vector<T> Piv(n);
        // Perform LU Decomp of Zi
        retCode = getrf(Zi,Piv);
        CHECK(retCode == 0); // Ensure we properly computed the LU
        // Now compute the inverse of Zi
        retCode = getri(Zi,Piv);
        CHECK(retCode == 0); // Ensure we properly computed the Inverse
        // This will store the final result for our first representativity check
        std::vector<T> firstCheck_; auto firstCheck = new_matrix( firstCheck_, n, n);
        // This will store the final result for our second representativity check
        std::vector<T> secondCheck_; auto secondCheck = new_matrix( secondCheck_, n, n);
        std::vector<T> tmp_; auto tmp = new_matrix( tmp_, n, n);


        // Now test VDV^{-1} - A
        gemm(Op::NoTrans, Op::NoTrans, real_t(1), Z, D, firstCheck);
        // We store this inside secondCheck mostly to prevent us from having to compute it again
        lacpy(Uplo::General, firstCheck, secondCheck);
        lacpy(Uplo::General, A, tmp);
        gemm(Op::NoTrans, Op::NoTrans, real_t(1), firstCheck, Zi, real_t(-1), tmp);
        // Compute the frobenius norm of the residual
        real_t normR = lange(tlapack::frob_norm, tmp);
        CHECK(normR <= tol * norm);

        // Now, we check VD = AV ie VD - AV
        // We already have VD from above, we create another temporary place to store AV
        gemm(Op::NoTrans, Op::NoTrans, real_t(1), A, Z, real_t(-1), secondCheck);

        normR = lange(tlapack::frob_norm, secondCheck);
        CHECK(normR <= tol * norm);

    }
}
