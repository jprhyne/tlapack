/// @file test_hqr_schurToEigen.cpp
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

#include <tlapack/lapack/hqr.hpp>
#include <tlapack/lapack/hqr_schurToEigen.hpp>

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

    idx_t n;

    n = GENERATE(5, 10, 30, 50, 100, 125, 150, 250,  300, 400, 500);
    const real_t eps = uroundoff<real_t>(); 
    const T tol = T( 10 * n *  eps);

    const T one = T(1);
    const T zeroT = T(0);

    // Function
    Create<matrix_t> new_matrix;


    // Create matrices
    std::vector<T> A_; auto A = new_matrix( A_, n, n);
    std::vector<T> U_; auto U = new_matrix( U_, n, n);
    std::vector<T> Z_; auto Z = new_matrix( Z_, n, n);
    std::vector<T> wr(n);
    std::vector<T> wi(n);

    // Generate A as a random symmetric matrix. This is to ensure it is diagonalizable
    for (idx_t i = 0; i < n; i++) {
        for (idx_t j = i; j < n; j++) {
           T val = rand_helper<T>(gen); 
           A(i,j) = val;
           if (i != j) 
               A(j,i) = val; // ensuring symmetry
        }
    }
    // Perform Hessenberg reduction on A.
    std::vector<T> tau(n);
    gehrd(0,n - 1, A, tau);

    // zero out the parts of A that represent Q
    for (idx_t i = 0; i < n; i++) 
        if (i != 0)
            for (idx_t j = 0; j < i - 1; j++) 
                A(i,j) = zeroT;
    // Copy A into U
    for (idx_t i = 0; i < n; i++)
        for (idx_t j = 0; j < n; j++)
            U(i,j) = A(i,j);

    // Start Z as eye(n)
    for (idx_t i = 0; i < n; i++)
        Z(i,i) = one;

    //Call hqr
    real_t norm = real_t(0.0);
    int retCode = hqr(U, 0, n - 1, wr, wi, true, Z, norm);
    CHECK(retCode == 0);

    // Getting here means that we have successfully ran all of 
    retCode = hqr_schurToEigen(U, 0, n - 1, wr, wi, Z, norm);
    CHECK(retCode == 0);
    // Zero out below quasi diagonal elements of T
    // First, zero out everything below the 1st subdiagonal
    for (idx_t i = 0; i < n; i++) 
        if (i != 0)
            for (idx_t j = 0; j < i - 1; j++) 
                U(i,j) = zeroT;
    // if wi[k]  = 0 then the sub diagonal elements need to be 0
    // If wi[k] != 0 then we have a schur block  
    idx_t k;
    for (k = 0; k < n-1; k++) {
        if (wi[k] == real_t(0)) {
            U(k+1,k) = zeroT;
        } else if (k < n-2){
            // This means we are in a schur block, so the next sub diagonal
            // element must be 0
            U(k+2,k+1) = zeroT;
            k++;
        }
    }
    /*
    std::vector<T> Zi_; auto Zi = new_matrix( Zi_, n, n);
    lacpy(Uplo::General, Zi, Z);
    std::vector<T> Piv(n);
    // Perform LU Decomp of Zi
    idx_t retCode = getrf(Zi,Piv);
    CHECK(retCode == 0);
    // Now compute the inverse of Zi
    idx_t getri(Zi,Piv);
    CHECK(retCode == 0);
    */
    std::vector<T> diffReal_; auto diffReal = new_matrix( diffReal_, n, n);
    std::vector<T> diffImag_; auto diffImag = new_matrix( diffImag_, n, n);
    real_t zero = real_t(0);
    for (idx_t k = 0; k < n; k++) {
        // Redeclare on each loop to not have to reset each value to 0
        // Vectors that will store Av = Ax + i * Ay
        std::vector<T> Ax(n);
        std::vector<T> Ay(n);
        // Vectors that will store \lambda v = (a + ib)(x+iy) = (ax - by) + i(ay + bx)
        std::vector<T> ax(n);
        std::vector<T> ay(n);
        // Check if the current eigenvalue is real or imaginary.
        // If it is real, we need only check Ax = ax
        if (wi[k] == zero) {
            // Real eigenvalue (eigenValsReal[k]) with eigenvector eigenMat(:,k)
            for (idx_t i = 0; i < n; i++) {
                for (idx_t j = 0; j < n; j++) {
                    Ax[i] += A(i,j) * Z(j,k);
                }
                ax[i] = wr[k] * Z(i,k);
            }
        } else if (wi[k] > real_t(0.0)) {
            //continue;
            // Imaginary with positive imaginary part so has associated eigenvector
            // eigenMat(:,k) + i * eigenMat(:,k+1)
            // compute Ax,Ay
            for (idx_t i = 0; i < n; i++) {
                for (idx_t j = 0; j < n; j++) {
                    Ax[i] += A(i,j) * Z(j,k);
                    Ay[i] += A(i,j) * Z(j, k + 1);
                }
                ax[i] = wr[k] * Z(i,k) - wi[k] * Z(i,k+1);
                ay[i] = wr[k] * Z(i,k + 1) + wi[k] * Z(i,k);
            }
        } else {
            //continue;
            // Imaginary with negative imaginary part so has associated eigenvector
            // eigenMat(:,k-1) - i * eigenMat(:,k)
            for (idx_t i = 0; i < n; i++) {
                for (idx_t j = 0; j < n; j++) {
                    Ax[i] += A(i,j) * Z(j,k - 1);
                    Ay[i] -= A(i,j) * Z(j, k);
                }
                ax[i] = wr[k] * Z(i, k - 1) + wi[k] * Z(i, k);
                ay[i] = wi[k] * Z(i, k - 1) - wr[k] * Z(i, k);
            }
        } 
        // Now we construct the kth column of diffReal and diffImag. Note that if we had a real eigenvalue the
        // kth column of diffImag will be exactly 0
        for (idx_t i = 0; i < n; i++) {
            // We want the kth column of diffReal to be Ax - ax
            // and the kth column of diffImag to be Ay - ay
            diffReal(i,k) = Ax[i] - ax[i];
            diffImag(i,k) = Ay[i] - ay[i];
        }
    }
    // Now we have computed AV - VD and stored the real and imaginary parts separately
    // Since we want to get a relative error we must compute \|AV - VD\| in order
    // for this to be as close to accurate as possible we will compute the sum 
    // of the absolute elements of AV - VD using the fact that |a + bi| = sqrt(a^2 + b^2)
    real_t normR = real_t(0.0);
    real_t normA = real_t(0.0);
    for (idx_t i = 0; i < n; i++) {
        for (idx_t j = 0; j < n; j++){
            real_t tmpA = diffReal(i,j);
            real_t tmpB = diffImag(i,j);
            normR += sqrt(tmpA * tmpA + tmpB * tmpB); 
            normA += tlapack::abs(A(i,j));
        }
    }
    CHECK(normR <= tol * normA);
}
