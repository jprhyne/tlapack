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

// Auxiliary routines

using namespace tlapack;

TEMPLATE_TEST_CASE("schur form is backwards stable", "[hqr][schur]", TLAPACK_TYPES_TO_TEST)
{
    srand(1);

    using matrix_t = TestType;
    using T = type_t<matrix_t>;
    using idx_t = size_type<matrix_t>;
    using range = std::pair<idx_t, idx_t>;
    typedef real_type<T> real_t;

    // Function
    Create<matrix_t> new_matrix;

    const T zero(0);
    const T one(1);

    idx_t m, n;

    n = GENERATE(10, 15, 20, 30);
    /*
    // First create a matrix A as upper hessenberg matrix.
    double *A = (double *) calloc(n*n,sizeof(double));
    double *T = (double *) calloc(n*n,sizeof(double));
    // These values are not needed for this particular file's
    // testing, however this is needed to run hqr.c
    double *eigValsReal = (double *) malloc(n*n*sizeof(double));
    double *eigValsImag = (double *) malloc(n*n*sizeof(double));
	// Generate A as a random matrix.
    for(int i = 0; i < n; i++) {
        int start = 0;
        if (i - 1 > 0)
            start = i - 1;
 	    for(int j = start; j < n; j++) {
            double val = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
	        a0(i,j) = val; 
	        t0(i,j) = val; 
        }
    }
    // Create a matrix to store the Schur Vectors
    double *schurMat = (double *) calloc(n*n,sizeof(double));
    for (int i = 0; i < n; i++)
        schurMat[i + i * n] = 1;
    // Now we call hqr. At the end schurMat will contain the schur vectors
    double norm = hqr(n,n,0,n-1,T,eigValsReal,eigValsImag,1,schurMat);
    if (norm < 0) {
        // This means that hqr did not converge to at some index,
        // so we print it out and terminate execution as our Schur
        // vectors will not be correct
        printf("Did not converge at index: %e\n",-norm);
        return 1;
    }
    // Getting here means that we have successfully ran all of 
    // hqr and got an answer, so now we check if our Schur vectors are correct
    //  check || Z' * Z - I ||_F
    double orthZ, tmp;
    orthZ = 0e+00;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            tmp = ( i == j ) ? 1.0e+00 : 0.00e+00;
            for (int k = 0; k < n; k++) {
                tmp -= schurMat0(k,i)*schurMat0(k,j);
            }
            orthZ += tmp * tmp;
        }
    }
    orthZ = sqrt( orthZ );

    // Zero out below quasi diagonal elements of T
    // First, zero out everything below the 1st subdiagonal
    for (int i = 0; i < n; i++) 
        for (int j = 0; j < i - 1; j++) 
            t0(i,j) = 0;
    // if eigValsImag[k]  = 0 then the sub diagonal elements need to be 0
    // If eigValsImag[k] != 0 then we have a schur block  
    int k;
    for (k = 0; k < n-1; k++) {
        if (eigValsImag[k] == 0) {
            t0(k+1,k) = 0;
        } else if (k < n-2){
            // This means we are in a schur block, so the next sub diagonal
            // element must be 0
            t0(k+2,k+1) = 0;
            k++;
        }
    }

    //  check || A - Z * T * Z^T ||_F / || A ||_F
    double normR, normA;
    normR = 0.0e+00;
    double *Zt = malloc(n * n * sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Zt[i + j * n] = schurMat0(j,i);
    // Yes, this is lazy and inefficient, however it 
    // accomplishes our goals in a reasonable time frame
    double *ZT = matmul(schurMat,n,n,T,n,n);
    double *rhs = matmul(ZT,n,n,Zt,n,n);
    double *ans = matsub(rhs,n,n,A,n,n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            normR += ans[i + j * n] * ans[i + j * n];
    normR = sqrt( normR );
    normA = 0.0e+00;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            normA += a0(i,j) * a0(i,j) ;
        }
    }
    normA = sqrt( normA );
    if (testFlag)
        printf("%% [ hqr2schur C ] n = %8d; seed = %8d; checks = [ %8.2e %8.2e ];\n", n, seed, orthZ, normR/normA);
    else
        printf( "%8d %8d %6.1e %6.1e\n", n, seed, orthZ, normR/normA);
    free(A);
    free(T);
    free(schurMat);
    free(Zt);
    free(ZT);
    free(rhs);
    free(ans);
    return 0;
    */    
}
