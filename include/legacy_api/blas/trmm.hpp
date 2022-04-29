// Copyright (c) 2017-2021, University of Tennessee. All rights reserved.
// Copyright (c) 2021-2022, University of Colorado Denver. All rights reserved.
//
// This file is part of <T>LAPACK.
// <T>LAPACK is free software: you can redistribute it and/or modify it under
// the terms of the BSD 3-Clause license. See the accompanying LICENSE file.

#ifndef __TLAPACK_LEGACY_TRMM_HH__
#define __TLAPACK_LEGACY_TRMM_HH__

#include "legacy_api/base/utils.hpp"
#include "legacy_api/base/types.hpp"
#include "blas/trmm.hpp"

namespace tlapack {

/**
 * Triangular matrix-matrix multiply:
 * \[
 *     B = \alpha op(A) B,
 * \]
 * or
 * \[
 *     B = \alpha B op(A),
 * \]
 * where $op(A)$ is one of
 *     $op(A) = A$,
 *     $op(A) = A^T$, or
 *     $op(A) = A^H$,
 * B is an m-by-n matrix, and A is an m-by-m or n-by-n, unit or non-unit,
 * upper or lower triangular matrix.
 *
 * Generic implementation for arbitrary data types.
 *
 * @param[in] layout
 *     Matrix storage, Layout::ColMajor or Layout::RowMajor.
 *
 * @param[in] side
 *     Whether $op(A)$ is on the left or right of B:
 *     - Side::Left:  $B = \alpha op(A) B$.
 *     - Side::Right: $B = \alpha B op(A)$.
 *
 * @param[in] uplo
 *     What part of the matrix A is referenced,
 *     the opposite triangle being assumed to be zero:
 *     - Uplo::Lower: A is lower triangular.
 *     - Uplo::Upper: A is upper triangular.
 *     - Uplo::General is illegal (see @ref gemm instead).
 *
 * @param[in] trans
 *     The form of $op(A)$:
 *     - Op::NoTrans:   $op(A) = A$.
 *     - Op::Trans:     $op(A) = A^T$.
 *     - Op::ConjTrans: $op(A) = A^H$.
 *
 * @param[in] diag
 *     Whether A has a unit or non-unit diagonal:
 *     - Diag::Unit:    A is assumed to be unit triangular.
 *     - Diag::NonUnit: A is not assumed to be unit triangular.
 *
 * @param[in] m
 *     Number of rows of matrix B. m >= 0.
 *
 * @param[in] n
 *     Number of columns of matrix B. n >= 0.
 *
 * @param[in] alpha
 *     Scalar alpha. If alpha is zero, A is not accessed.
 *
 * @param[in] A
 *     - If side = Left:
 *       the m-by-m matrix A, stored in an lda-by-m array [RowMajor: m-by-lda].
 *     - If side = Right:
 *       the n-by-n matrix A, stored in an lda-by-n array [RowMajor: n-by-lda].
 *
 * @param[in] lda
 *     Leading dimension of A.
 *     - If side = left:  lda >= max(1, m).
 *     - If side = right: lda >= max(1, n).
 *
 * @param[in,out] B
 *     The m-by-n matrix B, stored in an ldb-by-n array [RowMajor: m-by-ldb].
 *
 * @param[in] ldb
 *     Leading dimension of B. ldb >= max(1, m) [RowMajor: ldb >= max(1, n)].
 *
 * @ingroup trmm
 */
template< typename TA, typename TB >
void trmm(
    Layout layout,
    Side side,
    Uplo uplo,
    Op trans,
    Diag diag,
    idx_t m,
    idx_t n,
    scalar_type<TA, TB> alpha,
    TA const *A, idx_t lda,
    TB       *B, idx_t ldb )
{
    using internal::colmajor_matrix;

    // check arguments
    tblas_error_if( layout != Layout::ColMajor &&
                   layout != Layout::RowMajor );
    tblas_error_if( side != Side::Left &&
                   side != Side::Right );
    tblas_error_if( uplo != Uplo::Lower &&
                   uplo != Uplo::Upper );
    tblas_error_if( trans != Op::NoTrans &&
                   trans != Op::Trans &&
                   trans != Op::ConjTrans );
    tblas_error_if( diag != Diag::NonUnit &&
                   diag != Diag::Unit );
    tblas_error_if( m < 0 );
    tblas_error_if( n < 0 );
    tblas_error_if( lda < ((side == Side::Left) ? m : n) );
    tblas_error_if( ldb < ((layout == Layout::RowMajor) ? n : m) );

    // quick return
    if (m == 0 || n == 0)
        return;

    // adapt if row major
    if (layout == Layout::RowMajor) {
        side = (side == Side::Left)
            ? Side::Right
            : Side::Left;
        if (uplo == Uplo::Lower)
            uplo = Uplo::Upper;
        else if (uplo == Uplo::Upper)
            uplo = Uplo::Lower;
        std::swap( m , n );
    }

    // Matrix views
    const auto _A = (side == Side::Left)
                  ? colmajor_matrix<TA>( (TA*)A, m, m, lda )
                  : colmajor_matrix<TA>( (TA*)A, n, n, lda );
    auto _B = colmajor_matrix<TB>( B, m, n, ldb );

    trmm( side, uplo, trans, diag, alpha, _A, _B );
}

}  // namespace tlapack

#endif        //  #ifndef __TLAPACK_LEGACY_TRMM_HH__
