// Copyright (c) 2017-2021, University of Tennessee. All rights reserved.
// Copyright (c) 2021-2022, University of Colorado Denver. All rights reserved.
//
// This file is part of <T>LAPACK.
// <T>LAPACK is free software: you can redistribute it and/or modify it under
// the terms of the BSD 3-Clause license. See the accompanying LICENSE file.

#ifndef TLAPACK_BLAS_COPY_HH
#define TLAPACK_BLAS_COPY_HH

#include "tlapack/base/utils.hpp"

namespace tlapack {

/**
 * Copy vector, $y := x$.
 *
 * @param[in]  x A n-element vector x.
 * @param[out] y A vector of at least n elements.
 *
 * @ingroup copy
 */
template< class vectorX_t, class vectorY_t,
    class T = type_t<vectorY_t>,
    disable_if_allow_optblas_t<
        pair< vectorX_t, T >,
        pair< vectorY_t, T >
    > = 0
>
void copy( const vectorX_t& x, vectorY_t& y )
{
    using idx_t = size_type< vectorX_t >;

    // constants
    const idx_t n = size(x);

    // check arguments
    tlapack_check_false( size(y) < n );

    for (idx_t i = 0; i < n; ++i)
        y[i] = x[i];
}

#ifdef USE_LAPACKPP_WRAPPERS

    template< class vectorX_t, class vectorY_t,
        class T = type_t<vectorY_t>,
        enable_if_allow_optblas_t<
            pair< vectorX_t, T >,
            pair< vectorY_t, T >
        > = 0
    >
    inline
    void copy( const vectorX_t& x, vectorY_t& y )
    {
        using idx_t = size_type< vectorX_t >;

        // Legacy objects
        auto x_ = legacy_vector(x);
        auto y_ = legacy_vector(y);

        // Constants to forward
        const idx_t& n = x_.n;
        const idx_t incx = (x_.direction == Direction::Forward) ? idx_t(x_.inc) : idx_t(-x_.inc);
        const idx_t incy = (y_.direction == Direction::Forward) ? idx_t(y_.inc) : idx_t(-y_.inc);

        return ::blas::copy( n, x_.ptr, incx, y_.ptr, incy );
    }

#endif

}  // namespace tlapack

#endif        //  #ifndef TLAPACK_BLAS_COPY_HH
