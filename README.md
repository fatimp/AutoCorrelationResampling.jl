# AutoCorrelationResampling
[![CI](https://github.com/shamazmazum/AutoCorrelationResampling.jl/actions/workflows/test.yml/badge.svg)](https://github.com/shamazmazum/AutoCorrelationResampling.jl/actions/workflows/test.yml)

This package provides means to resample (currently, only upsampling is
supported) of autocorrelation functions.

Autocorrelation function of one variable is a Laurent polynomial on ℝ in the
form $s(x) = f(x)f(x^{-1})$ where $f(x)$ is a usual polynomial on ℝ.

This package provides resampling of autocorrelation function in the sense that
it changes the degree of $s(x)$ while maintaining the form $f(x)f(x^{-1})$.

Technically, this package works with autocorrelation function of three variables
$s(x,y,z)$ in the form of three-dimensional arrays and rescales them along the
third variable (or axis), i.e. an array of shape `(x, y, z)` becomes an array of
shape `(x, y, nz)` where `n` is a resampling factor.

## Why is it needed?

This package is an attempt to solve the problem of reconstruction of a porous
media when only a fraction of information about original media is available
(e.g. you have to reconstruct 3D cube from a stack of 2D slices taken along `z`
axis). If those slices are evenly sampled (i.e. you have each `n`-th slice), you
can do the following:

1. Stack the slices in 3D array. This array will have a length along `z` axis
   reduced by `n` times.
2. Calculate autocorrelation function (also known as two-point function) for the
   reduced array.
3. Upsample it by `n` times using this package
4. Reconstruct original 3D image (you can use `PhaseRec.jl` package for it).

## How to use?

Upsampling is done by `ac_upsample` function which takes an autocorrelation
array and a low-pass filter. A low-pass filter can be obtained with
`filter_coeffs` function. These filters are not unique, you can play with them
supplying different argument `initial` to `filter_coeffs`.
