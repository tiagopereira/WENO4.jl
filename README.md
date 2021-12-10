[![Build Status](https://github.com/tiagopereira/WENO4.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tiagopereira/WENO4.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/tiagopereira/WENO4.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/tiagopereira/WENO4.jl)


# WENO4.jl

A package to perform interpolation using the 4th order Weighted Essentially Non-Oscillatory (WENO) scheme of [Janett et al (2019)](https://ui.adsabs.harvard.edu/abs/2019A%26A...624A.104J/abstract). Based on [Weno4Interpolation](https://github.com/Goobley/Weno4Interpolation) by Chris Osborne.


## Installation

From the Julia REPL:

```julia
using Pkg
Pkg.add("WENO4")
```


## Example Usage

Create a grid `xp` and an array `fp` of values to be interpolated
```julia
xp = 1:0.2:5
f(x) = log(x)
fp = f.(xp)
```
Interpolate to a new set of points `xs`:
```julia
xs = [1.1, 2.1, 3.1, 4.1]
result = interpolate_weno4(xs, xp, fp)
```