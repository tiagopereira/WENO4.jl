"""
Interpolation routines for 4th order Weighted Essentially Non-Oscillatory (WENO)
schemes.
"""
module WENO4

export interpolate_weno4

"""
    function interpolate_weno4(
        xs::AbstractVector{<:Number},
        xp::AbstractVector{<:Number},
        fp::AbstractVector{<:Number};
        extrapolate=false
    )

Performs 1D interpolation using the fourth-order Weighted Essentially
Non-Oscillatory (WENO) scheme from [Janett et al (2019)](https://ui.adsabs.harvard.edu/abs/2019A%26A...624A.104J/abstract).
Based on https://github.com/Goobley/Weno4Interpolation by Chris Osborne.

# Arguments

- `xs::AbstractVector{<:Number}`: The x-coordinates at which to evaluate the
  interpolated values
- `xp::AbstractVector{<:Number}`: The x-coordinates of the data points. Must
  have at least four elements, be increasing, and have no NaNs.
- `fp::AbstractVector{<:Number}`: The y-coordinates of the data points. Must
  have same length as `xp`.

# Keywords

- `extrapolate::Bool`: whether or not to use extrapolation. If `true`,
  will perform quadratic extrapolation from the edge points. When
  `false`, it will not throw an error if `xs` is outside `xp`, but return
  the closest `fp` value.

"""
function interpolate_weno4(
    xs::AbstractVector{<:Number},
    xp::AbstractVector{<:Number},
    fp::AbstractVector{<:Number};
    extrapolate=false,
)
    Ngrid = length(xp)
    result = similar(xs)
    ε = 1.0f-6
    prevβ = -1
    β2 = zero(eltype(fp))
    β3 = zero(eltype(fp))
    @inbounds for (idx, x) in enumerate(xs)
        if x < xp[1]
            if !extrapolate
                result[idx] = fp[1]
                continue
            end
            # quadratic extrapolation from first three points
            i = 1
        elseif x > xp[end]
            if !extrapolate
                result[idx] = fp[end]
                continue
            end
            # quadratic extrapolation from last three points
            i = Ngrid
        else
            i = binary_search(xp, x)
        end
        if i == Ngrid
            i -= 1
        end
        # set stencil, pad with zeros when not relevant
        xi = xp[i]
        xip = xp[i+1]
        yi = fp[i]
        yip = fp[i+1]
        if i == 1
            xim = zero(eltype(xp))
            xipp = xp[i+2]
            yim = zero(eltype(fp))
            yipp = fp[i+2]
        elseif i == Ngrid - 1
            xim = xp[i-1]
            xipp = zero(eltype(xp))
            yim = fp[i-1]
            yipp = zero(eltype(fp))
        else
            xim = xp[i-1]
            xipp = xp[i+2]
            yim = fp[i-1]
            yipp = fp[i+2]
        end
        q2, q3 = weno4_q(x, xim, xi, xip, xipp, yim, yi, yip, yipp)
        if i == 1
            result[idx] = q3
        elseif i == Ngrid - 1
            result[idx] = q2
        else
            if i != prevβ  # reuse previously computed when possible
                β2, β3 = weno4_β(xim, xi, xip, xipp, yim, yi, yip, yipp)
                prevβ = i
            end
            γ2 = -(x - xipp) / (xipp - xim)
            γ3 = (x - xim) / (xipp - xim)
            α2 = γ2 / (ε + β2)
            α3 = γ3 / (ε + β3)
            ω2 = α2 / (α2 + α3)
            ω3 = α3 / (α2 + α3)
            result[idx] = ω2 * q2 + ω3 * q3
        end
    end
    return result
end

"""
    binary_search(a::AbstractVector{<:Number}, x)

Return the index of the last value in `a` less than or equal to `x`.
`a` is assumed to be sorted.
"""
function binary_search(a::AbstractVector{<:Number}, x)
    n = length(a)
    ileft = 2
    i = n - 1
    while ileft <= i
        mid = (ileft + i) >> 1
        if @inbounds x < a[mid]
            i = mid - 1
        else
            ileft = mid + 1
        end
    end
    return i
end

"""
Compute the smoothness indicators β2 and β3 for 4th order WENO interpolation.
Arguments are the x and y values for the 4-point stencil.
"""
function weno4_β(xim, xi, xip, xipp, yim, yi, yip, yipp)
    him = xi - xim
    hi = xip - xi
    hip = xipp - xip
    H = him + hi + hip
    # Derivatives
    yyim = -((2 * him + hi) * H + him * (him + hi)) / (him * (him + hi) * H) * yim
    yyim += ((him + hi) * H) / (him * hi * (hi + hip)) * yi
    yyim -= (him * H) / ((him + hi) * hi * hip) * yip
    yyim += (him * (him + hi)) / ((hi + hip) * hip * H) * yipp
    yyi = -(hi * (hi + hip)) / (him * (him + hi) * H) * yim
    yyi += (hi * (hi + hip) - him * (2 * hi + hip)) / (him * hi * (hi + hip)) * yi
    yyi += (him * (hi + hip)) / ((him + hi) * hi * hip) * yip
    yyi -= (him * hi) / ((hi + hip) * hip * H) * yipp
    yyip = (hi * hip) / (him * (him + hi) * H) * yim
    yyip -= (hip * (him + hi)) / (him * hi * (hi + hip)) * yi
    yyip += ((him + 2 * hi) * hip - (him + hi) * hi) / ((him + hi) * hi * hip) * yip
    yyip += ((him + hi) * hi) / ((hi + hip) * hip * H) * yipp
    yyipp = -((hi + hip) * hip) / (him * (him + hi) * H) * yim
    yyipp += (hip * H) / (him * hi * (hi + hip)) * yi
    yyipp -= ((hi + hip) * H) / ((him + hi) * hi * hip) * yip
    yyipp += ((2 * hip + hi) * H + hip * (hi + hip)) / ((hi + hip) * hip * H) * yipp
    # Smoothness indicators
    β2 = (hi + hip)^2 * (abs(yyip - yyi) / hi - abs(yyi - yyim) / him)^2
    β3 = (him + hi)^2 * (abs(yyipp - yyip) / hip - abs(yyip - yyi) / hi)^2
    return (β2, β3)
end

"""
Compute the q2 and q3 Lagrange interpolation factors for 4th order WENO interpolation.
Arguments are the value to interpolate plus the x and y values for the 4-point stencil.
"""
function weno4_q(x, xim, xi, xip, xipp, yim, yi, yip, yipp)
    him = xi - xim
    hi = xip - xi
    hip = xipp - xip
    q2 = yim * ((x - xi) * (x - xip)) / (him * (him + hi))
    q2 -= yi * ((x - xim) * (x - xip)) / (him * hi)
    q2 += yip * ((x - xim) * (x - xi)) / ((him + hi) * hi)
    q3 = yi * ((x - xip) * (x - xipp)) / (hi * (hi + hip))
    q3 -= yip * ((x - xi) * (x - xipp)) / (hi * hip)
    q3 += yipp * ((x - xi) * (x - xip)) / ((hi + hip) * hip)
    return (q2, q3)
end

end
