//
//  Cubic.swift
//  Numerical
//
//  Created by Adam Roberts on 4/11/19.
//

import Foundation

/// Inverse Cubic Interpolation
///
/// Given four points on a function, (a, f(a)), (b, f(b)), (c, f(c)), (d, f(d)),
/// this solves for the root of the inverse interpolated cubic polynomial. The
/// four y values must be distinct. The TOMS paper recommends this modified
/// version of the usual Aitken / Neville approach to avoid roundoff errors.
///
/// TOMS Algorithm 748 (1995), Section 3, Subroutine ipzero(a, b, e, d, x̄)
func inverseCubicInterpolation(a: Double,
                               b: Double,
                               c: Double,
                               d: Double,
                               fa: Double,
                               fb: Double,
                               fc: Double,
                               fd: Double) -> Double {
    let q11 = (c - d) * fc / (fd - fc)
    let q21 = (b - c) * fb / (fc - fb)
    let q31 = (a - b) * fa / (fb - fa)
    let d21 = (b - c) * fc / (fc - fb)
    let d31 = (a - b) * fb / (fb - fa)
    let q22 = (d21 - q11) * fb / (fd - fb)
    let q32 = (d31 - q21) * fa / (fc - fa)
    let d32 = (d31 - q21) * fc / (fc - fa)
    let q33 = (d32 - q22) * fa / (fd - fa)
    let r = a + q31 + q32 + q33
    return r
}
