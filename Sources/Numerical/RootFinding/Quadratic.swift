//
//  Quadratic.swift
//  Numerical
//
//  Created by Adam Roberts on 4/3/19.
//

import Foundation

/// Newton-Quadratic root finding step
///
/// Finds a root in the interval [a,b] of the quadratic polynomial
/// that goes through (a, f(a)) and (b, f(b) and near (d, f(d)):
///
/// P(a,b,d)(x) = f(a) + f[a,b] (x - a) + f[a,b,d] (x - a) (x - b)
///
/// P'(a,b,d)(x) = f[a,b] + f[a,b,d] (2x - a - b)
///
/// Then Newton's method lets us step to the root of P
///
/// TOMS Algorithm 748 (1995), Section 2, Subroutine Newton-Quadratic(a, b, d, r, k)
func newtonQuadtraticStep(a: Double,
                      b: Double,
                      d: Double,
                      fa: Double,
                      fb: Double,
                      fd: Double,
                      k: Int) -> Double {
    let fab = (fb - fa) / (b - a)
    let fbd = (fd - fb) / (d - b)
    let fabd = (fbd - fab) / (d - a)
    if fabd == 0 { return a - fa / fab } // note this is the same as secant
    let r0 = fa * fabd > 0 ? a : b
    let r = root(guess: r0, maxIter: k,
                 f: { (x: Double) in fa + (x - a) * (fab + (x - b) * fabd)},
                 f1: { (x: Double) in fab + fabd * (2.0 * x - a - b) }).value
    return r
}
