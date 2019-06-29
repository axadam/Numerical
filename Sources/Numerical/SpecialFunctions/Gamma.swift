//
//  Gamma.swift
//  Numerical
//
//  Created by Adam Roberts on 9/25/18.
//

import Foundation

/// Regularized Incomplete Gamma Function (lower), P(a,x)
///
/// P(a,x) = 𝛾(a,x) / 𝛤(a), 𝛾(a,x) = ∫0..x e^(-t) t^(a-1) dt, a > 0
///
/// For small x (less than a + 1) use the series approach, for large
/// use the continued fraction approach. NR also advises using quadrature
/// for large a but that is not yet implemented.
///
/// Numerical Receipes §6.2
public func p_gamma(_ a: Double, _ x: Double) -> Double {
    switch (a,x) {
    case (...0,_): return .nan
    case (_,..<0): return .nan
    case (_,0): return 0
    case (_,..<(a + 1)): return p_gamma_series(a: a, x: x)
    case (_,_): return 1 - q_gamma_frac(a: a, x: x)
    }
}

/// Regularized Incomplete Gamma Function (upper), Q(a,x)
///
/// Q(a,x) = 𝛤(a,x) / 𝛤(a), 𝛤(a,x) = ∫x..∞ e^(-t) t^(a-1) dt, a > 0
///
/// Implemented simply as complement of lower: Q(a,x) = 1 - P(a,x)
public func q_gamma(_ a: Double, _ x: Double) -> Double {
    return 1 - p_gamma(a, x)
}

/// Series approximation of P(a,x)
///
/// 𝛾(a,x) = e^(-x) x^a Σ0..∞ 𝛤(a) / 𝛤(a + 1 + n) x^n
///
/// Compute the denominator recursively: 𝛤(z + 1) = z 𝛤(z)
///
/// For initial term: 𝛤(a) / 𝛤(a + 1) = 1 / a
///
/// Numerical Receipes §6.2
func p_gamma_series(a: Double, x: Double) -> Double {
    let prefix = exp(a * log(x) - x - lgamma(a))
    let first = 1 / a
    let sum = recursiveSum(indices: 1..., sum0: first, state0: first, update: { i, state in
        let ap = a + Double(i)
        let state1 = state * x / ap
        return (state1, state1)
    }, until: { a, b in a.0 == b.0 })
    return prefix * sum
}

/// Continued fraction approximation of Q(a,x)
///
/// Q(a,x) = e^(-x) x^a 1 / (1 + x - a -) 1 (1 - a) / (3 + x - a -)  2 (2 - a) / (5 + x - a -)
///
/// This is the even part of the following (converges faster):
///
/// Q(a,x) = e^(-x) x^a 1 / (x +) (1 - a) / (1 +) 1 / (x +) (2 - a) / (2 +) 2 / (x +)
///
/// Numerical Receipes §6.2
func q_gamma_frac(a: Double, x: Double) -> Double {
    let prefix = exp(a * log(x) - x - lgamma(a))
    let b = 1 + x - a
    let c = 1 / Double.leastNormalMagnitude
    let d = 1 / b
    let frac = recursiveProduct(indices: 1..., product0: d, state0: (b: b, c: c, d: d), update: { i, state in
        let an = Double(-i) * (Double(i) - a)
        let (b0, c0, d0) = state
        let b1 = b0 + 2
        let c1 = max(Double.leastNormalMagnitude, b1 + an / c0)
        let d1 = 1 / max(Double.leastNormalMagnitude, b1 + an * d0)
        return (c1 * d1, (b: b1, c: c1, d: d1))
    }, until: { a, b in abs(b.1 - 1) < 1e-15 })
    return prefix * frac
}

/// Derivative of regularized lower incomplete gamma function, P
///
/// e^-x * x^(a-1) / Γ(a) = e^(-x + (a-1) * log(x) - logΓ(a))
public func p_gamma_deriv(a: Double, x: Double) -> Double {
    return exp(-x + (a - 1) * log(x) - lgamma(a))
}

/// Inverse of the regularized incomplete gamma P(a,x) function.
/// Gives x such that P(a,x) = p.
///
/// Start with approximation NR or from Handbook of Mathematical Functions vis NR.
/// Use Halley method to find root of P(a,x) - p. From Numerical Recipes §6.2.1
public func inv_p_gamma(_ a: Double, _ p: Double) -> Double {
    switch (a, p) {
        // handle domain edges
    case (...0,_): return .nan
    case (_,..<0): return .nan
    case (_,0): return 0
    case (_,1): return .infinity
    case (_,1...): return .nan
        
        // normal case
    case (_,_):
        let a1 = a - 1
        let lna1 = log(a1)
        let gln: Double = lgamma(a)
        let afac = exp(a1 * (lna1 - 1) - gln)
        
        // initial guess
        // a <= 1: P(a,1) ~ 0.253a + 0.12a²
        // a  > 1: normal quantile approx and then Handbook of Mathematical Functions, §26.4.17
        // Q(p) = 2 P⁻¹(ν / 2, p) = ν  ( 1 - 2/(9ν)  + xp √(2/(9ν))  )³
        //      = 2 P⁻¹(a    , p) = 2a ( 1 - 2/(18a) + xp √(2/(18a)) )³
        //          P⁻¹(a    , p) =  a ( 1 - 1/(9a)  + xp √(1/(9a))  )³
        let guess: Double = {
            switch a {
            case ...1:
                let t = 1 - a * (0.253 + a * 0.12)
                if p < t {
                    return pow(p / t, 1 / a)
                } else {
                    return 1 - log( 1 - (p-t) / (1-t))
                }
            case _:
                let xp = p < 0.5 ? qapprox(p: p) : -qapprox(p: 1 - p)
                return fmax(1e-3, a * (1 - 1/(9 * a) - xp / (3 * sqrt(a)))^^3)
            }
        }()
        
        // Halley method
        // p_gamma'(x) = e^-x * x^(a-1) / Γ(a)
        // p_gamma''(x) = e^-x (a - x - 1) x^(a-2) / Γ(a)
        // p_gamma''(x) / p_gamma'(x) = e^-x (a - x - 1) x^(a-2) / e^-x x^(a-1)
        //                            = (a - x - 1) / x = (a-1)/x - 1
        let x = rootSecondOrder(guess: guess,
                        xmin: 0,
                        maxIter: 11,
                        f: { x in p_gamma(a, x) - p },
                        f1: { x in
                            switch a {
                            case ...1:
                                return exp( -x + a1 * log(x) - gln)
                            case _:
                                return afac * exp( -(x - a1) + a1 * (log(x) - lna1))
                            }
                        },
                        f2f1: { x in a1 / x - 1 })
        return x
    }
}

