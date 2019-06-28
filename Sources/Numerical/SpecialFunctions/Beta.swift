//
//  Beta.swift
//  Numerical
//
//  Created by Adam Roberts on 9/25/18.
//

import Foundation

/// Beta function
///
/// B(a, b) = ∫0..1 t^(a - 1) (1 - t)^(b - 1) dt
///
/// = 𝛤(a) 𝛤(b) / 𝛤(a + b)
public func beta(a: Double, b: Double) -> Double {
    if (a + b > 100) { return exp(lbeta(a: a, b: b)) }
    return tgamma(a) * tgamma(b) / tgamma(a + b)
}

/// Log of Beta function
///
/// log B(a, b) = log 𝛤(a) + log 𝛤(b) - log 𝛤(a + b)
public func lbeta(a: Double, b: Double) -> Double {
    return lgamma(a) + lgamma(b) - lgamma(a + b)
}

/// Regularized Incomplete Beta function
///
/// I(x, a, b) = B(x, a, b) / B(a, b), B(x, a, b)
///
/// B(x, a, b) = ∫0..x t^(a - 1) (1 - t)^(b - 1) dt where a, b > 0
///
/// Evaluated using a continuous fraction. NR recommends for large a, b
/// to use quadrature but we haven't implemented yet.
///
/// Numerical Receipes §6.4
public func beta_reg(x: Double, a: Double, b: Double) -> Double {
    switch (x,a,b) {
        // handle domain edges
    case (_,...0,_): return .nan
    case (_,_,...0): return .nan
    case (..<0,_,_): return .nan
    case (0,_,_): return 0
    case (1,_,_): return 1
    case (1...,_,_): return .nan
        
        // FIXME: NR recommends quadrature due to slow convergence
//        case (_,3000...,3000...):
        
        // normal case
    case (..<((a + 1.0) / (a + b + 2.0)),_,_):
        return beta_reg_frac(x: x, a: a, b: b)
        
        // x above threshold due the flip
    case (_,_,_):
        return 1.0 - beta_reg_frac(x: 1.0 - x, a: b, b: a)
    }
}

/// Continued fraction approximation of I(x, a, b)
///
/// I(x, a, b) = prefix 1 / (1 +) d1 / (1 +) d2 / (1 +)
///
/// prefix = x^a (1 - x)^b / (a B(a, b))
///
/// d_(2m+1) = (a + m) (a + b + m) x / (a + 2m) / (a + 2m + 1)
///
/// d_(2m) = m (b - m) x / (a + 2m - 1) / (a + 2m)
///
/// Numerical Receipes §6.4
public func beta_reg_frac(x: Double, a: Double, b: Double) -> Double {
    let prefix = exp(a * log(x) + b * log(1 - x) - log(a) - lbeta(a: a, b: b))
    let qab = a + b
    let qap = a + 1
    let qam = a - 1
    let c = 1.0
    let d = 1 / absmax(1 - qab * x / qap)
    let frac = recursiveProduct(indices: 1..., product0: d, state0: (c: c, d: d), update: { i, state in
        let (c0, d0) = state
        let m = Double(i)
        let m2 = m * 2
        let aeven = m * (b - m) * x / ((qam + m2) * (a + m2))
        let d1 = 1 / absmax(1 + aeven * d0)
        let c1 = absmax(1 + aeven / c0)
        let aodd = -((a + m) * (qab + m) * x) / ((a + (m2)) * (qap + (m2)))
        let d2 = 1 / absmax(1 + aodd * d1)
        let c2 = absmax(1 + aodd / c1)
        return (d1 * c1 * d2 * c2, (c: c2, d: d2))
    }, until: { a, b in abs(b.1 - 1) < 1e-15 })
    let result = frac * prefix
    return result
}

/// Derivative of Regularized Incomplete Beta function
///
/// I'(x, a, b) = x^(a - 1) (1 - x)^(b - 1) / B(a, b)
public func beta_reg_deriv(x: Double, a: Double, b: Double) -> Double {
    switch (x,a,b) {
    case (_,...0,_): return .nan
    case (_,_,...0): return .nan
    case (..<0,_,_): return .nan
    case (_,_,_) where x > 1: return .nan
    case (0,..<1,_): return .nan
    case (0,1,_): return 1 / beta(a: a,b: b)
    case (0,_,_): return 0
    case (1,_,..<1): return .nan
    case (1,_,1): return 1 / beta(a: a, b: b)
    case (1,_,_): return 0
    case (_,_,_): return pow(x,a - 1) * pow(1 - x,b-1) / beta(a: a,b: b)
    }
}

/// Inverse of the Regularized Incomplete Beta function
///
/// xp such that I(xp, a, b) = p
///
/// Handles some special cases first, then general case uses an approximation
/// followed by Halley's method on I(x, a, b) - p
///
/// Numerical Receipes §6.4
public func inv_beta_reg(p: Double, a: Double, b: Double) -> Double {
    let q = 1 - p
    
    // Cases we can handle directly in closed form
    switch (p,a,b) {
    // handle domain edges
    case (_,...0,_): return .nan
    case (_,_,...0): return .nan
    case (..<0,_,_): return .nan
    case (0,_,_): return 0
    case (1,_,_): return 1
    case (1...,_,_): return .nan
        
    // a and b both 1 function is identity. I(x,1,1) = x
    case (_,1,1): return p
        
    // If only one is 1 we make it be b and I⁻¹(p,a,1) = p^(1/a)
    case (_,1,_): return 1 - inv_beta_reg(p: q, a: b, b: a)
    case (..<0.5,_,1): return pow(p, 1 / a)
    case (_,_,1): return exp(log1p(-q) / a)
        
    // Both one half, I⁻¹(p,1/2,1/2) = sin(p π / 2)²
    case (_,0.5,0.5): return sin(p * Double.pi/2)^^2

    // Otherwise we continue
    case (_,_,_): break
    }
    
    // Get our initial guess for root finding
    let guess: Double = {
        switch (a,b) {
        // Both a and b greater than 1
        // approximation from HMF §26.5.22
        case (1...,1...):
            // inverse normal approximation
            let yp = p < 0.5 ? qapprox(p: p) : -qapprox(p: 1 - p)
            
            let λ = (yp^^2 - 3) / 6
            let h = 2 / (1 / (2 * a - 1) + 1 / (2 * b - 1))
            let w = yp * sqrt(h + λ) / h - (1 / (2 * b - 1) - 1 / (2 * a - 1)) * (λ + 5 / 6 - 2 / (3 * h))
            return a / (a + b * exp(2 * w))
        // At least one of a and b < 1. Use NR approximation
        case (_,_):
            let lna = log(a / (a + b))
            let lnb = log(b / (a + b))
            let t = exp(a * lna) / a
            let u = exp(b * lnb) / b
            let w = t + u
            if p < t / w { return pow(a * w * p, 1 / a) }
            return 1 - pow(b * w * q, 1 / b)
        }
    }()
    
    // Halley's method
    // I'(x,a,b) = x^(a - 1) (1 - x)^(b - 1) / B(a,b)
    // I''(x,a,b) = -x^(a - 2) (1 - x)^(b - 2) (x (b - 2) + (x - 1) a + 1) / B(a,b)
    // I''(x,a,b) / I'(x,a,b) = (x (b - 2) + (x - 1) a + 1) / (x (1 - x))
    //            = (b - 2) / (1 - x) - a / x + 1 / (x (1 - x))
    //            = (b - 1) / (1 - x) - (a - 1) / x - 1 / (1 - x) - 1 / x + 1 / (x (1 - x))
    //            = (b - 1) / (1 - x) - (a - 1) / x - x / (x (1 - x)) - (1 - x) / (x (1 - x)) + 1 / (x (1 - x))
    //            = (b - 1) / (1 - x) - (a - 1) / x
    let afac = -lbeta(a: a, b: b)
    let a1 = a - 1
    let b1 = b - 1
    let x = rootSecondOrder(guess: guess,
                    xmin: 0,
                    xmax: 1,
                    maxIter: 10,
                    f: { x in beta_reg(x: x, a: a, b: b) - p },
                    f1: { x in exp(a1 * log(x) + b1 * log(1 - x) + afac) },
                    f2f1: { x in a1 / x - b1 / (1 - x) })
    return x
}

