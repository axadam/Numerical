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
/// We follow the approach advised by Temme in the Numerical Aspects part of
/// his chapter (§11.3.4), however the continued fraction is currently used everywhere.
///
/// Special Functions, N.M. Temme, 1996, §11.3
public func beta_reg(xy: Probability, a: Double, b: Double) -> Probability {
    let x = xy.p
    let y = xy.q
    
    // Inflection point x₀ = a / (a + b)
    let x₀ = a / (a + b)
    switch (x,y,a,b) {
        // handle domain edges
    case (_,_,...0,_): return .nan
    case (_,_,_,...0): return .nan
    case (..<0,_,_,_): return .nan
    case (1.nextUp...,_,_,_): return .nan
    case (0,_,_,_): return .p(0)
    case (_,0,_,_): return .q(0)
        
        // for x > x₀ we use I(1 - x,a,b) = 1 - I(x,b,a), Eq. 11.30
    case (x₀.nextUp...,_,_,_): return beta_reg(xy: xy.complement, a: b, b: a).complement
        
        // for most of the domain we can use the continued fraction
    case (_,_,_,_): return beta_reg_frac(xy: xy, a: a, b: b)
    }
}

public func beta_reg(x: Double, a: Double, b: Double) -> Probability {
    return beta_reg(xy: .p(x), a: a, b: b)
}

/// Continued Fraction estimate of Regularized Incomplete Beta Function
///
/// A continued fraction representation of the regularized incomplete beta function:
///
/// Iₓ(a,b) = x^a (1 - x)^b / (a B(a,b) ) (1 / (1 +) d₁ / (1 +) d₂ / (1 +) d₃ / (1 +) ), Eq. 8.17.22
///
/// d_{2m} = m (b - m) x / (a + 2m - 1) / (a + 2m)
///
/// d_{2m+1} = -(a + m) (a + b + m) x / (a + 2m) / (a + 2m - 1), Eq. 8.17.23
///
/// This continued fraction works well for a wide range of values. Avoid when x is close to x₀ = a / (a + b).
///
/// NIST Handbook of Mathematical Functions, 2010
public func beta_reg_frac(xy: Probability, a: Double, b: Double) -> Probability {
    let x = xy.p
    let y = xy.q
    let prefix = pow(x,a) * pow(y,b) / (a * beta(a: a, b: b))
    let cf = continued_fraction(b0: 0, a: { i in
        let mInt = (i - 1) / 2
        let m = Double(mInt)
        switch (i - 1) {
        case 0:
            return 1
        case 2 * mInt + 1:
            return -(a + m) * (a + b + m) / (a + 2 * m) / (a + 2 * m + 1) * x
        case _:
            return m * (b - m) / (a + 2 * m - 1) / (a + 2 * m) * x
        }
    }, b: { _ in 1 })
    return Probability(p: prefix * cf)
}

/// Asymptotic expansion of Incomplete Beta Function when `a` is large and a ≫ b
///
/// As a → ∞ with b and x fixed we have the following expansion (originally from Temme):
///
/// I(x, a, b) = 𝛤(a + b) / 𝛤(a) Σi=0... dᵢFᵢ, Eq. 2
///
/// Where the Fᵢ have the following recurrence relationship:
///
/// aFᵢ₊₁ = (i + b - aξ)Fᵢ + iξFᵢ₋₁, Eq. 3
///
/// F₀ = a^(-b) Q(b,aξ),
///
/// F₁ = (b - aξ) F₀ / a + ξ^b e^(-aξ) / (a𝛤(b))
///
/// The coefficients dᵢ have the generating function:
///
/// ((1 - e^(-t)) / t) ^ (b - 1) = Σi=0... dᵢ(t - ξ)ⁱ
///
/// "Uniform Asymptotic Expansion for the Incomplete Beta Function", Nemes & Olde Daalhius, 2016
///
/// Special Functions, N. M. Temme, 1996, §11.3.3.1
public func beta_reg_biga(x: Double, a: Double, b: Double) -> Double {
    let ξ = -log(x)
    
    /// 𝛤(a + b) / 𝛤(a), Eq. 2
    let prefix = exp(lgamma(a + b) - lgamma(a))
    
    let F₀ = pow(a, -b) * gamma_reg(b, a * ξ).q
    let F₁ = (b - a * ξ) / a * F₀ + pow(ξ, b) * exp(-a * ξ) / (a * tgamma(b))
    let Fᵢ₊₁ = { (n: Double, Fn: Double, Fnm1: Double) -> Double in
        return 1 / a * ((n + b - a * ξ) * Fn + n * ξ * Fnm1)
    }
    
    /// d₀ = ((1 - x) / ξ)^(b - 1), Eq. 4
    let d₀ = pow((1 - x) / ξ, b - 1)
    /// d₁ = (xξ + x - 1) (b - 1) d₀ / ((1 - x)ξ), Eq. 4
    let d₁ = (x * ξ + x - 1) * (b - 1) * d₀ / (1 - x) / ξ
    /// ξ(i + 1)(i + 2)d₀ dᵢ₊₂ = s₁ + s₂ + s₃, Eq. 5
    let dᵢ₊₂ = { (d: [Double]) -> Double in
        /// if b = 1 all terms after the first are zero
        if b == 1 { return 0 }
        
        /// infer i from number of previous terms. e.g. when computing d₂ we have |(d₀,d₁)| = 2 and i = 0.
        /// here we use the m and n terminology from Nemes.
        let nint = d.count - 2
        if nint < 0 { return .nan }
        let n = Double(nint)
        
        /// prefix = 1 / (ξ(i + 1)(i + 2)d₀)
        let prefix = 1 / ξ / (n + 1) / (n + 2) / d[0]
        
        /// s₁ = ξ Σm=0...i (m + 1)(n - 2m + 1 + (m - n - 1)/(b - 1)) d_(m+1) d_(n - m + 1)
        let s₁ = ξ * (0...nint).map { mint in
            let m = Double(mint)
            return (m + 1) * (n - 2 * m + 1 + (m - n - 1) / (b - 1)) * d[mint + 1] * d[nint - mint + 1]
        }.sum_naive()
        /// s₂ = Σm=0...i (m + 1)(n - 2m - 2 - ξ + (m - n) / (b - 1)) d_(m+1) d_(n - m)
        let s₂ = (0...nint).map { mint -> (Double) in
            let m = Double(mint)
            return (m + 1) * (n - 2 * m - 2 - ξ + (m - n) / (b - 1)) * d[mint + 1] * d[nint - mint]
        }.sum_naive()
        /// s₃ = Σm=0...i (1 - m - b) d_m d_(n - m)
        let s₃ = (0...nint).map { mint -> (Double) in
            let m = Double(mint)
            return (1 - m - b) * d[mint] * d[nint - mint]
        }.sum_naive()
        return prefix * (s₁ + s₂ + s₃)
    }
    
    let r = recursiveSum(indices: 1..., sum0: d₀ * F₀ + d₁ * F₁, state0: (Fᵢ₋₁: F₀, Fᵢ: F₁, d:[d₀,d₁]), update: { i, previous in
        let (Fᵢ₋₁,Fᵢ,d) = previous
        let F = Fᵢ₊₁(Double(i), Fᵢ, Fᵢ₋₁)
        let dnew = dᵢ₊₂(d)
        return (dnew * F, (Fᵢ₋₁: Fᵢ, Fᵢ: F, d:d + [dnew]))
    }, until: { a, b in print(prefix * b.0); return abs(b.0 - a.0) < 1e-12 * a.0 })
    
    return prefix * r
}

/// Derivative of Regularized Incomplete Beta function
///
/// I'(x, a, b) = x^(a - 1) (1 - x)^(b - 1) / B(a, b)
public func beta_reg_deriv(x: Double, a: Double, b: Double) -> Double {
    switch (x,a,b) {
    case (_,...0,_): return .nan
    case (_,_,...0): return .nan
    case (..<0,_,_): return .nan
    case (1.0.nextUp...,_,_): return .nan
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
public func inv_beta_reg(p: Probability, a: Double, b: Double) -> Double {
    // Cases we can handle directly in closed form
    switch (p.p,p.q,a,b) {
    // handle domain edges
    case (_,_,...0,_): return .nan
    case (_,_,_,...0): return .nan
    case (..<0,_,_,_): return .nan
    case (0,_,_,_): return 0
    case (_,0,_,_): return 1
    case (1...,_,_,_): return .nan
        
    // a and b both 1 function is identity. I(x,1,1) = x
    case (_,_,1,1): return p.p
        
    // If only one is 1 we make it be b and I⁻¹(p,a,1) = p^(1/a)
    case (_,_,1,_): return 1 - inv_beta_reg(p: .p(p.q), a: b, b: a)
    case (..<0.5,_,_,1): return pow(p.p, 1 / a)
    case (_,_,_,1): return exp(log1p(-p.q) / a)
        
    // Both one half, I⁻¹(p,1/2,1/2) = sin(p π / 2)²
    case (_,_,0.5,0.5): return sin(p.p * Double.pi/2)^^2

    // Otherwise we continue
    case (_,_,_,_): break
    }
    
    // Get our initial guess for root finding
    let guess: Double = {
        switch (a,b) {
        // Both a and b greater than 1
        // approximation from HMF §26.5.22
        case (1...,1...):
            // inverse normal approximation
            let yp = p.p < 0.5 ? qapprox(p: p.p) : -qapprox(p: p.q)
            
            let λ = (yp^^2 - 3.0) / 6.0
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
            if p.p < t / w { return pow(a * w * p.p, 1 / a) }
            return 1 - pow(b * w * p.q, 1 / b)
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
    let x = root(guess: guess,
                    xmin: 0,
                    xmax: 1,
                    maxIter: 10,
                    f: { x in beta_reg(x: x, a: a, b: b).difference(p) },
                    f1: { x in exp(a1 * log(x) + b1 * log(1 - x) + afac) },
                    f2f1: { x in a1 / x - b1 / (1 - x) })
    return x
}

