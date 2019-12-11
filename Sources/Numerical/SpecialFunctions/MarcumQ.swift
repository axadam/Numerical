//
//  MarcumQ.swift
//  Numerical
//
//  Created by Adam Roberts on 7/28/19.
//

import Foundation
import Scan

// MARK: Marcum Q

/// Generalized Marcum Q function
///
/// The Marcum Q function Qᵤ(x,y) is equivalent to the upper tail of the distribution
/// function of the non-central gamma distribution of order µ and non-centrality
/// parameter x evaluated at point y.
///
/// When x = 0 it is equivalent to the regularized
/// incomplete gamma function of order µ. That is,
///
/// Qᵤ(0,y) = Q(µ,y) = 𝛤(µ,y) / 𝛤(µ)
///
/// The function originally was developed in the field of communications. There the
/// parameter µ is the number of samples of the output of a square law detector.
///
/// Domains for different methods of calculation and the methods themselves
/// follow the cited paper. The domains are delineated in §6.
///
/// "Computation of the Marcum Q Function", Gil, Segura, Temme 2013
///
/// - Parameters:
///   - µ: the order of the non-central gamma function
///   - x: the non-centrality parameter
///   - y: the point at which to evaluate
///
/// - Returns: A tuple of the lower (p) and upper (q) CDF tails of the distribution
public func marcum(µ: Double, x: Double, y: Double) -> Probability {
    let ξ = 2 * sqrt(x * y)
    
    // these two functions draw a parabola outside of which quadrature
    // works well. Eq. 100
    let f₁ = x + µ - sqrt(4 * x + 2 * µ)
    let f₂ = x + µ + sqrt(4 * x + 2 * µ)
    
    // switch over all the domains
    switch (µ,x,y,ξ) {
        
    // 1. if x < 30 compute the series expansion
    case (     _,..<30,  ...(x + µ),                           _):
        let p = p_series(µ: µ, x: x, y: y)
        return Probability(p: p)
    case (     _,..<30,  (x + µ)...,                           _):
        let q = q_series(µ: µ, x: x, y: y)
        return Probability(q: q)
        
    // 2. if ξ > 30 and µ² < 2ξ use the large ξ asymptotic expansion
    case (     _,    _,           _,max(30,0.5 * µ^^2).nextUp...):
        return bigxy(µ: µ, x: x, y: y)
        
    // 3. if f₁(x,µ) < y < f₂(x,µ) and µ < 135 use recurrence
    case (..<135,    _,f₁...(x + µ),                           _):
        let p = p_recursion(µ: µ, x: x, y: y)
        return Probability(p: p)
    case (..<135,    _,(x + µ)...f₂,                           _):
        let q = q_recursion(µ: µ, x: x, y: y)
        return Probability(q: q)
        
    // 4. if f₁(x,µ) < y < f₂(x,µ) and µ ≥ 135 use the large µ aysmptotic expansion
    case (135...,    _,     f₁...f₂,                           _):
        return bigmu(µ: µ, x: x, y: y)
        
    // 5. in other cases use the quadrature method
    case (     _,    _,           _,                           _):
        return quadrature(µ: µ, x: x, y: y)
    }
}

// MARK: Marcum Q Derivative

/// Derivative of the Marcum Q function with respect to y
///
/// ∂Qᵤ(x,y)/∂y = Qᵤ₋₁(x,y) - Qᵤ(x,y), Eq. 16
///
/// Note that this is the derivative of the Q, or upper tail. If you want the derivative of the P
/// then you must negate.
///
/// "Computation of the Marcum Q Function", Gil, Segura, Temme 2013, §2.3
public func marcum_deriv(µ: Double, x: Double, y: Double) -> Double {
    // FIXME: use the recurrence relation instead of doing two full calculations
    return marcum(µ: µ - 1, x: x, y: y).q - marcum(µ: µ, x: x, y: y).q
}

// MARK: Inverse Marcum Q

/// Inverse of the Marcum Q function with respect to y
///
/// For a given µ and x and (p,q) pair find the y such that Qᵤ(x,y) = (p,q)
///
/// We first come up with a guess for y based on an approximation of ζ, then use
/// root finding to get the final result. We use just the first two terms of the approximation of ζ.
///
/// ζ ~ ζ₀ + Σi=1... ζ ᵢ / µⁱ, Eq. 5.23
///
/// "The asymptotic and numerical inversion of the Marcum Q−function", Gil, Segura, Temme, 2014, §5.2
public func inv_marcum(µ: Double, x µx: Double, p: Probability) -> Double {
    // this asymptotic method works on x normalized for µ
    let x = µx / µ
        
    // 1/2 erfc(ζ₀√(µ/2)) = q, Eq. 5.11
    // ζ₀ = erfc⁻¹(2q) √(2/µ)
    // to avoid rounding error we work on the smaller of p and q
    let (s,prob) = p.p < 0.5 ? (-1.0,p.p) : (1.0,p.q)
    let ζ₀ = s * invErfC(2 * prob) * sqrt(2 / µ)
    
    // find y₀ corresponding to ζ₀
    let y₀ = y(ζ₀,x)
    
    // find ζ₁ corresponding to ζ₀, x, and y₀
    let ζ₁ = zeta1(ζ₀,x,y₀)
    
    // two term approximation of ζ
    // ζ ~ ζ₀ + ζ₁ / µ, Eq. 5.23
    let ζ = ζ₀ + ζ₁ / µ
    
    // scale back up by µ for final guess
    let guess = µ * y(ζ,x)
    
    // root finding to get final answer
    // FIXME: default bracketing goes wayy too wide here, wasting work
    let r = root(guess: guess, tolerance: 1e-15 * prob, bracketFactor: 1.001) { marcum(µ: µ, x: µx, y: $0).difference(p) }
    return r
}

// MARK: Implementation

/// Find y given ζ (zeta) and x
///
/// Where x, y, and ζ are related by:
///
/// 1/2 ζ² = x + y - √(1 + 4xy) + log [(1 + √(1 + 4xy)) / (2y)], Eq. 7.1
///
/// When ζ is not small we find the root with Newton's method. Note that there will be two roots
/// and we need to find the one that satisfies sign(ζ) = sign(y - x - 1).
///
/// When ζ is small |y - x - 1| is small and we use an expansion of y - x - 1:
///
/// y = x + 1 + Σi=1... bᵢ(x) ζ ⁱ, Eq. 7.5
///
/// "The asymptotic and numerical inversion of the Marcum Q−function", Gil, Segura, Temme, 2014, §7
fileprivate func y(_ ζ: Double, _ x: Double) -> Double {
    switch abs(ζ) {
    
    // Small ζ use expansion in Eq. 7.5
    // b₁(x) = √(2x + 1), Eq. 7.6
    // b₂(x) = (3x + 1) / (3(2x + 1)),
    // b₃(x) = (6x + 1) / (36(2x + 1)⁵/²)
    case ..<0.5:
        let x2p1 = 2.0 * x + 1.0
        let b₁ = sqrt(x2p1)
        let b₂ = (1 + 3 * x) / (3 * x2p1)
        let b₃ = (1 + 6 * x) / (36 * pow(x2p1, 2.5))
        let sum = evaluate_polynomial(poly: [0.0,b₁,b₂,b₃], z: ζ)
        return x + 1 + sum
        
    // Normal size ζ, use Newton's method
    case    _:
        // 1/2 ζ²
        let hζ² = 0.5 * ζ^^2
        
        // There will be two roots. We need to choose the starting point
        // so that we get the right root to satisfy:
        // sign(ζ) = sign(y - x - 1)
        // so if ζ < 0 we want the root less than x + 1 and vice versa
        let guess = ζ > 0 ? x + 2 : 0.5
        
        // Newton's method
        // f(y) = x + y - √(1 + 4xy) + log [(1 + √(1 + 4xy)) / (2y)] - 1/2 ζ²
        // fʹ(y) = (y - 2xy - 1 + (y - 1) √(1 + 4xy)) / (y(1 + √(1 + 4xy)))
        // Eq. 7.12
        let r = root(guess: guess,
                     xmin: 1e-100,
                     f: { y in
                        let sq = sqrt(1 + 4 * x * y)
                        return x + y - sq + log((1 + sq) / (2 * y)) - hζ²
        },
                     f1: { y in
                        let sq = sqrt(1 + 4 * x * y)
                        return (y - 2 * x * y - 1 + (y - 1) * sq) / (y * (1 + sq))
        })
        return r
    }
}

/// Find ζ₁ given ζ₀, x, and y₀
///
/// ζ₁ = 1 / ζ₀ log f(ζ₀), Eq. 5.24
///
/// f(ζ₀) = ζ₀/ (y - x - 1) (1 + 2x + √(1 + 4xy)) / (2 (1 + 4xy)¹/⁴), Eq. 5.26
///
/// When ζ₀ is small we can use the expansion:
///
/// ζ₁ = Σi=0... dᵢ(x) ζ₀ ⁱ, Eq. 7.7
///
/// "The asymptotic and numerical inversion of the Marcum Q−function", Gil, Segura, Temme, 2014, §5, §7
fileprivate func zeta1(_ ζ₀: Double, _ x: Double, _ y: Double) -> Double {
    switch abs(ζ₀) {
    // Small |ζ₀|, use expansion
    // d₀(x) = -1/3 (3x + 1) (2x + 1)⁻³/², Eq. 7.9
    // d₁(x) = 1/36 (36x² + x + 1) (2x + 1)⁻³,
    // d₂(x) = -1/1620 (2160x³ - 594x² - 9x - 1) (2x + 1)⁻⁹/²
    case ..<0.5:
        let x2p1 = 2.0 * x + 1.0
        let xx = x2p1 * sqrt(x2p1)
        let d₀ = -(3.0 * x + 1.0) / (3.0 * xx)
        let d₁ = (36.0 * x^^2 + x + 1.0) / (36.0 * xx^^2)
        let d₂ = -(2160.0 * x^^3 - 594.0 * x^^2 - 9.0 * x - 1) / (1620.0 * xx^^3)
        let sum = evaluate_polynomial(poly: [d₀,d₁,d₂], z: ζ₀)
        return sum
        
    // Large |ζ₀| can directly evalute
    case      _:
        let sq1pξ² = sqrt(1 + 4 * x * y)
        let f = ζ₀ / (y - x - 1.0) * (1.0 + 2.0 * x + sq1pξ²) / (2.0 * sqrt(sq1pξ²))
        return log(f) / ζ₀
    }
}

/// Series method upper tail
///
/// Qᵤ(x,y) = e^(-x) Σi=0..∞ xⁱ / i! Qᵤ₊ᵢ(y), Eq. 7
///
/// We use the following recurrence relationships:
///
/// Qᵤ₊₁(y) = Qᵤ(y) + y^µ e^(-y) / 𝛤(µ+1), Eq. 18
///
/// 𝛤(µ+1) = µ 𝛤(µ)
///
/// For Qᵤ(x,y) recursion is stable in the forward direction so we start from
/// i = 0. There is a risk of underflow in the first terms.
///
/// "Computation of the Marcum Q Function", Gil, Segura, Temme 2013 §3
fileprivate func q_series(µ: Double, x: Double, y: Double) -> Double {
    // central gamma
    // FIXME: handle underflow of initial terms
    let Qᵤ = q_gamma(µ, y)
    
    // seed the modifier term: y^(µ-1) e^(-y) / 𝛤(µ)
    let d₀ = pow(y, µ - 1) * exp(-y) / tgamma(µ)

    // calculate the sum deriving the terms recursively
    let s = recursiveSum(indices: 1..., sum0: Qᵤ, state0: (Q: Qᵤ, d: d₀, p: 1.0), update: { iInt, state in
        let (Qᵤ₊ᵢ₋₁, dᵢ₋₁, pᵢ₋₁) = state
        let i = Double(iInt)
        
        // dᵢ = y^(µ + i - 1) e^(-y) / 𝛤(µ + i)
        //    = dᵢ₋₁ y / (µ + i - 1)
        let dᵢ = dᵢ₋₁ * y / (µ + i - 1)
        
        // pᵢ = xⁱ / i!
        //    = pᵢ₋₁ x / i
        let pᵢ = pᵢ₋₁ * x / i
        
        let Qᵤ₊ᵢ = Qᵤ₊ᵢ₋₁ + dᵢ
        let tᵢ = pᵢ * Qᵤ₊ᵢ
        
        return (tᵢ,(Qᵤ₊ᵢ,dᵢ,pᵢ))
    }, until: { a, b in abs(b.1 / b.0) < 1e-10 })
    return exp(-x) * s
}

/// Series method lower tail
///
/// Pᵤ(x,y) = e^(-x) Σi=0..∞ xⁱ / i! Pᵤ₊ᵢ(y), Eq. 7
///
/// We use the following recurrence relationships:
///
/// Pᵤ(y) = Pᵤ₊₁(y) + y^µ e^(-y) / 𝛤(µ+1), Eq. 18
///
/// Recursion is stable in the backwards direction so we must determine
/// a starting point such that truncation error, ε, is acceptable. This is
/// accomplished by finding the root f(n) = 0 of the following:
///
/// f(n) = (n + µ) log(n + µ) + n log(n) - 2n - n log(xy) - C, Eq. 25
///
/// C = log(𝛤(µ) / (2πε), Eq. 26
///
/// "Computation of the Marcum Q Function", Gil, Segura, Temme 2013 §3
fileprivate func p_series(µ: Double, x: Double, y: Double) -> Double {
    
    // first find the starting point for backwards recursion, n₀, by finding
    // the root of the following function:
    // f(n) = (n + µ) log(n + µ) + n log(n) - 2n - n log(xy) - C, Eq. 25
    // C = log(𝛤(µ) / (2πε), Eq. 26
    let ε = 1e-15
    let C = lgamma(µ) - log(2 * .pi * ε) + µ
    let f = { (n: Double) -> Double in
        return (n + µ) * log(n + µ) + n * log(n) - 2 * n - n * log(x * y) - C
    }
    
    // the first derivative of f is:
    // fʹ(n) = log(n (n + µ) / (xy)), Eq. 27
    let fʹ = { (n: Double) -> Double in
        return log(n) + log(n + µ) - log(x * y)
    }

    // a safe starting point is to the right of:
    // (-µ + √(µ²+ 4xy)) / 2, Eq. 27ff
    let guess = 1 + (-µ + sqrt(µ^^2 + 4 * x * y)) / 2
    
    // newton's method
    let n̂ = root(guess: guess, f: f, f1: fʹ)
    let n₀ = ceil(n̂)
    let n₀Int = Int(n₀)
    
    // first term of the sum is xⁱ / i! Pᵤ₊ᵢ(y), i = n₀
    let P₀ = p_gamma(µ + n₀, y)
    let p₀ = x^^n₀Int / tgamma(n₀ + 1)
    let d₀ = pow(y,µ + n₀) * exp(-y) / tgamma(µ + n₀ + 1)
    let t₀ = p₀ * P₀
    
    let s = recursiveSum(indices: (1...n₀Int).reversed(), sum0: t₀, state0: (P: P₀, d: d₀, p: p₀), update: { iInt, state in
        let (Pᵤ₊ᵢ₊₁, dᵢ₊₁, pᵢ₊₁) = state
        let i = Double(iInt)

        // dᵢ = y^(µ + i) e^(-y) / 𝛤(µ + i + 1)
        //    = dᵢ₊₁ (µ + i) / y
        let dᵢ = dᵢ₊₁ * (µ + i) / y
        
        // pᵢ = xⁱ / i!
        //    = pᵢ₊₁ x / i
        let pᵢ = pᵢ₊₁ * i / x
        
        let Pᵤ₊ᵢ = Pᵤ₊ᵢ₊₁ + dᵢ
        let tᵢ = pᵢ * Pᵤ₊ᵢ
        
        return (tᵢ,(Pᵤ₊ᵢ,dᵢ,pᵢ))
    }, until: { a, b in abs(b.1 / b.0) < 1e-10 })
    return exp(-x) * s
}

/// Use recursion to go from a point where quadrature works well to desired point
///
/// solve y = x + µ + √(4x + 2µ) for µ to get the point where quadrature works:
///
/// µʹ = y - x + 1 ± √(2x + 2y + 1)
///
/// Q is stable for forward recursion so choose the lower root and use the following
/// recursion rule to walk back up:
///
/// qᵤ₊₁ = (1 + cᵤ) qᵤ - cᵤ qᵤ₋₁, (Eq. 14)
///
/// cᵤ = √(y/x) Iᵤ(ξ) / Iᵤ₋₁(ξ),
///
/// Iν(x) / Iν₋₁(x) = 1 / 2ν x⁻¹ + 1 / 2(ν + 1) x⁻¹  + 1 / 2(ν + 2) x⁻¹ + ..., NIST 10.33.1
///
/// Note that the second term in the continued fraction for cᵤ is the continued
/// fraction for cᵤ₊₁, so we also have the following recursion for the coefficients:
///
/// cᵤ = √(y/x) / [2µ (2√xy)⁻¹ + √(x/y) cᵤ₊₁]
///
///  = y / [µ + x cᵤ₊₁]
///
/// Or going forwards:
///
/// cᵤ₊₁ = (y / cᵤ - µ) / x
///
/// "Computation of the Marcum Q Function", Gil, Segura, Temme 2013 §2.2
///
/// - Parameters:
///   - µ: the order of the non-central gamma function
///   - x: the non-centrality parameter
///   - y: the point at which to evaluate
///
/// - Returns: The upper tail CDF of the distribution
fileprivate func q_recursion(µ: Double, x: Double, y: Double) -> Double {
    
    // find the nearest µ that fits the criteria for using quadrature
    let µʹ = y - x + 1 - sqrt(2 * (x + y) + 1)
    
    // for recursion we need to move an integer distance
    let µ̃ = µ - ceil(µ - µʹ)
    let n = Int(ceil(µ - µʹ))
    
    // use quadrature twice to seed the three term recursion in Eq. 14
    let q₀  = quadrature(µ: µ̃    , x: x, y: y).q
    let q₋₁ = quadrature(µ: µ̃ - 1, x: x, y: y).q
    
    // find the coefficient for the first step
    let ξ = 2 * sqrt(x * y)
    let cµ̃ = sqrt(y / x) * continued_fraction(b0: 0, a: { _ in 1.0 }, b: { 2 * (µ̃ + Double($0)) / ξ })
    
    // recurse back up to the original µ
    let recurse = (0..<n).reduce((qᵢ₋₁: q₋₁, qᵢ: q₀, cᵢ: cµ̃)) { stateᵢ, iInt in
        let (qᵢ₋₁, qᵢ, cᵢ) = stateᵢ
        let i = Double(iInt)
        let cᵢ₊₁ = (y / cᵢ - (µ̃ + i)) / x
        let qᵢ₊₁ = (1 + cᵢ) * qᵢ - cᵢ * qᵢ₋₁
        return (qᵢ₋₁: qᵢ, qᵢ: qᵢ₊₁, cᵢ: cᵢ₊₁)
    }
    return recurse.qᵢ
}

/// Use recursion to go from a point where quadrature works well to desired point
///
/// solve y = x + µ + √(4x + 2µ) for µ to get the point where quadrature works:
///
/// µʹ = y - x + 1 ± √(2x + 2y + 1)
///
/// P is stable for forward recursion so choose the higher root and use the following
/// recursion rule to walk back up:
///
/// pᵤ₋₁ = ((1 + cᵤ) pᵤ - pᵤ₊₁) / cᵤ , (Eq. 14)
///
/// cᵤ = √(y/x) Iᵤ(ξ) / Iᵤ₋₁(ξ),
///
/// Iν(x) / Iν₋₁(x) = 1 / 2ν x⁻¹ + 1 / 2(ν + 1) x⁻¹  + 1 / 2(ν + 2) x⁻¹ + ..., NIST 10.33.1
///
/// Note that the second term in the continued fraction for cᵤ is the continued
/// fraction for cᵤ₊₁, so we also have the following recursion for the coefficients:
///
/// cᵤ = √(y/x) / [2µ (2√xy)⁻¹ + √(x/y) cᵤ₊₁]
///
///  = y / [µ + x cᵤ₊₁]
///
/// Or going forwards:
///
/// cᵤ₊₁ = (y / cᵤ - µ) / x
///
/// "Computation of the Marcum Q Function", Gil, Segura, Temme 2013 §2.2
///
/// - Parameters:
///   - µ: the order of the non-central gamma function
///   - x: the non-centrality parameter
///   - y: the point at which to evaluate
///
/// - Returns: The lower tail CDF of the distribution
fileprivate func p_recursion(µ: Double, x: Double, y: Double) -> Double {
    
    // find the nearest µ that fits the criteria for using quadrature
    let µʹ = y - x + 1 + sqrt(2 * (x + y) + 1)
    
    // for recursion we need to move an integer distance
    let µ̃ = µ + ceil(µʹ - µ)
    let n = Int(ceil(µʹ - µ))
    
    // use quadrature twice to seed the three term recursion in Eq. 14
    let p₀  = quadrature(µ: µ̃    , x: x, y: y).p
    let p₊₁ = quadrature(µ: µ̃ + 1, x: x, y: y).p
    
    // find the coefficient for the first step
    let ξ = 2 * sqrt(x * y)
    let cµ̃ = sqrt(y / x) * continued_fraction(b0: 0, a: { _ in 1.0 }, b: { 2 * (µ̃ + Double($0)) / ξ })
    
    // recurse back up to the original µ
    let recurse = (0..<n).reduce((pᵢ₊₁: p₊₁, pᵢ: p₀, cᵢ: cµ̃)) { stateᵢ, iInt in
        let (pᵢ₊₁, pᵢ, cᵢ) = stateᵢ
        let i = Double(iInt)
        let cᵢ₋₁ = y / (µ̃ - i - 1 + x * cᵢ)
        let pᵢ₋₁ = ((1 + cᵢ) * pᵢ - pᵢ₊₁) / cᵢ
        return (pᵢ₊₁: pᵢ, pᵢ: pᵢ₋₁, cᵢ: cᵢ₋₁)
    }
    return recurse.pᵢ
}

/// Quadrature method
///
/// We use a polar coordinates representation of the integral:
///
/// Qᵤ(µx,µy) = e^(-1/2µζ²)/2π ∫-π..π e^(µψ(θ)) f(θ) dθ, Eq. 3.22
///
/// -1/2ζ² = -x - y + φ(s₀), Eq. 3.21
///
/// φ(s) = x/s + ys - ln s, Eq. 3.11
///
/// s₀ = (1 + √(1 + ξ²)) / (2y) Eq. 3.13
///
/// For more detail on the integrand see the documentation of the
/// integrand function. Note that the integrand is symmetric around 0
/// so we may integrate over half the interval and double it.
///
/// "Recent software developments for special functions in the
/// Santander-Amsterdam project", Gil, Segura, Temme 2014
///
/// - Parameters:
///   - µ: the order of the non-central gamma function
///   - x: the non-centrality parameter. don't scale by µ, this is handled internally
///   - y: the point at which to evaluate. don't scale by µ, this is handled internally
///
/// - Returns: A tuple of the lower (p) and upper (q) CDF tails of the distribution
fileprivate func quadrature(µ: Double, x µx: Double, y µy: Double) -> Probability {
    // get unscaled parameters
    let x = µx / µ
    let y = µy / µ
    
    // precompute some constants in the integrand. we only use the square of ξ
    let ξ² = 4 * x * y
    let sq1pξ² = sqrt(1 + ξ²)
    
    // integrate over the half interval (0...π)
    let integral = trapezoidalQuadrature(range: 0...(Double.pi)) { integrand(θ: $0, µ: µ, y: y, ξ²: ξ², sq1pξ²: sq1pξ²)}
    
    // multiply by prefix to get result
    let ζ = zeta(x: x, y: y)
    let pq = exp(-0.5 * µ * ζ^^2) * integral / .pi
    
    // sign of ζ determines if we've got p or q
    return ζ < 0 ? Probability(q: pq) : Probability(p: -pq)
}

/// Integrand for the quadrature method of the Marcum Q method
///
/// This function represents the integrand in:
///
/// Qᵤ(µx,µy) = e^(-1/2µζ²)/2π ∫-π..π e^(µψ(θ)) f(θ) dθ, Eq. 3.22
///
/// f(θ) = [sinθ rʹ(θ) + (cosθ - r(θ)) r(θ)] / [r²(θ) - 2r(θ)cosθ + 1], Eq. 3.20
///
/// ψ(θ) = cosθ ρ(θ,ξ) - √(1 + ξ²) - ln[(θ/sinθ + ρ(θ,ξ)) / (1 + √(1 + ξ²))], Eq. 3.23
///
/// r(θ) = 1/(2y) (θ/sinθ + ρ(θ,ξ)), Eq. 3.16
///
/// rʹ(θ) = (sinθ -θcosθ) / (2ysin²θ) [1 + θ / (ρ(θ,ξ)sinθ)), Eq. 3.30
///
/// ρ(θ,ξ) = √((θ/sinθ)² + ξ²), Eq. 3.16
///
/// We decompose ψ(θ) into two terms and attack the computation of each according
/// to guidance in the cited paper. See the comments within the function for further
/// derivations.
///
/// We have special handling for the case of θ = 0 to avoid dividing by zero.
///
/// "Recent software developments for special functions in the
/// Santander-Amsterdam project", Gil, Segura, Temme 2014
fileprivate func integrand(θ: Double, µ: Double, y: Double, ξ²: Double, sq1pξ²: Double) -> Double {
    // Special handling of the θ = 0 case
    if θ == 0 {
        // ρ(0,ξ) = √(1 + ξ²),
        // r(0) = (1 + √(1 + ξ²)) / (2y),
        let r = (1 + sq1pξ²) / (2 * y)

        // f(0) = (1 - r) r / (r² - 2r + 1),
        //      = (1 - r) r / (1 - r)²,
        //      = r / (1 - r),
        // ψ(0) = √(1 + ξ²) - √(1 + ξ²) + ln[(1 + √(1 + ξ²)) / (1 + √(1 + ξ²))],
        //      = 0 + ln(1) = 0
        return r / (1 - r)
    }
    
    // precompute some quantities we'll use multiple times
    let θ² = θ^^2
    let sinθ = sin(θ)
    let sin²θ = sinθ^^2
    let cosθ = cos(θ)
    let θ_sinθ = θ / sinθ
    let θmsinθ = xmsin(θ)
    
    // compute the named quantities
    let ρ = sqrt(θ_sinθ^^2 + ξ²)
    let r = (θ_sinθ + ρ) / (2 * y)
    let rʹ = (-θmsinθ + 2.0 * θ * sin(0.5 * θ)^^2) / (2.0 * y * sin²θ) * (1.0 + θ_sinθ / ρ)
    let f = (sinθ * rʹ + (cosθ - r) * r) / (r * (r - 2 * cosθ) + 1)

    // ψ₁ = cos(θ)ρ(θ,ξ) - √(1 + ξ²),
    //    = [(θ / sin(θ))² - 1 - ξ² sin²θ] / [cos(θ)ρ(θ,ξ) + √(1 + ξ²)], Eq. 3.25
    // (θ / sin(θ))² - 1 = (θ - sin(θ)) (θ + sin(θ)) / sin²θ, Eq. 3.26
    let ψ₁ = (θmsinθ * (θ + sinθ) / sin²θ - θ² - ξ² * sin²θ) / (cosθ * ρ + sq1pξ²)

    // ψ₂ = -log( (θ / sin(θ) + ρ(θ,ξ)) / (1 +  √(1 + ξ²)) ),
    //    = -log(1 + z), Eq. 3.29
    // z  = (θ - sin(θ)) / sin(θ) / (1 + √(1 + ξ²)) [1 + (θ/sin(θ) + 1) / (ρ(θ,ξ) + √(1 + ξ²))]
    let ψ₂ = -log1p((θmsinθ / sinθ) / (1 + sq1pξ²) * (1 + (θ_sinθ + 1) / (ρ + sq1pξ²)))
    let ψ = ψ₁ + ψ₂

    // final result is e^(µψ(θ)) f(θ)
    return exp(µ * ψ) * f
}

/// Asymptotic for big ξ
///
/// We use the following expansions from the cited paper.
///
/// Qᵤ(x,y) = Σi=0..∞ ψᵢ, y > x
///
/// Pᵤ(x,y) = Σi=0..∞ ψ̃ᵢ, y < x, Eq. 37
///
/// ψᵢ = ρ^µ / (2√(2π)) (-1)ⁱ [Aᵢ(µ-1) + 1 / ρ Aᵢ(µ)] 𝜙ᵢ, Eq. 38
///
/// ψ̃ᵢ = -ψᵢ, i > 0,
///
/// ψ̃₀ = ρ^µ / (2√ρ) erfc(√(σξ)), Eq. 40
///
/// ξ = 2 √xy, σ = (√y - √x)² / ξ, ρ = √(y / x), Eq. 31
///
/// Where Aᵢ(µ) is defined as:
///
/// Aᵢ(µ) = 2^-i / i! Γ(1/2 + µ + i) / Γ(1/2 + µ - i), i = 0,1,2,...
///
/// And 𝜙ᵢ is defined as follows as an incomplete gamma function:
///
/// 𝜙ᵢ = σ^(i - 1/2) 𝛤(1/2 - i,σξ), Eq. 34
///
/// For further derivations see comments within the function.
///
/// "Computation of the Marcum Q Function", Gil, Segura, Temme 2013
///
/// - Parameters:
///   - µ: the order of the non-central gamma function
///   - x: the non-centrality parameter
///   - y: the point at which to evaluate
///
/// - Returns: A tuple of the lower (p) and upper (q) CDF tails of the distribution
fileprivate func bigxy(µ: Double, x: Double, y: Double) -> Probability {
    // ξ = 2 √xy           (xi)
    // σ = (√y - √x)² / ξ  (sigma)
    // ρ = √(y / x)        (rho)    Eq. 31
    let ξ = 2 * sqrt(x * y)
    let sqξ = sqrt(ξ)
    let sqσξ = abs(sqrt(y) - sqrt(x))
    let σξ = sqσξ^^2
    let ρ = sqrt(y / x)
    let ρµ = pow(ρ,µ)
    guard ρµ > 0 else {
        print("MarcumQ big xy underflow")
        return .nan
    }
    guard ρµ.isFinite else {
        print("MarcumQ big xy overflow")
        return .nan
    }
    
    // precompute some terms that are constant in the sum
    let fourµ² = 4.0 * µ^^2
    let iInvariantNum = (2 * µ - 1) * (ρ - 1)
    let iInvariantDenom = ρ * (2 * µ - 1)

    // Constant prefix for all terms. note that we have the constant e^(-σξ) here
    // because we factor it out of 𝜙ᵢ.
    // prefix = ρ^µ / 2√2π e^(-σξ)
    let prefix = ρµ / (2 * sqrt(2 * .pi)) * exp(-σξ)

    // A₀(µ) = 1 Γ(1/2 + µ + 0) / 1 Γ(1/2 + µ - 0)
    //       = 1
    let A₀ = 1.0
    
    // C₀(µ) = A₀(µ - 1) - 1 / ρ A₀(µ)
    //       = 1 - 1 / ρ
    let C₀ = 1 - 1 / ρ

    // We factor out from 𝜙ᵢ the constant term e^(-σξ) and the term ξ⁻ⁱ⁺¹/²
    // This makes the recursion step from 𝜙ᵢ₋₁ to 𝜙ᵢ simple
    // 𝜙₀ = √(π/σ) erfc(√y - √x)
    //    = √(π/σ) erfc(√(σξ))
    // 𝜙₀ e^(σξ) ξ⁰⁻¹/² = √(π/σξ) erfc(√(σξ)) e^(σξ)
    let 𝜙₀ = sqrt(.pi / σξ) * erfc(sqσξ) * exp(σξ)

    // seed the sign appropriately depending on whether we're doing P or Q
    let sgn₀ = y >= x ? 1.0 : -1.0

    // The first term depends on whether we are calculating P or Q
    // These have all the terms back in except the constant e^(-σξ)
    let ψ₀ = y >= x ?
        // ψ₀ = ρ^µ / 2√(2π) [1 + 1 / ρ] √(π/σ) erfc(√(σξ))
        //    = ρ^µ / 2√(2σ) [1 + 1 / ρ] erfc(√(σξ))
        prefix * sqξ * C₀ * 𝜙₀ :
        // ψ̃₀ = ρ^µ / 2√ρ erfc(√(σξ)), Eq. 40
        0.5 * ρµ / sqrt(ρ) * erfc(sqσξ) // * exp(σξ)
    
    // Now compute the sum ψ₀ + Σi=1..∞ ψᵢ
    let pq = recursiveSum(indices: 1..., sum0: ψ₀, state0: (Aᵢ: A₀, 𝜙ᵢ: 𝜙₀, ξ⁻ⁱsqξ: sqξ, sgnᵢ: sgn₀), update: { iInt, stateᵢ₋₁ in
        let (Aᵢ₋₁,𝜙ᵢ₋₁,ξ⁻ⁱ⁺¹sqξ,sgnᵢ₋₁) = stateᵢ₋₁
        let i = Double(iInt)

        // Aᵢ(µ) = 2^-i / i! Γ(1/2 + µ + i) / Γ(1/2 + µ - i), i = 0,1,2,...
        //       = Aᵢ₋₁(µ) * 1 / 2i * (1/2 + µ + i - 1) (1/2 + µ - i)
        //       = Aᵢ₋₁(µ) * 1 / 8i * (2µ + 2i - 1)     (2µ - 2i + 1)
        //       = Aᵢ₋₁(µ) * 1 / 8i * (4µ² - 4i² + 4i - 1)
        //       = Aᵢ₋₁(µ) * 1 / 8i * (4µ² - (2i - 1)²)
        let Aᵢ = Aᵢ₋₁ / (8 * i) * (fourµ² - (2 * i - 1)^^2)

        // Aᵢ(µ-1) = 2^-i / i! Γ(1/2 + µ - 1 + i) / Γ(1/2 + µ - 1 - i)
        //         = Aᵢ(µ) * (1/2 + µ - 1 - i) / (1/2 + µ - 1 - i)
        //         = Aᵢ(µ) * (2µ - 2i - 1) / (2µ + 2i - 1)
        // Cᵢ(µ) = Aᵢ(µ-1) - 1 / ρ Aᵢ(µ)
        //       = Aᵢ(µ) [(2µ - 2i - 1) / (2µ + 2i - 1) - 1 / ρ]
        //       = Aᵢ(µ) (2µρ - 2iρ - ρ - 2µ - 2i + 1) / ρ(2µ + 2i - 1)
        //       = Aᵢ(µ) [(2µ - 1) (ρ - 1) - 2i(ρ + 1)] / ρ(2µ + 2i - 1)
        let Cᵢ = Aᵢ * (iInvariantNum - i * 2 * (ρ + 1)) / (iInvariantDenom + i * 2 * ρ)

        // ξ⁻ⁱ⁺¹/² = ξ⁻ⁱ⁺¹⁺¹/² / ξ
        let ξ⁻ⁱsqξ = ξ⁻ⁱ⁺¹sqξ / ξ
        
        // Use the following recursive relation for 𝜙ᵢ:
        // (i - 1/2) 𝜙ᵢ = -σ 𝜙ᵢ₋₁ + e^(-σξ) ξ⁻ⁱ⁺¹/², Eq. 36
        // After factoring out e^(-σξ) and ξ⁻ⁱ⁺¹/²:
        // 𝜙ᵢ e^(σξ) ξⁱ⁻¹/² = (-σξ 𝜙ᵢ₋₁ + 1) / (i - 1/2)
        let 𝜙ᵢ = (-σξ * 𝜙ᵢ₋₁ + 1) / (i - 0.5)
        
        // (-1)ⁱ = -(-1)ⁱ⁻¹
        let sgnᵢ = -sgnᵢ₋₁
        
        // calculate final term, ψᵢ, with ξ⁻ⁱ⁺¹/² multiplied back in:
        // ψᵢ = ρ^µ / 2√(2π) (-1)ⁱ [Aᵢ(µ-1) + 1 / ρ Aᵢ(µ)] 𝜙ᵢ
        // ψᵢ e^(σξ) = ρ^µ / 2√(2π) (-1)ⁱ Cᵢ(µ) ξ⁻ⁱ⁺¹/² [𝜙ᵢ e^(σξ) ξⁱ⁻¹/²]
        let ψᵢ = prefix * sgnᵢ * ξ⁻ⁱsqξ * Cᵢ * 𝜙ᵢ
        return (ψᵢ, (Aᵢ,𝜙ᵢ,ξ⁻ⁱsqξ,sgnᵢ))
    }, until: { a, b in abs(b.1 / b.0) < 1e-10 })
    
    // We calculated either p or q depending on whether y > x
    return Probability(value: pq, isComplement: y >= x)
}

/// Asymptotic for big µ
///
/// Qᵤ₊₁(µx,µy) = √(µ / 2π) Σ i=0... Bᵢ, eq. 71
///
/// Bᵢ = Σ j=0...i fⱼ ᵢ₋ⱼ 𝛹ⱼ(ζ) / µⁱ⁻j, eq. 71
///
/// Pᵤ₊₁(µx,µy) = √(µ / 2π) Σ i=0... B*ᵢ, eq. 78
///
/// B*ᵢ = Σ j=0...i (-1)^j fⱼ ᵢ₋ⱼ 𝛹ⱼ(-ζ) / µⁱ⁻j, eq. 79
///
/// Compute P when y < x + µ
///
/// "Computation of the Marcum Q Function", Gil, Segura, Temme 2013
///
/// - Parameters:
///   - µ: the order of the non-central gamma function
///   - x: the non-centrality parameter. don't scale by µ, this is handled internally
///   - y: the point at which to evaluate. don't scale by µ, this is handled internally
///
/// - Returns: A tuple of the lower (p) and upper (q) CDF tails of the distribution
fileprivate func bigmu(µ µp1: Double, x µx: Double, y µy: Double) -> Probability {
    let µ = µp1 - 1
    let x = µx / µ
    let y = µy / µ
    
    // ζ = sign(x + 1 - y) √(2𝜙(ξ) - 2𝜙(z₀)), eq. 56
    let ζ = zeta(x: x, y: y)
    let sgn = ζ < 0 ? 1.0 : -1.0
    let sgnζ = sgn * ζ
    
    // ψ₀(ζ) = √(π / 2µ) erfc(-ζ √(µ/2)) eq. 67
    let ψ₀ = sqrt(.pi / (2 * µ)) * erfc(-sgnζ * sqrt(µ / 2))
    
    // e^(-1/2 µ ζ²)
    let ehµζ² = exp(-0.5 * µ * sgnζ^^2)
    guard ehµζ² > 0 else {
        print("Marcum Q big µ method underflow")
        return .nan
    }
    
    // ψᵢ(ζ) = (i - 1) / µ ψᵢ₋₂ + (-ζ)ⁱ⁻¹ / µ e^(-1/2 µ ζ²), eq. 68
    // both this and µⁿ below may waste effort if the main series converges
    // before we use some of the calculated items
    let ψⱼ = (1...3).scan((ψᵢ₋₂: 0.0, ψᵢ₋₁: ψ₀, ζⁱ⁻¹: 1.0)) { prev, i in
        let (ψᵢ₋₂,ψᵢ₋₁,ζⁱ⁻¹) = prev
        let ψᵢ = (Double(i - 1) * ψᵢ₋₂ + ζⁱ⁻¹ * ehµζ²) / µ
        let ζⁱ = -sgnζ * ζⁱ⁻¹
        return (ψᵢ₋₁,ψᵢ,ζⁱ)
        }.map { $0.ψᵢ₋₁ }
    
    // pre-calculate power of µ so we don't repeat ourselves below
    let µⁿ = (1...3).scan(1) { prev, n in prev * µ }
    
    // u = 1 / √(2x + 1) eq. 88
    let u = 1 / sqrt(2 * x + 1)
    
    let s = recursiveSum(indices: 1...3, sum0: ψ₀, state0: (), update: { i, stateᵢ₋₁ in
        // Bᵢ = Σ j=0...i fⱼ ᵢ₋ⱼ 𝛹ⱼ(ζ) / µⁱ⁻j, eq. 71
        let Bᵢ = (0...i).reduce(0.0) { accum, j in
            let fjij = f(j,i - j,u)
            let tⱼ = sgn^^j * fjij * ψⱼ[j] / µⁿ[i - j]
            return accum + tⱼ
        }
        return (Bᵢ, ())
    }, until: { a, b in abs(b.1 / b.0) < 1e-10 })
    let pq = sqrt(µ / (2 * .pi)) * s
    return Probability(value: pq, isComplement: ζ < 0)
}

fileprivate func f(_ j: Int, _ i: Int, _ u: Double) -> Double {
    let u² = u^^2
    let coef = fjic[ji(j,i)]!
    return u^^(j + 2 * i) * evaluate_polynomial(poly: coef, z: u²)
}

/// Key for fji coefficients dictionary
fileprivate struct ji: Hashable {
    let j: Int
    let i: Int
    init(_ j: Int, _ i: Int) {
        self.j = j
        self.i = i
    }
}

/// ζ (zeta) for large µ expansion
///
/// 1/2 ζ² = x + y - √(1 + 4xy) + log((1 + √(1 + 4xy)) / 2y), eq. 84
///
/// When y is very close to x + 1 we use:
///
/// ζ = (y - x - 1) / √(2x + 1) Σ i=0... cᵢ zⁱ, eq. 85
///
/// z = (y - x - 1) / (2x + 1)²
///
/// "Computation of the Marcum Q Function", Gil, Segura, Temme 2013
fileprivate func zeta(x: Double, y: Double) -> Double {
    let ymxm1 = y - x - 1
    switch abs(ymxm1) {
    case 0.5...:
        let xy4p1 = sqrt(1 + 4 * x * y)
        return -ymxm1.signum * sqrt(2 * (x + y - xy4p1 + log((1 + xy4p1) / (2 * y))))
    case      _:
        let z = ymxm1 / (2 * x + 1)^^2
        let coef = cic.map { evaluate_polynomial(poly: $0, z: x) }
        let s = evaluate_polynomial(poly: coef, z: z)
        return -ymxm1 / sqrt(2 * x + 1) * s
    }
}

// MARK: Coefficients

/// Coefficients for fji, eq. 90
///
/// The uᵢ for Bessel expansion are in oeis: numerators A144617, denomiators A144618
///
/// "Computation of the Marcum Q Function", Gil, Segura, Temme 2013
fileprivate let fjic: [ji: [Double]] = [
    ji(0,0) : [1],
    
    ji(0,1) : [3,0,-5].map { $0 / 24 },
    ji(1,0) : [3,1].map { $0 / 6 },
    
    ji(0,2) : [81,0,-462,0,385].map { $0 / 1152 },
    ji(1,1) : [-9,21,75,-95].map { $0 / 144 },
    ji(2,0) : [-3,0,5].map { $0 / 24 },
    
    ji(0,3) : [30375,0,369603,0,765765,0,-425425].map { $0 / 414720 },
    ji(1,2) : [-729,1053,9702,-11550,-12705,14245].map { $0 / 6912 },
    ji(2,1) : [27,-144,-402,1440,-925].map { $0 / 576 },
    ji(3,0) : [135,-117,-675,625].map { $0 / 2160 },
]

/// Coefficients for cᵢ, eq. 86
///
/// Related to Lambert W -1 branch expansion
///
/// "Computation of the Marcum Q Function", Gil, Segura, Temme 2013
fileprivate let cic: [[Double]] = [
    [1],
    [1,3].map { -$0 / 3 },
    [7,42,72].map { $0 / 36 },
    [73,657,2142,2700].map { -$0 / 540},
    [1331,15972,76356,177552,181440].map { $0 / 12960 }
]
