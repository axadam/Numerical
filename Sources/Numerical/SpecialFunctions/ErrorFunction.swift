//
//  ErrorFunction.swift
//  Numerical
//
//  Created by Adam Roberts on 9/25/18.
//

import Foundation

/// Inverse error function complement
///
/// We use the relation with the quantile of the Normal distribution and a fast
/// approximation thereof to get a guess. We then use Halley's method to get
/// our final value.
///
/// 𝚽⁻¹(p) = √2 erf⁻¹(2p - 1), where 𝚽⁻¹(p) is the inverse Normal CDF
///
/// Numerical Recipes 3rd edition, §6.2.2
public func invErfC(_ p: Double) -> Double {
    switch p {
        // handle domain limits
    case ..<0: return .nan
    case 0: return .infinity
    case 2: return -.infinity
    case 2...: return .nan
        
        // main case
    case _:
        let pp = p < 1 ? p : 2 - p

        // erfc⁻¹(p) = 1/√2 𝚽⁻¹(p/2)
        // 1 / √2 = 0.70711
        let guess = 0.70711 * qapprox(p: 0.5 * pp)

        // Halley method.
        // erf'(x) = -2/√π e^-x²
        // erf''(x) = 4/√π e^-x² x
        // erf''(x) / erf'(x) = -2 x
        // h = 2 f / (-4/√π e^-x² + f 2 x)
        //   = -f / (2/√π e^-x² - f x)
        let x = root(guess: guess,
                        maxIter: 2,
                        f: { x in erfc(x) - pp },
                        f1: { x in -1.12837916709551257 * exp(-x * x) },
                        f2f1: { x in -2 * x }).value
        return p < 1 ? x : -x
    }
}

/// Inverse error function
///
/// Implemented in terms of complement
public func invErf(_ p: Double) -> Double {
    return invErfC(1 - p)
}

/// Normal distribution quantile approximation
///
/// x = t - (a₀ + a₁t) / (1 + b₁t + b₂t²)
///
/// t = √ log p⁻²
///
/// a₀ = 2.30753, a₁ = 0.27061
///
/// b₁ = 0.99229, b₂ = 0.04481
///
/// 0 `<` p ≤ 1/2, accuracy ~ 1e-3
///
/// Handbook of Mathematical Functions, §26.2.22
func qapprox(p: Double) -> Double {
    guard 0 < p && p <= 0.5 else {
        return .nan
    }
    let t = sqrt(-2 * log(p))
    return t - (2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481))
}

