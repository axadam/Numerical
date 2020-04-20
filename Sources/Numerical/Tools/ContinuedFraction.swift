//
//  ContinuedFraction.swift
//  Numerical
//
//  Created by Adam Roberts on 9/4/19.
//

import Foundation

/// Evaluate a continued fraction of the form b₀ + a₁ / (b₁ +) a₂ / (b₂ +) ...
///
/// We use the modified Lentz method as described by Thompson and Barnett.
///
/// "Coulomb and Bessel Functions of Complex Arguments and Order", Thompson and Barnett, 1986
///
/// - Parameters:
///   - b0: The initial term b₀
///   - coeffs: A sequence of tuple pairs of (aᵢ,bᵢ)
///
/// - Returns: The convergent of the fraction (or the last step calculated if it didn't converge).
public func continued_fraction<S: Sequence>(b0: Double, coeffs: S) -> Double where S.Element == (a: Double, b: Double){
    let small = Double.leastNormalMagnitude * 10
    let h₀ = max(small, b0)
    let d₀ = 0.0
    let c₀ = h₀
    let cf = coeffs.lazy.scan((cᵢ₋₁: c₀, dᵢ₋₁: d₀, fracᵢ₋₁: h₀)) { accum, coeffᵢ in
        let (cᵢ₋₁,dᵢ₋₁,fracᵢ₋₁) = accum
        let (aᵢ,bᵢ) = coeffᵢ
        let cᵢ = absmax(small, bᵢ + aᵢ / cᵢ₋₁)
        let dᵢ = 1 / absmax(small, bᵢ + aᵢ * dᵢ₋₁)
        let delta = cᵢ * dᵢ
        let fracᵢ = fracᵢ₋₁ * delta
        return (cᵢ₋₁: cᵢ, dᵢ₋₁: dᵢ, fracᵢ₋₁: fracᵢ)
        }.until { b in abs(b.cᵢ₋₁ * b.dᵢ₋₁ - 1) < 1e-15 }
    guard let cfrac = cf?.result else {
        return .nan
    }
    return cfrac.fracᵢ₋₁
}

/// Evaluate a continued fraction of the form b₀ + a₁ / (b₁ +) a₂ / (b₂ +) ...
///
/// We use the modified Lentz method as described by Thompson and Barnett.
///
/// "Coulomb and Bessel Functions of Complex Arguments and Order", Thompson and Barnett, 1986
///
/// - Parameters:
///   - b0: The initial term b₀
///   - a: The ith numerator term, aᵢ, as a function of i = 1,2,3...
///   - b: The ith denominator term, bᵢ, as a function of i = 1,2,3...
///
/// - Returns: The convergent of the fraction (or the last step calculated if it didn't converge).
public func continued_fraction(b0: Double, a: @escaping (Int) -> (Double), b: @escaping (Int) -> Double) -> Double {
    let seq = (1...).lazy.map { return (a: a($0), b: b($0)) }
    return continued_fraction(b0: b0, coeffs: seq)
}

/// Returns the argument with the largest absolute value
///
/// absmax(x, y) = abs(x) < abs(y) ? y : x
///
/// Convenience function frequently used in Lentz form
/// of continued fractions.
func absmax(_ x: Double, _ y: Double = Double.leastNormalMagnitude) -> Double {
    return abs(x) < abs(y) ? y : x
}
