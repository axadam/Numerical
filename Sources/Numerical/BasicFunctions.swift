//
//  BasicFunctions.swift
//  Numerical
//
//  Created by Adam Roberts on 7/21/19.
//

import Foundation

/// Calculate e^x - 1 - x with precision even when x is small
///
/// Use built in sinh function and hyperbolic identies:
///
/// y = e^x               - 1 - x
///
///   = cosh(x) + sinh(x) - 1 - x
///
/// Subtracting 1 truncates so substitute for cosh(x) - 1:
///
///   = 2 sinh(x/2)² + sinh(x) - x
///
/// Would be nice to only evaluate sinh once (at x/2):
///
///   = 2 sinh(x/2)² + 2 sinh(x/2) cosh(x/2)         - x
///
///   = 2 sinh(x/2)² + 2 sinh(x/2) √(1 + sinh(x/2)²) - x
public func expm1mx(_ x: Double) -> Double {
    switch abs(x) {
    case 0.95...: return expm1(x) - x
    case _      :
        let shx2 = sinh(x/2)
        let sh²x2 = shx2^^2
        return (2 * sh²x2 + (2 * shx2 * sqrt(1 + sh²x2) - x))
    }
}

/// Calculate log(1 + x) - x with precision even when x is small
///
/// When x is small use alternative form:
///
/// y = -(e^log(1+x) - 1 - log(1+x))
///
///   = -(1 + x - 1 - log(1+x))
///
///   = log(1+x) - x
public func log1pmx(_ x: Double) -> Double {
    switch x {
    case ..<(-1): return .nan
    case -1     : return -.infinity
    case ..<0.95: return -expm1mx(log1p(x))
    case _      : return log1p(x) - x
    }
}

/// Calculate x - sin(x) with precision even when x is small
///
/// sin(x) = Σi=0... (-1)ⁱ x²ⁱ⁺¹ / (2i + 1)!
///
/// x - sin(x) = Σi=1... (-1)ⁱ⁺¹ x²ⁱ⁺¹ / (2i + 1)!
public func xmsin(_ x: Double) -> Double {
    switch x {
    case 1...:
        return x - sin(x)
    case _   :
        let x² = x^^2
        let s = recursiveSum(indices: 1..., sum0: 0.0, state0: -x, update: { i, prev in
            let j = Double(2 * i + 1)
            let t = -prev * x² / (j * (j - 1))
            return (t, t)
        }, until: { a, b in b.1.isApprox(.zero(scaleRelativeTo: b.0), threshold: .strict) })
        return s
    }
}

/// Returns the argument with the largest absolute value
///
/// absmax(x, y) = abs(x) < abs(y) ? y : x
func absmax(_ x: Double, _ y: Double = Double.leastNormalMagnitude) -> Double {
    return abs(x) < abs(y) ? y : x
}

public extension Sequence where Element: FloatingPoint {
    /// Sum a sequence of numbers
    ///
    /// Straightforward calculation of:
    ///
    /// Σ1...n aᵢ
    ///
    /// This method is subject to accumulation of rounding and truncation errors.
    func sum_naive() -> Element { return reduce(Element.zero, +) }
    
    /// Kahan Summation
    ///
    /// A compensated sum. After adding each term we check for error in low order digits
    /// and compensate. S2 is an estimate of the error of the last time S was rounded and
    /// truncated. The Neumaier variant, `sum_kbn`, is more accurate at the same cost
    /// and should be preferred.
    ///
    /// Sᵢ = Sᵢ₋₁ + (S2ᵢ₋₁ + Yᵢ),
    ///
    /// S2ᵢ = (Sᵢ₋₁ - Sᵢ) + (S2ᵢ₋₁ + Yᵢ)
    ///
    /// "Further remarks on reducing truncation errors", Kahan, 1965
    func sum_kahan() -> Element {
        return reduce((s: Element.zero, s2: Element.zero)) { accum, yᵢ in
            let (sᵢ₋₁,s2ᵢ₋₁) = accum
            let s2py = s2ᵢ₋₁ + yᵢ
            let sᵢ = sᵢ₋₁ + s2py
            let s2ᵢ = (sᵢ₋₁ - sᵢ) + s2py
            return (sᵢ,s2ᵢ)
        }.s
    }

    /// Kahan-Babuška-Neumaier Sum
    ///
    /// A variant of compensated sum that attends to the fact that sometimes the loss of
    /// precision is on the side of the sum and sometimes on the side of the term. Compensation
    /// is applied at the end.
    ///
    /// s₀ = 0, w₀ = 0
    ///
    /// sᵢ = aᵢ + sᵢ₋₁,
    ///
    /// wᵢ = wᵢ₋₁ + (aᵢ + (sᵢ₋₁ - sᵢ)), |aᵢ| ≤ |sᵢ₋₁|
    ///
    /// wᵢ = wᵢ₋₁ + (sᵢ₋₁ + (aᵢ - sᵢ)),  |aᵢ| > |sᵢ₋₁|
    ///
    /// Σ1...n aᵢ = s_n + w_n
    ///
    /// "Rundungsfehleranalyse einiger Verfahren zur Summation endlicher Summen", Neumaier, 1974, Eq. 2.IV
    func sum_kbn() -> Element {
        let (s,w) = reduce((s: Element.zero,w: Element.zero)) { accum, aᵢ in
            let (sᵢ₋₁,wᵢ₋₁) = accum
            let sᵢ = sᵢ₋₁ + aᵢ
            let Δ = abs(aᵢ) <= abs(sᵢ₋₁) ? (sᵢ₋₁ - sᵢ) + aᵢ : (aᵢ - sᵢ) + sᵢ₋₁
            let wᵢ = wᵢ₋₁ + Δ
            return (s: sᵢ, w: wᵢ)
        }
        return s + w
    }
}

public extension Collection where Index == Int, Element: FloatingPoint {
    /// Pairwise sum
    ///
    /// A divide and conquer algorithm where we divide the collection in two and then recursively
    /// apply the algorithm to each half. For collections below a certain size we revert to naive
    /// summation to avoid recurisve overhead.
    ///
    /// Pairwise summation is a special case of superblock accumulation where the number of
    /// blocks at each step, b, is equal to 2.
    func sum_pairwise(N: Int = 100) -> Element {
        switch count {
        case ..<N:
            return sum_naive()
        case    _:
            // Integer division to find index of mid-point
            let m = count / 2 + startIndex
            return self[startIndex...m].sum_pairwise(N: N) + self[(m+1)...(endIndex - 1)].sum_pairwise(N: N)
        }
    }
}
