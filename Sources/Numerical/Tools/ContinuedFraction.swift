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
public func continuedFraction<S: Sequence>(b0: Double, coeffs: S, maxIter: Int = 100) -> ContinuedFractionResult where S.Element == (a: Double, b: Double){
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
    }.until(maxIter: maxIter) { b in (b.cᵢ₋₁ * b.dᵢ₋₁).isApprox(.maybeZero(1, trusted: true), tolerance: EqualityTolerance(relative: 8 * Double.ulpOfOne, absolute: 0, absoluteForZero: 8 * Double.ulpOfOne)) }
    guard let cfrac = cf else {
        return .error
    }
    switch cfrac.exitState {
    case .exhaustedInput, .exceededMax: return .noConverge(terms: cfrac.iterations, estimate: cfrac.result.fracᵢ₋₁)
    case .converged: return .success(terms: cfrac.iterations, estimate: cfrac.result.fracᵢ₋₁)
    }
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
public func continuedFraction(b0: Double, a: @escaping (Int) -> (Double), b: @escaping (Int) -> Double, maxIter: Int = 100) -> ContinuedFractionResult {
    let seq = (1...).lazy.map { return (a: a($0), b: b($0)) }
    return continuedFraction(b0: b0, coeffs: seq, maxIter: maxIter)
}

public enum ContinuedFractionResult {
    case error
    case noConverge(terms: Int, estimate: Double)
    case success(terms: Int, estimate: Double)
}

public extension ContinuedFractionResult {
    var value: Double {
        switch self {
        case .error: return .nan
        case .noConverge(_, let e): return e
        case .success(_, let e): return e
        }
    }
    
    var iterations: Int {
        switch self {
        case .error: return 0
        case .noConverge(let n, _): return n
        case .success(let n, _): return n
        }
    }
    
    var converged: Bool {
        switch self {
        case .error, .noConverge: return false
        case .success: return true
        }
    }
}
