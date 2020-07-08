//
//  Trapezoidal.swift
//  Numerical
//
//  Created by Adam Roberts on 9/6/19.
//

import Foundation

/// Integration of a function by the trapezoidal rule
///
/// Integrates the function, f, over the specified range, [a,b].
///
/// The simplest version of the trapezoidal rule approximates f with a trapezoid
/// with edges at a and b. In that case the area is calculated as:
///
/// ∫a...b f(x) dx ≈ (b - a) (f(a) + f(b)) / 2
///
/// The composite trapezoidal rule works by splitting the range into subintervals. This
/// implementation makes the subintervals at evenly spaced points, x₀, x₁, x₂, ..., xN,
/// each ∆x apart. In this case the area is calculated as:
///
/// ∫a...b f(x) dx ≈ ∆x/2 Σi=1...N f(xᵢ₋₁) + f(xᵢ)
///
/// = ∆x/2 [ f(x₀) + 2 f(x₁) + 2 f(x₂) + ... + f(xN) ]
///
/// = ∆x [ (f(x₀) + f(xN)) / 2 + Σi=1...N-1 f(xᵢ) ]
///
/// This method keeps splitting the intervals in half until the sum converges, doing it
/// in such a way that we don't evaluate the function twice at any one point.
///
/// Note that the trapezoidal rule will overestimate concave up functions and underestimate
/// concave down ones.
///
/// https://en.wikipedia.org/wiki/Trapezoidal_rule
public func trapezoidal(range: ClosedRange<Double>, maxIter: Int = 10, f rawF: @escaping (Double) -> Double) -> QuadratureResult {
    let f = CountedFunction(f: rawF)

    // extract end points from Range
    let a = range.lowerBound
    let b = range.upperBound
    
    // initial interval is just the full interval
    let Δx₀ = (b - a)
    
    // initial estimate
    let I₀ = 0.5 * Δx₀ * (f(a) + f(b))
    
    // iteratively bisect the interval until we get desired convergence
    let q = sequence(first: (I: I₀, Δx: Δx₀, n: 1)) { accum in
        let (Iⱼ₋₁,Δxⱼ₋₁,nⱼ₋₁) = accum
        
        // bisect the interval
        let Δxⱼ = 0.5 * Δxⱼ₋₁
        
        // double the number of intervals
        let nⱼ = nⱼ₋₁ * 2
        
        // find the sequence of xᵢ (i odd), starting at half the previous interval
        // and going until one half the previous short of the end of the range
        let seq = sequence(first: a + Δxⱼ, next: { $0 + Δxⱼ₋₁ }).prefix(nⱼ₋₁)
        
        // sum up evaluations of the function at each point in the sequence
        let sum = seq.map { f($0) }.sumKBN()
        
        // new estimate is half the previous (since now intervals are half as wide)
        // plus ∆x/2 Σf(xᵢ), i odd
        let Iⱼ = 0.5 * Iⱼ₋₁ + sum * Δxⱼ

        return (Iⱼ,Δxⱼ,nⱼ)
    }.until(minIter: 3, maxIter: maxIter) { a, b in b.I.isApprox(.maybeZero(a.I), tolerance: .strict) }
    guard let quad = q else {
        return .error
    }
    switch quad.exitState {
    case .exhaustedInput: return .error
    case .exceededMax: return .noConverge(evals: f.count, estimate: quad.result.I)
    case .converged: return .success(evals: f.count, estimate: quad.result.I)
    }
}
