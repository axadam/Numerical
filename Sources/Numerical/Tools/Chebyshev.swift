//
//  Chebyshev.swift
//  Numerical
//
//  Created by Adam Roberts on 7/19/19.
//

/// Evaluate a Chebyshev polynomial
///
/// Evaluates a Chebyshev polynomial approximation of function f at a specified
/// point x.
///
/// f(x) ≈ p_n(x) = a₀ + a₁T₁(x) + a₂T₂(x) + ... + a_n T_n(x),
///
/// T₀(x) = 1, T₁(x) = x, Tᵢ₊₁(x) = 2x Tᵢ(x) - Tᵢ₋₁(x)
///
/// Using Clenshaw recursion:
///
/// p_n(x) = a₀ + x b₁ - b₂,
///
/// bᵢ = 2x bᵢ₊₁ - bᵢ₊₂ + aᵢ,
///
/// b_n = b_{n+1} = 0
///
/// Operates on interval [-1,1] but can be shifted to any other interval [a,b]:
///
/// y = (2x - (b + a)) / (b - a)
///
/// https://en.wikipedia.org/wiki/Clenshaw_algorithm
public func chebyshev(poly: [Double], z: Double, interval: ClosedRange<Double> = -1...1, m: Int? = nil) -> Double {
    // extract bounds of interval
    let a = interval.lowerBound
    let b = interval.upperBound
    
    // if we have an m argument and it is not more than our number of coefficients
    let mm = min(m ?? poly.count, poly.count)
    
    // Hold aside the first coefficient
    guard let a₀ = poly.first else {
        // we have no coefficients
        return 0
    }
    
    // Make sure we're evaluating at a point in our interval
    guard interval.contains(z) else {
        // evaluating outside of the allowed interval
        return .nan
    }
    
    // Map point in specified interval to a point in [-1,1]
    let x = a == -1 && b == 1 ? z : (2 * z - b - a) / (b - a)
    
    // Go through Clenshaw recursion
    let (b₁,b₂) = poly.prefix(mm).dropFirst().reversed().reduce((0.0,0.0)) { accum, aᵢ in
        let (bᵢ₊₁,bᵢ₊₂) = accum
        let bᵢ = 2 * x * bᵢ₊₁ - bᵢ₊₂ + aᵢ
        return (bᵢ,bᵢ₊₁)
    }
    
    // Final step to evalute polynomial
    let p = a₀ + x * b₁ - b₂
    return p
}
