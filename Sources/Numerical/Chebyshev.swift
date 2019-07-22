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
/// f(x) = c₀ / 2 + Σ i=1...N-1 cᵢTᵢ(x),
///
/// T₀(x) = 1, T₁(x) = x, Tᵢ₊₁(x) = 2x Tᵢ(x) - Tᵢ₋₁(x)
///
/// Using Clenshaw recursion:
///
/// f(x) = c₀ / 2 + x d₁ - d₂,
///
/// dᵢ = 2x dᵢ₊₁ - dᵢ₊₂ + cᵢ,
///
/// d_n = d_n+1 = 0
///
/// Operates on interval [-1,1] but can be shifted to any other interval [a,b]:
///
/// y = (2x - (b + a)) / (b - a)
public func chebyshev(poly: [Double], z: Double, interval: ClosedRange<Double> = -1...1, m: Int? = nil) -> Double {
    // extract bounds of interval
    let a = interval.lowerBound
    let b = interval.upperBound
    
    // if we have an m argument and it is not more than our number of coefficients
    let mm = min(m ?? poly.count, poly.count)
    
    // Hold aside the first coefficient
    guard let c₀ = poly.first else {
        // we have no coefficients
        return 0
    }
    
    // Make sure we're evaluating at a point in our interval
    guard interval.contains(z) else {
        // evaluating outside of the allowed interval
        return .nan
    }
    
    // Map point in specified interval to a point in [-1,1]
    let y = a == -1 && b == 1 ? z : (2 * z - b - a) / (b - a)
    
    // Go through Clenshaw recursion
    let (d₁,d₂) = poly.prefix(mm).dropFirst().reversed().reduce((0.0,0.0)) { accum, cᵢ in
        let (dᵢ₊₁,dᵢ₊₂) = accum
        let dᵢ = 2 * y * dᵢ₊₁ - dᵢ₊₂ + cᵢ
        return (dᵢ,dᵢ₊₁)
    }
    
    // Evaluate
    print("c₀: \(c₀), d₁: \(d₁), d₂: \(d₂)")
    let f = c₀ / 2 + y * d₁ - d₂
    return f
}
