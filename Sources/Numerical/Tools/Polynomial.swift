//
//  Polynomial.swift
//  Numeric
//
//  Created by Adam Roberts on 3/15/19.
//

import Foundation

/// Evaluate a rational expression of two polynomials at point z
///
/// The polynomials are each evaluated using Horner's Method.
///
/// In the case that z > 1 and both polynomials are the same length we use the
/// optimization of dividing both by zⁿ. This equates to evaluating the reverse
/// orders of the coefficients at the point 1 / z.
///
/// - Parameters:
///    - num: An array of coefficients for the polynomial in the numerator. Terms in increasing order.
///    - denom: An array of coefficients for the polynomial in the denominator. Terms in increasing order.
///    - z: the point at which to evalute
public func polynomialRatio(num: [Double], denom: [Double], z: Double) -> Double {
    let matched = num.count == denom.count
    switch (matched,z) {
    case (false,_): fallthrough
    case (true,...1):
        let s1 = polynomial(coeffs: num, z: z)
        let s2 = polynomial(coeffs: denom, z: z)
        return s1 / s2
    case (true,_):
        let s1 = polynomial(coeffs: num.reversed(), z: 1 / z)
        let s2 = polynomial(coeffs: denom.reversed(), z: 1 / z)
        return s1 / s2
    }
}

/// Evaluate the specified polynomial at point z
///
/// This implements Horner's Method which can evalute an n degree polynomial
/// with only n multiplications and n additions.
///
/// p(x) = a₀ + a₁ x + a₂ x² + a₃ x³ + ... + a_n xⁿ
///
/// = a₀ + x (a₁ + x(a₂ + x(a₃ + ... + x(a_n-1 + x a_n)...)))
///
/// - Parameters:
///    - poly: An array of coefficients for the polynomial. Terms in increasing order.
///    - z: The point at which to evalute the polynomial.
public func polynomial(coeffs: [Double], z: Double) -> Double {
    guard let sLast = coeffs.last else { return 0 }
    let rest = coeffs.dropLast()
    let ans = rest.reversed().reduce(sLast, { $0 * z + $1 } )
    return ans
}
