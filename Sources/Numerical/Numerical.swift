//
//  Numerical.swift
//  Numerical
//
//  Created by Adam Roberts on 9/25/18.
//

import Foundation

/// Restrict to some minimum absolute value
///
/// absmax(x, min) = abs(x) < min ? min : x
///
/// Convenience function frequently used in Lentz form
/// of continued fractions. Not currently exposed publically
func absmax(_ x: Double, min: Double = Double.leastNormalMagnitude) -> Double {
    return abs(x) < min ? min : x
}

public extension FloatingPoint {
    var signum: Self {
        if magnitude == 0 { return 0 }
        return Self(signOf: self, magnitudeOf: 1)
    }
}
