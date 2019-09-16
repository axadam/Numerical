//
//  Numerical.swift
//  Numerical
//
//  Created by Adam Roberts on 9/25/18.
//

import Foundation

public extension FloatingPoint {
    /// Unit magnitude with same sign as self
    var signum: Self {
        if magnitude == 0 { return 0 }
        return Self(signOf: self, magnitudeOf: 1)
    }
}

/// Round double to integer. Round half away from zero.
public func iround(_ x: Double) -> Int {
    let r = x.rounded(.toNearestOrAwayFromZero)
    let i = Int(r)
    return i
}

/// Round double to integer. Rounds towards zero (floor of magnitude)
public func itrunc(_ x: Double) -> Int {
    let r = x.rounded(.towardZero)
    let i = Int(r)
    return i
}
