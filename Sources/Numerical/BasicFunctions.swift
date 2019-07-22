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

