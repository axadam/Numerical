//
//  Roots.swift
//  Numerical
//
//  Created by Adam Roberts on 11/23/18.
//

import Foundation

/// Type signature for functions providing derivative free root finding within an interval
///
/// - Parameters:
///     - f: function whose root we seek.
///     - a: Lower end of the interval to search in.
///     - b: Upper end of the interval to search in.
///     - fa: Function f evaluated at a, f(a).
///     - fb: Function f evaluated at b, f(b).
///     - epsilon: Convergence criteria in terms of how close to zero we need to get.
public typealias bracketedRoot = (@escaping(Double) -> Double, Double, Double, Double, Double, Double) -> Double

/// Root finder for univariate functions
///
/// Starting from your guess it will attempt to bracket a root, and then
/// use the specified root finding method to find the root. Success can
/// depend on how good your guess is, especially if the function is not
/// monotonic.
///
/// If you can cheaply get one or two derivatives of your function then
/// use the appropriate variant of `root()`.
///
/// Default root finding method is Brent's.
///
/// - Parameters:
///     - guess: Initial guess for root.
///     - xmin: Minimum allowed value for root. Optional.
///     - xmax: Maximum allowed value for root. Optional.
///     - epsilon: Convergence criteria in terms of how close to zero we need to get. Default 1e-10.
///     - method: Root finding method. Defaults to Brent's method.
public func root(guess: Double, xmin: Double? = nil, xmax: Double? = nil, tolerance: Double = 10e-15, bracketFactor: Double = 1.6, method: bracketedRoot = brentRoot, f: (Double) -> Double) -> Double {
    guard let (a, b, fa, fb) = bracket(f: f, guess: guess, xmin: xmin, xmax: xmax, factor: bracketFactor) else { return .nan }
    let r = root(f: f, a: a, b: b, fa: fa, fb: fb, tolerance: tolerance, method: method)
    return r
}

func root(f: (Double) -> Double, a: Double, b: Double, fa: Double, fb: Double, tolerance: Double, method: bracketedRoot) -> Double {
    switch (fa,fb) {
    case (-tolerance...tolerance,_): return a
    case (_,-tolerance...tolerance): return b
    case (..<0,..<0): fallthrough
    case (0...,0...): return .nan
    default: return withoutActuallyEscaping(f) { g in method(g, a, b, fa, fb, tolerance) }
    }
}

/// Root finder for univariate function with one derivative
///
/// Starting from your guess uses Newton's method to find a root of your
/// function. Optionally accepts bounds. Success can depend on having a good
/// enough guess. Can be vulnerable to poorly behaved functions, especially
/// those with flat regions.
///
/// - Parameters:
///     - guess: Initial guess for root.
///     - xmin: Minimum allowed value for root. Optional.
///     - xmax: Maximum allowed value for root. Optional.
///     - maxIter: Maximum iterations allowed. Default 100.
///     - xtol: Convergence criteria in terms how small a step we've made. Default 1e-10.
///     - f: Function in which to find root
///     - f1: First derivative of f, f'
public func root(guess: Double, xmin: Double? = nil, xmax: Double? = nil, maxIter: Int = 10, xtol: Double = 1e-10, f: @escaping(Double) -> Double, f1: @escaping(Double) -> Double) -> Double {
    return newtonRoot(f: f, f1: f1, guess: guess, xmin: xmin, xmax: xmax, max_iter: maxIter, xtol: xtol)
}

/// Root finder for univariate function with two derivatives
///
/// Starting from your guess, uses the specified second order root finding
/// method to find the root of your function. Optionally accepts bounds. Success
/// can depend on having a good enough guess.
///
/// This version takes the ratio of the second derivative to the first as
/// an argument.
///
/// - Parameters:
///     - guess: Initial guess for root.
///     - xmin: Minimum allowed value for root. Optional.
///     - xmax: Maximum allowed value for root. Optional.
///     - maxIter: Maximum iterations allowed. Default 100.
///     - xtol: Convergence criteria in terms how small a step we've made. Default 1e-10.
///     - f: Function in which to find root
///     - f1: First derivative of f, f'
///     - f2f1: Ratio of second derivative of f to the first, f" / f'
public func root(guess: Double, xmin: Double? = nil, xmax: Double? = nil, maxIter: Int = 100, xtol: Double = 1e-10, f: @escaping(Double) -> Double, f1: @escaping(Double) -> Double, f2: @escaping(Double) -> Double) -> Double {
    return halleyRoot(guess: guess, xmin: xmin, xmax: xmax, maxIter: maxIter, xtol: xtol, f: f, f1: f1, f2: f2)
}

/// Root finder for univariate function with two derivatives
///
/// Starting from your guess, uses the specified second order root finding
/// method to find the root of your function. Optionally accepts bounds. Success
/// can depend on having a good enough guess.
///
/// - Parameters:
///     - guess: Initial guess for root.
///     - xmin: Minimum allowed value for root. Optional.
///     - xmax: Maximum allowed value for root. Optional.
///     - maxIter: Maximum iterations allowed. Default 100.
///     - xtol: Convergence criteria in terms how small a step we've made. Default 1e-10.
///     - f: Function in which to find root
///     - f1: First derivative of f, f'
///     - f2: Second derivative of f, f"
public func root(guess: Double, xmin: Double? = nil, xmax: Double? = nil, maxIter: Int = 100, xtol: Double = 1e-10, f: @escaping(Double) -> Double, f1: @escaping(Double) -> Double, f2f1: @escaping(Double) -> Double) -> Double {
    return halleyRoot(guess: guess, xmin: xmin, xmax: xmax, maxIter: maxIter, xtol: xtol, f: f, f1: f1, f2f1: f2f1)
}

/// Helper function for root finding
///
/// Provides bounds with bisection if your step goes out of bounds
func rootHelper(guess: Double, xmin: Double? = nil, xmax: Double? = nil, maxIter: Int = 100, xtol: Double = 1e-10, step: @escaping(Double) -> Double) -> Double {
    let r = sequence(first: guess) { x0 in
        let x1 = step(x0)
        if let xmin = xmin, x1 < xmin { return bisectionStep(a: x0, b: xmin) }
        if let xmax = xmax, x1 > xmax { return bisectionStep(a: x0, b: xmax) }
        return x1
    }.until(maxIter: maxIter) { s1, s2 in abs(s2 - s1) < xtol }
    guard let root = r?.result else { return .nan }
    return root
}
