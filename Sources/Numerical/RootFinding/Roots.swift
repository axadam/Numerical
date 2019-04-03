//
//  Roots.swift
//  Numerical
//
//  Created by Adam Roberts on 11/23/18.
//

import Foundation

/// Halley's Method to find a root. Builds on Newton's method by adding the
/// second term in the Taylor series.
///
/// Creates a closure to calculate the Halley step and sends to generic
/// root finding function.
///
/// This version uses the ratio of the second derivative to the first.
///
/// x_n+1 = x_n - 2 * f(x_n) * f'(x_n) / (2 * f'(x_n)^2 - f(x_n) * f''(x_n))
///
/// = x_n - 2 * f / (2 * f' - f * f'' / f')
///
/// = x_n - f / f' (1 - f'' f / (f')² / 2)⁻¹
///
/// = x_n - u (1 - t / 2)⁻¹, t = u * f'' / f', u = f / f'
///
/// - Parameters:
///     - guess: Initial guess for root.
///     - xmin: Minimum allowed value for root. Optional.
///     - xmax: Maximum allowed value for root. Optional.
///     - max_iter: Maximum iterations allowed. Optional.
///     - f: Function in which to find root
///     - f1: First derivative of f
///     - f2f1: Ratio of second derivative of f to first derivative of f
public func halley(guess: Double, xmin: Double? = nil, xmax: Double? = nil, max_iter: Int = 100, f: @escaping (Double) -> Double, f1: @escaping (Double) -> Double, f2f1: @escaping (Double) -> Double) -> Double {
    let step: (Double) -> Double = { x in
        let u = f(x) / f1(x)
        let tx = min(1, u * f2f1(x))
        return u / (1 - 0.5 * tx)
    }
    return stepper(prev: guess, max_iter: max_iter, min: xmin, max: xmax, step: step)
}

/// Halley's Method to find a root. Builds on Newton's method by adding the
/// second term in the Taylor series.
///
/// Creates a closure to calculate the Halley step and sends to generic
/// root finding function.
///
/// This version takes the second derivative as an argument.
///
/// x_n+1 = x_n - 2 * f(x_n) * f'(x_n) / (2 * f'(x_n)^2 - f(x_n) * f''(x_n))
///
/// = x_n - 2 * f / (2 * f' - f * f'' / f')
///
/// - Parameters:
///     - guess: Initial guess for root.
///     - xmin: Minimum allowed value for root. Optional.
///     - xmax: Maximum allowed value for root. Optional.
///     - max_iter: Maximum iterations allowed. Optional.
///     - f: Function in which to find root
///     - f1: First derivative of f
///     - f2: Second derivative of f
public func halley(guess: Double, xmin: Double? = nil, xmax: Double? = nil, max_iter: Int = 100, f: @escaping (Double) -> Double, f1: @escaping (Double) -> Double, f2: @escaping (Double) -> Double) -> Double {
    let step: (Double) -> Double = { x in
        let f = f(x)
        let f1 = f1(x)
        return 2 * f / (2 * f1 - f * f2(x) / f1)
    }
    return stepper(prev: guess, max_iter: max_iter, min: xmin, max: xmax, step: step)
}

/// Basic root finding function that lets you supply your own step function
///
/// This provides only the very basics and should only be used when you know
/// your function and your guess well enough. Checks for bounds and falls back on bisection.
///
/// This is a recursive implementation.
func stepper(prev: Double, i: Int = 0, max_iter: Int, min: Double? = nil, max: Double? = nil, step: (Double) -> Double) -> Double {
    if i >= max_iter { return prev }
    let h = step(prev)
    let next: Double = {
        if let min = min, h >= prev - min { return (prev + min) / 2 }
        if let max = max, h <= prev - max { return (prev + max) / 2 }
        return prev - h
    }()
    if abs(h) <= 1e-10 * abs(next) && i > 0 { return next }
    return stepper(prev: next, i: i + 1, max_iter: max_iter, min: min, max: max, step: step)
}

typealias bracketedRoot = (@escaping (Double) -> Double, Double, Double, Double, Double, Double) -> Double

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
func root(f: @escaping (Double) -> Double, guess: Double, method: bracketedRoot = brentRoot) -> Double {
    guard let (a, b, fa, fb) = bracket(f: f, guess: guess) else { return .nan }
    let r = root(f: f, a: a, b: b, fa: fa, fb: fb, method: method)
    return r
}

func root(f: @escaping (Double) -> Double, a: Double, b: Double, fa: Double, fb: Double, epsilon: Double = 1e-10, method: bracketedRoot) -> Double {
    switch (fa,fb) {
    case (-epsilon...epsilon,_): return a
    case (_,-epsilon...epsilon): return b
    case (..<0,..<0): fallthrough
    case (0...,0...): return .nan
    default: return method(f, a, b, fa, fb, epsilon)
    }
}
