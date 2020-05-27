//
//  Halley.swift
//  Numerical
//
//  Created by Adam Roberts on 4/4/19.
//

import Foundation

/// Halley's root finding step
///
/// Builds on Newton's method by adding the second term in the Taylor series.
///
/// xᵢ₊₁ = xᵢ - 2 f(xᵢ) f'(xᵢ) / (2 f'(xᵢ)² - f(xᵢ) f"(xᵢ))
///
/// = xᵢ - 2 f / (2 f' - f f" / f')
func halleyStep(x0: Double, f: (Double) -> Double, f1: (Double) -> Double, f2: (Double) -> Double) -> Double {
    let f = f(x0)
    let f1 = f1(x0)
    return x0 - 2 * f / (2 * f1 - f * f2(x0) / f1)
}

/// Halley's root finding step
///
/// Builds on Newton's method by adding the second term in the Taylor series.
/// This version takes the ratio of the second derivative to the first
/// as an argument.
///
/// xᵢ₊₁ = xᵢ - 2 f(xᵢ) f'(xᵢ) / (2 f'(xᵢ)² - f(xᵢ) f"(xᵢ))
///
/// = xᵢ - 2 f / (2 f' - f f" / f')
///
/// = xᵢ - f / f' (1 - f" f / (f')² / 2)⁻¹
///
/// = xᵢ - u (1 - t / 2)⁻¹, t = u f" / f', u = f / f'
func halleyStep(x0: Double, f: (Double) -> Double, f1: (Double) -> Double, f2f1: (Double) -> Double) -> Double {
    let u = f(x0) / f1(x0)
    let tx = min(1, u * f2f1(x0))
    return x0 - u / (1 - 0.5 * tx)
}

/// Halley's method of root finding
///
/// Builds on Newton's method by adding the second term of the Taylor
/// series.
///
/// https://en.wikipedia.org/wiki/Halley%27s_method
///
/// Numerical Recipes §9.4.2
public func halleyRoot(guess: Double, xmin: Double? = nil, xmax: Double? = nil, maxIter: Int = 100, xtol: Double = 1e-10, f: @escaping(Double) -> Double, f1: @escaping(Double) -> Double, f2: @escaping(Double) -> Double) -> RootResult {
    return rootHelper(guess: guess, xmin: xmin, xmax: xmax, maxIter: maxIter, xtol: xtol) { x0 in
        halleyStep(x0: x0, f: f, f1: f1, f2: f2)
    }
}

/// Halley's method of root finding
///
/// Builds on Newton's method by adding the second term of the Taylor
/// series. This version takes the ratio of the second derivative to the
/// first as an argument.
///
/// https://en.wikipedia.org/wiki/Halley%27s_method
///
/// Numerical Recipes §9.4.2
public func halleyRoot(guess: Double, xmin: Double? = nil, xmax: Double? = nil, maxIter: Int = 100, xtol: Double = 1e-10, f: @escaping(Double) -> Double, f1: @escaping(Double) -> Double, f2f1: @escaping(Double) -> Double) -> RootResult {
    return rootHelper(guess: guess, xmin: xmin, xmax: xmax, maxIter: maxIter, xtol: xtol) { x0 in
        halleyStep(x0: x0, f: f, f1: f1, f2f1: f2f1)
    }
}
