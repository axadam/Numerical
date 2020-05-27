//
//  Newton.swift
//  Numerical
//
//  Created by Adam Roberts on 4/4/19.
//

import Foundation

/// Newton-Raphson root finding step
///
/// xᵢ₊₁ = xᵢ - f(xᵢ) / f'(xᵢ)
func newtonStep(x0: Double, f: (Double) -> Double, f1: (Double) -> Double) -> Double {
    return x0 - f(x0) / f1(x0)
}

/// Newton-Raphson root finding method
///
/// Given a guess of a root, a function and its derivative finds a
/// root of the function. Optionally accepts bounds and falls back on
/// bisection when the Newton step goes outside the bounds.
///
/// https://en.wikipedia.org/wiki/Newton%27s_method
///
/// Numerical Recipes §9.4
public func newtonRoot(f: @escaping(Double) -> Double, f1: @escaping(Double) -> Double, guess: Double, xmin: Double? = nil, xmax: Double? = nil, max_iter: Int = 10, xtol: Double = 1e-10) -> RootResult {
    return rootHelper(guess: guess, xmin: xmin, xmax: xmax, maxIter: max_iter, xtol: xtol) { x0 in
        return newtonStep(x0: x0, f: f, f1: f1)
    }
}
