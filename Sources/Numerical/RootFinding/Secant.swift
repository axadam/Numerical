//
//  Secant.swift
//  Numerical
//
//  Created by Adam Roberts on 3/21/19.
//

import Foundation

/// Secant Method root finding step
///
/// Our next step is the intercept on the line through (x₀,y₀) and (x₁,y₁)
///
/// The slope of the line is (y₁ - y₀) / (x₁ - x₀) so the equation is:
///
/// y = y₀ + (x - x₀) (y₁ - y₀) / (x₁ - x₀)
///
/// Set y to 0 and solve for x to find the intercept:
///
/// x = x₀ - y₀ (x₁ - x₀) / (y₁ - y₀)
///
/// Numerical Recipes §9.2
func secantStep(x0: Double, x1: Double, y0: Double, y1: Double) -> Double {
    let x = x0 - y0 * (x1 - x0) / (y1 - y0)
    return x
}

/// Secant Method of root finding
///
/// Uses linear interpolation to find root in a <= x <= b. a and b must bracket a root.
/// This method can fail in some pathological cases where f is flat near root.
///
/// Numerical Recipes §9.2
func secantRoot(f: @escaping (Double) -> Double, a: Double, b: Double, fa: Double, fb: Double, epsilon: Double) -> Double {
    let maxIter = 30
    
    // for initial state make the bound that is a better estimate our last guess
    let (x0, x1, y0, y1) = abs(fa) < abs(fb) ? (b, a, fb, fa) : (a, b, fa, fb)
    
    let r = (0..<maxIter).lazy.scan( (state: (x0: x0, x1: x1, y0: y0, y1: y1), guess: x1) ) { arg0, i in
        let (x0, x1, y0, y1) = arg0.state
        let xnew = secantStep(x0: x0, x1: x1, y0: y0, y1: y1)
        let ynew = f(xnew)
        let state1 = (x0: x1, x1: xnew, y0: y1, y1: ynew)
        return (state1, xnew)
        }.converge { s1, s2 in abs(s2.guess - s1.guess) < epsilon || abs(s2.state.y1) < epsilon }
    guard let res = r else { return .nan }
    return res.guess
}

