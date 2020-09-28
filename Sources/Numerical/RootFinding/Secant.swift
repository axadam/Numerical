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
/// https://en.wikipedia.org/wiki/Secant_method
///
/// Numerical Recipes §9.2
func secantStep(x0: Double, x1: Double, y0: Double, y1: Double) -> Double {
    let x = x0 - y0 * (x1 - x0) / (y1 - y0)
    return x
}

/// Secant Method of root finding
///
/// Uses linear interpolation to find root in a ≤ x ≤ b. a and b must bracket a root.
/// This method can fail in some pathological cases where f is flat near root.
///
/// https://en.wikipedia.org/wiki/Secant_method
///
/// Numerical Recipes §9.2
public func secantRoot(bracket: BracketedRootEstimate, tolerance: EqualityTolerance<Double> = .strict, intercept: Double = 0, f rawF: @escaping(Double) -> Double) -> BracketedRootResult {
    let f = intercept == 0 ? CountedFunction(f: rawF) : CountedFunction { rawF($0) - intercept }
    let (a,b,fa,fb) = (bracket.a,bracket.b,bracket.fa,bracket.fb)

    // for initial state make the bound that is a better estimate our last guess
    let (x0, x1, y0, y1) = abs(fa) < abs(fb) ? (b, a, fb, fa) : (a, b, fa, fb)
    
    let r = sequence(first: (x0: x0, x1: x1, y0: y0, y1: y1)) { arg0 in
        let (x0, x1, y0, y1) = arg0
        let xnew = secantStep(x0: x0, x1: x1, y0: y0, y1: y1)
        let ynew = f(xnew)
        return (x0: x1, x1: xnew, y0: y1, y1: ynew)
    }.until(maxIter: 30) { s2 in s2.x0.isApprox(.maybeZero(s2.x1, trusted: true), tolerance: tolerance) || s2.y1.isApprox(.zero(scaleRelativeTo: intercept), tolerance: tolerance) }

    guard let res = r else { return .error } // shouldn't happen
    
    let e = BracketedRootEstimate(a: res.value.x0, b: res.value.x1, fa: res.value.y0, fb: res.value.y1)

    switch res {
    case .exhaustedInput: return .error // shouldn't happen
    case .exceededMax: return .noConverge(evals: f.count, estimate: e)
    case .success: return .success(evals: f.count, estimate: e)
    }
}

