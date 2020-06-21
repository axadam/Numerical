//
//  Bisection.swift
//  Numerical
//
//  Created by Adam Roberts on 3/21/19.
//

import Foundation

/// Bisect the distance between two points
func bisect(a: Double, b: Double) -> Double {
    let c = (a + b) / 2
    return c
}

extension BracketedRootEstimate {
    /// Bisection root finding step
    ///
    /// We find the midpoint and then pivot so that we still bracket the root
    func bisectionStep(f: CountedFunction<Double,Double>) -> BracketedRootEstimate {
        let xnew = bisect(a: a, b: b)
        let ynew = f(xnew)
        let (a1,b1,fa1,fb1) = ynew.sign == fa.sign ? (xnew, b, ynew, fb) : (a, xnew, fa, ynew)
        return BracketedRootEstimate(a: a1, b: b1, fa: fa1, fb: fb1)
    }
}

/// Bisection method of root finding
///
/// Starting with two points bracketing a root, keep stepping to the mid-point
/// while keeping the root bracketed. Slow but guaranteed to get there.
///
/// https://en.wikipedia.org/wiki/Bisection_method
public func bisectionRoot(bracket: BracketedRootEstimate, tolerance: EqualityTolerance<Double> = .strict, intercept: Double = 0, f rawF: @escaping(Double) -> Double) -> BracketedRootResult {
    let f = intercept == 0 ? CountedFunction(f: rawF) : CountedFunction { rawF($0) - intercept }
    let r = sequence(first: bracket) { state0 in
        return state0.bisectionStep(f: f)
    }.until(maxIter: 50) { s2 in s2.a.isApprox(.maybeZero(s2.b), threshold: tolerance) || s2.fb.isApprox(.zero(scaleRelativeTo: intercept), threshold: tolerance) }
    
    guard let res = r else { return .error } // shouldn't happen
    
    let e = res.result
    
    switch res.exitState {
    case .exhaustedInput: return .error // shouldn't happen
    case .exceededMax: return .noConverge(evals: f.count, estimate: e)
    case .converged: return .success(evals: f.count, estimate: e)
    }
}
