//
//  Dekker.swift
//  Numerical
//
//  Created by Adam Roberts on 3/22/19.
//

import Foundation

/// Dekker's method of root finding step
///
/// Find the next guess from both bisection and secant methods. Use the
/// secant method's only if it is between our last guess and bisection's.
func dekkerStep(a1: Double, b0: Double, b1: Double, fa1: Double, fb0: Double, fb1: Double) -> Double {
    let m = bisect(a: a1, b: b1)
    let s = fb0 == fb1 ? m : secantStep(x0: b0, x1: b1, y0: fb0, y1: fb1)
    let b2 = (s - b1).sign != (s - m).sign ? s : m
    return b2
}

/// Dekker's method of root finding
///
/// Find the next guess from both bisection and secant methods. Use the
/// secant method's only if it is between our last guess and bisection's.
///
/// "Finding a zero by means of successive linear interpolation", Th. J. Dekker, Constructive Aspects of the Fundamental Theorem of Algebra, 1969
///
/// https://en.wikipedia.org/wiki/Brent%27s_method#Dekker's_method
public func dekkerRoot(bracket: BracketedRootEstimate, tolerance: EqualityTolerance<Double>, intercept: Double, f rawF: @escaping(Double) -> Double) -> BracketedRootResult {
    let f = intercept == 0 ? CountedFunction(f: rawF) : CountedFunction { rawF($0) - intercept }
    let (a,b,fa,fb) = (bracket.a,bracket.b,bracket.fa,bracket.fb)
    
    // for initial state make the bound that is a better estimate our last guess
    let (a1,b1,fa1,fb1) = abs(fb) <= abs(fa) ? (a,b,fa,fb) : (b,a,fb,fa)
    
    let r = sequence(first: (a1: a1, b0: a1, b1: b1, fa1: fa1, fb0: fa1, fb1: fb1)) { arg0 in
        let (a1, b0, b1, fa1, fb0, fb1) = arg0
        let b2 = dekkerStep(a1: a1, b0: b0, b1: b1, fa1: fa1, fb0: fb0, fb1: fb1)
        let fb2 = f(b2)
        let (a2,fa2) = fb2.sign != fa1.sign ? (a1,fa1) : (b1,fb1)
        return abs(fa2) < abs(fb2) ? (a1: b2, b0: b1, b1: a2, fa1: fb2, fb0: fb1, fb1: fa2) : (a1: a2, b0: b1, b1: b2, fa1: fa2, fb0: fb1, fb1: fb2)
    }.until(maxIter: 50) { s2 in s2.a1.isApprox(.maybeZero(s2.b1, trusted: true), tolerance: tolerance) || s2.fb1.isApprox(.zero(scaleRelativeTo: intercept), tolerance: tolerance) }

    guard let res = r else { return .error } // shouldn't happen
    
    let e = BracketedRootEstimate(a: res.result.a1, b: res.result.b1, fa: res.result.fa1, fb: res.result.fb1)
    
    switch res.exitState {
    case .exhaustedInput: return .error // shouldn't happen
    case .exceededMax: return .noConverge(evals: f.count, estimate: e)
    case .converged: return .success(evals: f.count, estimate: e)
    }
}

