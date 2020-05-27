//
//  Bracket.swift
//  Numerical
//
//  Created by Adam Roberts on 4/2/19.
//

import Foundation

/// Bracket the root of a univariate function given an initial guess
///
/// This algorithm works by simply walking in whichever direction is
/// currently closer to zero. It can fail if your guess is not very
/// good, especially if your function is not monotonic.
///
/// There is optional support for bounding to your function's known
/// domain.
///
/// Numerical Recipes §9.1 with additions
public func bracket(f rawF: @escaping (Double) -> Double, guess a: Double, xmin: Double? = nil, xmax: Double? = nil, maxIter: Int = 30, factor: Double = 1.6) -> BracketResult {
    let f = CountedFunction(f: rawF)
    
    // Take our first step from initial guess. Distance determined by factor.
    // If we have upper or lower limits approach them asymptotically
    let b = { () -> Double in
        let factorb = a * factor
        if let xmin = xmin, factorb < a {
            return xmin - (xmin - a) / factor
        } else if let xmax = xmax, factorb > a {
            return xmax - (xmax - a) / factor
        } else {
            return factorb
        }
    }()
    
    // Evaluate at our initial guess and first step
    let fa = f(a)
    let fb = f(b)
    
    // If we already have different signs we're done
    if fa * fb < 0 { return .bracket(evals: f.count, estimate: BracketedRootEstimate(a: a, b: b, fa: fa, fb: fb)) }
    
    // Iteratively continue to take steps in which ever direction seems closer to zero
    // If we don't bracket before max iterations then returns nil
    let r = sequence(first: (a: a, b: b, fa: fa, fb: fb)) { state0 in
        let (a0, b0, fa0, fb0) = state0
        if abs(fa0) < abs(fb0) {
            let a1: Double
            if let xmin = xmin, a0 < b0 {
                a1 = xmin - (xmin - a0) / factor
            } else if let xmax = xmax, a0 > b0 {
                a1 = xmax - (xmax - a0) / factor
            } else {
                a1 = a0 + factor * (a0 - b0)
            }
            let fa1 = f(a1)
            return (a1,b0,fa1,fb0)
        } else {
            let b1: Double
            if let xmin = xmin, b0 < a0 {
                b1 = xmin - (xmin - b0) / factor
            } else if let xmax = xmax, b0 > a0 {
                b1 = xmax - (xmax - b0) / factor
            } else {
                b1 = b0 + factor * (b0 - a0)
            }
            let fb1 = f(b1)
            return (a0,b1,fa0,fb1)
        }
    }.prefix(maxIter).first { $0.fa * $0.fb < 0 }
    guard let rr = r else { return .noBracket(evals: f.count)}
    return .bracket(evals: f.count, estimate: BracketedRootEstimate(a: rr.a, b: rr.b, fa: rr.fa, fb: rr.fb))
}
