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
public func bracket(f: @escaping (Double) -> Double, guess: Double, xmin: Double? = nil, xmax: Double? = nil, maxIter: Int = 30) -> (a: Double, b: Double, fa: Double, fb: Double)? {
    let factor = 1.6
    let a = guess
    let fa = f(a)
    let b = guess * factor
    let fb = f(b)
    if fa * fb < 0 { return (a, b, fa, fb) }
    let r = recursiveSequence(indices: 0..<maxIter, initialState: (a: a, b: b, fa: fa, fb: fb), maxIter: maxIter, update: { i, state0 in
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
    }, until: { s1, s2 in print(s2); return s2.fa * s2.fb < 0 })
    return r
}
