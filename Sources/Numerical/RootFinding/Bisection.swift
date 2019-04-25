//
//  Bisection.swift
//  Numerical
//
//  Created by Adam Roberts on 3/21/19.
//

import Foundation

/// Bisection root finding step
///
/// The next step is simply the mid-point
func bisectionStep(a: Double, b: Double) -> Double {
    let c = (a + b) / 2
    return c
}

/// Bisection method of root finding
///
/// Starting with two points bracketing a root, keep stepping to the mid-point
/// while keeping the root bracketed. Slow but guaranteed to get there.
///
/// Numerical Recipes §9.1.1
public func bisectionRoot(f: (Double) -> Double, a: Double, b: Double, fa: Double, fb: Double, tolerance: Double) -> Double {
    let maxIter: Int = 50
    let q = recursiveSequence(indices: 0...,
                              initialState: (x0: a, x1: b, y0: fa, y1: fb),
                              maxIter: maxIter,
                              update: { i, state0 in
        let (x0, x1, y0, y1) = state0
        let xnew = bisectionStep(a: x0, b: x1)
        let ynew = f(xnew)
        let state1 = ynew.sign == y0.sign ? (xnew, x1, ynew, y1) : (x0, xnew, y0, ynew)
        return state1
    }, until: { s1, s2 in abs(s2.x1 - s2.x0) < tolerance || abs(s2.y1) < tolerance })
    
    guard let res = q else { return .nan }
    
    return res.x1    
}
