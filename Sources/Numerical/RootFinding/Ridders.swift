//
//  Ridders.swift
//  Numerical
//
//  Created by Adam Roberts on 3/22/19.
//

import Foundation

/// Ridders' root finding step
///
/// This is essentially false position but on f(a), f(m)e^Q, and f(b)e^2Q
/// where e^Q is the solution to quadratic equation f(a) - 2 f(m) e^Q + f(b) e^2Q = 0,
/// and m is the midpoint of a and b:
///
/// e^Q = (f(m) + sign[f(b)] √(f(m)² - f(a) f(b))) / f(b)
///
/// Yielding the updating formula
///
/// x' = m + (m - a) sign[f(a) - f(b)] f(m) / √(f(m)² - f(a) f(b))
///
/// Numerical Recipes §9.2.1
func riddersStep(f: (Double) -> Double, a: Double, b: Double, fa: Double, fb: Double) -> (x0: Double, x1: Double, y0: Double, y1: Double) {
    // find the mid point and evaluate there
    let m = (a + b) / 2
    let fm = f(m)
    
    // updating formula
    let s = sqrt(fm^^2 - fa * fb)
    let xnew = m + (m - a) * ((fa - fb).signum * fm / s)
    let fnew = f(xnew)

    // Bookkeeping to keep the root bracketed
    let (a1,b1,fa1,fb1): (Double, Double, Double, Double) = {
        if fnew.sign != fm.sign { return (m,xnew,fm,fnew) }
        if fnew.sign != fa.sign { return (a,xnew,fa,fnew) }
        return (xnew,b,fnew,fb)
    }()
    let (a1_,b1_,fa1_,fb1_) = abs(fa1) < abs(fb1) ? (b1,a1,fb1,fa1) : (a1,b1,fa1,fb1)
    return (a1_,b1_,fa1_,fb1_)
}

/// Ridders' root finding method
///
/// Approximates the functions curve with an exponential instead of a
/// polynomial. Robust.
public func riddersRoot(f: (Double) -> Double, a: Double, b: Double, fa: Double, fb: Double, epsilon: Double) -> Double {
    let maxIter = 30
    
    let r = recursiveSequence(indices: 0..<maxIter, initialState: (state: (x0: a, x1: b, y0: fa, y1: fb), guess: b), maxIter: maxIter, update: { i, arg0 in
        let (x0, x1, y0, y1) = arg0.state
        let state1 = riddersStep(f: f, a: x0, b: x1, fa: y0, fb: y1)
        return (state: state1, guess: state1.x1)
    }, until: { s1, s2 in abs(s2.state.y1) < epsilon })
    guard let res = r else { return .nan }
    return res.guess
}
