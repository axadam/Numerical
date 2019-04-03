//
//  Dekker.swift
//  Numerical
//
//  Created by Adam Roberts on 3/22/19.
//

import Foundation

/// Dekker's method of root finding
///
/// Find the next guess from both bisection and secant methods. Use the
/// secant method's only if it is between our last guess and bisection's.
func dekkerStep(a1: Double, b0: Double, b1: Double, fa1: Double, fb0: Double, fb1: Double) -> Double {
    let m = bisectionStep(a: a1, b: b1)
    let s = fb0 == fb1 ? m : secantStep(x0: b0, x1: b1, y0: fb0, y1: fb1)
    let b2 = (s - b1).sign != (s - m).sign ? s : m
    return b2
}

func dekkerRoot(f: @escaping (Double) -> Double, a: Double, b: Double, fa: Double, fb: Double, epsilon: Double) -> Double {
    let maxIter = 50
    
    // for initial state make the bound that is a better estimate our last guess
    let (a1,b1,fa1,fb1) = abs(fb) <= abs(fa) ? (a,b,fa,fb) : (b,a,fb,fa)
    let r = (0..<maxIter).lazy.scan( (state: (a1: a1, b0: a1, b1: b1, fa1: fa1, fb0: fa1, fb1: fb1), guess: fb1) ) { arg0, i in
        let (a1, b0, b1, fa1, fb0, fb1) = arg0.state
        let b2 = dekkerStep(a1: a1, b0: b0, b1: b1, fa1: fa1, fb0: fb0, fb1: fb1)
        let fb2 = f(b2)
        let (a2,fa2) = fb2.sign != fa1.sign ? (a1,fa1) : (b1,fb1)
        let state1 = abs(fa2) < abs(fb2) ? (a1: b2, b0: b1, b1: a2, fa1: fb2, fb0: fb1, fb1: fa2) : (a1: a2, b0: b1, b1: b2, fa1: fa2, fb0: fb1, fb1: fb2)
        return (state: state1, guess:  state1.b1)
        }.converge { s1, s2 in abs(s2.state.fb1) < epsilon }
    guard let res = r else { return .nan }
    return res.guess
}

