//
//  Toms748.swift
//  Numerical
//
//  Created by Adam Roberts on 4/11/19.
//

import Foundation

/// TOMS Algo 748 Root finding method
///
/// An advanced root finding algorithm making use of cubic and quadratic
/// interpolation and bracketing to keep the process safe.
///
/// TOMS Algorithm 748 (1995), Section 4, Algorithm 4.2
public func toms748Root(bracket: BracketedRootEstimate, tolerance: EqualityTolerance<Double> = .strict, intercept: Double = 0, f rawF: @escaping(Double) -> Double) -> BracketedRootResult {
    let f = intercept == 0 ? CountedFunction(f: rawF) : CountedFunction { rawF($0) - intercept }
    let (a,b,fa,fb) = (bracket.a,bracket.b,bracket.fa,bracket.fb)

    // factor by which we must shrink interval each iteration or we choose bisection
    // usually chosen as 0.5
    let µ = 0.5
    
    // 4.1.1 Initial guess using the secant method
    let c₁ = a - fa * (b - a) / (fb - fa)
    
    // 4.1.2 Bracket
    let (a₂, b₂, d₂, fa₂, fb₂, fd₂) = toms748Bracket(f: f, a: a, b: b, c: c₁, fa: fa, fb: fb, tolerance: tolerance.absolute)
    
    typealias Toms748State = (a: Double, b: Double, d: Double, e: Double, fa: Double, fb: Double, fd: Double, fe: Double)
    let r: UntilValue<Toms748State>? = (2...).lazy.scan((a: a₂, b: b₂, d: d₂, e: 1.0, fa: fa₂, fb: fb₂, fd: fd₂, fe: 1e5)) { (state: Toms748State, i: Int) -> Toms748State in
        let (aᵢ, bᵢ, dᵢ, eᵢ, faᵢ, fbᵢ, fdᵢ, feᵢ) = state
        
        // 4.2.3 First guess for the iteration is Inverse Cubic Interpolation if
        //       conditions are met. Otherwise Newton Quadratic. Also fall back
        //       if cubic takes us out of bounds
        let cᵢ: Double = {
            if i == 2 || (faᵢ - fbᵢ) * (faᵢ - fdᵢ) * (faᵢ - feᵢ) * (fbᵢ - fdᵢ) * (fbᵢ - feᵢ) * (fdᵢ - feᵢ) == 0 {
                return newtonQuadtraticStep(a: aᵢ, b: bᵢ, d: dᵢ, fa: faᵢ, fb: fbᵢ, fd: fdᵢ, k: 2)
            }
            let cCubic = inverseCubicInterpolation(a: aᵢ, b: bᵢ, c: dᵢ, d: eᵢ, fa: faᵢ, fb: fbᵢ, fc: fdᵢ, fd: feᵢ)
            if (cCubic - aᵢ) * (cCubic - bᵢ) >= 0 {
                return newtonQuadtraticStep(a: aᵢ, b: bᵢ, d: dᵢ, fa: faᵢ, fb: fbᵢ, fd: fdᵢ, k: 2)
            }
            return cCubic
        }()
        
        // 4.2.4 Keep track of previous d, then bracket. Check if we've found root
        let (ẽᵢ,fẽᵢ) = (dᵢ,fdᵢ)
        let (ãᵢ, b̃ᵢ, d̃ᵢ, fãᵢ, fb̃ᵢ, fd̃ᵢ) = toms748Bracket(f: f, a: aᵢ, b: bᵢ, c: cᵢ, fa: faᵢ, fb: fbᵢ, tolerance: tolerance.absolute)
        if fãᵢ == 0 || b̃ᵢ - ãᵢ < 2 * tole(a: ãᵢ, b: b̃ᵢ, fa: fãᵢ, fb: fb̃ᵢ, tolerance: tolerance.absolute) {
            return (ãᵢ, b̃ᵢ, d̃ᵢ, ẽᵢ, fãᵢ, fb̃ᵢ, fd̃ᵢ, fẽᵢ)
        }
        
        // 4.2.5 Next guess uses the same approach. Newton Quadratic is allowed
        //       three iterations this time
        let c̃ᵢ: Double = {
            if (fãᵢ - fb̃ᵢ) * (fãᵢ - fd̃ᵢ) * (fãᵢ - fẽᵢ) * (fb̃ᵢ - fd̃ᵢ) * (fb̃ᵢ - fẽᵢ) * (fd̃ᵢ - fẽᵢ) == 0 {
                return newtonQuadtraticStep(a: ãᵢ, b: b̃ᵢ, d: d̃ᵢ, fa: fãᵢ, fb: fb̃ᵢ, fd: fd̃ᵢ, k: 3)
            }
            let cCubic = inverseCubicInterpolation(a: ãᵢ, b: b̃ᵢ, c: d̃ᵢ, d: ẽᵢ, fa: fãᵢ, fb: fb̃ᵢ, fc: fd̃ᵢ, fd: fẽᵢ)
            if (cCubic - ãᵢ) * (cCubic - b̃ᵢ) >= 0 {
                return newtonQuadtraticStep(a: ãᵢ, b: b̃ᵢ, d: d̃ᵢ, fa: fãᵢ, fb: fb̃ᵢ, fd: fd̃ᵢ, k: 3)
            }
            return cCubic
        }()
        
        // 4.2.6 Bracket our latest guess and check if we've found root
        let (āᵢ, b̄ᵢ, d̄ᵢ, fāᵢ, fb̄ᵢ, fd̄ᵢ) = toms748Bracket(f: f, a: ãᵢ, b: b̃ᵢ, c: c̃ᵢ, fa: fãᵢ, fb: fb̃ᵢ, tolerance: tolerance.absolute)
        if fāᵢ == 0 || b̄ᵢ - āᵢ < tole(a: āᵢ, b: b̄ᵢ, fa: fāᵢ, fb: fb̄ᵢ, tolerance: tolerance.absolute) {
            return (āᵢ, b̄ᵢ, d̄ᵢ, ẽᵢ, fāᵢ, fb̄ᵢ, fd̄ᵢ, fẽᵢ)
        }
        
        // 4.1.5 Take the better estimate between a and b
        let (uᵢ,fuᵢ) = abs(fāᵢ) < abs(fb̄ᵢ) ? (āᵢ,fāᵢ) : (b̄ᵢ,fb̄ᵢ)
        
        // 4.1.6 Modified secant method for our next guess
        let c̄ᵢ = uᵢ - 2 * (b̄ᵢ - āᵢ) / (fb̄ᵢ - fāᵢ) * fuᵢ
        
        // 4.1.7 But don't trust secant further than half the width of the interval
        let ĉᵢ = abs(c̄ᵢ - uᵢ) > 0.5 * (b̄ᵢ - āᵢ) ? 0.5 * (b̄ᵢ + āᵢ) : c̄ᵢ
        
        // 4.1.8 Bracket our latest guess
        let (âᵢ,b̂ᵢ,d̂ᵢ,fâᵢ,fb̂ᵢ,fd̂ᵢ) = toms748Bracket(f: f, a: āᵢ, b: b̄ᵢ, c: ĉᵢ, fa: fāᵢ, fb: fb̄ᵢ, tolerance: tolerance.absolute)
        
        // 4.1.9 Accept our latest guess if we shrank the interval by µ, otherwise
        //       choose the midpoint and re-bracket using that
        if b̂ᵢ - âᵢ < µ * (bᵢ - aᵢ) {
            return (âᵢ,b̂ᵢ,d̂ᵢ,d̄ᵢ,fâᵢ,fb̂ᵢ,fd̂ᵢ,fd̄ᵢ)
        }
        let (eᵢ₊₁,feᵢ₊₁) = (d̂ᵢ,fd̂ᵢ)
        let (aᵢ₊₁,bᵢ₊₁,dᵢ₊₁,faᵢ₊₁,fbᵢ₊₁,fdᵢ₊₁) = toms748Bracket(f: f, a: âᵢ, b: b̂ᵢ, c: 0.5 * (âᵢ + b̂ᵢ), fa: fâᵢ, fb: fb̂ᵢ, tolerance: tolerance.absolute)
        return (aᵢ₊₁,bᵢ₊₁,dᵢ₊₁,eᵢ₊₁,faᵢ₊₁,fbᵢ₊₁,fdᵢ₊₁,feᵢ₊₁)
    }.until(maxIter: 30) { s2 in s2.fa == 0 || s2.b - s2.a <= 2 * tole(a: s2.a, b: s2.b, fa: s2.fa, fb: s2.fb, tolerance: tolerance.absolute) }

    guard let res = r else { return .error } // shouldn't happen
    
    let e = BracketedRootEstimate(a: res.value.a, b: res.value.b, fa: res.value.fa, fb: res.value.fb)
    
    switch res {
    case .exhaustedInput: return .error // shouldn't happen
    case .exceededMax: return .noConverge(evals: f.count, estimate: e)
    case .success: return .success(evals: f.count, estimate: e)
    }
}


/// Algo 748 Re-bracketing method
///
/// TOMS Algorithm 748 (1995), Section 6, Subroutine bracket(a, b, c, ā, b̄, d)
func toms748Bracket(f: CountedFunction<Double,Double>, a: Double, b: Double, c: Double, fa: Double, fb: Double, tolerance: Double) -> (a: Double, b: Double, d: Double, fa: Double, fb: Double, fd: Double) {
    let λ = 0.7
    let δ = λ * tole(a: a, b: b, fa: fa, fb: fb, tolerance: tolerance)
    
    // if c is outside interval or too close to edges move it towards the middle
    let c1: Double = {
        switch (b - a,c) {
        // a and b are too close together, bisect
        case (...(4 * δ),_): return (a + b) / 2
        // c is too close to a or outside a, place inside near a
        case (_,...(a + 2 * δ)): return a + 2 * δ
        // c is too close to b or outside b, place inside near b
        case (_,(b - 2 * δ)...): return b - 2 * δ
        // c is inside and not too close to edges
        default: return c
        }
    }()
    let fc1 = f(c1)

    if fc1 == 0 {
        return (a: c1, b: b, d: c1, fa: fc1, fb: fb, fd: fc1)
    }
    
    if fc1 * fa < 0 {
        return (a: a, b: c1, d: b, fa: fa, fb: fc1, fd: fb)
    }
    return (a: c1, b: b, d: a, fa: fc1, fb: fb, fd: fa)
}

/// Algo 748 Relative Tolerance method
///
/// TOMS Algorithm 748 (1995), Section 6, Equation 25
func tole(a: Double, b: Double, fa: Double, fb: Double, tolerance: Double) -> Double {
    let u = abs(fa) < abs(fb) ? a : b
    let tole = 2 * abs(u) * Double.ulpOfOne + tolerance
    return tole
}
