//
//  Roots.swift
//  Numerical
//
//  Created by Adam Roberts on 11/23/18.
//

import Foundation

/// Type signature for functions providing derivative free root finding within an interval
///
/// - Parameters:
///     - bracket: Interval that brackets at least one root of the function
///     - epsilon: Convergence criteria in terms of how close to zero we need to get.
///     - f: function whose root we seek.
public typealias bracketedRoot = (BracketedRootEstimate, Double, @escaping(Double) -> Double) -> BracketedRootResult

public enum BracketResult {
    case noBracket(evals: Int)
    case bracket(evals: Int, estimate: BracketedRootEstimate)
}

extension BracketResult {
    var evals: Int {
        switch self {
        case .noBracket(let n): return n
        case .bracket(let n, _): return n
        }
    }
}

public enum BracketedRootResult {
    case error
    case noConverge(evals: Int, estimate: BracketedRootEstimate)
    case success(evals: Int, estimate: BracketedRootEstimate)
}

public enum BracketAndRootResult {
    case error
    case noBracket(bracketEvals: Int)
    case noConverge(bracketEvals: Int, rootEvals: Int, estimate: BracketedRootEstimate)
    case success(bracketEvals: Int, rootEvals: Int, estimate: BracketedRootEstimate)
}

public extension BracketAndRootResult {
    var value: Double {
        switch self {
        case .error: return .nan
        case .noBracket: return .nan
        case .noConverge(_,_,let e): return e.value
        case .success(_,_,let e): return e.value
        }
    }
    
    var bracketEvals: Int {
        switch self {
        case .error: return 0
        case .noBracket(let n): return n
        case .noConverge(let n,_,_): return n
        case .success(let n,_,_): return n
        }
    }
    
    var rootEvals: Int {
        switch self {
        case .error: return 0
        case .noBracket: return 0
        case .noConverge(_, let n,_): return n
        case .success(_, let n,_): return n
        }
    }
    
    var evals: Int {
        return bracketEvals + rootEvals
    }
    
    var converged: Bool {
        switch self {
        case .success: return true
        default: return false
        }
    }
}

/// Estimate of the root of a function
///
/// Represents an interval, [a,b], and the function evaluted at the ends. a ≤ b and f(a) * f(b) ≤ 0
public struct BracketedRootEstimate {
    public let a: Double
    public let b: Double
    public let fa: Double
    public let fb: Double
}

public extension BracketedRootEstimate {
    init(x: Double, fx: Double) {
        a = x
        b = x
        fa = fx
        fb = fx
    }
    
    static var nan: BracketedRootEstimate {
        return self.init(x: .nan, fx: .nan)
    }
    
    var value: Double {
        switch (fa,fb) {
        case (0,_): return a
        case (_,0): return b
        case (_,_): return secantStep(x0: a, x1: b, y0: fa, y1: fb)
        }
    }
}

/// Root finder for univariate functions
///
/// Starting from your guess it will attempt to bracket a root, and then
/// use the specified root finding method to find the root. Success can
/// depend on how good your guess is, especially if the function is not
/// monotonic.
///
/// If you can cheaply get one or two derivatives of your function then
/// use the appropriate variant of `root()`.
///
/// Default root finding method is Brent's.
///
/// - Parameters:
///     - guess: Initial guess for root.
///     - xmin: Minimum allowed value for root. Optional.
///     - xmax: Maximum allowed value for root. Optional.
///     - epsilon: Convergence criteria in terms of how close to zero we need to get. Default 1e-10.
///     - method: Root finding method. Defaults to Brent's method.
public func root(guess: Double, xmin: Double? = nil, xmax: Double? = nil, tolerance: Double = 10e-15, bracketFactor: Double = 1.6, method: bracketedRoot = riddersRoot, f: @escaping (Double) -> Double) -> BracketAndRootResult {
    // attempt to bracket the root starting from the given guess
    let brac = bracket(f: f, guess: guess, xmin: xmin, xmax: xmax, factor: bracketFactor)
    
    // we couldn't bracket starting from this guess. give up
    guard case let .bracket(_, b) = brac else {
        return .noBracket(bracketEvals: brac.evals)
    }
    
    // if we get a successful bracket then find root with specified method
    let r = root(f: f, bracket: b, tolerance: tolerance, method: method)
    switch r {
    case .error: return .error
    case .noConverge(let n, let e): return .noConverge(bracketEvals: brac.evals, rootEvals: n, estimate: e)
    case .success(let n, let e): return .success(bracketEvals: brac.evals, rootEvals: n, estimate: e)
    }
}

func root(f: @escaping(Double) -> Double, bracket: BracketedRootEstimate, tolerance: Double, method: bracketedRoot) -> BracketedRootResult {
    switch (bracket.fa,bracket.fb) {
    case (-tolerance...tolerance,_): return .success(evals: 0, estimate: BracketedRootEstimate(x: bracket.a, fx: bracket.fa))
    case (_,-tolerance...tolerance): return .success(evals: 0, estimate: BracketedRootEstimate(x: bracket.b, fx: bracket.fb))
    case (..<0,..<0): fallthrough
    case (0...,0...): return .error
    default: return method(bracket,tolerance,f)
    }
}

public enum RootResult {
    case error
    case noConverge(evals: Int, estimate: Double)
    case success(evals: Int, estimate: Double)
}

public extension RootResult {
    var value : Double {
        switch self {
        case .error: return .nan
        case .noConverge(_,let e): return e
        case .success(_, let e): return e
        }
    }
    
    var evals: Int {
        switch self {
        case .error: return 0
        case .noConverge(let n,_): return n
        case .success(let n,_): return n
        }
    }
    
    var converged: Bool {
        switch self {
        case .success: return true
        default: return false
        }
    }
}

/// Root finder for univariate function with one derivative
///
/// Starting from your guess uses Newton's method to find a root of your
/// function. Optionally accepts bounds. Success can depend on having a good
/// enough guess. Can be vulnerable to poorly behaved functions, especially
/// those with flat regions.
///
/// - Parameters:
///     - guess: Initial guess for root.
///     - xmin: Minimum allowed value for root. Optional.
///     - xmax: Maximum allowed value for root. Optional.
///     - maxIter: Maximum iterations allowed. Default 100.
///     - xtol: Convergence criteria in terms how small a step we've made. Default 1e-10.
///     - f: Function in which to find root
///     - f1: First derivative of f, f'
public func root(guess: Double, xmin: Double? = nil, xmax: Double? = nil, maxIter: Int = 10, xtol: Double = 1e-10, f: @escaping(Double) -> Double, f1: @escaping(Double) -> Double) -> RootResult {
    return newtonRoot(f: f, f1: f1, guess: guess, xmin: xmin, xmax: xmax, max_iter: maxIter, xtol: xtol)
}

/// Root finder for univariate function with two derivatives
///
/// Starting from your guess, uses the specified second order root finding
/// method to find the root of your function. Optionally accepts bounds. Success
/// can depend on having a good enough guess.
///
/// This version takes the ratio of the second derivative to the first as
/// an argument.
///
/// - Parameters:
///     - guess: Initial guess for root.
///     - xmin: Minimum allowed value for root. Optional.
///     - xmax: Maximum allowed value for root. Optional.
///     - maxIter: Maximum iterations allowed. Default 100.
///     - xtol: Convergence criteria in terms how small a step we've made. Default 1e-10.
///     - f: Function in which to find root
///     - f1: First derivative of f, f'
///     - f2: Second derivative of f, f"
public func root(guess: Double, xmin: Double? = nil, xmax: Double? = nil, maxIter: Int = 100, xtol: Double = 1e-10, f: @escaping(Double) -> Double, f1: @escaping(Double) -> Double, f2: @escaping(Double) -> Double) -> RootResult {
    return halleyRoot(guess: guess, xmin: xmin, xmax: xmax, maxIter: maxIter, xtol: xtol, f: f, f1: f1, f2: f2)
}

/// Root finder for univariate function with two derivatives
///
/// Starting from your guess, uses the specified second order root finding
/// method to find the root of your function. Optionally accepts bounds. Success
/// can depend on having a good enough guess.
///
/// - Parameters:
///     - guess: Initial guess for root.
///     - xmin: Minimum allowed value for root. Optional.
///     - xmax: Maximum allowed value for root. Optional.
///     - maxIter: Maximum iterations allowed. Default 100.
///     - xtol: Convergence criteria in terms how small a step we've made. Default 1e-10.
///     - f: Function in which to find root
///     - f1: First derivative of f, f'
///     - f2f1: Ratio of second derivative of f to the first, f" / f'
public func root(guess: Double, xmin: Double? = nil, xmax: Double? = nil, maxIter: Int = 100, xtol: Double = 1e-10, f: @escaping(Double) -> Double, f1: @escaping(Double) -> Double, f2f1: @escaping(Double) -> Double) -> RootResult {
    return halleyRoot(guess: guess, xmin: xmin, xmax: xmax, maxIter: maxIter, xtol: xtol, f: f, f1: f1, f2f1: f2f1)
}

/// Helper function for root finding
///
/// Provides bounds with bisection if your step goes out of bounds
func rootHelper(guess: Double, xmin: Double? = nil, xmax: Double? = nil, maxIter: Int = 100, xtol: Double = 1e-10, step: @escaping(Double) -> Double) -> RootResult {
    let r = sequence(first: guess) { x0 in
        let x1 = step(x0)
        if let xmin = xmin, x1 < xmin { return bisect(a: x0, b: xmin) }
        if let xmax = xmax, x1 > xmax { return bisect(a: x0, b: xmax) }
        return x1
    }.until(maxIter: maxIter) { s1, s2 in abs(s2 - s1) < xtol }
    
    guard let res = r else { return .error }
    
    switch res.exitState {
    case .exhaustedInput: return .error
    case .exceededMax: return .noConverge(evals: res.iterations, estimate: res.result)
    case .converged: return .success(evals: res.iterations, estimate: res.result)
    }
}
