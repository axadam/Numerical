//
//  Quadrature.swift
//  Numerical
//
//  Created by Adam Roberts on 4/28/20.
//

import Foundation

public typealias UnivariateQuadrature = (ClosedRange<Double>,Int, @escaping (Double) -> Double) -> QuadratureResult

/// The result of a numerical integration
///
/// Contains the estimate from the numerical integration as well information about whether
/// the method converged and how many function evaluations were required.
public enum QuadratureResult {
    /// An unexpected error occurred
    case error
    
    /// Method did not converge. Includes numbers of evalutations and the last estimate
    case noConverge(evals: Int, estimate: Double)
    
    /// Method converged successfully. Includes number of evaluations and the estimate
    case success(evals: Int, estimate: Double)
}

public extension QuadratureResult {
    /// Numeric value
    ///
    /// Gives the best estimate available if we have one. Otherwise returns NaN.
    var value: Double {
        switch self {
        case .error: return .nan
        case .noConverge(_, let est): return est
        case .success(_, let est): return est
        }
    }
    
    /// Number of evaluations for this result
    ///
    /// Gives the number of evaluations or zero if there was an error
    var evals: Int {
        switch self {
        case .error: return 0
        case .noConverge(let n, _): return n
        case .success(let n, _): return n
        }
    }
    
    /// Whether the quadrature converged or not
    var converged: Bool {
        switch self {
        case .error, .noConverge: return false
        case .success: return true
        }
    }
}

/// Numerical integration of a function on a closed range
///
/// Using the specified numerical integration method it will attempt to perform numerical
/// integration of your function no the specified interval. The default method is Romberg's
/// Method which is a pretty good option in most cases. An implementation of the
/// Trapezoidal Rule is also available and is well suited to integrating periodic functions
/// over their periods or peak functions with fast decaying tails.
///
/// Returns a `QuadratureResult` which contains the integration estimate as well as
/// information about whether the method converged and how many function evaluations
/// were required.
public func integrate(range: ClosedRange<Double>, maxIter: Int = 10, method: UnivariateQuadrature = romberg, f: @escaping (Double) -> Double) -> QuadratureResult {
    let b = method(range,maxIter,f)
    return b
}
