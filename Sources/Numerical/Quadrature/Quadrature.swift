//
//  Quadrature.swift
//  Numerical
//
//  Created by Adam Roberts on 4/28/20.
//

import Foundation

public typealias UnivariateQuadrature = (CountedFunction<Double,Double>,ClosedRange<Double>,Int) -> QuadratureResult

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
    let countedF = CountedFunction(f: f)
    let b = method(countedF,range,maxIter)
    return b
}
