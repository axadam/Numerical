//
//  Quadrature.swift
//  Numerical
//
//  Created by Adam Roberts on 4/28/20.
//

import Foundation

public typealias UnivariateQuadrature = (ClosedRange<Double>,Int, @escaping (Double) -> Double) -> ConvergenceValue<Double>?

/// Numerical integration of a function on a closed range
///
/// Using the specified numerical integration method it will attempt to perform numerical
/// integration of your function no the specified interval. The default method is Romberg's
/// Method which is a pretty good option in most cases. An implementation of the
/// Trapezoidal Rule is also available and is well suited to integrating periodic functions
/// over their periods or peak functions with fast decaying tails.
///
/// Returns a `ConvergenceValue` which contains the integration estimate as well as
/// information about whether the method converged and how many function evaluations
/// were required.
public func integrate(range: ClosedRange<Double>, maxIter: Int = 10, method: UnivariateQuadrature = romberg, f: @escaping (Double) -> Double) -> ConvergenceValue<Double>? {
    let b = method(range,maxIter,f)
    return b
}
