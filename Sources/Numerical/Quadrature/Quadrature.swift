//
//  Quadrature.swift
//  Numerical
//
//  Created by Adam Roberts on 4/28/20.
//

import Foundation

public typealias UnivariateQuadrature = (CountedFunction<Double,Double>,ClosedRange<Double>,Int) -> QuadratureResult

public enum QuadratureResult {
    case error
    case noConverge(evals: Int, estimate: Double)
    case success(evals: Int, estimate: Double)
    
    var value: Double {
        switch self {
        case .error: return .nan
        case .noConverge(_, let est): return est
        case .success(_, let est): return est
        }
    }
}

public func integrate(range: ClosedRange<Double>, maxIter: Int = 10, method: UnivariateQuadrature = romberg, f: @escaping (Double) -> Double) -> Double {
    let countedF = CountedFunction(f: f)
    let b = method(countedF,range,maxIter)
    return b.value
}
