//
//  Series.swift
//  Numerical
//
//  Created by Adam Roberts on 12/9/18.
//

import Foundation
import Scan

public enum SeriesResult {
    case error
    case noConverge(terms: Int, estimate: Double)
    case success(terms: Int, estimate: Double)
}

public extension SeriesResult {
    var value: Double {
        switch self {
        case .error: return .nan
        case .noConverge(_, let e): return e
        case .success(_, let e): return e
        }
    }
    
    var iterations: Int {
        switch self {
        case .error: return 0
        case .noConverge(let n, _): return n
        case .success(let n, _): return n
        }
    }
    
    var converged: Bool {
        switch self {
        case .error, .noConverge: return false
        case .success: return true
        }
    }
}

public func indexedAccumulatingRecursiveSequence<IntSequence: Sequence, State>(indices: IntSequence, accum0: Double, state0: State, accumulate: @escaping (Double,Double) -> Double, update: @escaping (Int,State) -> (Double,State) ) -> LazyScanSequence<LazySequence<IntSequence>, (Double, Double, State)> where IntSequence.Element == Int {
    return indices.lazy.scan( (accum0,.nan,state0)) { arg0, i -> (Double,Double,State) in
        let (accumPrev,_,statePrev) = arg0
        let (term, state) = update(i, statePrev)
        let accum = accumulate(accumPrev,term)
        return (accum, term, state)
    }
}

public func series<IntSequence: Sequence, State>(indices: IntSequence, initialSum: Double = 0, initialState: State, maxIter: Int = 100, tolerance: EqualityTolerance<Double> = .strict, update: @escaping (Int,State) -> (Double,State)) -> SeriesResult where IntSequence.Element == Int {
    let r = indexedAccumulatingRecursiveSequence(
        indices: indices,
        accum0: initialSum,
        state0: initialState,
        accumulate: +,
        update: update
    ).until(maxIter: maxIter) { b in b.1.isApprox(.zero(scaleRelativeTo: b.0), tolerance: tolerance) }
    guard let res = r else { return .error }
    switch res.exitState {
    case .exhaustedInput, .exceededMax: return .noConverge(terms: res.iterations, estimate: res.result.0)
    case .converged: return .success(terms: res.iterations, estimate: res.result.0)
    }
}

public func product<IntSequence: Sequence, State>(indices: IntSequence, initialProduct: Double = 1, initialState: State, maxIter: Int = 100, tolerance: EqualityTolerance<Double> = .strict, update: @escaping (Int,State) -> (Double,State)) -> SeriesResult where IntSequence.Element == Int {
    let r = indexedAccumulatingRecursiveSequence(
        indices: indices,
        accum0: initialProduct,
        state0: initialState,
        accumulate: *,
        update: update
    ).until(maxIter: maxIter) { b in b.1.isApprox(.maybeZero(1, trusted: true), tolerance: tolerance) }
    guard let res = r else { return .error }
    switch res.exitState {
    case .exhaustedInput, .exceededMax: return .noConverge(terms: res.iterations, estimate: res.result.0)
    case .converged: return .success(terms: res.iterations, estimate: res.result.0)
    }
}
