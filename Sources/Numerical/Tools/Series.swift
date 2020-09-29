//
//  Series.swift
//  Numerical
//
//  Created by Adam Roberts on 12/9/18.
//

import Foundation
import Scan

public func indexedAccumulatingRecursiveSequence<IntSequence: Sequence, State>(indices: IntSequence, accum0: Double, state0: State, accumulate: @escaping (Double,Double) -> Double, update: @escaping (Int,State) -> (Double,State) ) -> LazyScanSequence<LazySequence<IntSequence>, (Double, Double, State)> where IntSequence.Element == Int {
    return indices.lazy.scan( (accum0,.nan,state0)) { arg0, i -> (Double,Double,State) in
        let (accumPrev,_,statePrev) = arg0
        let (term, state) = update(i, statePrev)
        let accum = accumulate(accumPrev,term)
        return (accum, term, state)
    }
}

/// Compute a truncated infinite series with terms defined in terms of an index and optionally some additional state
///
/// Computes S = Σ i ∈ `indices` { aᵢ } where `indices` is any sequence of Ints. Truncation happens at the earlier of
/// when the specified tolerance is reached for nearness between adjacent terms or when the number of terms reaches `maxiter`.
///
/// - Parameters:
///    - indices: Any sequence of integers. Can be an infinite sequence (e.g. 1...)
///    - initialSum: Initial value for the sum. Default is 0
///    - initialState: Initial values for any state that is used in calculating terms
///    - maxIter: Maximum number of terms to calculate before truncating. Default is 100
///    - tolerance: Allowable tolerance for deciding the series has converged and truncating. Default is fairly strict
///    - update: Closure defining the ith term of the sequence. Has access to the index, i, and any state as of the previous term
///
/// - Returns: A `SeriesResult` that encodes whether the series converged or not and how many terms were calculated
public func series<IntSequence: Sequence, State>(indices: IntSequence, initialSum: Double = 0, initialState: State, maxIter: Int = 100, tolerance: EqualityTolerance<Double> = .strict, update: @escaping (Int,State) -> (Double,State)) -> ConvergenceValue<Double>? where IntSequence.Element == Int {
    let r = indexedAccumulatingRecursiveSequence(
        indices: indices,
        accum0: initialSum,
        state0: initialState,
        accumulate: +,
        update: update
    ).until(maxIter: maxIter) { b in b.1.isApprox(.zero(scaleRelativeTo: b.0), tolerance: tolerance) }
    guard let res = r else { return nil }
    switch res {
    case .exhaustedInput, .exceededMax: return .didNotConverge(work: res.work, estimate: res.value.0)
    case .success: return .converged(work: res.work, estimate: res.value.0)
    }
}

/// Compute a truncated infinite product with terms defined in terms of an index and optionally some additional state
///
/// Computes S = ∏ i ∈ `indices` { aᵢ } where `indices` is any sequence of Ints. Truncation happens at the earlier of
/// when the specified tolerance is reached for nearness between adjacent terms or when the number of terms reaches `maxiter`.
///
/// - Parameters:
///    - indices: Any sequence of integers. Can be an infinite sequence (e.g. 1...)
///    - initialProduct: Initial value for the product. Default is 1
///    - initialState: Initial values for any state that is used in calculating terms
///    - maxIter: Maximum number of terms to calculate before truncating. Default is 100
///    - tolerance: Allowable tolerance for deciding the series has converged and truncating. Default is fairly strict
///    - update: Closure defining the ith term of the sequence. Has access to the index, i, and any state as of the previous term
///
/// - Returns: A `SeriesResult` that encodes whether the product converged or not and how many terms were calculated
public func product<IntSequence: Sequence, State>(indices: IntSequence, initialProduct: Double = 1, initialState: State, maxIter: Int = 100, tolerance: EqualityTolerance<Double> = .strict, update: @escaping (Int,State) -> (Double,State)) -> ConvergenceValue<Double>? where IntSequence.Element == Int {
    let r = indexedAccumulatingRecursiveSequence(
        indices: indices,
        accum0: initialProduct,
        state0: initialState,
        accumulate: *,
        update: update
    ).until(maxIter: maxIter) { b in b.1.isApprox(.maybeZero(1, trusted: true), tolerance: tolerance) }
    guard let res = r else { return nil }
    switch res {
    case .exhaustedInput, .exceededMax: return .didNotConverge(work: res.work, estimate: res.value.0)
    case .success: return .converged(work: res.work, estimate: res.value.0)
    }
}
