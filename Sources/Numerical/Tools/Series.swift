//
//  Series.swift
//  Numerical
//
//  Created by Adam Roberts on 12/9/18.
//

import Foundation
import Scan

public func recursiveSequence<Indices: Sequence, State>(indices: Indices, initialState: State, maxIter: Int, update: (Int, State) -> State, until: (State, State) -> Bool) -> State? where Indices.Element == Int {
    let result = withoutActuallyEscaping(update) { u in
        indices.lazy.scan( initialState ) { state0, index in
            return u(index, state0)
            }.converge(max_iter: maxIter, until: until)
    }
    return result?.result
}

public func recursiveSeries<IntSequence: Sequence, State>(indices: IntSequence, accum0: Double, state0: State, accumulate: (Double,Double) -> Double, update: (Int,State) -> (Double,State), until: ((Double, Double), (Double, Double)) -> Bool, max_iter: Int = 100 ) -> Double where IntSequence.Element == Int {
    let result = withoutActuallyEscaping(accumulate) { a in
        withoutActuallyEscaping(update) { u in
            indices.lazy.scan( (accum0,0,state0)) { arg0, i -> (Double,Double,State) in
                let (accumPrev,_,statePrev) = arg0
                let (term, state) = u(i, statePrev)
                let accum = a(accumPrev,term)
                return (accum, term, state)
                }.converge(max_iter: max_iter) { a, b in until((a.0,a.1),(b.0,b.1)) }
        }
    }
    guard let res = result?.result.0 else {
        return .nan
    }
    return res
}

public func recursiveSum<IntSequence: Sequence, State>(indices: IntSequence, sum0: Double, state0: State, update: (Int,State) -> (Double,State), until: ((Double, Double), (Double, Double)) -> Bool, max_iter: Int = 100 ) -> Double where IntSequence.Element == Int {
    return recursiveSeries(indices: indices,
                           accum0: sum0,
                           state0: state0,
                           accumulate: { sumPrev, term in sumPrev + term },
                           update: update,
                           until: until,
                           max_iter: max_iter)
}

public func recursiveProduct<IntSequence: Sequence, State>(indices: IntSequence, product0: Double, state0: State, update: (Int,State) -> (Double,State), until: ((Double, Double), (Double, Double)) -> Bool, max_iter: Int = 100 ) -> Double where IntSequence.Element == Int {
    return recursiveSeries(indices: indices,
                           accum0: product0,
                           state0: state0,
                           accumulate: { productPrev, term in productPrev * term },
                           update: update,
                           until: until,
                           max_iter: max_iter)
}
