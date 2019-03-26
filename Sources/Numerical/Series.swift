//
//  Series.swift
//  Numerical
//
//  Created by Adam Roberts on 12/9/18.
//

import Foundation
import Scan

public extension Sequence {
    func converge(max_iter: Int = 100, until f: (Element,Element) -> Bool) -> Element? {
        var g = makeIterator()
        var count = 0
        var last: Element? = nil
        while let e = g.next(), count < max_iter {
            if let lasty = last, f(lasty,e) {
                return e
            }
            last = e
            count += 1
        }
        if count >= max_iter {
            print("Did not converge in \(count) iterations")
        } else {
            print("Did not converge before end of sequence")
        }
        return last
    }
}

public func recursiveSeries<IntSequence: Sequence, State>(indices: IntSequence, accum0: Double, state0: State, accumulate: @escaping (Double,Double) -> Double, update: @escaping (Int,State) -> (Double,State), until: ((Double, Double), (Double, Double)) -> Bool, max_iter: Int = 100 ) -> Double where IntSequence.Element == Int {
    let result = indices.lazy.scan( (accum0,0,state0)) { arg0, i -> (Double,Double,State) in
        let (accumPrev,_,statePrev) = arg0
        let (term, state) = update(i, statePrev)
        let accum = accumulate(accumPrev,term)
        return (accum, term, state)
        }.converge(max_iter: max_iter) { a, b in until((a.0,a.1),(b.0,b.1)) }
    guard let res = result?.0 else {
        return 0
    }
    return res
}

public func recursiveSum<IntSequence: Sequence, State>(indices: IntSequence, sum0: Double, state0: State, update: @escaping (Int,State) -> (Double,State), until: ((Double, Double), (Double, Double)) -> Bool, max_iter: Int = 100 ) -> Double where IntSequence.Element == Int {
    return recursiveSeries(indices: indices,
                           accum0: sum0,
                           state0: state0,
                           accumulate: { sumPrev, term in sumPrev + term },
                           update: update,
                           until: until,
                           max_iter: max_iter)
}

public func recursiveProduct<IntSequence: Sequence, State>(indices: IntSequence, product0: Double, state0: State, update: @escaping (Int,State) -> (Double,State), until: ((Double, Double), (Double, Double)) -> Bool, max_iter: Int = 100 ) -> Double where IntSequence.Element == Int {
    return recursiveSeries(indices: indices,
                           accum0: product0,
                           state0: state0,
                           accumulate: { productPrev, term in productPrev * term },
                           update: update,
                           until: until,
                           max_iter: max_iter)
}
