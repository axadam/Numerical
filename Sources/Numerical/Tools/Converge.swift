//
//  Converge.swift
//  Numerical
//
//  Created by Adam Roberts on 11/28/19.
//

public extension Sequence {
    
    /// Returns the first member of the sequence that satisfies the condition,
    /// or the last member before stopping.
    ///
    /// It will stop before satisfying the condition if the maximum number of iterations is
    /// reached or the sequence is exhuasted. In either of those cases it will record its
    /// exit state in the result.
    ///
    /// Returns nil for the empty sequence.
    @inlinable
    func until(minIter: Int = 0, maxIter: Int = 100, _ predicate: (Element) -> Bool) -> IterativeResult<Element,ConvergenceState>? {
        var g = makeIterator()
        var count = 0
        var last: Element? = nil
        while let e = g.next() {
            if predicate(e) && count >= minIter {
                return IterativeResult(iterations: count, exitState: .converged, result: e)
            }
            last = e
            if count >= maxIter {
                return IterativeResult(iterations: count, exitState: .exceededMax, result: e)
            }
            count += 1
        }
        guard let e = last else { return nil }
        return IterativeResult(iterations: count, exitState: .exhaustedInput, result: e)
    }
    
    /// Returns the first member of the sequence that along with its preceding member
    /// satisfies the condition, or the last member before stopping.
    ///
    /// It will stop before satisfying the condition if the maximum number of iterations is
    /// reached or the sequence is exhuasted. In either of those cases it will record its
    /// exit state in the result.
    ///
    /// Returns nil for the empty sequence.
    @inlinable
    func until(minIter: Int = 0, maxIter: Int = 100, _ predicate: (Element,Element) -> Bool) -> IterativeResult<Element,ConvergenceState>? {
        var g = makeIterator()
        var count = 0
        var last: Element? = nil
        while let e = g.next() {
            if let lasty = last, predicate(lasty,e) && count >= minIter {
                return IterativeResult(iterations: count + 1, exitState: .converged, result: e)
            }
            last = e
            if count >= maxIter {
                return IterativeResult(iterations: count + 1, exitState: .exceededMax, result: e)
            }
            count += 1
        }
        guard let e = last else { return nil }
        return IterativeResult(iterations: count, exitState: .exhaustedInput, result: e)
    }
    
    /// Returns the first member, if any, of the sequence that along with its preceding member
    /// satisfies the predicate
    ///
    /// Returns nil if the predicate is never met
    @inlinable
    func first(where predicate: (Element,Element) throws -> Bool) rethrows -> Element? {
        var g = makeIterator()
        var last: Element? = nil
        while let e = g.next() {
            if let _last = last, try predicate(_last,e) {
                return e
            }
            last = e
        }
        return nil
    }
}

/// Value representing the exit state of the converge function on a sequence
public enum ConvergenceState {
    
    /// Success. Satisfied convergence condition
    case converged
    
    /// The sequence was used up before satisfying the convergence condition
    case exhaustedInput
    
    /// Reached the maximum number of iterations before satisfying the convergence condition
    case exceededMax
}
