//
//  Converge.swift
//  Numerical
//
//  Created by Adam Roberts on 11/28/19.
//

public extension Sequence {
    /// Returns the first member of the sequence that along with its preceding member
    /// satisfies the condition, or the last member before stopping.
    ///
    /// It will stop before satisfying the condition if the maximum number of iterations is
    /// reached or the sequence is exhuasted. In either of those cases it will record its
    /// exit state in the result.
    ///
    /// Returns nil for the empty sequence.
    @inlinable
    func converge(max_iter: Int = 100, until f: (Element,Element) -> Bool) -> IterativeResult<Element,ConvergenceState>? {
        var g = makeIterator()
        var count = 0
        var last: Element? = nil
        while let e = g.next() {
            if let lasty = last, f(lasty,e) {
                return IterativeResult(iterations: count + 1, exitState: .converged, result: e)
            }
            last = e
            if count >= max_iter {
                return IterativeResult(iterations: count + 1, exitState: .exceededMax, result: e)
            }
            count += 1
        }
        guard let e = last else { return nil }
        return IterativeResult(iterations: count, exitState: .exhaustedInput, result: e)
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
