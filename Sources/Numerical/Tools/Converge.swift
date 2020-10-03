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
    func until(minIter: Int = 0, maxIter: Int = 100, _ predicate: (Element) -> Bool) -> UntilValue<Element>? {
        var g = makeIterator()
        var count: UInt = 0
        var last: Element? = nil
        while let e = g.next() {
            if predicate(e) && count >= minIter {
                return .success(iterations: count, value: e)
            }
            last = e
            if count >= maxIter {
                return .exceededMax(iterations: count, value: e)
            }
            count += 1
        }
        guard let e = last else { return nil }
        return .exhaustedInput(iterations: count, value: e)
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
    func until(minIter: Int = 0, maxIter: Int = 100, _ predicate: (Element,Element) -> Bool) -> UntilValue<Element>? {
        var g = makeIterator()
        var count: UInt = 0
        var last: Element? = nil
        while let e = g.next() {
            if let lasty = last, predicate(lasty,e) && count >= minIter {
                return .success(iterations: count + 1, value: e)
            }
            last = e
            if count >= maxIter {
                return .exceededMax(iterations: count + 1, value: e)
            }
            count += 1
        }
        guard let e = last else { return nil }
        return .exhaustedInput(iterations: count, value: e)
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

/// The value of an `until()` operation on a `Sequence`
///
/// Stores the exit state of the operation, the last value at time of exit, and how many iterations
/// were completed.
public enum UntilValue<Element>: IterativeValue {
    /// The sequence was used up before satisfying the convergence condition
    case exceededMax(iterations: UInt, value: Element)
    /// Reached the maximum number of iterations before satisfying the convergence condition
    case exhaustedInput(iterations: UInt, value: Element)
    /// Satisfied convergence condition
    case success(iterations: UInt, value: Element)
}

public extension UntilValue {
    typealias Value = Element
    
    /// The last value at time of exit regardless of exit state
    var value: Element {
        switch self {
        case .exceededMax(_, let v): return v
        case .exhaustedInput(_, let v): return v
        case .success(_, let v): return v
        }
    }

    /// How many iterations were competed before exiting
    var work: UInt {
        switch self {
        case .exceededMax(let w, _): return w
        case .exhaustedInput(let w, _): return w
        case .success(let w, _): return w
        }
    }
}

