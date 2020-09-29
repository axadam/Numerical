//
//  IterativeResult.swift
//  Numerical
//
//  Created by Adam Roberts on 11/28/19.
//

/// A value that represents the product of an iterative process.
///
/// In addition to the last result of the iteration it stores the number of iterations that went into it.
public protocol IterativeValue {
    associatedtype Value
    /// The value produced by the iterative process
    var value: Value { get }
    
    /// A count of how many iterations or otherwise how much work went into computing the value
    var work: UInt { get }
}

/// When Optional wraps an IterativeValue that has a non-signalling exception value like NaN.
/// This is appropriate to use when you will only have a `nil` if no work has been done.
public extension Optional where Wrapped: IterativeValue, Wrapped.Value: FloatingPoint {
    /// The value produced by the iterative process or NaN if the Optional is empty
    var value: Wrapped.Value {
        switch self {
        case .none: return .nan
        case .some(let v): return v.value
        }
    }
    
    /// The work done in computing the value or 0 if the Optional is empty.
    var work: UInt {
        switch self {
        case .none: return 0
        case .some(let v): return v.work
        }
    }
}

public protocol Converging {
    var converged: Bool { get }
}

/// The value coming out of an iteratively converging process
public enum ConvergenceValue<Value>: IterativeValue, Converging {
    /// The process did not converge. Returns the last estimate and how much work was done
    case didNotConverge(work: UInt, estimate: Value)
    /// The process converged. Returns the converged estimate and how much work was done
    case converged(work: UInt, estimate: Value)
}

public extension ConvergenceValue {
    var value: Value {
        switch self {
        case .didNotConverge(_, let e): return e
        case .converged(_, let e): return e
        }
    }
    
    var work: UInt {
        switch self {
        case .didNotConverge(let w, _): return w
        case .converged(let w, _): return w
        }
    }
    
    var converged: Bool {
        switch self {
        case .didNotConverge: return false
        case .converged: return true
        }
    }
}

public extension Optional where Wrapped: Converging {
    var converged: Bool {
        switch self {
        case .none: return false
        case .some(let v): return v.converged
        }
    }
}
    var work: UInt { get }
}
