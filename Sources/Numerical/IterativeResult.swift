//
//  IterativeResult.swift
//  Numerical
//
//  Created by Adam Roberts on 11/28/19.
//

/// A value that represents the product of an iterative process.
///
/// In addition to the last result of the iteration it stores the number of iterations that went into it
/// and the exit state of the work. In some cases a process may return a result but exit in a state
/// other than complete success. The consumer can then decide what to do with the result.
///
/// Associated types are included for the product of the iteration and the exit state.
public struct IterativeResult<T,ExitState> {
    
    /// The product of the iterative process
    public let result: T
    
    /// The number of iterations
    public let iterations: Int
    
    /// The exit state of the iterative process
    public let exitState: ExitState
    
    @inlinable
    public init(iterations: Int, exitState: ExitState, result: T) {
        self.iterations = iterations
        self.exitState = exitState
        self.result = result
    }
}
