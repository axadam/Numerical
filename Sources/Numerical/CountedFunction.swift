//
//  CountedFunction.swift
//  Numerical
//
//  Created by Adam Roberts on 4/28/20.
//

/// Class to wrap a function and count how many times it is called
///
/// The counted execution function may be called by calling the instance itself in Swift 5.2 or later:
/// ```
///      func myFunc(x: Double) -> Double { return x + 1 }
///      let f = CountedFunction(myFunc)
///      let y = f(3.0) // 4.0
///      f.count // 1
/// ```
/// In older Swift versions you may call `f.eval(_ x: I)`
public class CountedFunction<I,O> {
    /// The number of times the function has been called through this instance
    public private(set) var count: Int
    private let f: (I) -> O
    public init(f: @escaping (I) -> O) {
        self.f = f
        count = 0
    }
    
    /// Call the function and increment the counter
    ///
    /// If the function takes more than one argument you must pass the arguments as a tuple
    public func callAsFunction(_ x: I) -> O { eval(x) }
    
    /// Call the function and increment the counter
    ///
    /// If the function takes more than one argument you must pass the arguments as a tuple
    public func eval(_ x: I) -> O {
        count += 1
        return f(x)
    }
}
