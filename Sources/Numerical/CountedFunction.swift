//
//  CountedFunction.swift
//  Numerical
//
//  Created by Adam Roberts on 4/28/20.
//

public class CountedFunction<I,O> {
    public private(set) var count: Int
    private let f: (I) -> O
    public init(f: @escaping (I) -> O) {
        self.f = f
        count = 0
    }
    public func callAsFunction(_ x: I) -> O { eval(x) }
    public func eval(_ x: I) -> O {
        count += 1
        return f(x)
    }
}
