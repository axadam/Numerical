//
//  ToolsTests.swift
//  NumericalTests
//
//  Created by Adam Roberts on 4/19/20.
//

import Foundation

import XCTest
import LogRelativeError
@testable import Numerical

final class ToolsTests: XCTestCase {
    override class func tearDown() {
        print(resultStore.md())
        super.tearDown()
    }

    static let resultStore: ResultStore = ResultStore()
    var rs: ResultStore { return type(of: self).resultStore }

    func testContinuedFraction() {
        let t = "Continued Fraction"
        let f = "LRE"
        
        let r2 = continued_fraction(b0: 1, a: { _ in 1 }, b: { _ in 2 })
        AssertLRE(r2, "1.414213562373095048801688724209698078569671875376948073176", resultStore: rs, table: t, testCase: "√2", field: f)
        
        let e  = continued_fraction(b0: 2, a: { $0 == 1 ? 1.0 : Double($0 - 1) }, b: { Double($0) })
        AssertLRE(e,  "2.718281828459045235360287471352662497757247093699959574966", resultStore: rs, table: t, testCase: "e", field: f)

        let pi = continued_fraction(b0: 3, a: { pow(2 * Double($0) - 1,2) }, b: { _ in 6 }, maxIter: 35000)
        AssertLRE(pi, "3.141592653589793238462643383279502884197169399375105820974", digits: 13.0, resultStore: rs, table: t, testCase: "π", field: f)
    }
}
