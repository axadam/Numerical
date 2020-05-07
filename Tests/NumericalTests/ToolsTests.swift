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
        AssertLRE(r2.value, "1.414213562373095048801688724209698078569671875376948073176", resultStore: rs, table: t, testCase: "√2", field: f, annotation: "\(r2.iterations)\(r2.converged ? "" : "*")")
        
        let e  = continued_fraction(b0: 2, a: { $0 == 1 ? 1.0 : Double($0 - 1) }, b: { Double($0) })
        AssertLRE(e.value,  "2.718281828459045235360287471352662497757247093699959574966", resultStore: rs, table: t, testCase: "e", field: f, annotation: "\(e.iterations)\(e.converged ? "" : "*")")

        let pi = continued_fraction(b0: 3, a: { pow(Double(2 * $0 - 1),2) }, b: { _ in 6 }, maxIter: 1000)
        AssertLRE(pi.value, "3.141592653589793238462643383279502884197169399375105820974", digits: 10.1, resultStore: rs, table: t, testCase: "π", field: f, annotation: "\(pi.iterations)\(pi.converged ? "" : "*")")
    }
    
    func testChebyshev() {
        let t = "Chebyshev, e^x test case from Approximation Theory Approximation Practice, Trefethen"
        let f = "exp(x)"
        let coeffsExp = [
        1.266065877752008,
        1.130318207984970,
        0.271495339534077,
        0.044336849848664,
        0.005474240442094,
        0.000542926311914,
        0.000044977322954,
        0.000003198436462,
        0.000000199212481,
        0.000000011036772,
        0.000000000550590,
        0.000000000024980,
        0.000000000001039,
        0.000000000000040,
        0.000000000000001
        ]
        let e = { (x: Double) -> Double in
            chebyshev(poly: coeffsExp, z: x)
        }
        AssertLRE(e(0.5), "\(exp(0.5))", digits: 14.9, resultStore: rs, table: t, testCase: "0.5", field: f)
        AssertLRE(e(-0.5), "\(exp(-0.5))", digits: 14.7, resultStore: rs, table: t, testCase: "-0.5", field: f)
    }
    
    func testHorners() {
        let t = "Horner's Method of polynomial evaluation"
        let coefs: [Double] = [-19,7,-4,6]
        let r = evaluate_polynomial(poly: coefs, z: 3)
        AssertLRE(r, "128", exact: true, resultStore: rs, table: t, testCase: "-19 + 7x - 4x² + 6x³; x=3", field: "LRE")
    }
}
