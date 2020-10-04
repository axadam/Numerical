//
//  BasicFunctionsTests.swift
//  NumericalTests
//
//  Created by Adam Roberts on 11/3/19.
//

import Foundation

import XCTest
import LogRelativeError
@testable import Numerical

final class BasicFunctionsTests: XCTestCase {
    override class func tearDown() {
        print(resultStore.md())
        super.tearDown()
    }

    static let resultStore: ResultStore = ResultStore()
    var rs: ResultStore { return type(of: self).resultStore }
    
    let t = "Basic Functions Near Zero"
    
    func testXMSinX() {
        let f = "x - sin(x)"
        AssertLRE(xmsin(9e-1) , "1.16673090372516611538617684286451376859e-1", resultStore: rs, table: t, testCase: "9e-1", field: f)
        AssertLRE(xmsin(3e-1) , "4.47979333866042489467925431497262632217e-3", resultStore: rs, table: t, testCase: "3e-1", field: f)
        AssertLRE(xmsin(1e-1) , "1.66583353171847693185801589377973010085e-4", resultStore: rs, table: t, testCase: "1e-1", field: f)
        AssertLRE(xmsin(3e-10), "4.499999999999999999979750000000000000000e-30", resultStore: rs, table: t, testCase: "3e-10", field: f)
        AssertLRE(xmsin(1e-10), "1.66666666666666666666583333333333333333e-31", resultStore: rs, table: t, testCase: "1e-10", field: f)
    }
    
    func testExpM1MX() {
        let f = "e^x - 1 - x"
        AssertLRE(expm1mx(9e-1) , "0.559603111156949663800126563602470695422", resultStore: rs, table: t, testCase: "9e-1", field: f)
        AssertLRE(expm1mx(3e-1) , "0.049858807576003103983744313328007330378", resultStore: rs, table: t, testCase: "3e-1", field: f)
        AssertLRE(expm1mx(1e-1) , "0.005170918075647624811707826490246668225", resultStore: rs, table: t, testCase: "1e-1", field: f)
        AssertLRE(expm1mx(3e-10), "4.5000000004500000000e-20", digits: 9.9, resultStore: rs, table: t, testCase: "3e-10", field: f)
        AssertLRE(expm1mx(1e-10), "5.000000000166666667e-21", digits: 10.4, resultStore: rs, table: t, testCase: "1e-10", field: f)
    }

    func testLog1PMX() {
        let f = "log(1 + x) - x"
        AssertLRE(log1pmx(9e-1) , "-0.258146113827605224008964022796510670364", resultStore: rs, table: t, testCase: "9e-1", field: f)
        AssertLRE(log1pmx(3e-1) , "-0.037635735532508947964504013119045602796", resultStore: rs, table: t, testCase: "3e-1", field: f)
        AssertLRE(log1pmx(1e-1) , "-0.004689820195675139956047876719234907779", digits: 14.5, resultStore: rs, table: t, testCase: "1e-1", field: f)
        AssertLRE(log1pmx(3e-10), "-4.4999999991000000002e-20", digits: 9.9, resultStore: rs, table: t, testCase: "3e-10", field: f)
        AssertLRE(log1pmx(1e-10), "-4.999999999666666667e-21", digits: 10.4, resultStore: rs, table: t, testCase: "1e-10", field: f)
    }

}
