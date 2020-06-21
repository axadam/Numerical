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

    let tSum = "Summation"
    
    let hardSum: [Double] = { () -> [Double] in
        let rand = randomArray(seed: 123, n: 10000, exponentRange: -10...10)
        return (rand.map { [$0,-$0] }.flatMap() { $0 } + [1.0]).shuffled(seed: 456)
    }()

    let easySum: [Double] = { () -> [Double] in
        let rand = randomArray(seed: 123, n: 10000, exponentRange: -0...0)
        return (rand.map { [$0,-$0] }.flatMap() { $0 } + [1.0]).shuffled(seed: 456)
    }()
    
    let petersSum: [Double] = [1.0, 1e100, 1.0, -1e100]

    func testNaive() {
        let f = "Naive"
        AssertLRE(easySum.sum_naive(), "1.0", exact: true, digits: 12.9, resultStore: rs, table: tSum, testCase: f, field: "Easy")
        AssertLRE(hardSum.sum_naive(), "1.0", exact: true, digits: 4.1, resultStore: rs, table: tSum, testCase: f, field: "Hard")
        AssertLRE(petersSum.sum_naive(), "2.0", exact: true, digits: 0.0, resultStore: rs, table: tSum, testCase: f, field: "Peters")
    }
    
    func testPairwise() {
        let f = "Pairwise"
        AssertLRE(easySum.sum_pairwise(), "1.0", exact: true, resultStore: rs, table: tSum, testCase: f, field: "Easy")
        AssertLRE(hardSum.sum_pairwise(), "1.0", exact: true, digits: 4.2, resultStore: rs, table: tSum, testCase: f, field: "Hard")
        AssertLRE(petersSum.sum_pairwise(), "2.0", exact: true, digits: 0.0, resultStore: rs, table: tSum, testCase: f, field: "Peters")
    }
    
    func testKahan() {
        let f = "Kahan"
        AssertLRE(easySum.sum_kahan(), "1.0", exact: true, resultStore: rs, table: tSum, testCase: f, field: "Easy")
        AssertLRE(hardSum.sum_kahan(), "1.0", exact: true, digits: 5.5, resultStore: rs, table: tSum, testCase: f, field: "Hard")
        AssertLRE(petersSum.sum_kahan(), "2.0", exact: true, digits: 0.0, resultStore: rs, table: tSum, testCase: f, field: "Peters")
    }

    func testKBN() {
        let f = "Kahan-Babuška-Neumaier"
        AssertLRE(easySum.sum_kbn(), "1.0", exact: true, resultStore: rs, table: tSum, testCase: f, field: "Easy")
        AssertLRE(hardSum.sum_kbn(), "1.0", exact: true, resultStore: rs, table: tSum, testCase: f, field: "Hard")
        AssertLRE(petersSum.sum_kbn(), "2.0", exact: true, resultStore: rs, table: tSum, testCase: f, field: "Peters")
    }
    
    func testXMSinX() {
        AssertLRE(xmsin(9e-1) , "1.16673090372516611538617684286451376859e-1", digits: 14.3)
        AssertLRE(xmsin(3e-1) , "4.47979333866042489467925431497262632217e-3", digits: 14.2)
        AssertLRE(xmsin(1e-1) , "1.66583353171847693185801589377973010085e-4", digits: 14.9)
        AssertLRE(xmsin(3e-10), "4.499999999999999999979750000000000000000e-30")
        AssertLRE(xmsin(1e-10), "1.66666666666666666666583333333333333333e-31")
    }
}
