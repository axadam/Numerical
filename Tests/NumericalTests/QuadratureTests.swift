//
//  QuadratureTests.swift
//  NumericalTests
//
//  Created by Adam Roberts on 4/24/20.
//

import Foundation


import XCTest
import LogRelativeError
@testable import Numerical

final class Quadrature: XCTestCase {
    override class func tearDown() {
        print(resultStore.md())
        super.tearDown()
    }

    static let resultStore: ResultStore = ResultStore()
    var rs: ResultStore { return type(of: self).resultStore }

    let cube = { (x: Double) in x * x * x }
    let recip = { (x: Double) in 1 / x }
    let iden = { (x: Double) in x }
    
    let t = "Quadrature"
    
    let fcube = "x³, [0,1]"
    let frecip = "1/x, [1,100]"
    let fiden = "x, [0,5000]"
    
    let vcube = "0.25"
    let vrecip = "4.605170185988091368035982909368728415202202977257545952066"
    let viden = "12500000"
    
    func testTrapezoidal() {
        let tc = "Trapezoidal"
        
        let c = trapezoidalQuadrature(range: 0...1, maxIter: 100, f: cube)
        let r = trapezoidalQuadrature(range: 1...100, maxIter: 1000, f: recip)
        let i = trapezoidalQuadrature(range: 0...5000, f: iden)

        AssertLRE(c, vcube, exact: true, digits: 10.8, resultStore: rs, table: t, testCase: tc, field: fcube)
        AssertLRE(r, vrecip, digits: 10.9, resultStore: rs, table: t, testCase: tc, field: frecip)
        AssertLRE(i, viden, exact: true, resultStore: rs, table: t, testCase: tc, field: fiden)
    }
}
