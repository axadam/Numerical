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
    // Poisson's example: perimeter of ellipse of 1/π by 0.6/π
    let pois = { (x: Double) in sqrt(1 - 0.36 * sin(x) * sin(x)) / (2 * Double.pi) }
    
    let t = "Quadrature"
    
    let fcube = "x³, [0,1]"
    let frecip = "1/x, [1,100]"
    let fiden = "x, [0,5000]"
    let fpois = "√(1 - 0.36sin²θ), [0,2π]"
    
    let vcube = "0.25"
    let vrecip = "4.605170185988091368035982909368728415202202977257545952066"
    let viden = "12500000"
    let vpois = "0.902779927772193884716139509461678950653101811833151562088"
    
    func testTrapezoidal() {
        let tc = "Trapezoidal"
        
        let c = trapezoidalQuadrature(range: 0...1, maxIter: 100, f: cube)
        let r = trapezoidalQuadrature(range: 1...100, maxIter: 1000, f: recip)
        let i = trapezoidalQuadrature(range: 0...5000, f: iden)
        let p = trapezoidalQuadrature(range: 0...(2 * Double.pi), f: pois)
        
        AssertLRE(c, vcube, exact: true, digits: 10.8, resultStore: rs, table: t, testCase: tc, field: fcube)
        AssertLRE(r, vrecip, digits: 10.9, resultStore: rs, table: t, testCase: tc, field: frecip)
        AssertLRE(i, viden, exact: true, resultStore: rs, table: t, testCase: tc, field: fiden)
        AssertLRE(p, vpois, resultStore: rs, table: t, testCase: tc, field: fpois)
    }
}
