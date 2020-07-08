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
    
    func testXMSinX() {
        AssertLRE(xmsin(9e-1) , "1.16673090372516611538617684286451376859e-1")
        AssertLRE(xmsin(3e-1) , "4.47979333866042489467925431497262632217e-3")
        AssertLRE(xmsin(1e-1) , "1.66583353171847693185801589377973010085e-4")
        AssertLRE(xmsin(3e-10), "4.499999999999999999979750000000000000000e-30")
        AssertLRE(xmsin(1e-10), "1.66666666666666666666583333333333333333e-31")
    }
}
