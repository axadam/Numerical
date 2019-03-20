import XCTest

import NumericalTests

var tests = [XCTestCaseEntry]()
tests += NumericalTests.allTests()
XCTMain(tests)