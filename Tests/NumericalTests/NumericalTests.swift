import XCTest
@testable import Numerical

final class NumericalTests: XCTestCase {
    func testInvErf() {
        let a = erf(1.2345)
        let b = invErf(a)
        XCTAssertEqual(b, 1.2345, accuracy: 1e-10)
    }
    
    func testIncompleteGamma() {
        let a = q_gamma(4, 0.7)
        // Q(4,7) = 0.9942465424077004166914...
        XCTAssertEqual(a, 0.9942465424077004166914, accuracy: 1e-10)
    }
    
    func testInverseIncompleteGamma() {
        let a = inv_p_gamma(4, 1 - 0.9942465424077004166914)
        XCTAssertEqual(a, 0.7, accuracy: 1e-10)
    }
    
    func testIncompleteBeta() {
        let a = beta_reg(x: 0.4, a: 3, b: 5)
        // I₀ٖ₄(3,5) = 0.580096
        XCTAssertEqual(a, 0.580096, accuracy: 1e-10)
    }
    
    func testInverseIncompleteBeta() {
        let a = inv_beta_reg(p: 0.580096, a: 3, b: 5)
        XCTAssertEqual(a, 0.4, accuracy: 1e-10)
    }

    static var allTests = [
        ("testInvErf", testInvErf),
    ]
}
