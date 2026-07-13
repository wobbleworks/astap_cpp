///----------------------------------------
///      @file live_stacking_test.cpp
///   @ingroup ASTAP++/tests
///     @brief Unit tests for the live-stacking colour-correction white-balance math.
/// @details Exercises the pure @c colour_correction_from_backgrounds (ports
///          @c auto_background_level). The full live-stacking loop has no headless
///          save hook, so the per-frame application is staff-reviewed; the factor
///          math is pinned here.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "stacking/live_colour_correction.h"

using astap::Background;
using astap::stacking::ColourCorrection;
using astap::stacking::colour_correction_from_backgrounds;

TEST_CASE("identical channels give an identity correction") {
    Background r; r.backgr = 10.0; r.star_level = 100.0;
    const Background g = r;
    const Background b = r;
    const auto cc = colour_correction_from_backgrounds(r, g, b);
    CHECK(cc.active);
    for (int i = 0; i < 3; ++i) {
        CHECK(cc.add[i]      == doctest::Approx(0.0));
        CHECK(cc.multiply[i] == doctest::Approx(1.0));
    }
    CHECK(cc.largest == doctest::Approx(1.0));
}

TEST_CASE("green and blue are white-balanced to red") {
    Background r; r.backgr = 10.0; r.star_level = 100.0;
    Background g; g.backgr = 20.0; g.star_level = 50.0;
    Background b; b.backgr = 5.0;  b.star_level = 200.0;
    const auto cc = colour_correction_from_backgrounds(r, g, b);

    // Red is the reference: identity.
    CHECK(cc.add[0]      == doctest::Approx(0.0));
    CHECK(cc.multiply[0] == doctest::Approx(1.0));

    // add_G = R.backgr·(G.star/R.star) − G.backgr = 10·0.5 − 20 = −15;  mul = 100/50 = 2.
    CHECK(cc.add[1]      == doctest::Approx(-15.0));
    CHECK(cc.multiply[1] == doctest::Approx(2.0));

    // add_B = 10·(200/100) − 5 = 15;  mul = 100/200 = 0.5.
    CHECK(cc.add[2]      == doctest::Approx(15.0));
    CHECK(cc.multiply[2] == doctest::Approx(0.5));

    // largest = max( sR=1, sG=(65535−15)·2/65535, sB=(65535+15)·0.5/65535 ) = sG.
    const double sG = (65535.0 - 15.0) * 2.0 / 65535.0;
    CHECK(cc.largest == doctest::Approx(sG));
    CHECK(cc.largest > 1.0);   // green scale dominates
}

TEST_CASE("a no-signal star_level==1 sentinel leaves that channel at identity") {
    Background r; r.backgr = 10.0; r.star_level = 100.0;
    Background g; g.backgr = 20.0; g.star_level = 1.0;   // sentinel: no green signal
    Background b; b.backgr = 5.0;  b.star_level = 200.0;
    const auto cc = colour_correction_from_backgrounds(r, g, b);

    CHECK(cc.add[1]      == doctest::Approx(0.0));   // green untouched
    CHECK(cc.multiply[1] == doctest::Approx(1.0));
    CHECK(cc.add[2]      == doctest::Approx(15.0));  // blue still balanced
    CHECK(cc.multiply[2] == doctest::Approx(0.5));
}

TEST_CASE("a no-signal red reference disables both green and blue corrections") {
    // R.star_level==1 makes both guards (`chan.star!=1 && R.star!=1`) false, so the
    // whole correction collapses to identity — but active is still true.
    Background r; r.backgr = 10.0; r.star_level = 1.0;   // sentinel: no red signal
    Background g; g.backgr = 20.0; g.star_level = 50.0;
    Background b; b.backgr = 5.0;  b.star_level = 200.0;
    const auto cc = colour_correction_from_backgrounds(r, g, b);
    CHECK(cc.active);
    for (int i = 0; i < 3; ++i) {
        CHECK(cc.add[i]      == doctest::Approx(0.0));
        CHECK(cc.multiply[i] == doctest::Approx(1.0));
    }
    CHECK(cc.largest == doctest::Approx(1.0));
}

TEST_CASE("largest selects the blue channel when blue's scale dominates") {
    // A small blue star_level → a large blue multiplier → sB is the max, so this
    // guards the channel-selection (a G/B index swap or wrong max would slip the
    // earlier green-dominant case).
    Background r; r.backgr = 10.0; r.star_level = 100.0;
    Background g; g.backgr = 20.0; g.star_level = 90.0;
    Background b; b.backgr = 5.0;  b.star_level = 25.0;
    const auto cc = colour_correction_from_backgrounds(r, g, b);
    CHECK(cc.add[2]      == doctest::Approx(-2.5));   // 10·(25/100) − 5
    CHECK(cc.multiply[2] == doctest::Approx(4.0));    // 100/25
    const double sB = (65535.0 + cc.add[2]) * cc.multiply[2] / 65535.0;
    CHECK(cc.largest == doctest::Approx(sB));
    CHECK(cc.largest > 3.0);   // blue dominates green (~1.11) and red (1)
}

TEST_CASE("a default-constructed ColourCorrection is an inactive identity") {
    // Pins the struct's member initializers — this is the value colour_correction_
    // factors returns for sub-3-colour input (that code path itself is not exercised
    // here, as it lives in the heavier live_stacking.cpp TU).
    const ColourCorrection cc;
    CHECK_FALSE(cc.active);
    for (int i = 0; i < 3; ++i) {
        CHECK(cc.add[i]      == doctest::Approx(0.0));
        CHECK(cc.multiply[i] == doctest::Approx(1.0));
    }
    CHECK(cc.largest == doctest::Approx(1.0));
}
