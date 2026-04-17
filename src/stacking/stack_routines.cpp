///----------------------------------------
///      @file stack_routines.cpp
///   @ingroup ASTAP++
///     @brief Core stacking algorithm implementations.
///    @author Ported from Han Kleijn's unit_stack_routines.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "stack_routines.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "stack.h"
#include "../core/demosaic.h"
#include "../core/fits.h"
#include "../core/globals.h"
#include "../core/photometry.h"
#include "../core/util.h"
#include "../core/wcs.h"
#include "../solving/star_align.h"

// Smedian isn't declared in any header; link_stubs.cpp provides the body.
namespace astap::core {
double Smedian(std::vector<double>& list, int len);
}  // namespace astap::core

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

using astap::bck;
using astap::solving::quad_star_distances1;
using astap::solving::quad_star_distances2;
using astap::solving::nr_references;
using astap::solving::nr_references2;
using astap::solving::find_quads;
using astap::solving::find_triples_using_quads;
using astap::solving::find_offset_and_rotation;
using astap::solving::reset_solution_vectors;
using astap::solving::bin_and_find_stars;
using astap::solving::report_binning;

///----------------------------------------
/// MARK: Forward declarations (TODO: replace with real headers)
///----------------------------------------

// Engine delegates pulled in from the real implementations.
using astap::core::load_fits;
using astap::core::save_fits;
using astap::core::get_background;
using astap::core::local_sd;
using astap::core::mode;
using astap::core::pixel_to_celestial;
using astap::core::celestial_to_pixel;
using astap::core::Smedian;
using astap::core::smedian;
using astap::core::strtofloat1;
using astap::core::strtofloat2;
using astap::core::strtoint2;

void demosaic_bayer(ImageArray& img);

// UI / progress (TODO: replace with ProgressReporter).
void memo2_message(const std::string& msg);
void progress_indicator(double value, const std::string& label);

// ListView subitem accessors.
// TODO: replace with a proper FileListModel.
std::string listview_subitem(int index, int column);
void        listview_set_subitem(int index, int column, const std::string& s);
void        listview_set_icon(int index, int column, int icon);
void        listview_mark_current(int index);

// Column indices.
constexpr auto L_X      = 0;
constexpr auto L_Y      = 1;
constexpr auto L_hfd    = 2;
constexpr auto L_result = 3;

///----------------------------------------
/// MARK: File-local helpers
///----------------------------------------

namespace {

///----------------------------------------
/// @brief Square a value.
///----------------------------------------

[[nodiscard]] constexpr auto sqr(auto v) noexcept { return v * v; }

///----------------------------------------
///  @brief Weighting factor for different exposure duration and gain.
/// @return Multiplicative weight relative to the reference frame.
///----------------------------------------

[[nodiscard]] double calc_weightF() {
    auto result = (head.exposure != 0.0)
                  ? head.exposure / head_ref.exposure
                  : 1.0;
                  
    if (head.egain != head_ref.egain) {
        // Check e-gain.
        auto gain1 = strtofloat1(head_ref.egain);
        auto gain2 = strtofloat1(head.egain);
        if (gain1 != 0.0) {
            result = result * gain2 / gain1;
        }
        if (std::abs(gain2 - gain1) > 0.01) {
            memo2_message(
                "Warning light with different EGAIN!! "
                + head.egain.substr(0, std::min<size_t>(5, head.egain.size()))
                + " instead of "
                + head_ref.egain.substr(0, std::min<size_t>(5, head_ref.egain.size()))
                + " [e-/ADU]. Will try to compensate accordingly.");
        }
    } else {
        if (head.gain != head_ref.gain) {
            memo2_message(
                "Warning light with different GAIN!! " + head.gain
                + " instead of " + head_ref.gain
                + ". Can not compensate unless EGAIN [e-/ADU] is added manually"
                  " to header.");
        }
    }
    return result;
}

///----------------------------------------
///      @brief Calculate the required output dimensions for mosaic stitching.
///   @param href Reference header.
///   @param h Current frame header.
/// @param[out] x_min Running minimum in reference-image X.
/// @param[out] x_max Running maximum in reference-image X.
/// @param[out] y_min Running minimum in reference-image Y.
/// @param[out] y_max Running maximum in reference-image Y.
///----------------------------------------

void calculate_required_dimensions(const Header& href, const Header& h,
                                   double& x_min, double& x_max,
                                   double& y_min, double& y_max) {
    auto ra  = 0.0;
    auto dec = 0.0;
    auto x   = 0.0;
    auto y   = 0.0;
    auto formalism = 0;  // TODO: from mainwindow.Polynomial1.itemindex
    
    // Left bottom.
    pixel_to_celestial(h, 1, 1, formalism, ra, dec);
    celestial_to_pixel(href, ra, dec, x, y);
    x_min = std::min(x_min, x); x_max = std::max(x_max, x);
    y_min = std::min(y_min, y); y_max = std::max(y_max, y);
    
    // Right bottom.
    pixel_to_celestial(h, h.width, 1, formalism, ra, dec);
    celestial_to_pixel(href, ra, dec, x, y);
    x_min = std::min(x_min, x); x_max = std::max(x_max, x);
    y_min = std::min(y_min, y); y_max = std::max(y_max, y);
    
    // Left top.
    pixel_to_celestial(h, 1, h.height, formalism, ra, dec);
    celestial_to_pixel(href, ra, dec, x, y);
    x_min = std::min(x_min, x); x_max = std::max(x_max, x);
    y_min = std::min(y_min, y); y_max = std::max(y_max, y);
    
    // Right top.
    pixel_to_celestial(h, h.width, h.height, formalism, ra, dec);
    celestial_to_pixel(href, ra, dec, x, y);
    x_min = std::min(x_min, x); x_max = std::max(x_max, x);
    y_min = std::min(y_min, y); y_max = std::max(y_max, y);
}

///----------------------------------------
///  @brief Minimum distance from a pixel to any image border.
///  @param fitsX Pixel X (FITS convention).
///  @param fitsY Pixel Y (FITS convention).
///  @param w Image width.
///  @param h Image height.
/// @return Distance in pixels.
///----------------------------------------

[[nodiscard]] float minimum_distance_borders(int fitsX, int fitsY, int w, int h) noexcept {
    auto r = static_cast<float>(std::min(fitsX, w - fitsX));
    r = std::min(static_cast<float>(fitsY), r);
    r = std::min(static_cast<float>(h - fitsY), r);
    return r;
}

///----------------------------------------
///  @brief Compute the drift vector for manual/ephemeris alignment.
/// @details For scale-1, 0..h/0..w range. @p c is the index into
///          @p files_to_process.
///  @param files_to_process File list.
///  @param c Index of the current frame.
///----------------------------------------

void calculate_manual_vector(std::span<FileToDo> files_to_process, int c) {
    auto ra1    = 0.0;
    auto dec1   = 0.0;
    auto x1     = 0.0;
    auto y1     = 0.0;
    auto shiftX = 0.0;
    auto shiftY = 0.0;
    
    if (head.cd1_1 == 0.0) {
        // Pure manual stacking.
        solution_vectorX[0] = 1.0;
        solution_vectorX[1] = 0.0;
        solution_vectorX[2] = referenceX
            - strtofloat2(listview_subitem(files_to_process[c].listview_index, L_X));
        solution_vectorY[0] = 0.0;
        solution_vectorY[1] = 1.0;
        solution_vectorY[2] = referenceY
            - strtofloat2(listview_subitem(files_to_process[c].listview_index, L_Y));
    } else {
        sincos(head.dec0, SIN_dec0, COS_dec0);
        astrometric_to_vector();
        
        // Touch the list view items (mirrors original).
        (void)listview_subitem(files_to_process[c].listview_index, L_X);
        (void)listview_subitem(files_to_process[c].listview_index, L_Y);
        
        // Convert the asteroid position to RA / DEC.
        pixel_to_celestial(head,
            strtofloat2(listview_subitem(files_to_process[c].listview_index, L_X)),
            strtofloat2(listview_subitem(files_to_process[c].listview_index, L_Y)),
            1, ra1, dec1);
            
        // Asteroid position in reference image pixels.
        celestial_to_pixel(head_ref, ra1, dec1, x1, y1);
        
        // Convert centre-based solution to origin at (0,0).
        if (solution_vectorX[0] < 0.0) {
            solution_vectorX[2] += head.width - 1;
        }
        if (solution_vectorY[1] < 0.0) {
            solution_vectorY[2] += head.height - 1;
        }
        
        shiftX = x1 - referenceX;
        shiftY = y1 - referenceY;
        
        solution_vectorX[2] -= shiftX;
        solution_vectorY[2] -= shiftY;
    }
}

///----------------------------------------
/// @brief Fill an image buffer with zeroes of the given shape.
/// @param img Image buffer to clear.
/// @param channels Number of colour channels.
/// @param height Image height.
/// @param width Image width.
///----------------------------------------

void clear_3d(ImageArray& img, int channels, int height, int width) {
    img.assign(channels,
               std::vector<std::vector<float>>(
                   height, std::vector<float>(width, 0.0f)));
}
  
}  // namespace

///----------------------------------------
/// MARK: astrometric_to_vector
///----------------------------------------

void astrometric_to_vector() {
    // SIP correction should be zero by definition.
    a_order = 0;
    
    // Only reliable for first order.
    calc_newx_newy(false, head.crpix1, head.crpix2);
    auto centerX = x_new_float;
    auto centerY = y_new_float;
    
    // One pixel in X.
    calc_newx_newy(false, head.crpix1 + 1, head.crpix2);
    solution_vectorX[0] = +(x_new_float - centerX);
    solution_vectorX[1] = -(y_new_float - centerY);
    
    // One pixel in Y.
    calc_newx_newy(false, head.crpix1, head.crpix2 + 1);
    solution_vectorY[0] = -(x_new_float - centerX);
    solution_vectorY[1] = +(y_new_float - centerY);
    
    // Range 0..width-1.
    solution_vectorX[2] = centerX - (head.crpix1 - 1);
    solution_vectorY[2] = centerY - (head.crpix2 - 1);
    
    // Check if the image is flipped (horizontal or vertical, but not both).
    auto flipped           = head.cd1_1 * head.cd2_2
                             - head.cd1_2 * head.cd2_1 > 0.0;
    auto flipped_reference = head_ref.cd1_1 * head_ref.cd2_2
                             - head_ref.cd1_2 * head_ref.cd2_1 > 0.0;
                             
    if (flipped != flipped_reference) {
        solution_vectorX[1] = -solution_vectorX[1];
        solution_vectorY[0] = -solution_vectorY[0];
    }
    
    // TODO: if solve_show_log1.checked -> log astrometric vector solution.
    memo2_message("Astrometric vector solution " + solution_str);
}

///----------------------------------------
/// MARK: initialise_calc_sincos_dec0
///----------------------------------------

void initialise_calc_sincos_dec0() {
    // Pre-compute sin/cos of dec0 to save work in the per-pixel loop.
    // For blink, head is used instead of head_ref.
    sincos(head.dec0, SIN_dec_ref, COS_dec_ref);
}

///----------------------------------------
/// MARK: test_bayer_matrix
///----------------------------------------

[[nodiscard]] bool test_bayer_matrix(const ImageArray& img) {
    constexpr auto steps = 100;
    
    const auto w = static_cast<int>(img[0][0].size());
    const auto h = static_cast<int>(img[0].size());
    
    auto middleY = h / 2;
    auto step_size = w / steps;
    if (step_size % 2 != 0) {
        step_size -= 1;  // even so it hits the 2x2 matrix
    }
    
    auto p11 = std::vector<double>(steps);
    auto p12 = std::vector<double>(steps);
    auto p21 = std::vector<double>(steps);
    auto p22 = std::vector<double>(steps);
    
    for (auto fitsX = 0; fitsX < steps; ++fitsX) {
        p11[fitsX] = img[0][middleY][step_size * fitsX];
        p12[fitsX] = img[0][middleY][step_size * fitsX + 1];
        p21[fitsX] = img[0][middleY + 1][step_size * fitsX];
        p22[fitsX] = img[0][middleY + 1][step_size * fitsX + 1];
    }
    
    auto m11 = Smedian(p11, steps);
    auto m12 = Smedian(p12, steps);
    auto m21 = Smedian(p21, steps);
    auto m22 = Smedian(p22, steps);
    auto lowest  = std::min(std::min(m11, m12), std::min(m21, m22));
    auto highest = std::max(std::max(m11, m12), std::max(m21, m22));
    
    return (highest - lowest) > 100.0;
}

///----------------------------------------
/// MARK: stack_LRGB
///----------------------------------------

void stack_LRGB(std::span<FileToDo> files_to_process, int& counter) {
    auto fitsX = 0;
    auto fitsY = 0;
    auto c = 0;
    auto width_max = 0;
    auto height_max = 0;
    auto x_new = 0;
    auto y_new = 0;
    auto binning = 0;
    auto max_stars = 0;
    auto col = 0;
    auto background_r = 0.0;
    auto background_g = 0.0;
    auto background_b = 0.0;
    auto background_l = 0.0;
    auto rgbsum = 0.0;
    auto red_f = 0.0;
    auto green_f = 0.0;
    auto blue_f = 0.0;
    auto value = 0.0;
    auto colr = 0.0;
    auto colg = 0.0;
    auto colb = 0.0;
    auto red_add = 0.0;
    auto green_add = 0.0;
    auto blue_add = 0.0;
    auto rr_factor = 0.0;
    auto rg_factor = 0.0;
    auto rb_factor = 0.0;
    auto gr_factor = 0.0;
    auto gg_factor = 0.0;
    auto gb_factor = 0.0;
    auto br_factor = 0.0;
    auto bg_factor = 0.0;
    auto bb_factor = 0.0;
    auto saturated_level = 0.0;
    auto hfd_min = 0.0;
    auto tempval = 0.0;
    auto aa = 0.0;
    auto bb2 = 0.0;
    auto cc2 = 0.0;
    auto dd = 0.0;
    auto ee = 0.0;
    auto ff = 0.0;
    auto init = false;
    auto solution = true;
    auto use_manual_align = false;
    auto use_ephemeris_alignment = false;
    auto use_astrometry_internal = false;
    auto use_sip = false;
    auto warning = std::string{};
    auto starlist1 = StarList{};
    auto starlist2 = StarList{};
    auto img_temp = ImageArray{};
    auto img_average = ImageArray{};
    
    use_manual_align        = astap::use_manual_align;
    use_ephemeris_alignment = astap::use_ephemeris_alignment;
    use_astrometry_internal = astap::use_astrometry_internal;
    hfd_min = std::max(0.8, strtofloat2(/*stackmenu1.min_star_size_stacking1.caption*/ ""));
    max_stars = strtoint2(/*stackmenu1.max_stars1.text*/ "", 500);
    use_sip   = /*stackmenu1.add_sip1.checked*/ false;
    
    counter = 0;
    jd_sum = 0.0;
    jd_start_first = 1e99;
    jd_end_last = 0.0;
    init = false;
    
    // Combining colours.
    memo2_message("Combining colours.");
    rr_factor = strtofloat2(/*rr1.text*/ "");
    rg_factor = strtofloat2(/*rg1.text*/ "");
    rb_factor = strtofloat2(/*rb1.text*/ "");
    gr_factor = strtofloat2(/*gr1.text*/ "");
    gg_factor = strtofloat2(/*gg1.text*/ "");
    gb_factor = strtofloat2(/*gb1.text*/ "");
    br_factor = strtofloat2(/*br1.text*/ "");
    bg_factor = strtofloat2(/*bg1.text*/ "");
    bb_factor = strtofloat2(/*bb1.text*/ "");
    
    red_add   = strtofloat2(/*red_filter_add1.text*/ "");
    green_add = strtofloat2(/*green_filter_add1.text*/ "");
    blue_add  = strtofloat2(/*blue_filter_add1.text*/ "");
    
    for (c = 0; c < static_cast<int>(files_to_process.size()); ++c) {
        // Should contain reference, R, G, B, RGB, L.
        if (c == 5) {
            // All colour files added: correct for number of samples per pixel.
            memo2_message("Correcting the number of pixels added together.");
            for (fitsY = 0; fitsY < height_max; ++fitsY) {
                for (fitsX = 0; fitsX < width_max;  ++fitsX) {
                    for (col = 0; col < 3; ++col) {
                        tempval = img_temp[col][fitsY][fitsX];
                        if (tempval > 0.0) {
                            // tempval > 1 is very rare (essentially only astrometric stacking).
                            img_average[col][fitsY][fitsX] =
                                500.0f + img_average[col][fitsY][fitsX] / tempval;
                        } else {
                            img_average[col][fitsY][fitsX] = 0.0f;
                        }
                    }
                }
            }
            memo2_message("Applying black spot filter on interim RGB image.");
            black_spot_filter(img_average);
        }
        
        if (files_to_process[c].name.empty()) {
            continue;
        }
        
        try {
            filename2 = files_to_process[c].name;
            if (c == 0) { memo2_message("Loading reference image: \"" + filename2 + "\"."); }
            if (c == 1) { memo2_message("Adding red file: \""   + filename2 + "\" to final image."); }
            if (c == 2) { memo2_message("Adding green file: \"" + filename2 + "\" to final image."); }
            if (c == 3) { memo2_message("Adding blue file: \""  + filename2 + "\" to final image."); }
            if (c == 4) { memo2_message("Adding RGB file: \""   + filename2 + "\" to final image."); }
            if (c == 5) { memo2_message("Using luminance file: \"" + filename2 + "\" for final image."); }
            
            // TODO: Application.ProcessMessages / esc handling.
            if (esc_pressed) {
                memo2_message("ESC pressed.");
                return;
            }
            if (!load_fits(filename2, true, true, !init, 0, astap::memo1_lines, head, img_loaded)) {
                memo2_message("Error loading " + filename2);
                return;
            }
            
            if (!init) {
                head_ref = head;
                initialise_calc_sincos_dec0();
            }
            
            if (!use_sip) {
                a_order = 0;
            }
            saturated_level = head.datamax_org * 0.97;
            
            if (c == 1) {
                get_background(0, img_loaded, true, false, bck);
                background_r = bck.backgr;
                counterR = head.light_count;  counterRdark = head.dark_count;
                counterRflat = head.flat_count; counterRbias = head.flatdark_count;
                exposureR = static_cast<int>(std::round(head.exposure));
                temperatureR = head.set_temperature;
            }
            if (c == 2) {
                get_background(0, img_loaded, true, false, bck);
                background_g = bck.backgr;
                counterG = head.light_count;  counterGdark = head.dark_count;
                counterGflat = head.flat_count; counterGbias = head.flatdark_count;
                exposureG = static_cast<int>(std::round(head.exposure));
                temperatureG = head.set_temperature;
            }
            if (c == 3) {
                get_background(0, img_loaded, true, false, bck);
                background_b = bck.backgr;
                counterB = head.light_count;  counterBdark = head.dark_count;
                counterBflat = head.flat_count; counterBbias = head.flatdark_count;
                exposureB = static_cast<int>(std::round(head.exposure));
                temperatureB = head.set_temperature;
            }
            if (c == 4) {
                get_background(0, img_loaded, true, false, bck);
                background_r = bck.backgr;
                background_g = background_r;
                background_b = background_r;
                counterRGB = head.light_count;  counterRGBdark = head.dark_count;
                counterRGBflat = head.flat_count; counterRGBbias = head.flatdark_count;
                exposureRGB = static_cast<int>(std::round(head.exposure));
                temperatureRGB = head.set_temperature;
            }
            if (c == 5) {
                get_background(0, img_loaded, true, false, bck);
                background_l = bck.backgr;
                counterL = head.light_count;  counterLdark = head.dark_count;
                counterLflat = head.flat_count; counterLbias = head.flatdark_count;
                exposureL = static_cast<int>(std::round(head.exposure));
                temperatureL = head.set_temperature;
            }
            
            if (use_astrometry_internal) {
                memo2_message("Preparing astrometric solution for interim file: " + filename2);
                if (head.cd1_1 == 0.0) {
                    solution = update_solution_and_save(img_loaded, head, astap::memo1_lines);
                } else {
                    solution = true;
                }
                if (!solution) {
                    memo2_message("Abort, No astrometric solution for " + filename2);
                    return;
                }
            } else if (!init) {
                if (use_manual_align || use_ephemeris_alignment) {
                    referenceX = strtofloat2(listview_subitem(files_to_process[c].listview_index, L_X));
                    referenceY = strtofloat2(listview_subitem(files_to_process[c].listview_index, L_Y));
                } else {
                    binning = report_binning(head.height);
                    bin_and_find_stars(img_loaded, binning, 1, hfd_min, max_stars,
                                       true, starlist1, warning);
                    find_quads(starlist1, quad_star_distances1);
                }
            }
            
            if (!init) {
                height_max = head.height;
                width_max  = head.width;
                clear_3d(img_average, 3, height_max, width_max);
                clear_3d(img_temp,    3, height_max, width_max);
            }
            
            solution = true;
            if (use_astrometry_internal) {
                sincos(head.dec0, SIN_dec0, COS_dec0);
            } else {
                if (init) {
                    // Second image.
                    if (use_manual_align || use_ephemeris_alignment) {
                        calculate_manual_vector(files_to_process, c);
                    } else {
                        bin_and_find_stars(img_loaded, binning, 1, hfd_min, max_stars,
                                           true, starlist2, warning);
                        find_quads(starlist2, quad_star_distances2);
                        if (find_offset_and_rotation(3, astap::quad_tolerance)) {
                            memo2_message(std::to_string(nr_references) + " of "
                                         + std::to_string(nr_references2)
                                         + " quads selected matching within tolerance. "
                                         + solution_str);
                        } else {
                            memo2_message("Not enough quad matches <3 or inconsistent solution, skipping this image.");
                            files_to_process[c].name.clear();
                            solution = false;
                            listview_set_icon(files_to_process[c].listview_index, L_result, 6);
                            listview_set_subitem(files_to_process[c].listview_index, L_result, "no solution");
                        }
                    }
                } else {
                    reset_solution_vectors(1);
                }
            }
            init = true;
            
            if (c != 0 && solution) {
                // Do not add reference channel (c=0, typically luminance).
                ++counter;
                date_to_jd(head.date_obs, head.date_avg, head.exposure);
                jd_start_first = std::min(jd_start, jd_start_first);
                jd_end_last    = std::max(jd_end,   jd_end_last);
                jd_sum        += jd_mid;
                
                if (use_astrometry_internal) {
                    astrometric_to_vector();
                }
                
                aa  = solution_vectorX[0];
                bb2 = solution_vectorX[1];
                cc2 = solution_vectorX[2];
                dd  = solution_vectorY[0];
                ee  = solution_vectorY[1];
                ff  = solution_vectorY[2];
                
                for (fitsY = 0; fitsY < head.height; ++fitsY) {
                    for (fitsX = 0; fitsX < head.width; ++fitsX) {
                        x_new = static_cast<int>(std::round(aa * fitsX + bb2 * fitsY + cc2));
                        y_new = static_cast<int>(std::round(dd * fitsX + ee * fitsY + ff));
                        
                        if (x_new < 0 || x_new > width_max - 1
                            || y_new < 0 || y_new > height_max - 1) {
                            continue;
                        }
                        
                        if (c == 1) {
                            // Red channel.
                            value = img_loaded[0][fitsY][fitsX];
                            if (value > saturated_level) {
                                for (col = 0; col < 3; ++col) {
                                    img_temp[col][y_new][x_new] = -9.0f;
                                }
                            } else {
                                value = value - background_r;
                                if (rr_factor > 1e-5) {
                                    img_average[0][y_new][x_new] += static_cast<float>(rr_factor * value);
                                    img_temp[0][y_new][x_new] += 1.0f;
                                }
                                if (rg_factor > 1e-5) {
                                    img_average[1][y_new][x_new] += static_cast<float>(rg_factor * value);
                                    img_temp[1][y_new][x_new] += 1.0f;
                                }
                                if (rb_factor > 1e-5) {
                                    img_average[2][y_new][x_new] += static_cast<float>(rb_factor * value);
                                    img_temp[2][y_new][x_new] += 1.0f;
                                }
                            }
                        } else if (c == 2) {
                            // Green channel.
                            value = img_loaded[0][fitsY][fitsX];
                            if (value > saturated_level) {
                                for (col = 0; col < 3; ++col) {
                                    img_temp[col][y_new][x_new] = -9.0f;
                                }
                            } else {
                                value = value - background_g;
                                if (gr_factor > 1e-5) {
                                    img_average[0][y_new][x_new] += static_cast<float>(gr_factor * value);
                                    img_temp[0][y_new][x_new] += 1.0f;
                                }
                                if (gg_factor > 1e-5) {
                                    img_average[1][y_new][x_new] += static_cast<float>(gg_factor * value);
                                    img_temp[1][y_new][x_new] += 1.0f;
                                }
                                if (gb_factor > 1e-5) {
                                    img_average[2][y_new][x_new] += static_cast<float>(gb_factor * value);
                                    img_temp[2][y_new][x_new] += 1.0f;
                                }
                            }
                        } else if (c == 3) {
                            // Blue channel.
                            value = img_loaded[0][fitsY][fitsX];
                            if (value > saturated_level) {
                                for (col = 0; col < 3; ++col) {
                                    img_temp[col][y_new][x_new] = -9.0f;
                                }
                            } else {
                                value = value - background_b;
                                if (br_factor > 1e-5) {
                                    img_average[0][y_new][x_new] += static_cast<float>(br_factor * value);
                                    img_temp[0][y_new][x_new] += 1.0f;
                                }
                                if (bg_factor > 1e-5) {
                                    img_average[1][y_new][x_new] += static_cast<float>(bg_factor * value);
                                    img_temp[1][y_new][x_new] += 1.0f;
                                }
                                if (bb_factor > 1e-5) {
                                    img_average[2][y_new][x_new] += static_cast<float>(bb_factor * value);
                                    img_temp[2][y_new][x_new] += 1.0f;
                                }
                            }
                        } else if (c == 4) {
                            // RGB image, naxis3=3.
                            img_average[0][y_new][x_new] += static_cast<float>(img_loaded[0][fitsY][fitsX] - background_r);
                            img_temp[0][y_new][x_new] += 1.0f;
                            img_average[1][y_new][x_new] += static_cast<float>(img_loaded[1][fitsY][fitsX] - background_g);
                            img_temp[1][y_new][x_new] += 1.0f;
                            img_average[2][y_new][x_new] += static_cast<float>(img_loaded[2][fitsY][fitsX] - background_b);
                            img_temp[2][y_new][x_new] += 1.0f;
                        } else if (c == 5) {
                            // Luminance: r := l*(0.33+r)/(r+g+b).
                            colr = img_average[0][y_new][x_new] - 475 + red_add;
                            colg = img_average[1][y_new][x_new] - 475 + green_add;
                            colb = img_average[2][y_new][x_new] - 475 + blue_add;
                            rgbsum = colr + colg + colb;
                            if (rgbsum < 0.1) {
                                rgbsum = 0.1;
                                red_f = rgbsum / 3;
                                green_f = red_f;
                                blue_f = red_f;
                            } else {
                                red_f   = std::clamp(colr / rgbsum, 0.0, 1.0);
                                green_f = std::clamp(colg / rgbsum, 0.0, 1.0);
                                blue_f  = std::clamp(colb / rgbsum, 0.0, 1.0);
                            }
                            img_average[0][y_new][x_new] = static_cast<float>(1000.0 + (img_loaded[0][fitsY][fitsX] - background_l) * red_f);
                            img_average[1][y_new][x_new] = static_cast<float>(1000.0 + (img_loaded[0][fitsY][fitsX] - background_l) * green_f);
                            img_average[2][y_new][x_new] = static_cast<float>(1000.0 + (img_loaded[0][fitsY][fitsX] - background_l) * blue_f);
                        }
                    }
                }
            }
            
            progress_indicator(94.0 + c, " LRGB");  // 95..99
        } catch (...) {
            // Silently continue on error.
        }
    }
    
    if (counter != 0) {
        // Restore reference solution (only works if not oversized).
        head = head_ref;
        head.naxis3 = 3;
        head.naxis  = 3;
        img_loaded  = img_average;
        head.width  = width_max;
        head.height = height_max;
        sum_exp = static_cast<double>(exposureR + exposureG + exposureG
                                      + exposureL + exposureRGB);
    }
}

///----------------------------------------
/// MARK: stack_average
///----------------------------------------

void stack_average(int process_as_osc,
                   std::span<FileToDo> files_to_process,
                   int& counter) {
    auto fitsX = 0;
    auto fitsY = 0;
    auto c = 0;
    auto width_max = 0;
    auto height_max = 0;
    auto old_width = 0;
    auto old_height = 0;
    auto old_naxis3 = 0;
    auto x_new = 0;
    auto y_new = 0;
    auto col = 0;
    auto binning = 0;
    auto max_stars = 0;
    auto background_correction = 0.0;
    auto weightF = 0.0;
    auto hfd_min = 0.0;
    auto aa = 0.0;
    auto bb = 0.0;
    auto cc = 0.0;
    auto dd = 0.0;
    auto ee = 0.0;
    auto ff = 0.0;
    auto init = false;
    auto solution = true;
    auto use_manual_align = false;
    auto use_ephemeris_alignment = false;
    auto use_astrometry_internal = false;
    auto use_sip = false;
    auto tempval = 0.0f;
    auto warning = std::string{};
    auto starlist1 = StarList{};
    auto starlist2 = StarList{};
    auto img_temp = ImageArray{};
    auto img_average = ImageArray{};
    
    use_manual_align        = astap::use_manual_align;
    use_ephemeris_alignment = astap::use_ephemeris_alignment;
    use_astrometry_internal = astap::use_astrometry_internal;
    hfd_min   = std::max(0.8, astap::hfd_min_setting);
    max_stars = astap::max_stars_setting;
    use_sip   = astap::add_sip;
    
    counter = 0;
    sum_exp = 0.0;
    sum_temp = 0.0;
    jd_sum = 0.0;
    jd_start_first = 1e99;
    jd_end_last = 0.0;
    
    for (c = 0; c < static_cast<int>(files_to_process.size()); ++c) {
        if (files_to_process[c].name.empty()) {
            continue;
        }
        try {
            listview_mark_current(files_to_process[c].listview_index);
            filename2 = files_to_process[c].name;
            
            if (esc_pressed) {
                memo2_message("ESC pressed.");
                return;
            }
            if (!load_fits(filename2, true, true, !init, 0, astap::memo1_lines, head, img_loaded)) {
                memo2_message("Error loading " + filename2);
                return;
            }
            
            if (!init) {
                old_width  = head.width;
                old_height = head.height;
                old_naxis3 = head.naxis3;
                
                // TODO: FITS comment write.
                head_ref = head;
                initialise_calc_sincos_dec0();
                if (bayerpat.empty() && process_as_osc == 2) {
                    if (/*bayer_pattern1.Text*/ std::string("auto") == "auto") {
                        memo2_message("Warning, Bayer colour pattern not in the header!");
                    } else if (!test_bayer_matrix(img_loaded)) {
                        memo2_message("Warning, grayscale image converted to colour!");
                    }
                }
            } else {
                if (old_width != head.width || old_height != head.height) {
                    memo2_message("Warning different size image!");
                }
                if (head.naxis3 > old_naxis3) {
                    memo2_message("Abort!! Can't combine colour to mono files.");
                    return;
                }
            }
            
            if (!use_sip) {
                a_order = 0;
            }
            (void)apply_dark_and_flat(img_loaded, head);
            
            memo2_message("Adding file: " + std::to_string(counter + 1)
                          + " \"" + filename2 + "\" to average.");
                          
            if (esc_pressed) {
                return;
            }
            
            if (process_as_osc > 0) {
                if (head.naxis3 > 1) {
                    memo2_message("Warning, light is already in colour! Will skip demosaic.");
                } else {
                    demosaic_bayer(img_loaded);
                }
            }
            
            if (!init) {
                jd_mid_reference = jd_mid;
                height_max = head.height;
                width_max  = head.width;
                binning    = report_binning(head.height);
                
                clear_3d(img_average, head.naxis3, height_max, width_max);
                clear_3d(img_temp,    1,           height_max, width_max);
                
                if (use_manual_align || use_ephemeris_alignment) {
                    referenceX = strtofloat2(listview_subitem(files_to_process[c].listview_index, L_X));
                    referenceY = strtofloat2(listview_subitem(files_to_process[c].listview_index, L_Y));
                } else if (!use_astrometry_internal) {
                    bin_and_find_stars(img_loaded, binning, 1, hfd_min, max_stars,
                                       true, starlist1, warning);
                    find_quads(starlist1, quad_star_distances1);
                    pedestal_s = bck.backgr;
                    if (pedestal_s < 500.0) {
                        pedestal_s = 500.0;
                    }
                    background_correction = pedestal_s - bck.backgr;
                    head.datamax_org += background_correction;
                    if (head.datamax_org > 0xFFFF) {
                        head.datamax_org = 0xFFFF;
                    }
                    head.pedestal = background_correction;
                }
            }
            
            solution = true;
            if (use_astrometry_internal) {
                sincos(head.dec0, SIN_dec0, COS_dec0);
            } else {
                if (init) {
                    if (use_manual_align || use_ephemeris_alignment) {
                        calculate_manual_vector(files_to_process, c);
                    } else {
                        bin_and_find_stars(img_loaded, binning, 1, hfd_min, max_stars,
                                           true, starlist2, warning);
                        background_correction = pedestal_s - bck.backgr;
                        head.datamax_org += background_correction;
                        if (head.datamax_org > 0xFFFF) {
                            head.datamax_org = 0xFFFF;
                        }
                        head.pedestal = background_correction;
                        find_quads(starlist2, quad_star_distances2);
                        if (find_offset_and_rotation(3, astap::quad_tolerance)) {
                            memo2_message(std::to_string(nr_references) + " of "
                                          + std::to_string(nr_references2)
                                          + " quads selected. " + solution_str);
                        } else {
                            memo2_message("Not enough quad matches <3 or inconsistent solution, skipping this image.");
                            files_to_process[c].name.clear();
                            solution = false;
                            listview_set_icon(files_to_process[c].listview_index, L_result, 6);
                            listview_set_subitem(files_to_process[c].listview_index, 2, "no solution");
                        }
                    }
                } else {
                    reset_solution_vectors(1);
                }
            }
            init = true;
            
            if (solution) {
                ++counter;
                sum_exp  += head.exposure;
                sum_temp += head.set_temperature;
                weightF   = calc_weightF();
                
                date_to_jd(head.date_obs, head.date_avg, head.exposure);
                jd_start_first = std::min(jd_start, jd_start_first);
                jd_end_last    = std::max(jd_end, jd_end_last);
                jd_sum        += jd_mid;
                airmass_sum   += airmass;
                
                if (use_astrometry_internal) {
                    astrometric_to_vector();
                }
                
                aa = solution_vectorX[0];
                bb = solution_vectorX[1];
                cc = solution_vectorX[2];
                dd = solution_vectorY[0];
                ee = solution_vectorY[1];
                ff = solution_vectorY[2];
                
                for (fitsY = 0; fitsY < head.height; ++fitsY) {
                    for (fitsX = 0; fitsX < head.width; ++fitsX) {
                        x_new = static_cast<int>(std::round(aa * fitsX + bb * fitsY + cc));
                        y_new = static_cast<int>(std::round(dd * fitsX + ee * fitsY + ff));
                        if (x_new < 0 || x_new > width_max - 1
                            || y_new < 0 || y_new > height_max - 1) {
                            continue;
                        }
                        for (col = 0; col < head.naxis3; ++col) {
                            img_average[col][y_new][x_new] += static_cast<float>(img_loaded[col][fitsY][fitsX] * weightF);
                        }
                        img_temp[0][y_new][x_new] += static_cast<float>(weightF);
                    }
                }
            }
            
            progress_indicator(10.0 + 89.0 * counter / images_selected, " Stacking");
        } catch (...) {
            // Silently continue on error.
        }
    }
    
    if (counter != 0) {
        head_ref.naxis3 = head.naxis3;
        head_ref.naxis  = head.naxis;
        head_ref.datamax_org = head.datamax_org;
        head = head_ref;
        head.height = height_max;
        head.width  = width_max;
        img_loaded.assign(head.naxis3,
                          std::vector<std::vector<float>>(
                              head.height, std::vector<float>(head.width, 0.0f)));
                              
        for (fitsY = 0; fitsY < head.height; ++fitsY) {
            for (fitsX = 0; fitsX < head.width; ++fitsX) {
                tempval = img_temp[0][fitsY][fitsX];
                for (col = 0; col < head.naxis3; ++col) {
                    if (tempval != 0.0f) {
                        img_loaded[col][fitsY][fitsX] =
                            static_cast<float>(background_correction + img_average[col][fitsY][fitsX] / tempval);
                    } else {
                        // Black-spot / missing-value filter.
                        if (fitsX > 0 && img_temp[0][fitsY][fitsX - 1] != 0.0f) {
                            img_loaded[col][fitsY][fitsX] = static_cast<float>(background_correction + img_loaded[col][fitsY][fitsX - 1]);
                        } else if (fitsY > 0 && img_temp[0][fitsY - 1][fitsX] != 0.0f) {
                            img_loaded[col][fitsY][fitsX] = static_cast<float>(background_correction + img_loaded[col][fitsY - 1][fitsX]);
                        } else {
                            img_loaded[col][fitsY][fitsX] = 0.0f;
                        }
                    }
                }
            }
        }
    }
}

///----------------------------------------
/// MARK: stack_mosaic
///----------------------------------------

void stack_mosaic(int process_as_osc,
                  std::span<FileToDo> files_to_process,
                  double max_dev_backgr,
                  int& counter) {
    auto fitsX = 0;
    auto fitsY = 0;
    auto c = 0;
    auto width_max = 0;
    auto height_max = 0;
    auto x_new = 0;
    auto y_new = 0;
    auto col = 0;
    auto cropW = 0;
    auto cropH = 0;
    auto iterations = 0;
    auto greylevels = 0;
    auto count = 0;
    auto formalism = 0;  // TODO: mainwindow.Polynomial1.itemindex
    auto value = 0.0;
    auto dummy = 0.0;
    auto median = 0.0;
    auto median2 = 0.0;
    auto delta_median = 0.0;
    auto correction = 0.0;
    auto maxlevel = 0.0;
    auto mean = 0.0;
    auto noise = 0.0;
    auto hotpixels = 0.0;
    auto coverage = 0.0;
    auto raMiddle = 0.0;
    auto decMiddle = 0.0;
    auto x_min = 0.0;
    auto x_max = 0.0;
    auto y_min = 0.0;
    auto y_max = 0.0;
    auto total_fov = 0.0;
    auto fw = 0.0;
    auto fh = 0.0;
    auto tempval = 0.0f;
    auto init = false;
    auto vector_based = false;
    auto merge_overlap = false;
    auto equalise_background = false;
    auto use_sip = false;
    auto background_correction = std::array<double, 3>{};
    auto background_correction_center = std::array<double, 3>{};
    auto background = std::array<double, 3>{};
    auto counter_overlap = std::array<int, 3>{};
    auto bckArr = std::array<double, 4>{};
    auto oldsip = false;
    auto img_temp = ImageArray{};
    auto img_average = ImageArray{};
    
    // Find dimensions of this package.
    memo2_message("Analysing and calculating celestial field-of-view dimensions.");
    
    count = 0;
    total_fov = 0.0;
    init = false;
    oldsip = sip;
    sip = false;  // prevent large error due to SIP outside image
    
    for (c = 0; c < static_cast<int>(files_to_process.size()); ++c) {
        if (files_to_process[c].name.empty()) {
            continue;
        }
        if (!load_fits(files_to_process[c].name, true, false, false, 0, astap::memo1_lines, head, img_loaded)) {
            memo2_message("Error loading " + filename2);
            return;
        }
        if (!init) {
            head_ref = head;
            init = true;
        }
        calculate_required_dimensions(head_ref, head, x_min, x_max, y_min, y_max);
        total_fov += head.cdelt1 * head.cdelt2 * head.width * head.height;
        ++count;
    }
    sip = oldsip;
    if (std::abs(x_max - x_min) < 1.0) {
        memo2_message("Abort. Failed to calculate mosaic dimensions!");
        return;
    }
    
    merge_overlap       = /*merge_overlap1.checked*/ false;
    equalise_background = /*Equalise_background1.checked*/ false;
    counter = 0;
    sum_exp = 0.0;
    sum_temp = 0.0;
    jd_sum = 0.0;
    jd_start_first = 1e99;
    jd_end_last = 0.0;
    init = false;
    use_sip = /*add_sip1.checked*/ false;
    dummy = 0.0;
    
    // TODO: classify-object warning.
    
    for (c = 0; c < static_cast<int>(files_to_process.size()); ++c) {
        if (files_to_process[c].name.empty()) {
            continue;
        }
        try {
            listview_mark_current(files_to_process[c].listview_index);
            filename2 = files_to_process[c].name;
            
            if (esc_pressed) {
                memo2_message("ESC pressed.");
                return;
            }
            if (!load_fits(filename2, true, true, !init, 0, astap::memo1_lines, head, img_loaded)) {
                memo2_message("Error loading " + filename2);
                return;
            }
            
            if (init) {
                if (head.naxis3 > static_cast<int>(img_average.size())) {
                    memo2_message("Abort!! Can't combine mono and colour files.");
                    return;
                }
            }
            
            if (!init) {
                head_ref = head;
                fw = head.cdelt1 * std::abs(x_max - x_min);
                fh = head.cdelt2 * std::abs(y_max - y_min);
                coverage = total_fov / (fw * fh);
                if (coverage < 0.5) {
                    memo2_message("Abort!! Too many missing tiles. Coverage too low.");
                    return;
                }
                pixel_to_celestial(head, (x_min + x_max) / 2.0,
                                   (y_min + y_max) / 2.0, formalism, raMiddle, decMiddle);
                sincos(decMiddle, SIN_dec_ref, COS_dec_ref);
                head_ref.ra0    = raMiddle;
                head_ref.crpix1 = std::abs(x_max - x_min) / 2.0;
                head_ref.crpix2 = std::abs(y_max - y_min) / 2.0;
            }
            
            if (!use_sip) {
                a_order = 0;
            }
            
            memo2_message("Adding file: " + std::to_string(counter + 1)
                          + " \"" + filename2 + "\" to mosaic.");
            if (a_order == 0) {
                memo2_message("Warning. Image distortion correction not working.");
            }
            
            if (esc_pressed) {
                return;
            }
            
            if (process_as_osc > 0) {
                if (head.naxis3 > 1) {
                    memo2_message("Warning, light is already in colour! Will skip demosaic.");
                } else {
                    demosaic_bayer(img_loaded);
                }
            }
            
            if (!init) {
                width_max  = std::abs(static_cast<int>(std::round(x_max - x_min)));
                height_max = std::abs(static_cast<int>(std::round(y_max - y_min)));
                clear_3d(img_average, head.naxis3, height_max, width_max);
                clear_3d(img_temp,    1,           height_max, width_max);
            }
            
            for (col = 0; col < head.naxis3; ++col) {
                if (equalise_background) {
                    bckArr[0] = mode(img_loaded, false, col,
                                     0, static_cast<int>(std::round(0.2 * head.width)),
                                     0, static_cast<int>(std::round(0.2 * head.height)),
                                     32000, greylevels);
                    bckArr[1] = mode(img_loaded, false, col,
                                     0, static_cast<int>(std::round(0.2 * head.width)),
                                     static_cast<int>(std::round(0.8 * head.height)),
                                     head.height - 1, 32000, greylevels);
                    bckArr[2] = mode(img_loaded, false, col,
                                     static_cast<int>(std::round(0.8 * head.width)),
                                     head.width - 1,
                                     0, static_cast<int>(std::round(0.2 * head.height)),
                                     32000, greylevels);
                    bckArr[3] = mode(img_loaded, false, col,
                                     static_cast<int>(std::round(0.8 * head.width)),
                                     head.width - 1,
                                     static_cast<int>(std::round(0.8 * head.height)),
                                     head.height - 1, 32000, greylevels);
                    background[col] = smedian(bckArr, 4);
                    background_correction_center[col] = 1000.0 - background[col];
                } else {
                    background[col] = 0.0;
                    background_correction_center[col] = 0.0;
                }
            }
            
            // Always astrometric for mosaic.
            sincos(head.dec0, SIN_dec0, COS_dec0);
            
            ++counter;
            sum_exp  += head.exposure;
            sum_temp += head.set_temperature;
            date_to_jd(head.date_obs, head.date_avg, head.exposure);
            jd_start_first = std::min(jd_start, jd_start_first);
            jd_end_last    = std::max(jd_end, jd_end_last);
            jd_sum        += jd_mid;
            
            vector_based = false;
            if (a_order == 0) {
                astrometric_to_vector();
                vector_based = true;
            }
            ap_order = 0;  // don't correct for RA to XY for mosaic
            
            cropW = static_cast<int>(/*mosaic_crop1.Position*/ 0 * head.width  / 200);
            cropH = static_cast<int>(/*mosaic_crop1.Position*/ 0 * head.height / 200);
            
            background_correction = {0.0, 0.0, 0.0};
            
            if (init) {
                counter_overlap = {0, 0, 0};
                for (fitsY = 1 + cropH; fitsY <= head.height - (1 + 1 + cropH); ++fitsY) {
                    for (fitsX = 1 + cropW; fitsX <= head.width - (1 + 1 + cropW); ++fitsX) {
                        calc_newx_newy(vector_based, fitsX, fitsY);
                        x_new = static_cast<int>(std::round(x_new_float));
                        y_new = static_cast<int>(std::round(y_new_float));
                        if (x_new < 0 || x_new > width_max - 1
                            || y_new < 0 || y_new > height_max - 1) {
                            continue;
                        }
                        if (img_loaded[0][fitsY][fitsX] <= 0.0001) {
                            continue;  // black border
                        }
                        if (img_average[0][y_new][x_new] != 0.0f) {
                            for (col = 0; col < head.naxis3; ++col) {
                                correction = std::round(img_average[col][y_new][x_new]
                                    - (img_loaded[col][fitsY][fitsX]
                                       + background_correction_center[col]));
                                if (std::abs(correction) < max_dev_backgr * 1.5) {
                                    background_correction[col] += correction;
                                    counter_overlap[col]       += 1;
                                }
                            }
                        }
                    }
                }
                for (auto k = 0; k < 3; ++k) {
                    if (counter_overlap[k] > 0) {
                        background_correction[k] /= counter_overlap[k];
                    }
                }
            }
            init = true;
            
            for (fitsY = 1 + cropH; fitsY <= head.height - (1 + 1 + cropH); ++fitsY) {
                for (fitsX = 1 + cropW; fitsX <= head.width - (1 + 1 + cropW); ++fitsX) {
                    calc_newx_newy(vector_based, fitsX, fitsY);
                    x_new = static_cast<int>(std::round(x_new_float));
                    y_new = static_cast<int>(std::round(y_new_float));
                    if (x_new < 0 || x_new > width_max - 1
                        || y_new < 0 || y_new > height_max - 1) {
                        continue;
                    }
                    if (img_loaded[0][fitsY][fitsX] <= 0.0001) {
                        continue;
                    }
                    
                    dummy = 1.0 + minimum_distance_borders(fitsX, fitsY, head.width, head.height);
                    if (img_temp[0][y_new][x_new] == 0.0f) {
                        for (col = 0; col < head.naxis3; ++col) {
                            img_average[col][y_new][x_new] = static_cast<float>(
                                img_loaded[col][fitsY][fitsX]
                                + background_correction_center[col]
                                + background_correction[col]);
                        }
                        img_temp[0][y_new][x_new] = static_cast<float>(dummy);
                    } else {
                        for (col = 0; col < head.naxis3; ++col) {
                            median = background_correction_center[col]
                                   + background_correction[col]
                                   + median_background(img_loaded, col, 15, 15, fitsX, fitsY);
                            if (!merge_overlap) {
                                // Method 2.
                                median2 = median_background(img_average, col, 15, 15, x_new, y_new);
                                delta_median = median - median2;
                                img_average[col][y_new][x_new] = static_cast<float>(
                                    img_average[col][y_new][x_new]
                                    + delta_median * (1.0
                                      - img_temp[0][y_new][x_new]
                                        / (dummy + img_temp[0][y_new][x_new])));
                            } else {
                                // Method 1.
                                value = img_loaded[col][fitsY][fitsX]
                                      + background_correction_center[col];
                                local_sd(fitsX - 15, fitsY - 15, fitsX + 15, fitsY + 15,
                                         col, img_loaded, noise, mean, iterations);
                                maxlevel = median + noise * 5;
                                if (value < maxlevel
                                    && img_loaded[col][fitsY][fitsX - 1] < maxlevel
                                    && img_loaded[col][fitsY][fitsX + 1] < maxlevel
                                    && img_loaded[col][fitsY - 1][fitsX] < maxlevel
                                    && img_loaded[col][fitsY + 1][fitsX] < maxlevel) {
                                    img_average[col][y_new][x_new] = static_cast<float>(
                                        +img_average[col][y_new][x_new]
                                            * img_temp[0][y_new][x_new]
                                            / (dummy + img_temp[0][y_new][x_new])
                                        + (value + background_correction[col]) * dummy
                                            / (dummy + img_temp[0][y_new][x_new]));
                                }
                            }
                        }
                        img_temp[0][y_new][x_new] = static_cast<float>(dummy);
                    }
                }
            }
            
            progress_indicator(10.0 + 89.0 * counter / images_selected, " Stacking");
        } catch (...) {
            // Silently continue on error.
        }
    }
    
    if (counter != 0) {
        head_ref.naxis3 = head.naxis3;
        head_ref.naxis  = head.naxis;
        head = head_ref;
        head.height = height_max;
        head.width  = width_max;
        img_loaded.assign(head.naxis3,
                          std::vector<std::vector<float>>(
                              head.height, std::vector<float>(head.width, 0.0f)));
                              
        for (fitsY = 0; fitsY < head.height; ++fitsY) {
            for (fitsX = 0; fitsX < head.width; ++fitsX) {
                tempval = img_temp[0][fitsY][fitsX];
                for (col = 0; col < head.naxis3; ++col) {
                    if (tempval != 0.0f) {
                        img_loaded[col][fitsY][fitsX] = img_average[col][fitsY][fitsX];
                    } else {
                        if (fitsX > 0 && img_temp[0][fitsY][fitsX - 1] != 0.0f) {
                            img_loaded[col][fitsY][fitsX] = img_loaded[col][fitsY][fitsX - 1];
                        } else if (fitsY > 0 && img_temp[0][fitsY - 1][fitsX] != 0.0f) {
                            img_loaded[col][fitsY][fitsX] = img_loaded[col][fitsY - 1][fitsX];
                        } else {
                            img_loaded[col][fitsY][fitsX] = 0.0f;
                        }
                    }
                }
            }
        }
    }
    
    // Disable SIP.
    // TODO: mainwindow.Polynomial1.itemindex = 0;  // switch to WCS
    a_order = 0;
}

///----------------------------------------
/// MARK: stack_sigmaclip
///----------------------------------------

struct SigmaClipSolution {
    double solution_vectorX[3]{};
    double solution_vectorY[3]{};
    double cblack = 0.0;
};

void stack_sigmaclip(int process_as_osc,
                     std::span<FileToDo> files_to_process,
                     int& counter) {
    auto solutions = std::vector<SigmaClipSolution>(files_to_process.size());
    
    auto fitsX = 0;
    auto fitsY = 0;
    auto c = 0;
    auto width_max = 0;
    auto height_max = 0;
    auto old_width = 0;
    auto old_height = 0;
    auto old_naxis3 = 0;
    auto x_new = 0;
    auto y_new = 0;
    auto col = 0;
    auto binning = 0;
    auto max_stars = 0;
    auto variance_factor = 0.0;
    auto value = 0.0;
    auto weightF = 0.0;
    auto hfd_min = 0.0;
    auto aa = 0.0;
    auto bb = 0.0;
    auto cc = 0.0;
    auto dd = 0.0;
    auto ee = 0.0;
    auto ff = 0.0;
    auto init = false;
    auto solution = true;
    auto use_manual_align = false;
    auto use_ephemeris_alignment = false;
    auto use_astrometry_internal = false;
    auto use_sip = false;
    auto tempval = 0.0f;
    auto target_background = 500.0f;
    auto background_correction = 0.0f;
    auto warning = std::string{};
    auto starlist1 = StarList{};
    auto starlist2 = StarList{};
    auto img_temp = ImageArray{};
    auto img_average = ImageArray{};
    auto img_final = ImageArray{};
    auto img_variance = ImageArray{};
    
    variance_factor = sqr(astap::sigma_clip_factor);
    hfd_min   = std::max(0.8, astap::hfd_min_setting);
    max_stars = astap::max_stars_setting;
    use_sip   = astap::add_sip;
    use_manual_align        = astap::use_manual_align;
    use_ephemeris_alignment = astap::use_ephemeris_alignment;
    use_astrometry_internal = astap::use_astrometry_internal;
    
    counter = 0;
    sum_exp = 0.0;
    sum_temp = 0.0;
    jd_sum = 0.0;
    jd_start_first = 1e99;
    jd_end_last = 0.0;
    
    init = false;
    background_correction = 0.0f;
    
    // --------- Pass 1: weighted average ---------
    init = false;
    for (c = 0; c < static_cast<int>(files_to_process.size()); ++c) {
        if (files_to_process[c].name.empty()) {
            continue;
        }
        try {
            listview_mark_current(files_to_process[c].listview_index);
            filename2 = files_to_process[c].name;
            
            if (esc_pressed) {
                memo2_message("ESC pressed.");
                return;
            }
            if (!load_fits(filename2, true, true, !init, 0, astap::memo1_lines, head, img_loaded)) {
                memo2_message("Error loading " + filename2);
                return;
            }
            if (!init) {
                old_width  = head.width;
                old_height = head.height;
                old_naxis3 = head.naxis3;
                head_ref = head;
                initialise_calc_sincos_dec0();
                if (bayerpat.empty() && process_as_osc == 2) {
                    // TODO: bayer_pattern1.Text check
                    if (!test_bayer_matrix(img_loaded)) {
                        memo2_message("Warning, grayscale image converted to colour!");
                    }
                }
            } else {
                if (old_width != head.width || old_height != head.height) {
                    memo2_message("Warning different size image!");
                }
                if (head.naxis3 > old_naxis3) {
                    memo2_message("Abort!! Can't combine colour to mono files.");
                    return;
                }
            }
            if (!use_sip) {
                a_order = 0;
            }
            
            (void)apply_dark_and_flat(img_loaded, head);
            
            memo2_message("Adding light file: " + std::to_string(counter + 1)
                          + " \"" + filename2 + "\" dark compensated.");
                          
            if (esc_pressed) {
                return;
            }
            
            if (process_as_osc > 0) {
                if (head.naxis3 > 1) {
                    memo2_message("Warning, light is already in colour! Will skip demosaic.");
                } else {
                    demosaic_bayer(img_loaded);
                }
            }
            
            if (use_astrometry_internal) {
                get_background(0, img_loaded, true, false, bck);
                solutions[c].cblack = bck.backgr;
            }
            
            if (!init) {
                binning = report_binning(head.height);
                if (!use_astrometry_internal) {
                    if (use_manual_align || use_ephemeris_alignment) {
                        referenceX = strtofloat2(listview_subitem(files_to_process[c].listview_index, L_X));
                        referenceY = strtofloat2(listview_subitem(files_to_process[c].listview_index, L_Y));
                    } else {
                        bin_and_find_stars(img_loaded, binning, 1, hfd_min, max_stars,
                                           true, starlist1, warning);
                        find_quads(starlist1, quad_star_distances1);
                    }
                }
                height_max = head.height;
                width_max  = head.width;
                clear_3d(img_average, head.naxis3, height_max, width_max);
                clear_3d(img_temp,    head.naxis3, height_max, width_max);
                target_background = static_cast<float>(std::max(500.0, bck.backgr));
                memo2_message("Target background for all images is "
                              + std::to_string(target_background));
            }
            
            solution = true;
            if (use_astrometry_internal) {
                sincos(head.dec0, SIN_dec0, COS_dec0);
            } else {
                if (init) {
                    if (use_manual_align || use_ephemeris_alignment) {
                        calculate_manual_vector(files_to_process, c);
                    } else {
                        bin_and_find_stars(img_loaded, binning, 1, hfd_min, max_stars,
                                           true, starlist2, warning);
                        find_quads(starlist2, quad_star_distances2);
                        if (find_offset_and_rotation(3, astap::quad_tolerance)) {
                            memo2_message(std::to_string(nr_references) + " of "
                                          + std::to_string(nr_references2) + " quads. "
                                          + solution_str);
                            std::copy(std::begin(solution_vectorX), std::end(solution_vectorX),
                                      std::begin(solutions[c].solution_vectorX));
                            std::copy(std::begin(solution_vectorY), std::end(solution_vectorY),
                                      std::begin(solutions[c].solution_vectorY));
                            solutions[c].cblack = bck.backgr;
                        } else {
                            memo2_message("Not enough quad matches <3, skipping this image.");
                            files_to_process[c].name.clear();
                            solution = false;
                            listview_set_icon(files_to_process[c].listview_index, L_result, 6);
                            listview_set_subitem(files_to_process[c].listview_index, L_result, "no solution");
                        }
                    }
                } else {
                    // First image.
                    reset_solution_vectors(1);
                    std::copy(std::begin(solution_vectorX), std::end(solution_vectorX),
                              std::begin(solutions[c].solution_vectorX));
                    std::copy(std::begin(solution_vectorY), std::end(solution_vectorY),
                              std::begin(solutions[c].solution_vectorY));
                    solutions[c].cblack = bck.backgr;
                }
            }
            init = true;
            
            if (solution) {
                ++counter;
                sum_exp  += head.exposure;
                sum_temp += head.set_temperature;
                weightF   = calc_weightF();
                background_correction = static_cast<float>(solutions[c].cblack - target_background);
                head.datamax_org = std::min(65535.0, head.datamax_org - background_correction);
                
                date_to_jd(head.date_obs, head.date_avg, head.exposure);
                jd_start_first = std::min(jd_start, jd_start_first);
                jd_end_last    = std::max(jd_end, jd_end_last);
                jd_sum        += jd_mid;
                airmass_sum   += airmass;
                
                if (use_astrometry_internal) {
                    astrometric_to_vector();
                }
                
                aa = solution_vectorX[0];
                bb = solution_vectorX[1];
                cc = solution_vectorX[2];
                dd = solution_vectorY[0];
                ee = solution_vectorY[1];
                ff = solution_vectorY[2];
                
                for (fitsY = 0; fitsY < head.height; ++fitsY) {
                    for (fitsX = 0; fitsX < head.width; ++fitsX) {
                        x_new = static_cast<int>(std::round(aa * fitsX + bb * fitsY + cc));
                        y_new = static_cast<int>(std::round(dd * fitsX + ee * fitsY + ff));
                        if (x_new < 0 || x_new > width_max - 1
                            || y_new < 0 || y_new > height_max - 1) {
                            continue;
                        }
                        for (col = 0; col < head.naxis3; ++col) {
                            img_average[col][y_new][x_new] += static_cast<float>(
                                (img_loaded[col][fitsY][fitsX] - background_correction) * weightF);
                            img_temp[col][y_new][x_new] += static_cast<float>(weightF);
                        }
                    }
                }
            }
            progress_indicator(10.0 + std::round(0.3333 * 90.0 * counter / images_selected), " []..");
        } catch (...) {
            // Silently continue on error.
        }
    }
    if (counter != 0) {
        for (fitsY = 0; fitsY < height_max; ++fitsY) {
            for (fitsX = 0; fitsX < width_max; ++fitsX) {
                for (col = 0; col < head.naxis3; ++col) {
                    if (img_temp[col][fitsY][fitsX] != 0.0f) {
                        img_average[col][fitsY][fitsX] /= img_temp[col][fitsY][fitsX];
                    }
                }
            }
        }
    }
    
    // --------- Pass 2: standard deviation ---------
    counter = 0;
    init = false;
    for (c = 0; c < static_cast<int>(files_to_process.size()); ++c) {
        if (files_to_process[c].name.empty()) {
            continue;
        }
        try {
            listview_mark_current(files_to_process[c].listview_index);
            filename2 = files_to_process[c].name;
            
            if (esc_pressed) {
                memo2_message("ESC pressed.");
                return;
            }
            if (!load_fits(filename2, true, true, !init, 0, astap::memo1_lines, head, img_loaded)) {
                memo2_message("Error loading " + filename2);
                return;
            }
            
            (void)apply_dark_and_flat(img_loaded, head);
            
            memo2_message("Calculating pixels sigma of light file "
                          + std::to_string(counter + 1) + " " + filename2);
            if (esc_pressed) {
                return;
            }
            
            if (process_as_osc > 0) {
                if (head.naxis3 > 1) {
                    memo2_message("Warning, light is already in colour! Will skip demosaic.");
                } else {
                    demosaic_bayer(img_loaded);
                }
            }
            
            if (!init) {
                clear_3d(img_variance, head.naxis3, height_max, width_max);
            }
            
            ++counter;
            
            if (use_astrometry_internal) {
                sincos(head.dec0, SIN_dec0, COS_dec0);
            } else {
                if (use_manual_align || use_ephemeris_alignment) {
                    if (!init) {
                        reset_solution_vectors(1);
                    } else {
                        calculate_manual_vector(files_to_process, c);
                    }
                } else {
                    // Reuse solution from first pass.
                    std::copy(std::begin(solutions[c].solution_vectorX),
                              std::end(solutions[c].solution_vectorX),
                              std::begin(solution_vectorX));
                    std::copy(std::begin(solutions[c].solution_vectorY),
                              std::end(solutions[c].solution_vectorY),
                              std::begin(solution_vectorY));
                    bck.backgr = solutions[c].cblack;
                }
            }
            init = true;
            
            weightF = calc_weightF();
            background_correction = static_cast<float>(solutions[c].cblack - target_background);
            head.datamax_org = std::min(65535.0, head.datamax_org - background_correction);
            
            if (use_astrometry_internal) {
                astrometric_to_vector();
            }
            
            aa = solution_vectorX[0];
            bb = solution_vectorX[1];
            cc = solution_vectorX[2];
            dd = solution_vectorY[0];
            ee = solution_vectorY[1];
            ff = solution_vectorY[2];
            
            for (fitsY = 0; fitsY < head.height; ++fitsY) {
                for (fitsX = 0; fitsX < head.width; ++fitsX) {
                    x_new = static_cast<int>(std::round(aa * fitsX + bb * fitsY + cc));
                    y_new = static_cast<int>(std::round(dd * fitsX + ee * fitsY + ff));
                    if (x_new < 0 || x_new > width_max - 1
                        || y_new < 0 || y_new > height_max - 1) {
                        continue;
                    }
                    for (col = 0; col < head.naxis3; ++col) {
                        img_variance[col][y_new][x_new] += static_cast<float>(sqr(
                            (img_loaded[col][fitsY][fitsX] - background_correction) * weightF
                            - img_average[col][y_new][x_new]));
                    }
                }
            }
            progress_indicator(10.0 + 30.0 + std::round(0.33333 * 90.0 * counter / images_selected), " [][] .");
        } catch (...) {
            // Silently continue on error.
        }
    }
    if (counter != 0) {
        for (fitsY = 0; fitsY < height_max; ++fitsY) {
            for (fitsX = 0; fitsX < width_max; ++fitsX) {
                for (col = 0; col < head.naxis3; ++col) {
                    if (img_temp[col][fitsY][fitsX] != 0.0f) {
                        img_variance[col][fitsY][fitsX] = static_cast<float>(
                            1.0f + img_variance[col][fitsY][fitsX] / img_temp[col][fitsY][fitsX]);
                    }
                }
            }
        }
    }
    
    // --------- Pass 3: throw out outliers ---------
    counter = 0;
    init = false;
    for (c = 0; c < static_cast<int>(files_to_process.size()); ++c) {
        if (files_to_process[c].name.empty()) {
            continue;
        }
        try {
            listview_mark_current(files_to_process[c].listview_index);
            filename2 = files_to_process[c].name;
            
            if (esc_pressed) {
                memo2_message("ESC pressed.");
                return;
            }
            if (!load_fits(filename2, true, true, !init, 0, astap::memo1_lines, head, img_loaded)) {
                memo2_message("Error loading " + filename2);
                return;
            }
            (void)apply_dark_and_flat(img_loaded, head);
            
            memo2_message("Combining " + std::to_string(counter + 1) + " \""
                          + filename2 + "\", ignoring outliers.");
            if (esc_pressed) {
                return;
            }
            
            if (process_as_osc > 0) {
                if (head.naxis3 > 1) {
                    memo2_message("Warning, light is already in colour! Will skip demosaic.");
                } else {
                    demosaic_bayer(img_loaded);
                }
            }
            
            if (!init) {
                clear_3d(img_temp,  head.naxis3, height_max, width_max);
                clear_3d(img_final, head.naxis3, height_max, width_max);
            }
            ++counter;
            
            if (use_astrometry_internal) {
                sincos(head.dec0, SIN_dec0, COS_dec0);
            } else {
                if (use_manual_align || use_ephemeris_alignment) {
                    if (!init) {
                        reset_solution_vectors(1);
                    } else {
                        calculate_manual_vector(files_to_process, c);
                    }
                } else {
                    std::copy(std::begin(solutions[c].solution_vectorX),
                              std::end(solutions[c].solution_vectorX),
                              std::begin(solution_vectorX));
                    std::copy(std::begin(solutions[c].solution_vectorY),
                              std::end(solutions[c].solution_vectorY),
                              std::begin(solution_vectorY));
                    bck.backgr = solutions[c].cblack;
                }
            }
            init = true;
            
            weightF = calc_weightF();
            background_correction = static_cast<float>(solutions[c].cblack - target_background);
            head.datamax_org = std::min(65535.0, head.datamax_org - background_correction);
            
            if (use_astrometry_internal) {
                astrometric_to_vector();
            }
            
            aa = solution_vectorX[0];
            bb = solution_vectorX[1];
            cc = solution_vectorX[2];
            dd = solution_vectorY[0];
            ee = solution_vectorY[1];
            ff = solution_vectorY[2];
            
            for (fitsY = 0; fitsY < head.height; ++fitsY) {
                for (fitsX = 0; fitsX < head.width; ++fitsX) {
                    x_new = static_cast<int>(std::round(aa * fitsX + bb * fitsY + cc));
                    y_new = static_cast<int>(std::round(dd * fitsX + ee * fitsY + ff));
                    if (x_new < 0 || x_new > width_max - 1
                        || y_new < 0 || y_new > height_max - 1) {
                        continue;
                    }
                    for (col = 0; col < head.naxis3; ++col) {
                        value = (img_loaded[col][fitsY][fitsX] - background_correction) * weightF;
                        if (sqr(value - img_average[col][y_new][x_new])
                            < variance_factor * img_variance[col][y_new][x_new]) {
                            img_final[col][y_new][x_new] += static_cast<float>(value);
                            img_temp[col][y_new][x_new]  += static_cast<float>(weightF);
                        }
                    }
                }
            }
            progress_indicator(10.0 + 60.0 + std::round(0.33333 * 90.0 * counter / images_selected), " [][][]");
        } catch (...) {
            // Silently continue on error.
        }
    }
    
    if (counter != 0) {
        head_ref.naxis3 = head.naxis3;
        head_ref.naxis  = head.naxis;
        head_ref.datamax_org = head.datamax_org;
        head = head_ref;
        head.height = height_max;
        head.width  = width_max;
        img_loaded.assign(head.naxis3,
                          std::vector<std::vector<float>>(
                              head.height, std::vector<float>(head.width, 0.0f)));
                              
        for (col = 0; col < head.naxis3; ++col) {
            for (fitsY = 0; fitsY < head.height; ++fitsY) {
                for (fitsX = 0; fitsX < head.width; ++fitsX) {
                    tempval = img_temp[col][fitsY][fitsX];
                    if (tempval != 0.0f) {
                        img_loaded[col][fitsY][fitsX] = img_final[col][fitsY][fitsX] / tempval;
                    } else {
                        if (fitsX > 0 && fitsY > 0) {
                            if (img_temp[col][fitsY][fitsX - 1] != 0.0f) {
                                img_loaded[col][fitsY][fitsX] = img_loaded[col][fitsY][fitsX - 1];
                            } else if (img_temp[col][fitsY - 1][fitsX] != 0.0f) {
                                img_loaded[col][fitsY][fitsX] = img_loaded[col][fitsY - 1][fitsX];
                            } else {
                                img_loaded[col][fitsY][fitsX] = 0.0f;
                            }
                        } else {
                            img_loaded[col][fitsY][fitsX] = 0.0f;
                        }
                    }
                }
            }
        }
    }
}

///----------------------------------------
/// MARK: stack_comet
///----------------------------------------

struct CometSolution {
    double solution_vectorX[3]{};
    double solution_vectorY[3]{};
    float  cblack[3]{};
};

void stack_comet(int process_as_osc,
                 std::span<FileToDo> files_to_process,
                 int& counter) {
    auto solutions = std::vector<CometSolution>(files_to_process.size());
    
    auto fitsX = 0;
    auto fitsY = 0;
    auto c = 0;
    auto width_max = 0;
    auto height_max = 0;
    auto old_width = 0;
    auto old_height = 0;
    auto old_naxis3 = 0;
    auto x_new = 0;
    auto y_new = 0;
    auto col = 0;
    auto value = 0.0;
    auto weightF = 0.0;
    auto hfd_min = 0.0;
    auto aa = 0.0;
    auto bb = 0.0;
    auto cc = 0.0;
    auto dd = 0.0;
    auto ee = 0.0;
    auto ff = 0.0;
    auto delta_JD_required = 0.0;
    auto target_background = 500.0;
    auto JD_reference = 0.0;
    auto init = false;
    auto solution = true;
    auto use_manual_align = false;
    auto use_ephemeris_alignment = false;
    auto use_astrometry_internal = false;
    auto use_sip = false;
    auto tempval = 0.0f;
    auto jd_fraction = 0.0f;
    float background_correction[3] = {0.0f, 0.0f, 0.0f};
    auto img_temp = ImageArray{};
    auto img_final = ImageArray{};
    auto img_variance = ImageArray{};
    
    hfd_min = std::max(0.8, astap::hfd_min_setting);
    use_sip = astap::add_sip;
    use_manual_align = astap::use_manual_align;
    use_ephemeris_alignment = astap::use_ephemeris_alignment;
    use_astrometry_internal = astap::use_astrometry_internal;
    
    counter = 0;
    sum_exp = 0.0;
    sum_temp = 0.0;
    jd_sum = 0.0;
    jd_start_first = 1e99;
    jd_end_last = 0.0;
    init = false;
    
    // --------- Pass 1: find JD when each pixel is at max value ---------
    init = false;
    for (c = 0; c < static_cast<int>(files_to_process.size()); ++c) {
        if (files_to_process[c].name.empty()) {
            continue;
        }
        try {
            listview_mark_current(files_to_process[c].listview_index);
            filename2 = files_to_process[c].name;
            
            if (esc_pressed) {
                memo2_message("ESC pressed.");
                return;
            }
            if (!load_fits(filename2, true, true, !init, 0, astap::memo1_lines, head, img_loaded)) {
                memo2_message("Error loading " + filename2);
                return;
            }
            
            if (!init) {
                old_width  = head.width;
                old_height = head.height;
                old_naxis3 = head.naxis3;
                head_ref = head;
                initialise_calc_sincos_dec0();
                if (bayerpat.empty() && process_as_osc == 2) {
                    if (!test_bayer_matrix(img_loaded)) {
                        memo2_message("Warning, grayscale image converted to colour!");
                    }
                }
            } else {
                if (old_width != head.width || old_height != head.height) {
                    memo2_message("Warning different size image!");
                }
                if (head.naxis3 > old_naxis3) {
                    memo2_message("Abort!! Can't combine colour to mono files.");
                    return;
                }
            }
            if (!use_sip) {
                a_order = 0;
            }
            
            (void)apply_dark_and_flat(img_loaded, head);
            memo2_message("Registrating drifting stars movements: "
                          + std::to_string(counter + 1) + " \"" + filename2 + "\"");
            if (esc_pressed) {
                return;
            }
            
            if (process_as_osc > 0) {
                if (head.naxis3 > 1) {
                    memo2_message("Warning, light is already in colour! Will skip demosaic.");
                } else {
                    demosaic_bayer(img_loaded);
                }
            }
            
            // Measure background for each colour.
            memo2_message("Measuring background for all colours");
            for (col = 0; col < head.naxis3; ++col) {
                get_background(col, img_loaded, true, false, bck);
                solutions[c].cblack[col] = static_cast<float>(bck.backgr);
            }
            
            if (!init) {
                referenceX = strtofloat2(listview_subitem(files_to_process[c].listview_index, L_X));
                referenceY = strtofloat2(listview_subitem(files_to_process[c].listview_index, L_Y));
                
                height_max = head.height;
                width_max  = head.width;
                clear_3d(img_variance, 2, height_max, width_max);
                target_background = std::max(500.0f, solutions[c].cblack[0]);
                memo2_message("Target background for all images is "
                              + std::to_string(target_background));
            }
            
            solution = true;
            if (init) {
                calculate_manual_vector(files_to_process, c);
            } else {
                reset_solution_vectors(1);
                std::copy(std::begin(solution_vectorX), std::end(solution_vectorX),
                          std::begin(solutions[c].solution_vectorX));
                std::copy(std::begin(solution_vectorY), std::end(solution_vectorY),
                          std::begin(solutions[c].solution_vectorY));
            }
            init = true;
            
            if (solution) {
                ++counter;
                sum_exp  += head.exposure;
                sum_temp += head.set_temperature;
                weightF   = calc_weightF();
                for (col = 0; col < head.naxis3; ++col) {
                    background_correction[col] = static_cast<float>(solutions[c].cblack[col] - target_background);
                }
                head.datamax_org = std::min(65535.0, head.datamax_org - background_correction[0]);
                
                date_to_jd(head.date_obs, head.date_avg, head.exposure);
                jd_start_first = std::min(jd_start, jd_start_first);
                jd_end_last    = std::max(jd_end, jd_end_last);
                jd_sum        += jd_mid;
                airmass_sum   += airmass;
                
                jd_fraction = static_cast<float>(jd_mid - std::floor(jd_mid));
                
                if (counter == 1) {
                    JD_reference = jd_start;
                } else if (counter == 2) {
                    // Drift vs reference image.
                    auto hfd_row = strtofloat2(listview_subitem(files_to_process[c].listview_index, L_hfd));
                    delta_JD_required = std::abs(jd_start - JD_reference)
                        * 3.0 * hfd_row
                        / std::sqrt(sqr(solution_vectorX[2]) + sqr(solution_vectorY[2]));
                    memo2_message("For stars 3*HFD drift takes "
                                  + std::to_string(delta_JD_required * 24.0 * 3600.0) + "sec");
                }
                
                aa = solution_vectorX[0];
                bb = solution_vectorX[1];
                cc = solution_vectorX[2];
                dd = solution_vectorY[0];
                ee = solution_vectorY[1];
                ff = solution_vectorY[2];
                
                for (fitsY = 0; fitsY < head.height; ++fitsY) {
                    for (fitsX = 0; fitsX < head.width; ++fitsX) {
                        x_new = static_cast<int>(std::round(aa * fitsX + bb * fitsY + cc));
                        y_new = static_cast<int>(std::round(dd * fitsX + ee * fitsY + ff));
                        if (x_new < 0 || x_new > width_max - 1
                            || y_new < 0 || y_new > height_max - 1) {
                            continue;
                        }
                        
                        value = 0.0;
                        for (col = 0; col < head.naxis3; ++col) {
                            value += (img_loaded[col][fitsY][fitsX] - background_correction[col]) * weightF;
                        }
                        if (value > img_variance[0][y_new][x_new]) {
                            img_variance[0][y_new][x_new] = static_cast<float>(value);
                            img_variance[1][y_new][x_new] = jd_fraction;
                        }
                    }
                }
            }
            progress_indicator(10.0 + std::round(0.5 * 90.0 * counter / images_selected), " []..");
        } catch (...) {
            // Silently continue on error.
        }
    }
    
    // --------- Pass 2: combine, reject frames where star is passing ---------
    counter = 0;
    init = false;
    for (c = 0; c < static_cast<int>(files_to_process.size()); ++c) {
        if (files_to_process[c].name.empty()) {
            continue;
        }
        try {
            listview_mark_current(files_to_process[c].listview_index);
            filename2 = files_to_process[c].name;
            
            if (esc_pressed) {
                memo2_message("ESC pressed.");
                return;
            }
            if (!load_fits(filename2, true, true, !init, 0, astap::memo1_lines, head, img_loaded)) {
                memo2_message("Error loading " + filename2);
                return;
            }
            (void)apply_dark_and_flat(img_loaded, head);
            
            date_to_jd(head.date_obs, head.date_avg, head.exposure);
            jd_fraction = static_cast<float>(jd_mid - std::floor(jd_mid));
            
            memo2_message("Combining " + std::to_string(counter + 1) + " \""
                          + filename2 + "\", ignoring moving stars.");
            if (esc_pressed) {
                return;
            }
            
            if (process_as_osc > 0) {
                if (head.naxis3 > 1) {
                    memo2_message("Warning, light is already in colour! Will skip demosaic.");
                } else {
                    demosaic_bayer(img_loaded);
                }
            }
            
            if (!init) {
                clear_3d(img_temp,  1,           height_max, width_max);
                clear_3d(img_final, head.naxis3, height_max, width_max);
            }
            ++counter;
            
            if (use_astrometry_internal) {
                sincos(head.dec0, SIN_dec0, COS_dec0);
            } else {
                if (use_manual_align || use_ephemeris_alignment) {
                    if (!init) {
                        reset_solution_vectors(1);
                    } else {
                        calculate_manual_vector(files_to_process, c);
                    }
                } else {
                    std::copy(std::begin(solutions[c].solution_vectorX),
                              std::end(solutions[c].solution_vectorX),
                              std::begin(solution_vectorX));
                    std::copy(std::begin(solutions[c].solution_vectorY),
                              std::end(solutions[c].solution_vectorY),
                              std::begin(solution_vectorY));
                }
            }
            
            weightF = calc_weightF();
            for (col = 0; col < head.naxis3; ++col) {
                background_correction[col] = static_cast<float>(solutions[c].cblack[col] - target_background);
            }
            head.datamax_org = std::min(65535.0, head.datamax_org - background_correction[0]);
            
            aa = solution_vectorX[0];
            bb = solution_vectorX[1];
            cc = solution_vectorX[2];
            dd = solution_vectorY[0];
            ee = solution_vectorY[1];
            ff = solution_vectorY[2];
            
            for (fitsY = 0; fitsY < head.height; ++fitsY) {
                for (fitsX = 0; fitsX < head.width; ++fitsX) {
                    x_new = static_cast<int>(std::round(aa * fitsX + bb * fitsY + cc));
                    y_new = static_cast<int>(std::round(dd * fitsX + ee * fitsY + ff));
                    if (x_new < 0 || x_new > width_max - 1
                        || y_new < 0 || y_new > height_max - 1) {
                        continue;
                    }
                    
                    auto keep = !init
                        || std::abs(jd_fraction - img_variance[1][y_new][x_new]) > delta_JD_required;
                    if (!keep) {
                        continue;
                    }
                    
                    for (col = 0; col < head.naxis3; ++col) {
                        value = (img_loaded[col][fitsY][fitsX] - background_correction[col]) * weightF;
                        img_final[col][y_new][x_new] += static_cast<float>(value);
                        img_temp[0][y_new][x_new]    += static_cast<float>(weightF);
                    }
                }
            }
            
            init = true;
            progress_indicator(10.0 + 45.0 + std::round(0.5 * 90.0 * counter / images_selected), " [][] ");
        } catch (...) {
            // Silently continue on error.
        }
    }
    
    if (counter != 0) {
        head_ref.naxis3 = head.naxis3;
        head_ref.naxis  = head.naxis;
        head_ref.datamax_org = head.datamax_org;
        head = head_ref;
        head.height = height_max;
        head.width  = width_max;
        img_loaded.assign(head.naxis3,
                          std::vector<std::vector<float>>(
                              head.height, std::vector<float>(head.width, 0.0f)));
                              
        for (col = 0; col < head.naxis3; ++col) {
            for (fitsY = 0; fitsY < head.height; ++fitsY) {
                for (fitsX = 0; fitsX < head.width; ++fitsX) {
                    tempval = img_temp[0][fitsY][fitsX];
                    if (tempval != 0.0f) {
                        img_loaded[col][fitsY][fitsX] = img_final[col][fitsY][fitsX] / tempval;
                    } else {
                        if (fitsX > 0 && fitsY > 0) {
                            if (img_temp[0][fitsY][fitsX - 1] != 0.0f) {
                                img_loaded[col][fitsY][fitsX] = img_loaded[col][fitsY][fitsX - 1];
                            } else if (img_temp[0][fitsY - 1][fitsX] != 0.0f) {
                                img_loaded[col][fitsY][fitsX] = img_loaded[col][fitsY - 1][fitsX];
                            } else {
                                img_loaded[col][fitsY][fitsX] = 0.0f;
                            }
                        } else {
                            img_loaded[col][fitsY][fitsX] = 0.0f;
                        }
                    }
                }
            }
        }
    }
}

///----------------------------------------
/// MARK: calibration_and_alignment
///----------------------------------------

void calibration_and_alignment(int process_as_osc,
                               std::span<FileToDo> files_to_process,
                               int& counter) {
    auto fitsX = 0;
    auto fitsY = 0;
    auto c = 0;
    auto width_max = 0;
    auto height_max = 0;
    auto old_width = 0;
    auto old_height = 0;
    auto old_naxis3 = 0;
    auto x_new = 0;
    auto y_new = 0;
    auto col = 0;
    auto binning = 0;
    auto max_stars = 0;
    auto background_correction = 0.0;
    auto hfd_min = 0.0;
    auto aa = 0.0;
    auto bb = 0.0;
    auto cc = 0.0;
    auto dd = 0.0;
    auto ee = 0.0;
    auto ff = 0.0;
    auto init = false;
    auto solution = true;
    auto use_manual_align = false;
    auto use_ephemeris_alignment = false;
    auto use_astrometry_internal = false;
    auto use_sip = false;
    auto warning = std::string{};
    auto starlist1 = StarList{};
    auto starlist2 = StarList{};
    auto img_temp = ImageArray{};
    auto img_average = ImageArray{};
    
    hfd_min   = std::max(0.8, astap::hfd_min_setting);
    max_stars = astap::max_stars_setting;
    use_sip   = astap::add_sip;
    use_manual_align        = astap::use_manual_align;
    use_ephemeris_alignment = astap::use_ephemeris_alignment;
    use_astrometry_internal = astap::use_astrometry_internal;
    
    counter = 0;
    sum_exp = 0.0;
    
    for (c = 0; c < static_cast<int>(files_to_process.size()); ++c) {
        if (files_to_process[c].name.empty()) {
            continue;
        }
        try {
            listview_mark_current(files_to_process[c].listview_index);
            filename2 = files_to_process[c].name;
            
            if (esc_pressed) {
                memo2_message("ESC pressed.");
                return;
            }
            if (!load_fits(filename2, true, true, true, 0, astap::memo1_lines, head, img_loaded)) {
                memo2_message("Error loading " + filename2);
                return;
            }
            if (!init) {
                old_width  = head.width;
                old_height = head.height;
                old_naxis3 = head.naxis3;
                head_ref = head;
                initialise_calc_sincos_dec0();
                if (bayerpat.empty() && process_as_osc == 2) {
                    if (!test_bayer_matrix(img_loaded)) {
                        memo2_message("Warning, grayscale image converted to colour!");
                    }
                }
            } else {
                if (old_width != head.width || old_height != head.height) {
                    memo2_message("Warning different size image!");
                }
                if (head.naxis3 > old_naxis3) {
                    memo2_message("Abort!! Can't combine colour to mono files.");
                    return;
                }
            }
            
            if (!use_sip) {
                a_order = 0;
            }
            (void)apply_dark_and_flat(img_loaded, head);
            
            memo2_message("Calibrating and aligning file: "
                          + std::to_string(counter + 1) + " \"" + filename2 + "\"");
            if (esc_pressed) {
                return;
            }
            
            if (process_as_osc > 0) {
                if (head.naxis3 > 1) {
                    memo2_message("Warning, light is already in colour! Will skip demosaic.");
                } else {
                    demosaic_bayer(img_loaded);
                }
            } else if (!bayerpat.empty()) {
                memo2_message("Warning, alignment will ruin Bayer pattern!");
            }
            
            if (!init) {
                binning = report_binning(head.height);
            }
            if (!init && !use_astrometry_internal) {
                if (use_manual_align || use_ephemeris_alignment) {
                    referenceX = strtofloat2(listview_subitem(files_to_process[c].listview_index, L_X));
                    referenceY = strtofloat2(listview_subitem(files_to_process[c].listview_index, L_Y));
                } else {
                    bin_and_find_stars(img_loaded, binning, 1, hfd_min, max_stars,
                                       true, starlist1, warning);
                    find_quads(starlist1, quad_star_distances1);
                    pedestal_s = bck.backgr;
                    if (pedestal_s < 500.0) {
                        pedestal_s = 500.0;
                    }
                    background_correction = pedestal_s - bck.backgr;
                    head.datamax_org += background_correction;
                    if (head.datamax_org > 0xFFFF) {
                        head.datamax_org = 0xFFFF;
                    }
                    head.pedestal = background_correction;
                }
            }
            
            if (!init) {
                height_max = head.height;
                width_max  = head.width;
                clear_3d(img_average, head.naxis3, height_max, width_max);
                clear_3d(img_temp,    head.naxis3, height_max, width_max);
                if (use_manual_align || use_ephemeris_alignment) {
                    referenceX = strtofloat2(listview_subitem(files_to_process[c].listview_index, L_X));
                    referenceY = strtofloat2(listview_subitem(files_to_process[c].listview_index, L_Y));
                }
            }
            
            // Clear image_average / img_temp per image.
            for (fitsY = 0; fitsY < height_max; ++fitsY) {
                for (fitsX = 0; fitsX < width_max; ++fitsX) {
                    for (col = 0; col < head.naxis3; ++col) {
                        img_average[col][fitsY][fitsX] = 0.0f;
                        img_temp[col][fitsY][fitsX]    = 0.0f;
                    }
                }
            }
            
            solution = true;
            if (use_astrometry_internal) {
                sincos(head.dec0, SIN_dec0, COS_dec0);
            } else {
                if (init) {
                    if (use_manual_align || use_ephemeris_alignment) {
                        calculate_manual_vector(files_to_process, c);
                    } else {
                        bin_and_find_stars(img_loaded, binning, 1, hfd_min, max_stars,
                                           true, starlist2, warning);
                        background_correction = pedestal_s - bck.backgr;
                        head.datamax_org += background_correction;
                        if (head.datamax_org > 0xFFFF) {
                            head.datamax_org = 0xFFFF;
                        }
                        head.pedestal = background_correction;
                        find_quads(starlist2, quad_star_distances2);
                        if (find_offset_and_rotation(3, astap::quad_tolerance)) {
                            memo2_message(std::to_string(nr_references) + " of "
                                          + std::to_string(nr_references2) + " quads. "
                                          + solution_str);
                        } else {
                            memo2_message("Not enough quad matches, skipping.");
                            files_to_process[c].name.clear();
                            solution = false;
                            listview_set_icon(files_to_process[c].listview_index, L_result, 6);
                            listview_set_subitem(files_to_process[c].listview_index, 2, "no solution");
                        }
                    }
                } else {
                    reset_solution_vectors(1);
                }
            }
            init = true;
            
            if (solution) {
                ++counter;
                if (use_astrometry_internal) {
                    astrometric_to_vector();
                }
                
                aa = solution_vectorX[0];
                bb = solution_vectorX[1];
                cc = solution_vectorX[2];
                dd = solution_vectorY[0];
                ee = solution_vectorY[1];
                ff = solution_vectorY[2];
                
                for (fitsY = 0; fitsY < head.height; ++fitsY) {
                    for (fitsX = 0; fitsX < head.width; ++fitsX) {
                        x_new = static_cast<int>(std::round(aa * fitsX + bb * fitsY + cc));
                        y_new = static_cast<int>(std::round(dd * fitsX + ee * fitsY + ff));
                        if (x_new < 0 || x_new > width_max - 1
                            || y_new < 0 || y_new > height_max - 1) {
                            continue;
                        }
                        for (col = 0; col < head.naxis3; ++col) {
                            img_average[col][y_new][x_new] += static_cast<float>(
                                img_loaded[col][fitsY][fitsX] + background_correction);
                            img_temp[col][y_new][x_new] += 1.0f;
                        }
                    }
                }
            }
            
            // Scale to number of pixels.
            head.height = height_max;
            head.width  = width_max;
            img_loaded.assign(head.naxis3,
                              std::vector<std::vector<float>>(
                                  head.height, std::vector<float>(head.width, 0.0f)));
                                  
            for (col = 0; col < head.naxis3; ++col) {
                for (fitsY = 0; fitsY < head.height; ++fitsY) {
                    for (fitsX = 0; fitsX < head.width; ++fitsX) {
                        if (img_temp[col][fitsY][fitsX] != 0.0f) {
                            img_loaded[col][fitsY][fitsX] =
                                img_average[col][fitsY][fitsX] / img_temp[col][fitsY][fitsX];
                        } else {
                            if (fitsX > 0 && fitsY > 0) {
                                if (img_temp[col][fitsY][fitsX - 1] != 0.0f) {
                                    img_loaded[col][fitsY][fitsX] = img_loaded[col][fitsY][fitsX - 1];
                                } else if (img_temp[col][fitsY - 1][fitsX] != 0.0f) {
                                    img_loaded[col][fitsY][fitsX] = img_loaded[col][fitsY - 1][fitsX];
                                } else {
                                    img_loaded[col][fitsY][fitsX] = 0.0f;
                                }
                            } else {
                                img_loaded[col][fitsY][fitsX] = 0.0f;
                            }
                        }
                    }
                }
            }
            
            // Save _aligned.fit. TODO: path manipulation.
            auto dot = filename2.find_last_of('.');
            auto aligned = (dot == std::string::npos)
                ? filename2 + "_aligned.fit"
                : filename2.substr(0, dot) + "_aligned.fit";
            filename2 = aligned;
            
            // TODO: FITS header edits (CRPIX1/2 correction, COMMENTs).
            if (head.cd1_1 != 0.0) {
                head.crpix1 = solution_vectorX[0] * (head.crpix1 - 1)
                            + solution_vectorX[1] * (head.crpix2 - 1)
                            + solution_vectorX[2];
                head.crpix2 = solution_vectorY[0] * (head.crpix1 - 1)
                            + solution_vectorY[1] * (head.crpix2 - 1)
                            + solution_vectorY[2];
            }
            
            // TODO: update header comment / pedestal / counters.
            
            if (nrbits == 16) {
                if (!save_fits(img_loaded, astap::memo1_lines, filename2, 16, true)) {
                    return;
                }
            } else {
                if (!save_fits(img_loaded, astap::memo1_lines, filename2, -32, true)) {
                    return;
                }
            }
            memo2_message("New aligned image created: " + filename2);
            // TODO: report_results and plot_fits.
            progress_indicator(10.0 + std::round(90.0 * counter / images_selected), "Cal");
        } catch (...) {
            // Silently continue on error.
        }
    }
    
    // TODO: plot_fits(mainwindow.image1, true, true) -> UI refresh.
}

///----------------------------------------
/// MARK: Stub definitions
///----------------------------------------

// Thin UI-layer adapters. The engine computation lives elsewhere; these
// wrappers resolve global state and forward. UI sinks (memo / progress /
// listview) remain no-ops until the GUI installs real callbacks — see
// set_memo2_sink() etc.

void demosaic_bayer(ImageArray& img) {
    // The legacy Pascal code pulled pattern + method from stackmenu1 combos.
    // Phase 5b defaults to the RGGB hint with bilinear interpolation; pattern
    // parity is resolved from the live globals (BAYERPAT, XBAYROFF, YBAYROFF,
    // ROWORDER) so the result matches the CLI.
    const auto pattern = astap::core::get_demosaic_pattern(
        2, astap::xbayroff, astap::ybayroff, astap::roworder);
    astap::core::demosaic_bayer(img, head, pattern,
                                astap::core::DemosaicMethod::Bilinear);
}

namespace {
MemoSink g_memo_sink;
ProgressSink g_progress_sink;
}  // namespace

void set_memo2_sink(MemoSink sink) {
    g_memo_sink = std::move(sink);
}

void set_progress_sink(ProgressSink sink) {
    g_progress_sink = std::move(sink);
}

void memo2_message(const std::string& msg) {
    if (g_memo_sink) {
        g_memo_sink(msg);
    }
}

void progress_indicator(double value, const std::string& label) {
    if (g_progress_sink) {
        g_progress_sink(value, label);
    }
}

std::string listview_subitem(int /*index*/, int /*col*/) { return {}; }
void        listview_set_subitem(int, int, const std::string&) {}
void        listview_set_icon(int, int, int) {}
void        listview_mark_current(int) {}

void sincos(double a, double& s, double& c) {
    s = std::sin(a);
    c = std::cos(a);
}
  
} // namespace
