// stb_image-backed IImageDecoder implementation.
//
// This translation unit carries the STB_IMAGE_IMPLEMENTATION definition so
// the library's single header compiles into object code exactly once. The
// header must be on the include path — download from:
//
//   https://raw.githubusercontent.com/nothings/stb/master/stb_image.h
//
// and drop into a `third_party/` directory (or anywhere the compiler can
// find it). If stb_image.h is NOT available, the adapter compiles but
// install_stb_image_decoder() returns a decoder that fails every request
// with a helpful error.

#include "stb_image_decoder.h"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Conditional stb_image.h inclusion. If the header is missing, the adapter
// compiles into a "decoder not available" stub so the build doesn't break.
// -----------------------------------------------------------------------------
#if __has_include(<stb_image.h>)
    #define STB_IMAGE_IMPLEMENTATION
    #define STB_IMAGE_STATIC
    // We handle 16-bit PNGs explicitly; tell stb to give us 16-bit samples
    // when the source file carries them.
    #include <stb_image.h>
    #define ASTAP_HAVE_STB_IMAGE 1
#else
    #define ASTAP_HAVE_STB_IMAGE 0
#endif

namespace astap::cli {

namespace {

#if ASTAP_HAVE_STB_IMAGE

// Decode via stbi_load_16 first (preserves 16-bit PNGs etc.); fall back to
// stbi_load for 8-bit formats. Both paths produce an interleaved H*W*C
// buffer that we split into astap::ImageArray's [channel][row][col] layout.
bool decode_stb(const std::filesystem::path&              path,
                astap::core::IImageDecoder::DecodedImage& out,
                std::string&                              error_out) {
    const std::string path_str = path.string();

    int w = 0, h = 0, ch = 0;

    // Always request the native channel count from stb (last arg = 0).
    // stbi_is_16_bit tells us whether the file carries >8-bit samples.
    const bool is_16 = stbi_is_16_bit(path_str.c_str()) != 0;

    std::vector<std::vector<std::vector<float>>> channels;

    if (is_16) {
        stbi_us* data = stbi_load_16(path_str.c_str(), &w, &h, &ch, 0);
        if (!data) {
            error_out = stbi_failure_reason() ? stbi_failure_reason() : "stbi_load_16 failed";
            return false;
        }
        out.bits_per_sample = 16;
        channels.assign(ch, std::vector<std::vector<float>>(
            h, std::vector<float>(w, 0.0f)));
        // stb returns floats via division-by-65535 — we keep raw counts.
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                const stbi_us* pixel = data + (y * w + x) * ch;
                for (int c = 0; c < ch; ++c) {
                    channels[c][y][x] = static_cast<float>(pixel[c]);
                }
            }
        }
        stbi_image_free(data);
    } else {
        stbi_uc* data = stbi_load(path_str.c_str(), &w, &h, &ch, 0);
        if (!data) {
            error_out = stbi_failure_reason() ? stbi_failure_reason() : "stbi_load failed";
            return false;
        }
        out.bits_per_sample = 8;
        channels.assign(ch, std::vector<std::vector<float>>(
            h, std::vector<float>(w, 0.0f)));
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                const stbi_uc* pixel = data + (y * w + x) * ch;
                for (int c = 0; c < ch; ++c) {
                    channels[c][y][x] = static_cast<float>(pixel[c]);
                }
            }
        }
        stbi_image_free(data);
    }

    out.width     = w;
    out.height    = h;
    out.channels  = ch;
    out.pixels    = std::move(channels);
    return true;
}

struct StbDecoder : astap::core::IImageDecoder {
    bool decode_raster(const std::filesystem::path& path,
                       std::string_view             ext_upper,
                       DecodedImage&                out,
                       std::string&                 error_out) override {
        // stb_image doesn't do TIFF; the PC TIFFs from cameras / planetary
        // captures won't load here. Delegate errors up the stack with a
        // clear message so users know to plug in a libtiff adapter.
        if (ext_upper == ".TIF" || ext_upper == ".TIFF") {
            error_out = "stb_image does not support TIFF; install a libtiff "
                        "or libtiff-cxx adapter in addition to this one";
            return false;
        }
        return decode_stb(path, out, error_out);
    }

    bool decode_raw(const std::filesystem::path& /*path*/,
                    bool                         /*save_intermediate*/,
                    DecodedImage&                /*out*/,
                    std::string&                 error_out) override {
        error_out = "stb_image does not support camera RAW; use LibRaw or "
                    "the convert_raw() shell-out path";
        return false;
    }
};

#else  // !ASTAP_HAVE_STB_IMAGE

struct StbDecoder : astap::core::IImageDecoder {
    bool decode_raster(const std::filesystem::path&, std::string_view,
                       DecodedImage&, std::string& e) override {
        e = "stb_image.h not found at build time; add it to the include path "
            "and rebuild cli/stb_image_decoder.cpp";
        return false;
    }
    bool decode_raw(const std::filesystem::path&, bool,
                    DecodedImage&, std::string& e) override {
        e = "stb_image.h not found at build time";
        return false;
    }
};

#endif

}  // anonymous namespace

std::shared_ptr<astap::core::IImageDecoder> install_stb_image_decoder() {
    auto dec = std::make_shared<StbDecoder>();
    astap::core::set_image_decoder(dec);
    return dec;
}

}  // namespace astap::cli
