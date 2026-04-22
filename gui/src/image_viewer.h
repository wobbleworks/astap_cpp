///----------------------------------------
///      @file image_viewer.h
///   @ingroup ASTAP++
///     @brief Central image canvas for the ASTAP++ Qt GUI.
///   @details Owns the loaded @c astap::ImageArray, caches a stretched
///            @c QImage, and renders it pannable / zoomable.
///
///            Phase 1: fit-to-window, mouse-wheel zoom centred on cursor,
///            click-drag pan, 1:1 actual-size, default linear stretch from
///            header @c datamin_org / @c datamax_org.
///
///            Phase 2: explicit stretch range (black/white points), optional
///            logarithmic stretch, horizontal/vertical flip, exposed 256-bin
///            histogram per channel.
///
///            Phase 3+ adds annotation overlays.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include "../../src/types.h"
#include "annotation_scanner.h"
#include "star_detector.h"

#include <QImage>
#include <QPoint>
#include <QPointF>
#include <QTimer>
#include <QWidget>

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @class ImageViewer
/// @brief Renders an astap::ImageArray onto a pannable, zoomable canvas.
///----------------------------------------

class ImageViewer final : public QWidget {
	Q_OBJECT

public:
	using HistogramBins = std::array<std::uint32_t, 256>;

	///----------------------------------------
	/// @struct ChannelStats
	/// @brief Per-channel summary statistics computed at image load.
	///----------------------------------------
	struct ChannelStats {
		float min = 0.0f;
		float max = 0.0f;
		float mean = 0.0f;
		float median = 0.0f;
		float stddev = 0.0f;      // global σ (all pixels, useful for stats)
		float noise = 0.0f;       // sigma-clipped sky noise (from engine)
		float background = 0.0f;  // histogram mode (from engine)
		float starLevel = 0.0f;   // detection threshold above background
		std::uint64_t count = 0;
	};

	explicit ImageViewer(QWidget* parent = nullptr);
	~ImageViewer() override = default;

	/// @brief Replace the displayed image. Resets zoom to fit-to-window and
	///        recomputes histogram + auto-stretch defaults.
	void setImage(astap::ImageArray image, astap::Header header);

	/// @brief Drop the displayed image and repaint empty.
	void clear();

	[[nodiscard]] bool hasImage() const noexcept { return !_image.empty(); }
	[[nodiscard]] int imageWidth() const noexcept { return _header.width; }
	[[nodiscard]] int imageHeight() const noexcept { return _header.height; }
	[[nodiscard]] int channels() const noexcept { return static_cast<int>(_image.size()); }

	/// @name Stretch / display state
	///@{

	/// @brief Sample value mapped to display black.
	[[nodiscard]] float stretchLo() const noexcept { return _stretchLo; }

	/// @brief Sample value mapped to display white.
	[[nodiscard]] float stretchHi() const noexcept { return _stretchHi; }

	/// @brief Update both stretch endpoints. Re-renders.
	void setStretchRange(float lo, float hi);

	/// @brief Whether logarithmic stretch is applied between lo/hi.
	[[nodiscard]] bool logarithmicStretch() const noexcept { return _logStretch; }

	/// @brief Switch between linear and logarithmic stretch. Re-renders.
	void setLogarithmicStretch(bool on);

	[[nodiscard]] bool flipHorizontal() const noexcept { return _flipH; }
	[[nodiscard]] bool flipVertical() const noexcept { return _flipV; }

	void setFlipHorizontal(bool on);
	void setFlipVertical(bool on);

	/// @brief 256-bin histogram per channel. Empty when no image is loaded.
	[[nodiscard]] const std::vector<HistogramBins>& histogram() const noexcept { return _histogram; }

	/// @brief Per-channel summary statistics. Empty when no image is loaded.
	[[nodiscard]] const std::vector<ChannelStats>& stats() const noexcept { return _stats; }

	/// @brief Initial / "auto" stretch endpoints computed at load time.
	[[nodiscard]] std::pair<float, float> dataRange() const noexcept { return {_dataLo, _dataHi}; }

	/// @brief Lowest sample seen in the image.
	[[nodiscard]] float sampleMin() const noexcept { return _sampleMin; }

	/// @brief Highest sample seen in the image.
	[[nodiscard]] float sampleMax() const noexcept { return _sampleMax; }

	/// @brief Re-apply the auto-stretch heuristic and emit @c stretchChanged.
	void autoStretch();

	/// @brief Current stretched / flipped 8-bit render, suitable for JPG/PNG-8 export.
	/// @details Empty when no image is loaded. Mirrors exactly what the viewer paints.
	[[nodiscard]] const QImage& renderedImage() const noexcept { return _rendered; }
	///@}

	/// @name Zoom controls
	///@{
	void zoomIn();
	void zoomOut();
	void fitToWindow();
	void actualSize();
	///@}

	/// @brief Convert a viewport-local point to image-pixel coordinates
	///        (FITS convention: 1-based, origin bottom-left, y goes up).
	/// @param viewportPoint Cursor position relative to the viewport.
	/// @param[out] inImage True if the cursor is over the image.
	/// @return Image pixel position; meaningful only when @p inImage is true.
	[[nodiscard]] QPointF viewportToImage(QPointF viewportPoint, bool& inImage) const;

	/// @name Star overlay
	///@{

	/// @brief Replace the set of overlaid stars. Triggers a repaint.
	void setStarMarkers(std::vector<DetectedStar> stars);

	/// @brief Remove any overlaid stars and repaint.
	void clearStarMarkers();

	[[nodiscard]] const std::vector<DetectedStar>& starMarkers() const noexcept { return _stars; }
	///@}

	/// @name Annotation overlay
	///@{
	void setAnnotations(std::vector<AnnotationMarker> annotations);
	void clearAnnotations();
	[[nodiscard]] const std::vector<AnnotationMarker>& annotations() const noexcept { return _annotations; }

	void setConstellations(ConstellationOverlay overlay);
	void clearConstellations();
	///@}

	/// @name Saturation
	///@{
	[[nodiscard]] float saturation() const noexcept { return _saturation; }
	void setSaturation(float sat);
	///@}

signals:
	/// @brief Emitted after @ref setImage has updated histograms / defaults.
	void imageLoaded();

	/// @brief Emitted whenever stretch range or stretch mode changes.
	void stretchChanged();

	/// @brief Emitted whenever flip H or flip V changes.
	void flipChanged();

	/// @brief Emitted whenever the cursor moves over the image (or off it).
	/// @param imagePos FITS pixel position (valid when @p inImage is true).
	/// @param inImage True if the cursor is over the image, false off-image.
	void cursorMoved(QPointF imagePos, bool inImage);

protected:
	void paintEvent(QPaintEvent* event) override;
	void wheelEvent(QWheelEvent* event) override;
	void mousePressEvent(QMouseEvent* event) override;
	void mouseMoveEvent(QMouseEvent* event) override;
	void mouseReleaseEvent(QMouseEvent* event) override;
	void resizeEvent(QResizeEvent* event) override;

private:
	void renderImage();
	void scheduleRender();
	void computeStatistics();
	void clampPan();
	void applyZoomAtViewport(double newZoom, QPointF viewportPoint);

	QTimer _renderTimer;

	astap::ImageArray _image;
	astap::Header _header{};
	QImage _rendered;

	// Histogram + per-image statistics, computed once on setImage.
	std::vector<HistogramBins> _histogram;
	std::vector<ChannelStats> _stats;
	float _sampleMin = 0.0f;
	float _sampleMax = 0.0f;
	float _dataLo = 0.0f;
	float _dataHi = 65535.0f;

	// User-controlled display state.
	float _stretchLo = 0.0f;
	float _stretchHi = 65535.0f;
	bool _logStretch = false;
	bool _flipH = false;
	bool _flipV = false;

	// Pan / zoom state.
	double _zoom = 1.0;
	bool _fitMode = true;
	QPointF _pan{0.0, 0.0};

	bool _panning = false;
	QPoint _panAnchor;
	QPointF _panAnchorOffset;

	float _saturation = 1.0f;

	std::vector<DetectedStar> _stars;
	std::vector<AnnotationMarker> _annotations;
	ConstellationOverlay _constellations;
};

} // namespace astap::gui
