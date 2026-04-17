///----------------------------------------
///      @file image_inspector_dialog.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the Image Inspector dialog.
///    @author Created by John Stephen on 4/17/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "image_inspector_dialog.h"

#include "../../src/core/globals.h"

#include <QBrush>
#include <QColor>
#include <QDialogButtonBox>
#include <QHeaderView>
#include <QLabel>
#include <QTableWidget>
#include <QVBoxLayout>
#include <QtConcurrent>

#include <algorithm>
#include <array>
#include <cmath>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

namespace {

constexpr auto kCols = 3;
constexpr auto kRows = 3;

struct CellStats {
	int count = 0;
	double medianHfd = 0.0;
	double medianFwhm = 0.0;
};

// Median of a vector; leaves the vector reordered.
double median_of(std::vector<double>& v) {
	if (v.empty()) {
		return 0.0;
	}
	const auto n = v.size();
	std::nth_element(v.begin(), v.begin() + n / 2, v.end());
	auto mid = v[n / 2];
	if ((n % 2) == 0 && n > 1) {
		// Average with the next-lower value for an even-sized sample.
		std::nth_element(v.begin(), v.begin() + n / 2 - 1,
		                 v.begin() + n / 2);
		mid = (mid + v[n / 2 - 1]) * 0.5;
	}
	return mid;
}

std::array<CellStats, kRows * kCols> bin_stars(
		const std::vector<DetectedStar>& stars, int width, int height) {
	std::array<std::vector<double>, kRows * kCols> hfds;
	std::array<std::vector<double>, kRows * kCols> fwhms;

	const auto cellW = std::max(1, width / kCols);
	const auto cellH = std::max(1, height / kRows);

	for (const auto& s : stars) {
		// Stars are in FITS convention (1-based); convert to 0-based index.
		const auto col = std::clamp(static_cast<int>((s.x - 1) / cellW),
		                            0, kCols - 1);
		// FITS Y is bottom-up, but cell (0,0) is visually top-left. Flip.
		const auto rowFromTop = std::clamp(
			static_cast<int>((height - (s.y - 1)) / cellH),
			0, kRows - 1);
		const auto idx = rowFromTop * kCols + col;
		hfds[idx].push_back(s.hfd);
		fwhms[idx].push_back(s.fwhm);
	}

	std::array<CellStats, kRows * kCols> out{};
	for (auto i = 0; i < kRows * kCols; ++i) {
		out[i].count = static_cast<int>(hfds[i].size());
		out[i].medianHfd = median_of(hfds[i]);
		out[i].medianFwhm = median_of(fwhms[i]);
	}
	return out;
}

// Blend green (sharp) → yellow → red (soft) based on cell HFD vs best.
QColor shade_for_hfd(double hfd, double bestHfd, double worstHfd) {
	if (bestHfd <= 0.0 || worstHfd <= bestHfd + 0.01) {
		return QColor(64, 160, 64, 80);  // uniform: green
	}
	const auto t = std::clamp((hfd - bestHfd) / (worstHfd - bestHfd),
	                          0.0, 1.0);
	// Green (64,160,64) → yellow (200,180,48) → red (200,64,64)
	const auto r = 64 + static_cast<int>(std::round((200 - 64) * t));
	const auto g = 160 + static_cast<int>(std::round((180 - 160) * (1.0 - std::abs(t - 0.5) * 2.0)));
	const auto b = 64;
	return QColor(r, g, b, 80);
}

}  // namespace

ImageInspectorDialog::ImageInspectorDialog(QWidget* parent) :
	QDialog(parent) {
	setWindowTitle(tr("Image Inspector"));
	resize(520, 460);
	buildUi();
}

ImageInspectorDialog::~ImageInspectorDialog() = default;

void ImageInspectorDialog::buildUi() {
	auto* root = new QVBoxLayout(this);

	_status = new QLabel(tr("Analysing image…"), this);
	_status->setStyleSheet("color: gray;");
	root->addWidget(_status);

	_table = new QTableWidget(kRows, kCols, this);
	_table->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
	_table->verticalHeader()->setSectionResizeMode(QHeaderView::Stretch);
	_table->horizontalHeader()->setVisible(false);
	_table->verticalHeader()->setVisible(false);
	_table->setSelectionMode(QAbstractItemView::NoSelection);
	_table->setFocusPolicy(Qt::NoFocus);
	_table->setEditTriggers(QAbstractItemView::NoEditTriggers);
	root->addWidget(_table, 1);

	_summary = new QLabel(this);
	_summary->setWordWrap(true);
	root->addWidget(_summary);

	auto* buttons = new QDialogButtonBox(QDialogButtonBox::Close, this);
	connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::close);
	root->addWidget(buttons);
}

void ImageInspectorDialog::analyseCurrentImage() {
	if (astap::img_loaded.empty()) {
		_status->setText(tr("No image loaded."));
		return;
	}

	_imageWidth = astap::head.width;
	_imageHeight = astap::head.height;

	// Copy the image for the worker thread; detect_stars promises thread
	// safety only against its own scratch, and we don't want to race a
	// viewer refresh.
	const auto imgCopy = astap::img_loaded;

	_watcher = std::make_unique<QFutureWatcher<DetectionResult>>();
	connect(_watcher.get(), &QFutureWatcher<DetectionResult>::finished,
	        this, &ImageInspectorDialog::onDetectionFinished);

	_status->setText(tr("Detecting stars…"));
	_watcher->setFuture(QtConcurrent::run([img = std::move(imgCopy)]() {
		return detect_stars(img, /*snr_min=*/5.0);
	}));
}

void ImageInspectorDialog::onDetectionFinished() {
	if (!_watcher) {
		return;
	}
	const auto result = _watcher->result();
	_watcher.reset();

	populateTable(result);
}

void ImageInspectorDialog::populateTable(const DetectionResult& result) {
	const auto cells = bin_stars(result.stars, _imageWidth, _imageHeight);

	// Find best / worst HFD across cells that have any stars.
	auto bestHfd = 1e9;
	auto worstHfd = 0.0;
	for (const auto& c : cells) {
		if (c.count > 0) {
			bestHfd = std::min(bestHfd, c.medianHfd);
			worstHfd = std::max(worstHfd, c.medianHfd);
		}
	}
	if (bestHfd == 1e9) {
		_status->setText(tr("No stars detected — check exposure / focus."));
		_summary->clear();
		return;
	}

	for (auto r = 0; r < kRows; ++r) {
		for (auto c = 0; c < kCols; ++c) {
			const auto& cell = cells[r * kCols + c];
			auto text = (cell.count == 0)
				? tr("(empty)")
				: tr("%1 stars\nHFD %2\nFWHM %3")
					.arg(cell.count)
					.arg(cell.medianHfd, 0, 'f', 2)
					.arg(cell.medianFwhm, 0, 'f', 2);

			auto* item = new QTableWidgetItem(text);
			item->setTextAlignment(Qt::AlignCenter);
			if (cell.count > 0) {
				item->setBackground(shade_for_hfd(cell.medianHfd,
				                                  bestHfd, worstHfd));
			}
			_table->setItem(r, c, item);
		}
	}

	const auto spread = (worstHfd - bestHfd) / bestHfd * 100.0;
	auto verdict = QString{};
	if (spread < 10.0) {
		verdict = tr("Flat field — no significant tilt.");
	} else if (spread < 25.0) {
		verdict = tr("Mild corner softening — check collimation / "
		             "backspacing.");
	} else {
		verdict = tr("Significant tilt — HFD varies by %1% across field.")
			.arg(spread, 0, 'f', 0);
	}

	_status->setText(tr("%1 stars detected, median HFD %2.")
		.arg(static_cast<int>(result.stars.size()))
		.arg(result.medianHfd, 0, 'f', 2));
	_summary->setText(tr("Best-cell HFD %1 · Worst-cell HFD %2 · "
	                     "Spread %3%.\n%4")
		.arg(bestHfd, 0, 'f', 2)
		.arg(worstHfd, 0, 'f', 2)
		.arg(spread, 0, 'f', 0)
		.arg(verdict));
}

} // namespace astap::gui
