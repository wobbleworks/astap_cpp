///----------------------------------------
///      @file controls_panel.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the right-hand display controls panel.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "controls_panel.h"
#include "histogram_widget.h"
#include "image_viewer.h"

#include <QCheckBox>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QPushButton>
#include <QSlider>
#include <QSpinBox>
#include <QTableWidget>
#include <QVBoxLayout>

#include <algorithm>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

namespace {

// Stretch sliders / spin boxes operate in the working 0..65535 domain that
// load_image normalises sources into.
constexpr int kSliderMin = 0;
constexpr int kSliderMax = 65535;

} // namespace

///----------------------------------------
/// MARK: ControlsPanel
///----------------------------------------

ControlsPanel::ControlsPanel(QWidget* parent) :
	QWidget(parent) {

	auto* root = new QVBoxLayout(this);
	root->setContentsMargins(8, 8, 8, 8);
	root->setSpacing(8);

	// Histogram
	_histogram = new HistogramWidget(this);
	root->addWidget(_histogram);

	// Stretch group
	auto* stretchGroup = new QGroupBox(tr("Stretch"), this);
	auto* stretchLayout = new QVBoxLayout(stretchGroup);
	stretchLayout->setSpacing(4);

	auto makeRow = [stretchGroup](const QString& label, QSlider*& slider, QSpinBox*& spin) {
		auto* row = new QHBoxLayout();
		row->setSpacing(6);
		auto* lbl = new QLabel(label, stretchGroup);
		lbl->setMinimumWidth(36);
		slider = new QSlider(Qt::Horizontal, stretchGroup);
		slider->setRange(kSliderMin, kSliderMax);
		slider->setSingleStep(64);
		slider->setPageStep(2048);
		spin = new QSpinBox(stretchGroup);
		spin->setRange(kSliderMin, kSliderMax);
		spin->setButtonSymbols(QAbstractSpinBox::NoButtons);
		spin->setMaximumWidth(72);
		row->addWidget(lbl);
		row->addWidget(slider, 1);
		row->addWidget(spin);
		return row;
	};

	// Black-point row
	stretchLayout->addLayout(makeRow(tr("Black"), _loSlider, _loSpin));

	// White-point row
	stretchLayout->addLayout(makeRow(tr("White"), _hiSlider, _hiSpin));

	// Auto / log row
	auto* btnRow = new QHBoxLayout();
	_autoButton = new QPushButton(tr("Auto"), stretchGroup);
	_logCheck = new QCheckBox(tr("Logarithmic"), stretchGroup);
	btnRow->addWidget(_autoButton);
	btnRow->addStretch(1);
	btnRow->addWidget(_logCheck);
	stretchLayout->addLayout(btnRow);

	root->addWidget(stretchGroup);

	// Display group
	auto* displayGroup = new QGroupBox(tr("Display"), this);
	auto* displayLayout = new QVBoxLayout(displayGroup);
	auto* flipRow = new QHBoxLayout();
	_flipHCheck = new QCheckBox(tr("Flip H"), displayGroup);
	_flipVCheck = new QCheckBox(tr("Flip V"), displayGroup);
	flipRow->addWidget(_flipHCheck);
	flipRow->addWidget(_flipVCheck);
	flipRow->addStretch(1);
	displayLayout->addLayout(flipRow);

	auto* satRow = new QHBoxLayout();
	auto* satLabel = new QLabel(tr("Saturation"), displayGroup);
	_satSlider = new QSlider(Qt::Horizontal, displayGroup);
	_satSlider->setRange(0, 300);
	_satSlider->setValue(100);
	_satSlider->setTickInterval(50);
	satRow->addWidget(satLabel);
	satRow->addWidget(_satSlider, 1);
	displayLayout->addLayout(satRow);

	root->addWidget(displayGroup);

	// Statistics group
	auto* statsGroup = new QGroupBox(tr("Statistics"), this);
	auto* statsLayout = new QVBoxLayout(statsGroup);
	_statsTable = new QTableWidget(0, 7, statsGroup);
	_statsTable->setHorizontalHeaderLabels(
		{tr("Ch"), tr("Min"), tr("Max"), tr("Mean"), tr("Bgnd"), tr("Noise"), tr("Stars")});
	_statsTable->verticalHeader()->setVisible(false);
	_statsTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
	_statsTable->setSelectionMode(QAbstractItemView::NoSelection);
	_statsTable->setFocusPolicy(Qt::NoFocus);
	_statsTable->setShowGrid(false);
	_statsTable->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
	_statsTable->horizontalHeader()->setStretchLastSection(true);
	_statsTable->setSizeAdjustPolicy(QAbstractScrollArea::AdjustToContents);
	statsLayout->addWidget(_statsTable);
	root->addWidget(statsGroup);

	// Push everything to the top
	root->addStretch(1);

	// Wire control signals
	connect(_loSlider, &QSlider::valueChanged, this, &ControlsPanel::onLoSliderChanged);
	connect(_hiSlider, &QSlider::valueChanged, this, &ControlsPanel::onHiSliderChanged);
	connect(_loSpin,  qOverload<int>(&QSpinBox::valueChanged), this, &ControlsPanel::onLoSpinChanged);
	connect(_hiSpin,  qOverload<int>(&QSpinBox::valueChanged), this, &ControlsPanel::onHiSpinChanged);
	connect(_logCheck, &QCheckBox::toggled, this, &ControlsPanel::onLogToggled);
	connect(_autoButton, &QPushButton::clicked, this, &ControlsPanel::onAutoClicked);
	connect(_flipHCheck, &QCheckBox::toggled, this, &ControlsPanel::onFlipHToggled);
	connect(_flipVCheck, &QCheckBox::toggled, this, &ControlsPanel::onFlipVToggled);
	connect(_satSlider, &QSlider::valueChanged, this, [this](int v) {
		if (!_suppress && _viewer) {
			_viewer->setSaturation(v / 100.0f);
		}
	});
}

void ControlsPanel::attachViewer(ImageViewer* viewer) {
	if (_viewer == viewer) {
		return;
	}

	if (_viewer) {
		disconnect(_viewer, nullptr, this, nullptr);
	}

	_viewer = viewer;
	_histogram->attachViewer(viewer);

	if (_viewer) {
		connect(_viewer, &ImageViewer::imageLoaded,
			this, &ControlsPanel::onViewerLoaded, Qt::UniqueConnection);
		connect(_viewer, &ImageViewer::stretchChanged,
			this, &ControlsPanel::onViewerStretchChanged, Qt::UniqueConnection);
		connect(_viewer, &ImageViewer::flipChanged,
			this, &ControlsPanel::onViewerFlipChanged, Qt::UniqueConnection);
	}

	syncFromViewer();
}

void ControlsPanel::onLoSliderChanged(int v) {
	if (_suppress) {
		return;
	}
	// Keep lo < hi by nudging the high slider up if needed.
	if (v >= _hiSlider->value()) {
		_suppress = true;
		_hiSlider->setValue(std::min(kSliderMax, v + 1));
		_hiSpin->setValue(_hiSlider->value());
		_suppress = false;
	}
	_suppress = true;
	_loSpin->setValue(v);
	_suppress = false;
	pushStretchToViewer();
}

void ControlsPanel::onHiSliderChanged(int v) {
	if (_suppress) {
		return;
	}
	// Keep hi > lo by nudging the low slider down if needed.
	if (v <= _loSlider->value()) {
		_suppress = true;
		_loSlider->setValue(std::max(kSliderMin, v - 1));
		_loSpin->setValue(_loSlider->value());
		_suppress = false;
	}
	_suppress = true;
	_hiSpin->setValue(v);
	_suppress = false;
	pushStretchToViewer();
}

void ControlsPanel::onLoSpinChanged(int v) {
	if (_suppress) {
		return;
	}
	_suppress = true;
	_loSlider->setValue(v);
	_suppress = false;
	pushStretchToViewer();
}

void ControlsPanel::onHiSpinChanged(int v) {
	if (_suppress) {
		return;
	}
	_suppress = true;
	_hiSlider->setValue(v);
	_suppress = false;
	pushStretchToViewer();
}

void ControlsPanel::onLogToggled(bool on) {
	if (_viewer) {
		_viewer->setLogarithmicStretch(on);
	}
}

void ControlsPanel::onAutoClicked() {
	if (_viewer) {
		_viewer->autoStretch();
	}
}

void ControlsPanel::onFlipHToggled(bool on) {
	if (_viewer) {
		_viewer->setFlipHorizontal(on);
	}
}

void ControlsPanel::onFlipVToggled(bool on) {
	if (_viewer) {
		_viewer->setFlipVertical(on);
	}
}

void ControlsPanel::onViewerLoaded() {
	syncFromViewer();
	refreshStats();
}

void ControlsPanel::onViewerStretchChanged() {
	syncFromViewer();
}

void ControlsPanel::onViewerFlipChanged() {
	syncFromViewer();
}

void ControlsPanel::syncFromViewer() {
	if (!_viewer) {
		return;
	}

	_suppress = true;
	const auto lo = std::clamp(static_cast<int>(_viewer->stretchLo() + 0.5f),
		kSliderMin, kSliderMax);
	const auto hi = std::clamp(static_cast<int>(_viewer->stretchHi() + 0.5f),
		kSliderMin, kSliderMax);
	_loSlider->setValue(lo);
	_loSpin->setValue(lo);
	_hiSlider->setValue(hi);
	_hiSpin->setValue(hi);
	_logCheck->setChecked(_viewer->logarithmicStretch());
	_flipHCheck->setChecked(_viewer->flipHorizontal());
	_flipVCheck->setChecked(_viewer->flipVertical());
	_satSlider->setValue(static_cast<int>(_viewer->saturation() * 100.0f + 0.5f));
	_suppress = false;
}

void ControlsPanel::refreshStats() {
	if (!_viewer) {
		_statsTable->setRowCount(0);
		return;
	}

	const auto& s = _viewer->stats();
	_statsTable->setRowCount(static_cast<int>(s.size()));

	static constexpr const char* kMonoLabels[] = {"Y"};
	static constexpr const char* kRgbLabels[] = {"R", "G", "B"};
	const bool mono = (s.size() == 1);

	for (int row = 0; row < static_cast<int>(s.size()); ++row) {
		const auto& st = s[static_cast<std::size_t>(row)];
		const auto label = mono
			? QString::fromLatin1(kMonoLabels[0])
			: (row < 3
				? QString::fromLatin1(kRgbLabels[row])
				: QString::number(row));

		auto setCell = [this, row](int col, const QString& text) {
			auto* item = new QTableWidgetItem(text);
			item->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
			_statsTable->setItem(row, col, item);
		};

		auto fmt = [](float v) {
			return QString::number(v, 'f', (std::abs(v) < 10.0f) ? 2 : 0);
		};

		setCell(0, label);
		setCell(1, fmt(st.min));
		setCell(2, fmt(st.max));
		setCell(3, fmt(st.mean));
		setCell(4, fmt(st.background));
		setCell(5, fmt(st.noise));
		setCell(6, fmt(st.starLevel));
	}

	_statsTable->resizeColumnsToContents();
}

void ControlsPanel::pushStretchToViewer() {
	if (!_viewer || _suppress) {
		return;
	}
	_viewer->setStretchRange(static_cast<float>(_loSlider->value()),
	                         static_cast<float>(_hiSlider->value()));
}

} // namespace astap::gui
