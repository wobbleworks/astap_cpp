///----------------------------------------
///      @file light_curve_view.cpp
///   @ingroup ASTAP++
///     @brief Implementation of @ref astap::gui::LightCurveView.
///    @author Created by John Stephen on 4/25/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "light_curve_view.h"

#include "../../src/core/aavso_report.h"

#include <QFont>
#include <QFontMetrics>
#include <QPaintEvent>
#include <QPainter>
#include <QPen>
#include <QRectF>

#include <algorithm>
#include <cmath>
#include <limits>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

LightCurveView::LightCurveView(QWidget* parent) : QWidget(parent) {
	setMinimumHeight(160);
	setBackgroundRole(QPalette::Base);
	setAutoFillBackground(true);
}

void LightCurveView::setRows(
		const std::vector<astap::core::AavsoMeasurement>& rows) {
	_jd.clear();
	_var.clear();
	_chk.clear();
	_cmp.clear();
	_jd.reserve(rows.size());
	_var.reserve(rows.size());
	_chk.reserve(rows.size());
	_cmp.reserve(rows.size());
	for (const auto& m : rows) {
		_jd.push_back(m.jd);
		_var.push_back(m.var_magnitude);
		_chk.push_back(m.check_magnitude);
		_cmp.push_back(m.comp_magnitude);
	}
	update();
}

void LightCurveView::setShowComp(bool show) {
	_showComp = show;
	update();
}

QSize LightCurveView::sizeHint() const {
	return {480, 200};
}

void LightCurveView::paintEvent(QPaintEvent*) {
	auto p = QPainter{this};
	p.setRenderHint(QPainter::Antialiasing, true);
	p.fillRect(rect(), palette().base());

	if (_jd.empty()) {
		p.setPen(palette().mid().color());
		p.drawText(rect(), Qt::AlignCenter,
			tr("No measurements yet — run “Measure all frames” to plot."));
		return;
	}

	// Domain of the plot. We accept the comp series only when it has at
	// least one non-zero value (rows from ensemble mode set comp = 0).
	auto jd_min = *std::min_element(_jd.begin(), _jd.end());
	auto jd_max = *std::max_element(_jd.begin(), _jd.end());
	if (jd_max == jd_min) jd_max = jd_min + 1.0;  // single point ⇒ avoid /0

	auto mag_min = std::numeric_limits<double>::infinity();
	auto mag_max = -std::numeric_limits<double>::infinity();
	const auto track = [&](double v) {
		if (v == 0.0) return;
		mag_min = std::min(mag_min, v);
		mag_max = std::max(mag_max, v);
	};
	for (auto v : _var) track(v);
	for (auto v : _chk) track(v);
	if (_showComp) {
		for (auto v : _cmp) track(v);
	}
	if (!std::isfinite(mag_min) || !std::isfinite(mag_max)) {
		p.setPen(palette().mid().color());
		p.drawText(rect(), Qt::AlignCenter,
			tr("All magnitudes are zero — calibration may be missing."));
		return;
	}

	// Pad the magnitude axis so the points don't sit on the frame.
	const auto range = std::max(mag_max - mag_min, 0.05);
	mag_min -= range * 0.05;
	mag_max += range * 0.05;

	const auto fm     = QFontMetrics{font()};
	const auto leftPx = fm.horizontalAdvance("12.345") + 10;
	const auto rightPx = 8;
	const auto topPx  = fm.height() + 6;
	const auto botPx  = fm.height() + 8;

	const auto plot = QRectF(leftPx, topPx,
	                         std::max(40, width()  - leftPx - rightPx),
	                         std::max(40, height() - topPx  - botPx));

	// Frame.
	p.setPen(palette().mid().color());
	p.drawRect(plot);

	// Axis labels.
	p.setPen(palette().text().color());
	p.drawText(QPointF(4, fm.ascent() + 2), tr("Mag"));

	const auto jdLabel = tr("JD");
	const auto jdLW = fm.horizontalAdvance(jdLabel);
	p.drawText(QPointF(plot.right() - jdLW,
	                   plot.bottom() + fm.ascent() + 4), jdLabel);

	// Y ticks (5 — every 25%).
	for (auto i = 0; i <= 4; ++i) {
		const auto t = i / 4.0;
		const auto y = plot.top() + plot.height() * t;
		p.drawLine(QPointF(plot.left() - 3, y), QPointF(plot.left(), y));
		const auto m = mag_min + (mag_max - mag_min) * t;
		p.drawText(QPointF(2, y + fm.ascent() / 2.0 - 1),
			QString::number(m, 'f', 2));
	}
	// X ticks (3).
	for (auto i = 0; i <= 3; ++i) {
		const auto t = i / 3.0;
		const auto x = plot.left() + plot.width() * t;
		p.drawLine(QPointF(x, plot.bottom()),
		           QPointF(x, plot.bottom() + 3));
		const auto j = jd_min + (jd_max - jd_min) * t;
		const auto label = QString::number(j, 'f', 4);
		const auto labW  = fm.horizontalAdvance(label);
		auto labX = x - labW / 2.0;
		labX = std::clamp<double>(labX, plot.left(),
		                          plot.right() - labW);
		p.drawText(QPointF(labX, plot.bottom() + fm.ascent() + 3), label);
	}

	const auto plotPoint = [&](double jd, double mag, const QColor& c) {
		if (mag == 0.0) return;
		const auto tx = (jd - jd_min) / (jd_max - jd_min);
		const auto ty = (mag - mag_min) / (mag_max - mag_min);
		// Magnitude axis is inverted: brighter (smaller mag) at the top.
		const auto px = plot.left() + plot.width() * tx;
		const auto py = plot.top()  + plot.height() * ty;
		p.setPen(QPen(c, 1.5));
		p.setBrush(c);
		p.drawEllipse(QPointF(px, py), 3.0, 3.0);
	};

	// Series legend at the top of the plot (above the frame).
	auto legendX = double(leftPx);
	const auto legendY = fm.ascent();
	const auto drawLegend = [&](const QString& label, const QColor& c) {
		p.setPen(QPen(c, 2));
		p.setBrush(c);
		p.drawEllipse(QPointF(legendX + 4, legendY - fm.ascent() / 2.0 + 2),
			3.0, 3.0);
		p.setPen(palette().text().color());
		p.drawText(QPointF(legendX + 14, legendY + 2), label);
		legendX += 20 + fm.horizontalAdvance(label);
	};
	drawLegend(tr("Variable"), QColor(192, 47, 47));
	drawLegend(tr("Check"),    QColor(47, 142, 47));
	if (_showComp) drawLegend(tr("Comp"), QColor(53, 96, 173));

	// Series.
	const auto cVar = QColor(192, 47, 47);
	const auto cChk = QColor(47, 142, 47);
	const auto cCmp = QColor(53, 96, 173);
	for (auto i = std::size_t{0}; i < _jd.size(); ++i) {
		plotPoint(_jd[i], _var[i], cVar);
		plotPoint(_jd[i], _chk[i], cChk);
		if (_showComp) plotPoint(_jd[i], _cmp[i], cCmp);
	}
}

} // namespace astap::gui
