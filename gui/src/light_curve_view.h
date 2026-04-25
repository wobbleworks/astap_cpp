///----------------------------------------
///      @file light_curve_view.h
///   @ingroup ASTAP++
///     @brief Compact JD-vs-magnitude plot widget for time-series photometry.
///   @details Mirrors the @c plot_graph routine in Pascal's @c unit_aavso.pas.
///            Renders the variable, check, and (optional) comp series for a
///            collected @c astap::core::AavsoCollectResult so the user can
///            see the light curve without leaving the dialog.
///
///            Magnitude axis is inverted (brighter at top) per astronomical
///            convention. Empty data yields a blank plot with explanatory
///            text.
///    @author Created by John Stephen on 4/25/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <QWidget>

#include <vector>

namespace astap::core { struct AavsoMeasurement; }

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @class LightCurveView
/// @brief QWidget that paints a JD-vs-magnitude scatter plot.
///----------------------------------------

class LightCurveView final : public QWidget {
	Q_OBJECT

public:
	explicit LightCurveView(QWidget* parent = nullptr);

	/// @brief Set the rows to plot. Pass an empty vector to clear.
	void setRows(const std::vector<astap::core::AavsoMeasurement>& rows);

	/// @brief Toggle whether the comp-magnitude series is included.
	void setShowComp(bool show);

protected:
	void paintEvent(QPaintEvent* event) override;
	[[nodiscard]] QSize sizeHint() const override;

private:
	std::vector<double> _jd;
	std::vector<double> _var;
	std::vector<double> _chk;
	std::vector<double> _cmp;   // 0.0 means "no comp value"
	bool _showComp = true;
};

} // namespace astap::gui
