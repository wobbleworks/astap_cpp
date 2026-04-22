///----------------------------------------
///      @file preferences_dialog.h
///   @ingroup ASTAP++
///     @brief Modal Preferences dialog (Catalog / Site / Files tabs).
///    @author Created by John Stephen on 4/22/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <QDialog>

class QComboBox;
class QDoubleSpinBox;
class QLineEdit;
class QPushButton;

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @class PreferencesDialog
/// @brief Centralised editor for persistent user preferences.
/// @details Pulls existing values from @c QSettings (same keys the rest of
///          the app uses) and writes them back on OK. Tabs:
///          - Catalog — star database directory + preferred database name.
///          - Site — default observing location used by the SQM dialog.
///          - Files — clear the recent-files list, reset window layout.
///----------------------------------------

class PreferencesDialog final : public QDialog {
	Q_OBJECT

public:
	explicit PreferencesDialog(QWidget* parent = nullptr);
	~PreferencesDialog() override = default;

protected:
	void accept() override;

private slots:
	void browseDatabase();
	void clearRecents();
	void resetWindowLayout();

private:
	void buildLayout();
	void loadFromSettings();

	// Catalog tab
	QLineEdit*   _databasePath = nullptr;
	QPushButton* _browseButton = nullptr;
	QComboBox*   _databaseName = nullptr;

	// Site tab
	QLineEdit*      _latitude    = nullptr;
	QLineEdit*      _longitude   = nullptr;
	QDoubleSpinBox* _temperature = nullptr;
	QDoubleSpinBox* _pressure    = nullptr;

	// Files tab
	QPushButton* _clearRecentsButton = nullptr;
	QPushButton* _resetLayoutButton  = nullptr;
};

} // namespace astap::gui
