///----------------------------------------
///      @file stack_window.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the stacking window.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "stack_window.h"
#include "image_viewer.h"

#include "../../src/core/fits.h"
#include "../../src/core/globals.h"
#include "../../src/core/image_io.h"
#include "../../src/stacking/stack.h"
#include "../../src/stacking/stack_driver.h"   // create_master
#include "../../src/stacking/stack_routines.h"

#include <QCheckBox>
#include <QComboBox>
#include <QCoreApplication>
#include <QDebug>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QProgressBar>
#include <QPushButton>
#include <QSettings>
#include <QSpinBox>
#include <QRegularExpression>
#include <QSet>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QTabWidget>
#include <QVBoxLayout>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <functional>
#include <map>
#include <vector>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

namespace {

enum MethodIndex { kMethodAverage = 0, kMethodSigmaClip = 1, kMethodLRGB = 2 };
enum AlignmentIndex {
	kAlignStarMatch = 0,
	kAlignAstrometric = 1,
	kAlignManual = 2,
	kAlignNone = 3,
	kAlignWcs = 4
};

// Channel role for LRGB stacking. "Light" is untagged / treated as a mono
// frame by Average and Sigma-clip methods.
enum Channel {
	kChanLight = 0,
	kChanL = 1,
	kChanR = 2,
	kChanG = 3,
	kChanB = 4,
	kChanRGB = 5,
};

QStringList channelLabels() {
	return {
		QObject::tr("Light"),
		QStringLiteral("L"),
		QStringLiteral("R"),
		QStringLiteral("G"),
		QStringLiteral("B"),
		QStringLiteral("RGB"),
	};
}

// Read FITS header only (skip pixel data) to get width × height. Returns
// area = w * h, or 0 on failure.
long long probe_area(const QString& path) {
	astap::Header h;
	astap::ImageArray img;
	auto memo = std::vector<std::string>{};
	if (!astap::core::load_fits(std::filesystem::path(path.toStdString()),
	                            /*light=*/true, /*load_data=*/false,
	                            /*update_memo=*/false, /*get_ext=*/0,
	                            memo, h, img)) {
		return 0;
	}
	return static_cast<long long>(h.width) * h.height;
}

// Guess a channel from a filename. Matches "_L_", "_R_", "_G_", "_B_",
// "_RGB_", "_Lum_", etc. Falls back to Light.
Channel guess_channel(const QString& path) {
	const auto name = QFileInfo(path).completeBaseName().toLower();
	// Check RGB first so "_rgb_" doesn't match as just R.
	const auto tokens = name.split(QRegularExpression("[^a-z0-9]+"),
	                               Qt::SkipEmptyParts);
	for (const auto& t : tokens) {
		if (t == "rgb") return kChanRGB;
		if (t == "l" || t == "lum" || t == "luminance") return kChanL;
		if (t == "r" || t == "red")   return kChanR;
		if (t == "g" || t == "green") return kChanG;
		if (t == "b" || t == "blue")  return kChanB;
	}
	return kChanLight;
}

}  // namespace

StackWindow::StackWindow(QWidget* parent) :
	QWidget(parent, Qt::Window) {

	setWindowTitle(tr("Stack"));
	resize(520, 600);

	auto* root = new QVBoxLayout(this);

	_tabs = new QTabWidget(this);
	root->addWidget(_tabs, 1);

	buildLightsTab();
	buildDarksTab();
	buildFlatsTab();
	buildFlatDarksTab();
	buildSettingsTab();

	// Shared "Classify by" bar, visible under the tabs (Pascal classify_groupbox1).
	root->addWidget(buildClassifyBar());

	_phaseLabel = new QLabel(this);
	_phaseLabel->setStyleSheet("color: gray;");
	_phaseLabel->setVisible(false);
	root->addWidget(_phaseLabel);

	_progress = new QProgressBar(this);
	_progress->setRange(0, 100);
	_progress->setValue(0);
	_progress->setTextVisible(true);
	root->addWidget(_progress);

	auto* stackRow = new QHBoxLayout();
	_stackButton = new QPushButton(tr("Stack"), this);
	_stackButton->setDefault(true);
	stackRow->addStretch(1);
	stackRow->addWidget(_stackButton);
	root->addLayout(stackRow);

	connect(_stackButton, &QPushButton::clicked, this, &StackWindow::startStack);
}

namespace {

// Column layout shared by the dark and flat candidate tables (mirrors the Pascal
// listview columns the port's Header can supply; Type/Background/σ are omitted —
// no IMAGETYP field / no pixel pass here).
enum MasterCol {
	kMFile = 0, kMExp, kMTemp, kMBin, kMSize, kMFilter, kMGain, kMDate, kMJd,
	kMCalstat, kMCompat, kMColCount
};

// Read the first Lights-tab frame's header (no pixels) as the reference light for
// the Compatibility column. Returns false when there are no lights.
[[nodiscard]] bool first_light_header(QTableWidget* lights, astap::Header& out) {
	if (!lights || lights->rowCount() == 0) {
		return false;
	}
	auto* item = lights->item(0, 0);
	if (!item) {
		return false;
	}
	const auto path = item->data(Qt::UserRole).toString();
	auto memo = std::vector<std::string>{};
	auto img  = astap::ImageArray{};
	return astap::core::load_fits(std::filesystem::path(path.toStdString()),
		/*light=*/true, /*load_data=*/false, /*update_memo=*/false, /*get_ext=*/0,
		memo, out, img);
}

// Assemble the match options from the shared "Classify by" controls.
[[nodiscard]] astap::stacking::MasterMatchOptions match_options_from(
		QCheckBox* exp, QCheckBox* temp, QCheckBox* gain, QCheckBox* filter,
		QSpinBox* deltaTemp) {
	astap::stacking::MasterMatchOptions o;
	o.classify_exposure    = exp && exp->isChecked();
	o.classify_temperature = temp && temp->isChecked();
	o.classify_gain        = gain && gain->isChecked();
	o.classify_filter      = filter && filter->isChecked();
	o.delta_temp           = deltaTemp ? deltaTemp->value() : 3;
	return o;
}

// Format an exposure the way the Pascal listview does (unit_stack.pas:4415-4417):
// integer seconds for >= 10 s, otherwise up to 6 significant figures. The flat
// path keys its grouping / flat-dark pairing / master filename off this string
// (Pascal compares the raw F_exposure display string, not a rounded value), so a
// 0.8 s and a 1.2 s flat stay in separate masters instead of both rounding to 1.
[[nodiscard]] QString format_exposure(double seconds) {
	if (seconds >= 10.0) {
		return QString::number(static_cast<long long>(std::lrint(seconds)));
	}
	return QString::number(seconds, 'g', 6);
}

// Create an empty candidate table with the shared columns.
[[nodiscard]] QTableWidget* make_master_table(QWidget* parent) {
	auto* t = new QTableWidget(0, kMColCount, parent);
	t->setHorizontalHeaderLabels({
		QObject::tr("File"), QObject::tr("Exposure"), QObject::tr("Temperature"),
		QObject::tr("Binning"), QObject::tr("Size"), QObject::tr("Filter"),
		QObject::tr("Gain"), QObject::tr("Date"), QObject::tr("JD"),
		QObject::tr("Calibration"), QObject::tr("Compatibility")});
	t->horizontalHeader()->setSectionResizeMode(kMFile, QHeaderView::Stretch);
	t->horizontalHeader()->setSectionResizeMode(kMCompat, QHeaderView::Stretch);
	t->verticalHeader()->setVisible(false);
	t->setSelectionBehavior(QAbstractItemView::SelectRows);
	t->setSelectionMode(QAbstractItemView::ExtendedSelection);
	t->setEditTriggers(QAbstractItemView::NoEditTriggers);
	t->setWordWrap(false);
	// Click-to-sort with indicators (Pascal AutoSortIndicator); default by Date.
	t->setSortingEnabled(true);
	t->sortByColumn(kMDate, Qt::AscendingOrder);
	return t;
}

}  // namespace

void StackWindow::addMasterFiles(QTableWidget* table, const QString& title) {
	QSettings settings;
	const auto lastDir = settings.value("files/lastCalDir").toString();
	const auto paths = QFileDialog::getOpenFileNames(this,
		tr("Add %1 frames").arg(title), lastDir,
		tr("FITS images (*.fit *.fits *.fts);;All files (*)"));
	if (paths.isEmpty()) {
		return;
	}
	settings.setValue("files/lastCalDir", QFileInfo(paths.first()).absolutePath());

	auto existing = QSet<QString>{};
	for (int r = 0; r < table->rowCount(); ++r) {
		existing.insert(table->item(r, kMFile)->data(Qt::UserRole).toString());
	}

	// Populate with sorting off so a mid-insert re-sort can't shuffle half-filled
	// rows; restore + re-sort afterward.
	const bool wasSorting = table->isSortingEnabled();
	table->setSortingEnabled(false);
	for (const auto& p : paths) {
		if (!existing.contains(p)) {
			insertMasterRow(table, p);   // unreadable headers are skipped
		}
	}
	table->setSortingEnabled(wasSorting);
	refreshCompatibility();
}

// Analyse one candidate FITS (header only) and append a fully-populated row.
// Shared by addMasterFiles and the Replace-by-master flow. Returns false (and
// adds nothing) when the header cannot be read.
bool StackWindow::insertMasterRow(QTableWidget* table, const QString& path) {
	astap::stacking::MasterMetadata m;
	if (!astap::stacking::analyse_master(
	        std::filesystem::path(path.toStdString()), m)) {
		return false;  // unreadable header — skip
	}
	_masterMeta.insert(path, m);   // cache for refreshCompatibility

	const int row = table->rowCount();
	table->insertRow(row);

	auto* fileItem = new QTableWidgetItem(QFileInfo(path).fileName());
	fileItem->setData(Qt::UserRole, path);
	fileItem->setToolTip(path);
	fileItem->setFlags(fileItem->flags() | Qt::ItemIsUserCheckable);
	fileItem->setCheckState(Qt::Checked);
	table->setItem(row, kMFile, fileItem);

	auto cell = [&](int col, const QString& text) {
		auto* it = new QTableWidgetItem(text);
		it->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
		table->setItem(row, col, it);
	};
	cell(kMExp,     m.exposure_seconds > 0.0 ? format_exposure(m.exposure_seconds)
	                                         : QStringLiteral("—"));
	cell(kMTemp,    QString::number(m.set_temperature));
	cell(kMBin,     QString::number(m.xbinning, 'g', 3));
	cell(kMSize,    QString("%1×%2").arg(m.width).arg(m.height));
	cell(kMFilter,  QString::fromStdString(m.filter_name));
	cell(kMGain,    QString::fromStdString(m.gain));
	// Calendar date (the YYYY-MM-DD portion of DATE-OBS) — Pascal's default sort.
	cell(kMDate,    QString::fromStdString(m.date_obs).left(10));
	cell(kMJd,      QString::number(m.jd, 'f', 3));
	// CALSTAT: a non-empty value marks a flat that has already been calibrated
	// (flat-darks subtracted) — the Replace-by-master-flat flow skips those.
	cell(kMCalstat, QString::fromStdString(m.calstat));
	cell(kMCompat,  QString());
	return true;
}

void StackWindow::removeSelectedRows(QTableWidget* table) {
	auto rows = QList<int>{};
	for (const auto* item : table->selectedItems()) {
		if (item->column() == kMFile) {
			rows.append(item->row());
		}
	}
	std::sort(rows.begin(), rows.end(), std::greater<int>());
	for (const auto r : rows) {
		table->removeRow(r);
	}
}

void StackWindow::refreshCompatibility() {
	astap::Header light;
	const bool haveLight = first_light_header(_fileTable, light);
	const auto opts = match_options_from(_classDarkExp, _classDarkTemp,
		_classDarkGain, _classFlatFilter, _deltaTemp);
	updateCompatibility(_darkTable, /*isFlat=*/false, light, opts, haveLight);
	updateCompatibility(_flatTable, /*isFlat=*/true,  light, opts, haveLight);
}

// Kept as a plain member (not a nested lambda) and free of iterator-deref /
// ternary-of-calls: those tripped an MSVC front-end ICE (C1001) here.
void StackWindow::updateCompatibility(QTableWidget* table, bool isFlat,
		const astap::Header& light,
		const astap::stacking::MasterMatchOptions& opts, bool haveLight) {
	for (int r = 0; r < table->rowCount(); ++r) {
		auto* fileItem = table->item(r, kMFile);
		if (!fileItem) {
			continue;
		}
		QString text = QStringLiteral("—");
		const QString path = fileItem->data(Qt::UserRole).toString();
		if (haveLight && _masterMeta.contains(path)) {
			const astap::stacking::MasterMetadata m = _masterMeta.value(path);
			std::string issue;
			if (isFlat) {
				issue = astap::stacking::master_flat_issue(light, m, opts);
			} else {
				issue = astap::stacking::master_dark_issue(light, m, opts);
			}
			if (issue.empty()) {
				text = tr("OK");
			} else {
				text = QString::fromStdString(issue);
			}
		}
		auto* it = table->item(r, kMCompat);
		if (!it) {
			it = new QTableWidgetItem();
			table->setItem(r, kMCompat, it);
		}
		it->setText(text);
	}
}

void StackWindow::applyLibrariesToEngine() {
	auto gather = [](QTableWidget* table) {
		auto paths = std::vector<std::filesystem::path>{};
		for (int r = 0; r < table->rowCount(); ++r) {
			auto* it = table->item(r, kMFile);
			if (it && it->checkState() == Qt::Checked) {
				paths.emplace_back(it->data(Qt::UserRole).toString().toStdString());
			}
		}
		return paths;
	};
	astap::stacking::set_master_match_options(match_options_from(
		_classDarkExp, _classDarkTemp, _classDarkGain, _classFlatFilter, _deltaTemp));
	const auto darks = gather(_darkTable);
	const auto flats = gather(_flatTable);
	astap::stacking::set_dark_library(darks);
	astap::stacking::set_flat_library(flats);
}

namespace {

// One check-marked candidate row plus its cached (header-only) metadata.
struct MasterCand {
	QString path;
	astap::stacking::MasterMetadata m;
};

// Average of the group's set temperatures (Pascal temperatureRound = round(avg)).
[[nodiscard]] int average_set_temperature(const std::vector<MasterCand>& group) {
	if (group.empty()) {
		return 0;
	}
	long long sum = 0;
	for (const auto& c : group) {
		sum += c.m.set_temperature;
	}
	return static_cast<int>(std::lrint(static_cast<double>(sum) / group.size()));
}

// Pascal extract_letters_and_numbers_only: keep A–Z, a–z, 0–9 and '-'.
[[nodiscard]] QString sanitize_filter(const std::string& in) {
	QString out;
	for (const char ch : in) {
		if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '-') {
			out.append(QChar::fromLatin1(ch));
		}
	}
	return out;
}

// Gather the check-marked rows of @p table, using the metadata cache (no header
// re-read). Rows whose filename contains @p masterMarker are skipped (already a
// master); pass nullptr to keep every check-marked row. Rows without cached
// metadata are skipped.
[[nodiscard]] std::vector<MasterCand> checked_raws(
		QTableWidget* table, int fileCol,
		const QHash<QString, astap::stacking::MasterMetadata>& cache,
		const char* masterMarker) {
	std::vector<MasterCand> pool;
	for (int r = 0; r < table->rowCount(); ++r) {
		auto* it = table->item(r, fileCol);
		if (!it || it->checkState() != Qt::Checked) {
			continue;
		}
		const QString path = it->data(Qt::UserRole).toString();
		if (masterMarker &&
		    QFileInfo(path).fileName().contains(masterMarker, Qt::CaseInsensitive)) {
			continue;  // already a master file
		}
		if (!cache.contains(path)) {
			continue;
		}
		pool.push_back({path, cache.value(path)});
	}
	return pool;
}

}  // namespace

void StackWindow::replaceDarksByMaster() {
	auto pool = checked_raws(_darkTable, kMFile, _masterMeta, "master_dark");
	if (pool.empty()) {
		QMessageBox::information(this, tr("Master dark"),
			tr("Check-mark one or more raw dark frames to combine."));
		return;
	}

	const bool byExp  = _classDarkExp  && _classDarkExp->isChecked();
	const bool byTemp = _classDarkTemp && _classDarkTemp->isChecked();
	const bool byGain = _classDarkGain && _classDarkGain->isChecked();
	const bool byDate = _classDarkDate && _classDarkDate->isChecked();
	const int  deltaT = _deltaTemp ? _deltaTemp->value() : 3;

	std::vector<QString> consumed;
	std::vector<QString> created;

	// Peel off successive groups matching the first remaining frame's spec, one
	// master per group (Pascal replace_by_master_dark repeat-until loop). Working on
	// a snapshot vector (not the live table) keeps a failed create from re-matching
	// the same rows forever.
	while (!pool.empty()) {
		const astap::stacking::MasterMetadata spec = pool.front().m;
		std::vector<MasterCand> group;
		std::vector<MasterCand> remain;
		for (const auto& cand : pool) {
			const astap::stacking::MasterMetadata& m = cand.m;
			bool match = (m.width == spec.width);   // width always
			if (match && byExp)  { match = (m.exposure == spec.exposure); }
			if (match && byTemp) { match = (std::abs(m.set_temperature - spec.set_temperature) <= deltaT); }
			if (match && byGain) { match = (m.gain == spec.gain); }
			if (match && byDate) { match = (std::abs(m.jd - spec.jd) <= 0.5); }
			if (match) { group.push_back(cand); }
			else       { remain.push_back(cand); }
		}
		pool.swap(remain);
		if (group.empty()) {
			break;   // the spec always matches itself, so this is unreachable
		}

		std::vector<std::filesystem::path> frames;
		frames.reserve(group.size());
		for (const auto& c : group) {
			frames.emplace_back(c.path.toStdString());
		}

		const std::filesystem::path folder =
			std::filesystem::path(group.front().path.toStdString()).parent_path();
		const int tempAvg = average_set_temperature(group);
		const QString name = QString("master_dark_%1x%2s_at_%3C_%4.fit")
			.arg(static_cast<int>(group.size()))
			.arg(spec.exposure)
			.arg(tempAvg)
			.arg(QString::fromStdString(spec.date_obs).left(10));
		const std::filesystem::path out = folder / name.toStdString();

		const auto res = astap::stacking::create_master(
			frames, astap::stacking::MasterKind::Dark, out);
		if (!res.ok) {
			QMessageBox::warning(this, tr("Master dark"),
				tr("Could not create %1: %2")
					.arg(name, QString::fromStdString(res.message)));
			continue;   // group already removed from pool → no infinite loop
		}
		for (const auto& c : group) {
			consumed.push_back(c.path);
		}
		created.push_back(QString::fromStdString(out.string()));
	}

	replaceRows(_darkTable, consumed, created);
	if (!created.empty()) {
		const int nm = static_cast<int>(created.size());
		_phaseLabel->setText(tr("Created %1 master dark%2 from %3 frames")
			.arg(nm)
			.arg(nm == 1 ? QString() : QStringLiteral("s"))
			.arg(static_cast<int>(consumed.size())));
		_phaseLabel->setVisible(true);
	}
}

void StackWindow::replaceFlatsByMaster() {
	// Gather check-marked raw flats, skipping master files and any flat that is
	// already calibrated (non-empty CALSTAT) — Pascal replace_by_master_flat gate.
	std::vector<MasterCand> pool;
	for (auto& cand : checked_raws(_flatTable, kMFile, _masterMeta, "master_flat")) {
		if (cand.m.calstat.empty()) {
			pool.push_back(std::move(cand));
		}
	}
	if (pool.empty()) {
		QMessageBox::information(this, tr("Master flat"),
			tr("Check-mark one or more uncalibrated raw flat frames to combine."));
		return;
	}

	// Check-marked flat-darks available for subtraction (all of them; filtered per
	// group below when classifying by exposure — Pascal average_flatdarks).
	const auto flatDarks = checked_raws(_flatDarkTable, kMFile, _masterMeta, nullptr);

	const bool byFilter   = _classFlatFilter   && _classFlatFilter->isChecked();
	const bool byDuration = _classFlatDuration && _classFlatDuration->isChecked();
	const bool byDate     = _classFlatDate     && _classFlatDate->isChecked();

	std::vector<QString> consumed;
	std::vector<QString> created;
	bool flatDarkUsed = false;

	while (!pool.empty()) {
		const astap::stacking::MasterMetadata spec = pool.front().m;
		const QString specExp = format_exposure(spec.exposure_seconds);
		std::vector<MasterCand> group;
		std::vector<MasterCand> remain;
		for (const auto& cand : pool) {
			const astap::stacking::MasterMetadata& m = cand.m;
			bool match = (m.width == spec.width) && (m.height == spec.height);  // dims always
			if (match && byFilter)   { match = (m.filter_name == spec.filter_name); }
			if (match && byDuration) { match = (format_exposure(m.exposure_seconds) == specExp); }
			if (match && byDate)     { match = (std::abs(m.jd - spec.jd) <= 0.5); }
			if (match) { group.push_back(cand); }
			else       { remain.push_back(cand); }
		}
		pool.swap(remain);
		if (group.empty()) {
			break;
		}

		// Flat-darks paired with this group: when classifying by duration, only
		// those whose exposure matches; otherwise every check-marked flat-dark.
		std::vector<std::filesystem::path> fdFrames;
		for (const auto& fd : flatDarks) {
			if (byDuration && format_exposure(fd.m.exposure_seconds) != specExp) {
				continue;
			}
			fdFrames.emplace_back(fd.path.toStdString());
		}

		std::vector<std::filesystem::path> frames;
		frames.reserve(group.size());
		for (const auto& c : group) {
			frames.emplace_back(c.path.toStdString());
		}

		const std::filesystem::path folder =
			std::filesystem::path(group.front().path.toStdString()).parent_path();
		const QString filter = sanitize_filter(spec.filter_name);
		const QString name =
			QString("master_flat_corrected_with_flat_darks_%1_%2xF_%3xFD_%4sec_%5.fit")
				.arg(filter)
				.arg(static_cast<int>(group.size()))
				.arg(static_cast<int>(fdFrames.size()))
				.arg(specExp)
				.arg(QString::fromStdString(spec.date_obs).left(10));
		const std::filesystem::path out = folder / name.toStdString();

		const auto res = astap::stacking::create_master(
			frames, astap::stacking::MasterKind::Flat, out, fdFrames);
		if (!res.ok) {
			QMessageBox::warning(this, tr("Master flat"),
				tr("Could not create %1: %2")
					.arg(name, QString::fromStdString(res.message)));
			continue;
		}
		if (res.flatdarks_combined > 0) {
			flatDarkUsed = true;
		}
		for (const auto& c : group) {
			consumed.push_back(c.path);
		}
		created.push_back(QString::fromStdString(out.string()));
	}

	replaceRows(_flatTable, consumed, created);
	if (!created.empty()) {
		const int nm = static_cast<int>(created.size());
		_phaseLabel->setText(tr("Created %1 master flat%2 from %3 frames")
			.arg(nm)
			.arg(nm == 1 ? QString() : QStringLiteral("s"))
			.arg(static_cast<int>(consumed.size())));
		_phaseLabel->setVisible(true);
	}

	// Pascal clears the flat-darks list once any were consumed into a master.
	if (flatDarkUsed) {
		_flatDarkTable->setRowCount(0);
	}
}

// Remove the @p consumed raw rows from @p table and append the @p created master
// rows in their place, preserving the table's sort state.
void StackWindow::replaceRows(QTableWidget* table,
		const std::vector<QString>& consumed,
		const std::vector<QString>& created) {
	if (consumed.empty() && created.empty()) {
		return;
	}
	const bool wasSorting = table->isSortingEnabled();
	table->setSortingEnabled(false);
	const QSet<QString> gone(consumed.begin(), consumed.end());
	for (int r = table->rowCount() - 1; r >= 0; --r) {
		auto* it = table->item(r, kMFile);
		if (it && gone.contains(it->data(Qt::UserRole).toString())) {
			table->removeRow(r);
		}
	}
	for (const auto& p : created) {
		insertMasterRow(table, p);
	}
	table->setSortingEnabled(wasSorting);
	refreshCompatibility();
}

void StackWindow::addDarkFiles()     { addMasterFiles(_darkTable, tr("dark")); }
void StackWindow::addFlatFiles()     { addMasterFiles(_flatTable, tr("flat")); }
void StackWindow::addFlatDarkFiles() { addMasterFiles(_flatDarkTable, tr("flat dark")); }
void StackWindow::removeDarks()      { removeSelectedRows(_darkTable); }
void StackWindow::removeFlats()      { removeSelectedRows(_flatTable); }
void StackWindow::removeFlatDarks()  { removeSelectedRows(_flatDarkTable); }
void StackWindow::clearDarks()       { _darkTable->setRowCount(0); refreshCompatibility(); }
void StackWindow::clearFlats()       { _flatTable->setRowCount(0); refreshCompatibility(); }
void StackWindow::clearFlatDarks()   { _flatDarkTable->setRowCount(0); }

void StackWindow::buildLightsTab() {
	auto* page = new QWidget(_tabs);
	auto* layout = new QVBoxLayout(page);

	_fileTable = new QTableWidget(0, 2, page);
	_fileTable->setHorizontalHeaderLabels({tr("File"), tr("Channel")});
	_fileTable->horizontalHeader()->setSectionResizeMode(
		0, QHeaderView::Stretch);
	_fileTable->horizontalHeader()->setSectionResizeMode(
		1, QHeaderView::ResizeToContents);
	_fileTable->verticalHeader()->setVisible(false);
	_fileTable->setWordWrap(false);
	_fileTable->setSelectionBehavior(QAbstractItemView::SelectRows);
	_fileTable->setSelectionMode(QAbstractItemView::ExtendedSelection);
	_fileTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
	layout->addWidget(_fileTable, 1);

	auto* buttonRow = new QHBoxLayout();
	_addButton = new QPushButton(tr("Add…"), page);
	_removeButton = new QPushButton(tr("Remove"), page);
	_clearButton = new QPushButton(tr("Clear"), page);
	buttonRow->addWidget(_addButton);
	buttonRow->addWidget(_removeButton);
	buttonRow->addWidget(_clearButton);
	buttonRow->addStretch(1);
	layout->addLayout(buttonRow);

	connect(_addButton, &QPushButton::clicked, this, &StackWindow::addFiles);
	connect(_removeButton, &QPushButton::clicked, this, &StackWindow::removeSelected);
	connect(_clearButton, &QPushButton::clicked, this, &StackWindow::clearList);

	_tabs->addTab(page, tr("Lights"));
}

void StackWindow::buildDarksTab() {
	auto* page = new QWidget(_tabs);
	auto* layout = new QVBoxLayout(page);
	_darkTable = make_master_table(page);
	_darkTable->setColumnHidden(kMFilter, true);   // darks have no filter (Pascal listview2)
	_darkTable->setColumnHidden(kMCalstat, true);  // …nor a calibration status
	layout->addWidget(_darkTable, 1);

	auto* row = new QHBoxLayout();
	auto* add = new QPushButton(tr("Add…"), page);
	auto* rem = new QPushButton(tr("Remove"), page);
	auto* clr = new QPushButton(tr("Clear"), page);
	row->addWidget(add);
	row->addWidget(rem);
	row->addWidget(clr);
	row->addStretch(1);
	layout->addLayout(row);

	// Master-creation row (Pascal replace_by_master_dark1 + classify_dark_date1).
	auto* makeRow = new QHBoxLayout();
	auto* replace = new QPushButton(tr("Replace check-marked by one or more master dark"), page);
	replace->setToolTip(tr("Combine the check-marked raw darks into one or more master "
	                       "darks and swap them into this list."));
	_classDarkDate = new QCheckBox(tr("Classify by date for master creation"), page);
	_classDarkDate->setToolTip(tr("Group frames taken within 12 hours into separate masters."));
	makeRow->addWidget(replace);
	makeRow->addWidget(_classDarkDate);
	makeRow->addStretch(1);
	layout->addLayout(makeRow);

	connect(add, &QPushButton::clicked, this, &StackWindow::addDarkFiles);
	connect(rem, &QPushButton::clicked, this, &StackWindow::removeDarks);
	connect(clr, &QPushButton::clicked, this, &StackWindow::clearDarks);
	connect(replace, &QPushButton::clicked, this, &StackWindow::replaceDarksByMaster);

	_tabs->addTab(page, tr("Darks"));
}

void StackWindow::buildFlatsTab() {
	auto* page = new QWidget(_tabs);
	auto* layout = new QVBoxLayout(page);
	_flatTable = make_master_table(page);
	layout->addWidget(_flatTable, 1);

	auto* row = new QHBoxLayout();
	auto* add = new QPushButton(tr("Add…"), page);
	auto* rem = new QPushButton(tr("Remove"), page);
	auto* clr = new QPushButton(tr("Clear"), page);
	row->addWidget(add);
	row->addWidget(rem);
	row->addWidget(clr);
	row->addStretch(1);
	layout->addLayout(row);

	// Master-creation row (Pascal replace_by_master_flat1 + classify_flat_duration1 /
	// classify_flat_date1). The check-marked frames on the Flat darks tab are averaged
	// and subtracted per the formula below.
	auto* makeRow = new QHBoxLayout();
	auto* replace = new QPushButton(
		tr("Replace check-marked by master flat (flat darks included)"), page);
	replace->setToolTip(tr("Combine the check-marked raw flats into one or more master "
	                       "flats, subtracting the check-marked flat darks."));
	_classFlatDuration =
		new QCheckBox(tr("Classify by exposure duration for master creation"), page);
	_classFlatDuration->setToolTip(
		tr("Group flats (and flat darks) by exposure duration into separate masters."));
	_classFlatDate = new QCheckBox(tr("Classify by date for master creation"), page);
	_classFlatDate->setToolTip(tr("Group flats taken within 12 hours into separate masters."));
	makeRow->addWidget(replace);
	// Stack the two flat classify checkboxes vertically so the long button + both
	// captions do not fight for width in the 520-px window.
	auto* flatClassifyCol = new QVBoxLayout();
	flatClassifyCol->addWidget(_classFlatDuration);
	flatClassifyCol->addWidget(_classFlatDate);
	makeRow->addLayout(flatClassifyCol);
	makeRow->addStretch(1);
	layout->addLayout(makeRow);

	// Pascal formula label (unit_stack.lfm: 'Master flat := (1/n ∑ [flats] - 1/n ∑ [flat darks])').
	auto* formula = new QLabel(
		tr("Master flat  :=  (1/n ∑ [flats]  −  1/n ∑ [flat darks])"), page);
	formula->setStyleSheet("color: gray;");
	layout->addWidget(formula);

	connect(add, &QPushButton::clicked, this, &StackWindow::addFlatFiles);
	connect(rem, &QPushButton::clicked, this, &StackWindow::removeFlats);
	connect(clr, &QPushButton::clicked, this, &StackWindow::clearFlats);
	connect(replace, &QPushButton::clicked, this, &StackWindow::replaceFlatsByMaster);

	_tabs->addTab(page, tr("Flats"));
}

void StackWindow::buildFlatDarksTab() {
	auto* page = new QWidget(_tabs);
	auto* layout = new QVBoxLayout(page);
	_flatDarkTable = make_master_table(page);
	// Raw flat-dark/bias frames only — no filter, no calibration status, and they
	// are never matched against a light so no compatibility column either.
	_flatDarkTable->setColumnHidden(kMFilter, true);
	_flatDarkTable->setColumnHidden(kMCalstat, true);
	_flatDarkTable->setColumnHidden(kMCompat, true);
	layout->addWidget(_flatDarkTable, 1);

	auto* row = new QHBoxLayout();
	auto* add = new QPushButton(tr("Add…"), page);
	auto* rem = new QPushButton(tr("Remove"), page);
	auto* clr = new QPushButton(tr("Clear"), page);
	row->addWidget(add);
	row->addWidget(rem);
	row->addWidget(clr);
	row->addStretch(1);
	layout->addLayout(row);

	connect(add, &QPushButton::clicked, this, &StackWindow::addFlatDarkFiles);
	connect(rem, &QPushButton::clicked, this, &StackWindow::removeFlatDarks);
	connect(clr, &QPushButton::clicked, this, &StackWindow::clearFlatDarks);

	_tabs->addTab(page, tr("Flat darks"));
}

QWidget* StackWindow::buildClassifyBar() {
	auto* group = new QGroupBox(tr("Classify by"), this);
	auto* row = new QHBoxLayout(group);

	// Light classification (Pascal classify_filter_light1 / classify_object1).
	_lightFilter = new QCheckBox(tr("Light filter"), group);
	_lightFilter->setToolTip(
		tr("Combine as LRGB/RGB using each light's channel tag (Lights tab)."));
	_lightObject = new QCheckBox(tr("Light object"), group);
	_lightObject->setToolTip(
		tr("Stack images of different objects in one run "
		   "(multi-object stacking is not yet implemented)."));
	_lightObject->setEnabled(false);   // no per-object classification in the engine yet

	_classDarkExp    = new QCheckBox(tr("Dark exposure duration"), group);
	_classDarkTemp   = new QCheckBox(tr("Dark temperature"), group);
	_classDarkGain   = new QCheckBox(tr("Dark gain"), group);
	_classFlatFilter = new QCheckBox(tr("Flat filter"), group);

	_deltaTemp = new QSpinBox(group);
	_deltaTemp->setRange(1, 5);          // Pascal delta_temp_updown1 (Min 1, Max 5)
	_deltaTemp->setValue(3);             // default 3
	_deltaTemp->setPrefix(tr("± "));
	_deltaTemp->setSuffix(tr(" °C"));
	_deltaTemp->setToolTip(tr("Temperature tolerance for dark matching"));

	row->addWidget(_lightFilter);
	row->addWidget(_lightObject);
	row->addWidget(_classDarkExp);
	row->addWidget(_classDarkTemp);
	row->addWidget(_deltaTemp);
	row->addWidget(_classDarkGain);
	row->addWidget(_classFlatFilter);
	row->addStretch(1);

	// The Compatibility column depends on these — recompute on any change.
	connect(_classDarkExp,    &QCheckBox::toggled, this, &StackWindow::refreshCompatibility);
	connect(_classDarkTemp,   &QCheckBox::toggled, this, &StackWindow::refreshCompatibility);
	connect(_classDarkGain,   &QCheckBox::toggled, this, &StackWindow::refreshCompatibility);
	connect(_classFlatFilter, &QCheckBox::toggled, this, &StackWindow::refreshCompatibility);
	connect(_deltaTemp, qOverload<int>(&QSpinBox::valueChanged),
	        this, &StackWindow::refreshCompatibility);

	return group;
}

void StackWindow::buildSettingsTab() {
	auto* page = new QWidget(_tabs);
	auto* form = new QFormLayout(page);

	// Combine algorithm. LRGB is not a method here — it is activated by the
	// "Light filter" checkbox (Pascal: method = Average/Sigma-clip, and LRGB
	// classification is switched on separately by classify_filter_light1).
	_methodCombo = new QComboBox(page);
	_methodCombo->addItem(tr("Average (weighted)"), kMethodAverage);
	_methodCombo->addItem(tr("Sigma-clip average"), kMethodSigmaClip);
	form->addRow(tr("Method:"), _methodCombo);

	_alignmentCombo = new QComboBox(page);
	_alignmentCombo->addItem(tr("Star-match (quads)"), kAlignStarMatch);
	_alignmentCombo->addItem(tr("Astrometric (per-frame solve)"), kAlignAstrometric);
	_alignmentCombo->addItem(tr("Manual reference star"), kAlignManual);
	_alignmentCombo->addItem(tr("None (pre-aligned)"), kAlignNone);
	_alignmentCombo->addItem(tr("WCS (use each header's solution)"), kAlignWcs);
	form->addRow(tr("Alignment:"), _alignmentCombo);

	_sigmaFactor = new QDoubleSpinBox(page);
	_sigmaFactor->setRange(0.5, 10.0);
	_sigmaFactor->setSingleStep(0.1);
	_sigmaFactor->setValue(2.0);
	_sigmaFactor->setDecimals(1);
	_sigmaFactor->setSuffix(tr(" σ"));
	form->addRow(tr("Sigma-clip threshold:"), _sigmaFactor);

	_maxStars = new QSpinBox(page);
	_maxStars->setRange(50, 5000);
	_maxStars->setSingleStep(50);
	_maxStars->setValue(500);
	form->addRow(tr("Max stars for matching:"), _maxStars);

	// LRGB colour-mix matrix (Pascal rr1..bb1 calibration text boxes). Rows are the
	// input channel, columns the output channel; identity = a straight RGB combine.
	// Applied only when LRGB is active (the Light-filter checkbox).
	auto* mixGroup = new QGroupBox(tr("Colour-mix matrix (LRGB)"), page);
	mixGroup->setToolTip(
		tr("Weight of each colour input channel in each output channel. The identity "
		   "(diagonal 1) is a straight RGB combine; e.g. set R-in→B and B-in→R to swap "
		   "red and blue. Weights of zero or below are ignored (as in the original). "
		   "Applies when Light filter (LRGB) is enabled."));
	auto* grid = new QGridLayout(mixGroup);
	grid->addWidget(new QLabel(tr("out →"), mixGroup), 0, 0);
	const char* chan[3] = {"R", "G", "B"};
	for (int j = 0; j < 3; ++j) {
		auto* h = new QLabel(QString::fromLatin1(chan[j]), mixGroup);
		h->setAlignment(Qt::AlignCenter);
		grid->addWidget(h, 0, j + 1);
	}
	for (int i = 0; i < 3; ++i) {
		grid->addWidget(new QLabel(tr("%1 in").arg(QString::fromLatin1(chan[i])), mixGroup), i + 1, 0);
		for (int j = 0; j < 3; ++j) {
			auto* sp = new QDoubleSpinBox(mixGroup);
			sp->setRange(-9.999, 9.999);
			sp->setSingleStep(0.1);
			sp->setDecimals(3);
			sp->setValue(i == j ? 1.0 : 0.0);   // identity
			_colorMix[i][j] = sp;
			grid->addWidget(sp, i + 1, j + 1);
		}
	}
	form->addRow(mixGroup);

	_tabs->addTab(page, tr("Settings"));
}

void StackWindow::applySettingsToEngine() {
	const auto align = _alignmentCombo->currentData().toInt();
	astap::use_manual_align        = (align == kAlignManual);
	astap::use_ephemeris_alignment = false;  // no ephemeris UI yet
	astap::use_astrometry_internal = (align == kAlignAstrometric);
	astap::skip_alignment          = (align == kAlignNone);
	astap::use_wcs_alignment       = (align == kAlignWcs);

	astap::sigma_clip_factor = _sigmaFactor->value();
	astap::max_stars_setting = _maxStars->value();

	// LRGB colour-mix matrix (rows = input R/G/B, cols = output R/G/B). Harmless for
	// the non-LRGB methods, which never read it.
	astap::lrgb_colour_mix = {
		_colorMix[0][0]->value(), _colorMix[0][1]->value(), _colorMix[0][2]->value(),
		_colorMix[1][0]->value(), _colorMix[1][1]->value(), _colorMix[1][2]->value(),
		_colorMix[2][0]->value(), _colorMix[2][1]->value(), _colorMix[2][2]->value(),
	};
}

void StackWindow::addFiles() {
	QSettings settings;
	const auto lastDir = settings.value("files/lastStackDir").toString();

	const auto paths = QFileDialog::getOpenFileNames(
		this,
		tr("Add light frames"),
		lastDir,
		tr("FITS images (*.fit *.fits *.fts *.new);;"
		   "All images (*.fit *.fits *.fts *.new "
		               "*.png *.jpg *.jpeg *.bmp *.tif *.tiff);;"
		   "All files (*)"));
	if (paths.isEmpty()) {
		return;
	}

	settings.setValue("files/lastStackDir",
		QFileInfo(paths.first()).absolutePath());

	// Find already-present paths to skip duplicates. The full path is stashed
	// in Qt::UserRole so the visible cell text can be a short basename.
	auto existing = QSet<QString>{};
	for (int r = 0; r < _fileTable->rowCount(); ++r) {
		existing.insert(_fileTable->item(r, 0)->data(Qt::UserRole).toString());
	}

	const auto labels = channelLabels();
	for (const auto& p : paths) {
		if (existing.contains(p)) {
			continue;
		}
		const auto row = _fileTable->rowCount();
		_fileTable->insertRow(row);
		auto* pathItem = new QTableWidgetItem(QFileInfo(p).fileName());
		pathItem->setData(Qt::UserRole, p);
		pathItem->setToolTip(p);
		_fileTable->setItem(row, 0, pathItem);

		auto* combo = new QComboBox(_fileTable);
		combo->addItems(labels);
		combo->setCurrentIndex(guess_channel(p));
		_fileTable->setCellWidget(row, 1, combo);
	}
	refreshCompatibility();  // the first light is the compatibility reference
}

void StackWindow::removeSelected() {
	auto rows = QList<int>{};
	for (const auto* item : _fileTable->selectedItems()) {
		if (item->column() == 0) {
			rows.append(item->row());
		}
	}
	std::sort(rows.begin(), rows.end(), std::greater<int>());
	for (const auto r : rows) {
		_fileTable->removeRow(r);
	}
	refreshCompatibility();
}

void StackWindow::clearList() {
	_fileTable->setRowCount(0);
	refreshCompatibility();
}

namespace {

struct Row { QString path; int channel; long long area; };

/// @brief Return the friendly name for a channel enum value.
[[nodiscard]] const char* channel_name(int c) {
	switch (c) {
		case kChanL:   return "L";
		case kChanR:   return "R";
		case kChanG:   return "G";
		case kChanB:   return "B";
		case kChanRGB: return "RGB";
		default:       return "?";
	}
}

/// @brief Run a plain Average stack on the supplied rows and save the result
///        to @p outPath as FITS.
/// @details Uses stack_average on a file list built from @p files, largest-
///          area first (so the engine's reference frame matches the canonical
///          behaviour). Returns the number of frames combined on success, or
///          0 if fewer than 2 frames could be stacked.
[[nodiscard]] int run_average_to_fits(const std::vector<Row>& files,
                                       const QString& outPath,
                                       int osc) {
	// Largest first so stack_average's internal ref picks the biggest frame.
	auto ordered = files;
	std::stable_sort(ordered.begin(), ordered.end(),
		[](const Row& a, const Row& b) { return a.area > b.area; });

	auto todo = std::vector<astap::FileToDo>{};
	todo.reserve(ordered.size());
	for (int i = 0; i < static_cast<int>(ordered.size()); ++i) {
		todo.push_back({ordered[i].path.toStdString(), i});
	}

	auto counter = 0;
	astap::stacking::stack_average(osc, std::span<astap::FileToDo>(todo), counter);

	if (counter < 2 || astap::img_loaded.empty()) {
		return 0;
	}

	// Write the pre-stacked master to the temp FITS. memo1_lines already
	// carries a valid header (copied from the reference frame during the
	// stack), so save_fits can drop it in as-is.
	auto memo = astap::memo1_lines;
	const auto ok = astap::core::save_fits(
		astap::img_loaded,
		memo,
		std::filesystem::path(outPath.toStdString()),
		/*type1=*/-32,   // float-32, preserves the combined pedestal+noise floor
		/*override2=*/true);
	return ok ? counter : 0;
}

}  // namespace

void StackWindow::startStack() {
	const auto count = _fileTable->rowCount();
	if (count < 1) {
		QMessageBox::information(this, tr("Stack"),
			tr("Add at least one frame to stack."));
		return;
	}

	applySettingsToEngine();
	applyLibrariesToEngine();   // install the checked dark/flat candidate libraries

	_stackButton->setEnabled(false);
	_progress->setValue(0);
	_phaseLabel->setVisible(false);

	// Gather (path, channel, area) rows from the table.
	auto rows = std::vector<Row>{};
	rows.reserve(count);
	for (int i = 0; i < count; ++i) {
		const auto path = _fileTable->item(i, 0)->data(Qt::UserRole).toString();
		auto* combo = qobject_cast<QComboBox*>(
			_fileTable->cellWidget(i, 1));
		const auto chan = combo ? combo->currentIndex() : kChanLight;
		rows.push_back({path, chan, probe_area(path)});
	}

	astap::stacking::set_progress_sink(
		[this](double value, const std::string& /*label*/) {
			_progress->setValue(static_cast<int>(std::round(value)));
			QCoreApplication::processEvents();
		});
	astap::stacking::set_memo2_sink([](const std::string& msg) {
		qDebug().noquote() << "[stack]" << QString::fromStdString(msg);
	});

	// Centralised cleanup: run on every exit path so progress sinks don't
	// leak into later runs.
	auto tear_down = [this]() {
		astap::stacking::set_progress_sink(nullptr);
		astap::stacking::set_memo2_sink(nullptr);
		_stackButton->setEnabled(true);
		_phaseLabel->setVisible(false);
	};

	auto counter = 0;
	const auto method = _methodCombo->currentData().toInt();
	const auto osc = 0;  // TODO: OSC/Bayer toggle
	// LRGB is activated by the "Light filter" classify checkbox (Pascal
	// classify_filter_light1), not by a Method entry.
	const bool lrgb = _lightFilter && _lightFilter->isChecked();

	if (lrgb) {
		// Group rows by channel. Skip untagged ("Light") rows in LRGB mode —
		// those files aren't meant for the combine.
		auto byChannel = std::map<int, std::vector<Row>>{};
		for (const auto& r : rows) {
			if (r.channel != kChanLight) {
				byChannel[r.channel].push_back(r);
			}
		}

		auto missing = QStringList{};
		if (!byChannel.count(kChanR)) missing << "R";
		if (!byChannel.count(kChanG)) missing << "G";
		if (!byChannel.count(kChanB)) missing << "B";
		if (!missing.isEmpty()) {
			QMessageBox::warning(this, tr("LRGB"),
				tr("LRGB combine needs at least one file tagged for each of "
				   "R, G, B. Missing: %1.").arg(missing.join(", ")));
			tear_down();
			return;
		}

		// Count channels that need pre-stacking so the phase label is
		// meaningful. One phase per multi-file channel, plus one for the
		// final combine.
		auto autoChainPhases = 0;
		for (const auto& [ch, files] : byChannel) {
			if (files.size() > 1) {
				++autoChainPhases;
			}
		}

		// Resolve each channel to a single file path — either the lone tagged
		// file, or a pre-stacked master written to the temp dir.
		if (autoChainPhases > 0 && !_tempDir) {
			_tempDir = std::make_unique<QTemporaryDir>();
			if (!_tempDir->isValid()) {
				QMessageBox::warning(this, tr("LRGB"),
					tr("Could not create a temporary directory for channel "
					   "masters: %1").arg(_tempDir->errorString()));
				tear_down();
				return;
			}
		}

		auto channelMaster = std::map<int, QString>{};
		auto currentPhase = 0;
		for (const auto& [ch, files] : byChannel) {
			if (files.size() == 1) {
				channelMaster[ch] = files[0].path;
				continue;
			}

			++currentPhase;
			_phaseLabel->setVisible(true);
			_phaseLabel->setText(
				tr("Stacking channel %1 (%2 of %3) — %4 frames")
					.arg(channel_name(ch))
					.arg(currentPhase)
					.arg(autoChainPhases + 1)        // +1 for the final combine
					.arg(files.size()));
			_progress->setValue(0);
			QCoreApplication::processEvents();

			const auto masterPath =
				_tempDir->filePath(QString("master_%1.fits").arg(channel_name(ch)));

			const auto n = run_average_to_fits(files, masterPath, osc);
			if (n == 0) {
				QMessageBox::warning(this, tr("LRGB auto-chain"),
					tr("Could not pre-stack the %1 channel (%2 frames).")
						.arg(channel_name(ch))
						.arg(files.size()));
				tear_down();
				return;
			}
			channelMaster[ch] = masterPath;
		}

		// Final phase: LRGB combine. Re-probe the master paths for area so the
		// reference-picking logic handles cross-channel size differences.
		if (autoChainPhases > 0) {
			_phaseLabel->setText(
				tr("Combining LRGB (%1 of %1)").arg(autoChainPhases + 1));
			_progress->setValue(0);
			QCoreApplication::processEvents();
		}

		auto pathFor = [&](int ch) -> QString {
			auto it = channelMaster.find(ch);
			return (it != channelMaster.end()) ? it->second : QString{};
		};

		const auto lPath   = pathFor(kChanL);
		const auto rPath   = pathFor(kChanR);
		const auto gPath   = pathFor(kChanG);
		const auto bPath   = pathFor(kChanB);
		const auto rgbPath = pathFor(kChanRGB);

		auto areaOf = [&](const QString& p) -> long long {
			if (p.isEmpty()) return 0;
			return probe_area(p);
		};

		auto refPath = !lPath.isEmpty() ? lPath : rPath;
		auto refArea = areaOf(refPath);
		for (const auto* p : {&rPath, &gPath, &bPath, &lPath}) {
			if (!p->isEmpty() && areaOf(*p) > refArea) {
				refPath = *p;
				refArea = areaOf(*p);
			}
		}

		auto files = std::vector<astap::FileToDo>{};
		files.push_back({refPath.toStdString(), 0});
		files.push_back({rPath.toStdString(), 0});
		files.push_back({gPath.toStdString(), 0});
		files.push_back({bPath.toStdString(), 0});
		files.push_back({rgbPath.toStdString(), 0});  // may be empty
		files.push_back({lPath.toStdString(), 0});    // may be empty

		astap::stacking::stack_LRGB(std::span<astap::FileToDo>(files), counter);
	} else {
		// Average / Sigma-clip: feed everything untagged by channel. Put
		// the largest-dim frame first so it becomes the engine's reference.
		auto ordered = rows;
		std::stable_sort(ordered.begin(), ordered.end(),
			[](const Row& a, const Row& b) { return a.area > b.area; });
		auto files = std::vector<astap::FileToDo>{};
		files.reserve(ordered.size());
		for (int i = 0; i < static_cast<int>(ordered.size()); ++i) {
			files.push_back({ordered[i].path.toStdString(), i});
		}
		const auto span = std::span<astap::FileToDo>(files);
		if (method == kMethodSigmaClip) {
			astap::stacking::stack_sigmaclip(osc, span, counter);
		} else {
			astap::stacking::stack_average(osc, span, counter);
		}
	}

	astap::stacking::set_progress_sink(nullptr);
	astap::stacking::set_memo2_sink(nullptr);

	if (counter < 2 || astap::img_loaded.empty()) {
		QMessageBox::warning(this, tr("Stack"),
			tr("Could not stack enough frames (%1 of %2).")
				.arg(counter).arg(count));
		_progress->setValue(0);
		_stackButton->setEnabled(true);
		_phaseLabel->setVisible(false);
		return;
	}

	astap::head.light_count = counter;
	astap::filename2 = "Stacked_" + std::to_string(counter) + "_frames";

	if (_viewer) {
		_viewer->setImage(astap::img_loaded, astap::head);
	}

	_progress->setValue(100);
	_stackButton->setEnabled(true);
	_phaseLabel->setVisible(false);

	emit stackCompleted(counter);
}

} // namespace astap::gui
