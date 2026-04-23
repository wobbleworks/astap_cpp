///----------------------------------------
///      @file qt_http_client.cpp
///   @ingroup ASTAP++
///     @brief Implementation of @ref QtHttpClient.
///    @author Created by John Stephen on 4/23/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "qt_http_client.h"

#include <QByteArray>
#include <QEventLoop>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QString>
#include <QTimer>
#include <QUrl>

#include <memory>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

QtHttpClient::QtHttpClient(std::chrono::milliseconds timeout)
	: _timeout(timeout) {
}

std::expected<std::string, std::string> QtHttpClient::get(const std::string& url) {
	const auto qurl = QUrl(QString::fromStdString(url));
	if (!qurl.isValid()) {
		return std::unexpected(std::string("Invalid URL: ") + url);
	}

	QNetworkRequest request(qurl);
	request.setAttribute(QNetworkRequest::RedirectPolicyAttribute,
		QNetworkRequest::NoLessSafeRedirectPolicy);
	request.setRawHeader("User-Agent", "ASTAP++/1.0");

	auto reply = std::unique_ptr<QNetworkReply>(_nam.get(request));
	if (!reply) {
		return std::unexpected(std::string("Failed to dispatch request"));
	}

	QEventLoop loop;
	QObject::connect(reply.get(), &QNetworkReply::finished, &loop, &QEventLoop::quit);

	QTimer timeoutTimer;
	timeoutTimer.setSingleShot(true);
	auto timedOut = false;
	QObject::connect(&timeoutTimer, &QTimer::timeout, &loop, [&]() {
		timedOut = true;
		reply->abort();
		loop.quit();
	});
	timeoutTimer.start(static_cast<int>(_timeout.count()));

	loop.exec();

	if (timedOut) {
		return std::unexpected(std::string("Request timed out: ") + url);
	}

	if (reply->error() != QNetworkReply::NoError) {
		return std::unexpected(reply->errorString().toStdString());
	}

	const auto body = reply->readAll();
	return std::string(body.constData(), static_cast<std::size_t>(body.size()));
}

} // namespace astap::gui
