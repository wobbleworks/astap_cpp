///----------------------------------------
///      @file qt_http_client.h
///   @ingroup ASTAP++
///     @brief Qt-backed implementation of @ref astap::IHttpClient.
///   @details Wraps @c QNetworkAccessManager with a blocking @c get(...) for
///            use by the online-catalog modules (VSP/VSX, Simbad, Vizier).
///            Construct one per worker thread; do NOT share across threads.
///    @author Created by John Stephen on 4/23/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include "../../src/types.h"

#include <QNetworkAccessManager>

#include <chrono>
#include <expected>
#include <string>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @class QtHttpClient
/// @brief Synchronous HTTP GET via QNetworkAccessManager + QEventLoop.
/// @details Each instance owns its own @c QNetworkAccessManager. The
///          blocking @c get(...) is intended to run from a worker thread
///          (e.g. inside @c QtConcurrent::run); calling it from the GUI
///          thread will spin a nested event loop and freeze input until
///          the request completes or times out.
///----------------------------------------

class QtHttpClient final : public astap::IHttpClient {
public:
	explicit QtHttpClient(std::chrono::milliseconds timeout = std::chrono::seconds{30});
	~QtHttpClient() override = default;

	std::expected<std::string, std::string> get(const std::string& url) override;

private:
	QNetworkAccessManager _nam;
	std::chrono::milliseconds _timeout;
};

} // namespace astap::gui
