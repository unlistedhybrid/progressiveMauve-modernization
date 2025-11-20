/*******************************************************************************
 * $Id: Files.h,v 1.23 2004/04/19 23:10:13 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 ******************************************************************************/

#ifndef __libMems_Files_h__
#define __libMems_Files_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <climits>
#endif

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdlib> // getenv
#include <cstdio>  // close
#include <mutex>

#include "boost/filesystem.hpp"
#include "boost/algorithm/string.hpp"

// --- Thread-safe file deletion registration ---
inline std::vector<std::string>& registerFileToDelete(const std::string& fname = "") {
    static std::vector<std::string>* files = new std::vector<std::string>();
    static std::mutex files_mutex;
    if (!fname.empty()) {
        std::lock_guard<std::mutex> lock(files_mutex);
        files->push_back(fname);
    }
    return *files;
}

inline void deleteRegisteredFiles() {
    std::vector<std::string>& del_files = registerFileToDelete();
    for (size_t fileI = 0; fileI < del_files.size(); ++fileI) {
        boost::system::error_code ec;
        boost::filesystem::remove(del_files[fileI], ec); // ignore errors
    }
    del_files.clear();
}

// --- Cross-platform temporary file name creation (does NOT create file) ---
inline std::string CreateTempFileName(const std::string& prefix) {
    std::string dir, name, ret_path;
    boost::filesystem::path path(prefix);

#if BOOST_FILESYSTEM_VERSION >= 3
    dir = path.parent_path().string();
    name = path.filename().string();
#else
    dir = path.branch_path().string();
    name = path.leaf();
#endif

    if (name == "/") {
        dir += name;
        name.clear();
    }

#ifdef _WIN32
    char buf[MAX_PATH + 1] = {0};

    if (dir.empty()) {
        DWORD ret = GetTempPathA(MAX_PATH, buf);
        if (!ret) {
            std::cerr << "GetTempPath failed\n";
            dir = ".";
        } else {
            dir = buf;
        }
    } else {
        boost::algorithm::replace_all(dir, "/", "\\");
    }

    UINT uRet = GetTempFileNameA(dir.c_str(), name.empty() ? "tmp" : name.c_str(), 0, buf);
    if (!uRet) {
        std::cerr << "GetTempFileName failed\n";
        return std::string();
    }
    ret_path = buf;

#else // POSIX
    if (dir.empty()) {
        const char* env_val = getenv("TMP");
        dir = env_val ? env_val : "";
        if (dir.empty()) {
            env_val = getenv("TMPDIR");
            dir = env_val ? env_val : "";
        }
        if (dir.empty()) {
            env_val = getenv("TEMP");
            dir = env_val ? env_val : "";
        }
        if (dir.empty()) {
#ifdef __DOS__
            dir = ".";
#else
            dir = "/tmp";
#endif
        }
    }

    // Compose a path for mkstemp
    boost::filesystem::path tmppath = boost::filesystem::path(dir) / (name.empty() ? "file" : name);
    std::string path_str = tmppath.string() + "XXXXXX";
    char buf[PATH_MAX + 1];
    std::snprintf(buf, sizeof(buf), "%s", path_str.c_str());
    buf[sizeof(buf) - 1] = '\0';

    // Use mkstemp if available, with safe fallback for unique name generation
    int fdTemp = mkstemp(buf);
    if (fdTemp == -1) {
        // mkstemp failed, fall back to safe unique name generation
        unsigned my_pid = 0;
#ifndef __DOS__
        my_pid = getpid();
#endif
        static const size_t numTries = 10000;
        for (size_t n = 0; n < numTries; n++) {
            std::ostringstream oss;
            oss << tmppath.string() << my_pid << "." << std::setfill('0') << std::setw(5) << n;
            boost::filesystem::path try_path(oss.str());
            if (!boost::filesystem::exists(try_path)) {
                ret_path = oss.str();
                break;
            }
        }
        if (ret_path.empty()) {
            // All attempts failed
            return std::string();
        }
    } else {
        ret_path = buf;
        close(fdTemp);
    }

#endif // _WIN32
    return ret_path;
}

#endif // __libMems_Files_h__
