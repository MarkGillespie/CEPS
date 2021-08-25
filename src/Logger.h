#pragma once

#include <fstream>
#include <iomanip> // setw
#include <iostream>
#include <limits>
#include <tuple>
#include <vector>

enum class LogType { STRING, DOUBLE };
std::string to_string(LogType lot);

std::ostream& operator<<(std::ostream& output, const LogType& lo);

class Logger {
  public:
    Logger();

    // use templates to prevent casting of input types
    template <typename T>
    void log(std::string name, T val);

    bool writeLog(std::string filename);
    void writeLog(std::ostream& out);

  protected:
    std::vector<std::tuple<std::string, LogType, size_t>> logs;
    std::vector<std::string> stringLogs;
    std::vector<double> doubleLogs;
};

// explicit specializations
// clang-format off
template <> void Logger::log(std::string name, const char* val);
template <> void Logger::log(std::string name, std::string val);
template <> void Logger::log(std::string name, bool val);
template <> void Logger::log(std::string name, double val);
template <> void Logger::log(std::string name, int val);
template <> void Logger::log(std::string name, size_t val);
//clang-format on
