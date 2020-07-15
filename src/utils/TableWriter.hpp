#pragma once

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <stddef.h>
#include <string>
#include <utility>
#include <vector>

struct Column {
  std::string name;
  int         width;
  int         precision = 6;

  explicit Column(std::string const &name);

  Column(std::string const &name, int width);

  Column(std::string const &name, int width, int precision);
};

class Table {
public:
  std::vector<Column> cols;
  std::string         sepChar = "|";
  char                padding = ' ';
  std::ostream &      out     = std::cout;

  Table();

  Table(std::ostream &out);

  /// Adds a column of given name, width and float precision
  template <class... T>
  void addColumn(T &&... arg)
  {
    cols.emplace_back(std::forward<T>(arg)...);
  }

  /// Prints the formatted header
  void printHeader();

  /// Prints a line, accepting arbitrary arguments
  template <class... Ts>
  void printRow(Ts... args)
  {
    printRow(static_cast<size_t>(0), args...);
  }

  /// Prints a ostream convertible type
  template <class T, class... Ts>
  void printRow(size_t index, T a, Ts... args)
  {
    out << padding << std::setw(cols[index].width) << std::setprecision(cols[index].precision)
        << a << padding << sepChar;
    printRow(index + 1, args...);
  }

  /// Prints a duration as milliseconds
  template <class Rep, class Period, class... Ts>
  void printRow(size_t index, std::chrono::duration<Rep, Period> duration, Ts... args)
  {
    double ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
    out << padding << std::setw(cols[index].width) << std::setprecision(cols[index].precision)
        << ms << padding << sepChar;
    printRow(index + 1, args...);
  }

  /// Recursion anchor, prints the last entry and the endl
  template <class T>
  void printRow(size_t index, T a)
  {
    out << padding << std::setw(cols[index].width) << std::setprecision(cols[index].precision)
        << a << padding << sepChar << '\n';
  }
};
