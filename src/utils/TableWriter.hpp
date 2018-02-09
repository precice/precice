#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <string>
#include <vector>

struct Column
{
  std::string name;
  int width;

  explicit Column(std::string const & colName);
  
  Column(std::string const & colName, int colWidth);
};


class Table
{
public:

  std::vector<Column> cols;
  std::string sepChar = "|";
  std::string padding = " ";
  std::ostream* out = &std::cout;

  Table() {};

  /// Initialize the Table with a list of table headers
  explicit Table(std::initializer_list<std::string> headers);

  /// Initialize the Table with a list of pairs of table headers and width
  explicit Table(std::initializer_list<std::pair<int, std::string>> headers);

  /// Prints the formatted header
  void printHeader();

  /// Prints a line, accepting arbitrary arguments
  template<class ... Ts>
  void printLine(Ts... args)
  {
    printLine(static_cast<size_t>(0), args...);    
  }

  /// Prints a ostream convertible type
  template<class T, class ... Ts>
  void printLine(size_t index, T a, Ts... args)
  {
    int width = index < cols.size() ? cols[index].width : 0;
    *out << padding << std::setw(width) << a << padding << sepChar;
    printLine(index+1, args...);    
  }

  /// Prints a duration as milliseconds
  template<class Rep, class Period, class ... Ts>
  void printLine(size_t index, std::chrono::duration<Rep, Period> duration, Ts... args)
  {
    int width = index < cols.size() ? cols[index].width : 0;
    double ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
    *out << padding << std::setw(width) << ms << padding << sepChar;
    printLine(index+1, args...);    
  }
    
  /// Recursion anchor, prints the last entry and the endl
  template<class T>
  void printLine(size_t index, T a)
  {
    int width = index < cols.size() ? cols[index].width : 0;
    *out << padding << std::setw(width) << a << padding << sepChar
         << std::endl;
  }

    
};
