#include "TableWriter.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

Column::Column(std::string const &name)
    : name(name),
      width(name.size())
{
}

Column::Column(std::string const &name, int width)
    : name(name),
      width(std::max(width, static_cast<int>(name.size())))
{
  precision = std::min(this->precision, this->width - 1);
}

Column::Column(std::string const &name, int width, int precision)
    : name(name),
      width(std::max(width, static_cast<int>(name.size())))
{
  this->precision = std::min(precision, this->width - 1);
}

Table::Table()
    : out(std::cout)
{
}

Table::Table(std::ostream &out)
    : out(out)
{
}

void Table::printHeader()
{
  using namespace std;

  for (auto &h : cols) {
    out << padding << setw(h.width) << h.name << padding << sepChar;
  }

  out << '\n';
  int         headerLength = std::accumulate(cols.begin(), cols.end(), 0, [this](int count, Column col) {
    return count + col.width + sepChar.size() + 2;
  });
  std::string sepLine(headerLength, '-');
  out << sepLine << '\n';
}
