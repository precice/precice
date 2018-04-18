#include "TableWriter.hpp"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <string>
#include <vector>


Column::Column(std::string const & colName)
  : name(colName),
    width(colName.size())
{}

Column::Column(std::string const & colName, int colWidth) :
  name(colName),
  width(colWidth)
{}

Table::Table(std::initializer_list<std::string> headers)
{
  for (auto & h : headers) {
    cols.emplace_back(h);
  }    
}

Table::Table(std::initializer_list<std::pair<int, std::string>> headers)
{
  for (auto & h : headers) {
    cols.emplace_back(std::get<1>(h), std::get<0>(h));
  }    
}

void Table::printHeader()
{
  using namespace std;

  for (auto & h : cols) {
    *out << padding << setw(h.width) << h.name << padding << sepChar;
  }
    
  *out << endl;
  int headerLength = std::accumulate(cols.begin(), cols.end(), 0, [](int count, Column col){
      return count + col.width + 3;
    });
  std::string sepLine(headerLength, '-');
  *out << sepLine << endl;
}
