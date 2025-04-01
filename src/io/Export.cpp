#include "io/Export.hpp"
#include <filesystem>
#include "utils/fmt.hpp"

namespace precice::io {

bool Export::isParallel() const
{
  return _size > 1;
};

std::string Export::formatIndex(int index) const
{
  if (index == 0) {
    return "init";
  }
  using std::string_literals::operator""s;
  return ((_kind == ExportKind::TimeWindows) ? "dt"s : "it"s).append(std::to_string(index));
}

void Export::writeSeriesFile(std::string_view filename) const
{
  if (_records.empty())
    return;

  namespace fs = std::filesystem;
  fs::path outfile(_location);
  if (not _location.empty()) {
    fs::create_directories(outfile);
  }
  outfile /= filename;

  // Prepare filestream
  std::ofstream outFile(outfile.string(), std::ios::trunc);

  outFile << R"({ "file-series-version" : "1.0", "files" : [)";

  for (std::size_t i = 0; i < _records.size() - 1; ++i) {
    outFile << fmt::format(R"( {{ "name" : "{}", "time" : {} }},)", _records[i].filename, _records[i].time);
  }
  outFile << fmt::format(R"( {{ "name" : "{}", "time" : {} }} ] }})", _records.back().filename, _records.back().time);
}

void Export::recordExport(std::string filename, double time)
{
  _records.push_back(Record{filename, time});
}

} // namespace precice::io
