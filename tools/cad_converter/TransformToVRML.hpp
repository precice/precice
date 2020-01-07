#ifndef TOOLS_DATACONVERTER_HPP_
#define TOOLS_DATACONVERTER_HPP_
#include <string>

class TransformToVRML {
public:
  enum FileType {
    UNDEFINED,
    IGES,
    STEP,
    VRML
  };

  TransformToVRML(
      const std::string &inputFileName,
      const std::string &outputFileName);

  virtual ~TransformToVRML(){};

  bool parseFileType();

  bool doTransform();

private:
  std::string _inputFileName;

  std::string _outputFileName;

  FileType _inputFileType;

  void eliminateRedundancy();

  bool transformIGES();

  bool transformSTEP();
};

#endif /* TOOLS_DATACONVERTER_HPP_ */
