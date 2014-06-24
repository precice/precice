#ifndef XMLWriter_h
#define XMLWriter_h

#ifdef Parallel
#include <mpi.h>
#endif
#include <string>
#include <vector>
#include <stack>
#include <cstdio>
#include "tarch/logging/Log.h"

namespace tarch{
  namespace xmlwriter {
    class XMLWriter;
  }}

class tarch::xmlwriter::XMLWriter{
public:
  /**
   * Constructor
   *
   * @sTmp
   */
  XMLWriter(std::string sTmp);
  XMLWriter(const std::string& sTmp, const std::string& headerInfo);
  ~XMLWriter();
  /**
   * creates a new child to an open XML Tag
   */
  void createChild(std::string sTag, std::string sValue);

  /**
   * creates a new child to an open XML Tag
   */
  template<class type>
  void createChild(std::string sTag, std::string sIdentifier, type value)
  {
    fprintf(fp,"\n");
    //Indent properly
    for(int iTmp =0;iTmp<iLevel;iTmp++)
      fprintf(fp,"\t");
    fprintf(fp,"<%s",sTag.c_str());
    //Add Attributes
    while(0 < vectAttrData.size()/2)
    {
      std::string sTmp = vectAttrData.back();
      fprintf(fp," %s=", sTmp.c_str());
      vectAttrData.pop_back();
      sTmp = vectAttrData.back();
      fprintf(fp,"\"%s\"", sTmp.c_str());
      vectAttrData.pop_back();
    }
    vectAttrData.clear();
    //add value and close tag
    std::stringstream ss;
    ss <<value;
    fprintf(fp," %s=\"%s\"/>",sIdentifier.c_str(),ss.str().c_str());
  }

  /**
   * creates a new XML Tag
   */
  void createTag(std::string sTag);
  void createTagWithInformation(std::string sTag, std::string sInformation);
  void closeLasttag();
  const std::string& getFileName() const {
    return sXmlFile;
  }
  void closeAlltags();
  void addAttribute(std::string sAttrName, std::string sAttrvalue);
  void addAttribute(std::string sAttrName, double value);
  void addAttribute(std::string sAttrName, int value);
  void addComment(std::string sComment);
private:
  static tarch::logging::Log _log;
  std::string sXmlFile;
  std::vector<std::string> vectAttrData;
  FILE *fp;
  int iLevel;
  std::stack<std::string>  sTagStack;
};

#endif // XMLWriter_h
