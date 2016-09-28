#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/xmlwriter/XMLWriter.h"
#include <cstdarg>
#include <cstdlib>
#include <sstream>

using namespace tarch::xmlwriter;

logging::Logger XMLWriter::_log("tarch::logging::XMLWriter");

XMLWriter::XMLWriter(std::string sTmp)
{
  sXmlFile = sTmp;
  fp = NULL;
  iLevel = 0;
  fp = fopen(sXmlFile.c_str(),"w");
  if(fp == NULL) {
    std::stringstream ss;
    ss << "Unable to open output file " << sTmp;
    preciceWarning("XMLWriter",ss.str());
  }
}

XMLWriter::XMLWriter(const std::string& sTmp, const std::string& headerInfo)
{
  sXmlFile = sTmp;
  fp = NULL;
  iLevel = 0;
  fp = fopen(sXmlFile.c_str(),"w");
  if(fp == NULL) {
    std::stringstream ss;
    ss << "Unable to open output file " << sTmp;
    preciceWarning("XMLWriter",ss.str());
  } else {
    std::stringstream ss;
    ss << "<?xml version=\"1.0\" encoding=\"UTF-8\" \?>" <<std::endl
        <<"<?" << headerInfo <<"\?>";
    fprintf(fp, "%s" , ss.str().c_str() );
  }
}

XMLWriter::~XMLWriter()
{
  if(fp != NULL) {
    fclose(fp);
    fp = NULL;
  }
  vectAttrData.clear();
}

void XMLWriter::createTag(std::string sTag)
{
  fprintf(fp,"\n");
  //Indent properly
  for(int iTmp =0;iTmp<iLevel;iTmp++)
    fprintf(fp,"\t");
  fprintf(fp,"<%s",sTag.c_str());
  //Add Attributes
  while(0 < vectAttrData.size()/2) {
    std::string sTmp = vectAttrData.back();
    fprintf(fp," %s=", sTmp.c_str());
    vectAttrData.pop_back();
    sTmp = vectAttrData.back();
    fprintf(fp,"\"%s\"", sTmp.c_str());
    vectAttrData.pop_back();
  }
  vectAttrData.clear();
  fprintf(fp,">");
  sTagStack.push(sTag);
  iLevel++;

}

void XMLWriter::createTagWithInformation(std::string sTag, std::string sInformation)
{
  fprintf(fp,"\n");
  //Indent properly
  for(int iTmp =0;iTmp<iLevel;iTmp++)
    fprintf(fp,"\t");
  fprintf(fp,"<%s %s ",sTag.c_str(), sInformation.c_str());
  //Add Attributes
  while(0 < vectAttrData.size()/2) {
    std::string sTmp = vectAttrData.back();
    fprintf(fp," %s=", sTmp.c_str());
    vectAttrData.pop_back();
    sTmp = vectAttrData.back();
    fprintf(fp,"\"%s\"", sTmp.c_str());
    vectAttrData.pop_back();
  }
  vectAttrData.clear();
  fprintf(fp,">");
  sTagStack.push(sTag);
  iLevel++;

}

void XMLWriter::closeLasttag()
{
  fprintf(fp,"\n");
  iLevel--;
  //Indent properlyxml
  for(int iTmp =0;iTmp<iLevel;iTmp++)
    fprintf(fp,"\t");
  fprintf(fp,"</%s>",sTagStack.top().c_str());
  sTagStack.pop();//pop out the last tag
  return;
}

void XMLWriter::closeAlltags()
{
  while(sTagStack.size() != 0) {
    fprintf(fp,"\n");
    iLevel--;
    //Indent properly
    for(int iTmp =0;iTmp<iLevel;iTmp++)
      fprintf(fp,"\t");
    fprintf(fp,"</%s>",sTagStack.top().c_str());
    sTagStack.pop();//pop out the last tag
  }
  return;
}
void XMLWriter::createChild(std::string sTag,std::string sValue)
{
  fprintf(fp,"\n");
  //Indent properly
  for(int iTmp =0;iTmp<iLevel;iTmp++)
    fprintf(fp,"\t");
  fprintf(fp,"<%s",sTag.c_str());
  //Add Attributes
  while(0 < vectAttrData.size()/2) {
    std::string sTmp = vectAttrData.back();
    fprintf(fp," %s=", sTmp.c_str());
    vectAttrData.pop_back();
    sTmp = vectAttrData.back();
    fprintf(fp,"\"%s\"", sTmp.c_str());
    vectAttrData.pop_back();
  }
  vectAttrData.clear();
  //add value and close tag
  fprintf(fp,">%s</%s>",sValue.c_str(),sTag.c_str());
}

void XMLWriter::addAttribute(std::string sKey, std::string sVal)
{
  vectAttrData.push_back(sVal);
  vectAttrData.push_back(sKey);
}

void XMLWriter::addAttribute(std::string sKey, double value)
{
  std::ostringstream strs;
  strs << value;
  vectAttrData.push_back(strs.str());
  vectAttrData.push_back(sKey);
}

void XMLWriter::addAttribute(std::string sKey, int value)
{
  std::ostringstream strs;
  strs << value;
  vectAttrData.push_back(strs.str());
  vectAttrData.push_back(sKey);
}


void XMLWriter::addComment(std::string sComment)
{
  fprintf(fp,"\n");
  //Indent properly
  for(int iTmp =0;iTmp<iLevel;iTmp++)
    fprintf(fp,"\t");
  fprintf(fp,"<!--%s-->",sComment.c_str());
}
