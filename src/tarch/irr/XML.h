// Copyright (C) 2002-2005 Nikolaus Gebhardt
// This file is part of the "Irrlicht Engine" and the "irrXML" project.
// For conditions of distribution and use, see copyright notice in irrlicht.h and/or irrXML.h

#ifndef __IRR_XML_H_INCLUDED__
#define __IRR_XML_H_INCLUDED__

#ifdef Parallel
#include <mpi.h>
#endif
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include "tarch/la/Scalar.h"
#include "tarch/la/Vector.h"
#include "tarch/la/DynamicVector.h"
#include "tarch/logging/Log.h"

//
// Added by the Peano project
//
namespace tarch {
  namespace irr
  {
    namespace io
    {
      //! Enumeration of all supported source text file formats
      enum ETEXT_FORMAT
      {
        //! ASCII, file without byte order mark, or not a text file
        ETF_ASCII,

        //! UTF-8 format
        ETF_UTF8,

        //! UTF-16 format, big endian
        ETF_UTF16_BE,

        //! UTF-16 format, little endian
        ETF_UTF16_LE,

        //! UTF-32 format, big endian
        ETF_UTF32_BE,

        //! UTF-32 format, little endian
        ETF_UTF32_LE
      };


      //! Enumeration for all xml nodes which are parsed by IrrXMLReader
      enum EXML_NODE
      {
        //! No xml node. This is usually the node if you did not read anything yet.
        EXN_NONE = 0,

        //! A xml element, like "smaller_than"foo>
        EXN_ELEMENT = 1,

        //! End of an xml element, like "smaller_than"/foo>
        EXN_ELEMENT_END = 2,

        //! Text within a xml element: "smaller_than" this is the text. "smaller_than"/foo>
        EXN_TEXT = 3,

        //! An xml comment like &lt;!-- I am a comment --&gt; or a DTD definition.
        EXN_COMMENT = 4,

        //! An xml cdata section like &lt;![CDATA[ this is some CDATA ]]&gt;
        EXN_CDATA = 5,

        //! Unknown element.
        EXN_UNKNOWN = 6
      };

      //! Callback class for file read abstraction.
      /** With this, it is possible to make the xml parser read in other things
	than just files. The Irrlicht engine is using this for example to
	read xml from compressed .zip files. To make the parser read in
	any other data, derive a class from this interface, implement the
	two methods to read your data and give a pointer to an instance of
	your implementation when calling createIrrXMLReader(),
	createIrrXMLReaderUTF16() or createIrrXMLReaderUTF32() */
      class IFileReadCallBack
      {

      public:

        //! virtual destructor
        virtual ~IFileReadCallBack() {};

        //! Reads an amount of bytes from the file.
        /** \param buffer: Pointer to buffer where to read bytes will be written to.
		\param sizeToRead: Amount of bytes to read from the file.
		\return Returns how much bytes were read. */
        virtual int read(void* buffer, int sizeToRead) = 0;

        //! Returns size of file in bytes
        virtual int getSize() = 0;
      };

      //! Empty class to be used as parent class for IrrXMLReader.
      /** If you need another class as base class for the xml reader, you can do this by creating
	the reader using for example new CXMLReaderImpl<char, YourBaseClass>(yourcallback);
	The Irrlicht Engine for example needs IUnknown as base class for every object to
	let it automaticly reference countend, hence it replaces IXMLBase with IUnknown.
	See irrXML.cpp on how this can be done in detail. */
      class IXMLBase
      {
      };

      template<class char_type, class super_class>
      class IIrrXMLReader : public super_class
      {

      private:
        static tarch::logging::Log _log;
      public:

        //! Destructor
        virtual ~IIrrXMLReader() {};

        //! Reads forward to the next xml node.
        /** \return Returns false, if there was no further node.  */
        virtual bool read() = 0;

        //! Returns the type of the current XML node.
        virtual EXML_NODE getNodeType() const = 0;

        //! Returns attribute count of the current XML node.
        /** This is usually
		non null if the current node is EXN_ELEMENT, and the element has attributes.
		\return Returns amount of attributes of this xml node. */
        virtual int getAttributeCount() const = 0;

        //! Returns name of an attribute.
        /** \param idx: Zero based index, should be something between 0 and getAttributeCount()-1.
		\return Name of the attribute, 0 if an attribute with this index does not exist. */
        virtual const char_type* getAttributeName(int idx) const = 0;

        //! Returns the value of an attribute.
        /** \param idx: Zero based index, should be something between 0 and getAttributeCount()-1.
		\return Value of the attribute, 0 if an attribute with this index does not exist. */
        virtual const char_type* getAttributeValue(int idx) const = 0;

        //! Returns the value of an attribute.
        /** \param name: Name of the attribute.
		\return Value of the attribute, 0 if an attribute with this name does not exist. */
        virtual const char_type* getAttributeValue(const char_type* name) const = 0;

        //! Returns the value of an attribute.
        /** \param name: Name of the attribute.
            \return Value of the attribute, 0 if an attribute with this name does not exist. */
        const char_type* getAttributeValue(const std::string& name) const {
          return getAttributeValue(name.c_str());
        }

        //! Returns the value of an attribute in a safe way.
        /** Like getAttributeValue(), but does not
		return 0 if the attribute does not exist. An empty string ("") is returned then.
		\param name: Name of the attribute.
		\return Value of the attribute, and "" if an attribute with this name does not exist */
        virtual const char_type* getAttributeValueSafe(const char_type* name) const = 0;

        //! Returns the value of an attribute as integer.
        /** \param name Name of the attribute.
		\return Value of the attribute as integer, and 0 if an attribute with this name does not exist or
		the value could not be interpreted as integer. */
        virtual int getAttributeValueAsInt(const char_type* name) const = 0;

        /* Returns the value of an attribute as integer.
         * param name Name of the attribute.
         * return Value of the attribute as integer, and 0 if an attribute with this name does not exist or
         * the value could not be interpreted as integer.
         **/
        inline int getAttributeValueAsInt(
            const std::string& name
            ) const {
          return getAttributeValueAsInt(name.c_str());
        }

        //! Returns the value of an attribute as integer.
        /** \param idx: Zero based index, should be something between 0 and getAttributeCount()-1.
		\return Value of the attribute as integer, and 0 if an attribute with this index does not exist or
		the value could not be interpreted as integer. */
        virtual int getAttributeValueAsInt(int idx) const = 0;

        //! Returns the value of an attribute as float.
        /** \param name: Name of the attribute.
		\return Value of the attribute as float, and 0 if an attribute with this name does not exist or
		the value could not be interpreted as float. */
        virtual float getAttributeValueAsFloat(const char_type* name) const = 0;

        //! Returns the value of an attribute as float.
        /** \param name: Name of the attribute.
            \return Value of the attribute as float, and 0 if an attribute with this name does not exist or
            the value could not be interpreted as float. */
        inline float getAttributeValueAsFloat(const std::string& name) const {
          return getAttributeValueAsFloat(name.c_str());
        }

        //! Returns the value of an attribute as float.
        /** \param idx: Zero based index, should be something between 0 and getAttributeCount()-1.
		\return Value of the attribute as float, and 0 if an attribute with this index does not exist or
		the value could not be interpreted as float. */
        virtual float getAttributeValueAsFloat(int idx) const = 0;

        /**
         * Added by Tobias Weinzierl 30.8.2006
         */
        virtual double getAttributeValueAsDouble(const char_type* name) const = 0;

        /**
         * Added by Michael Lieb 06.08.2012
         */
        inline double getAttributeValueAsDouble(const std::string& name) const {
          return getAttributeValueAsDouble(name.c_str());
        }

        /**
         * Added by Tobias Weinzierl 30.8.2006
         */
        virtual double getAttributeValueAsDouble(int idx) const = 0;


        /**
         * Added by Tobias Neckel 27.04.2007. Additional support of fraction
         * specification (originally programmed by Tobias Weinzierl): One may now use
         * fractions such as "1/3" etc. and the real value of pi (more accurate than
         * "0.3333333" and "3.141592").
         *
         * This method encapsulates everything necessary to get a double value out of
         * a string.
         *
         * This method returns a double value of an attribute. It is possible to
         * specify fractions (such as "1/9" instead of "0.1111111") and pi. If some
         * problems exist, numeric_limits::signaling_NaN() is returned. It is assumed
         * that numeric_limits::has_signaling_NaN() is true in order to do so.
         *
         * @param value String containing the desired double value (possibly as a
         *              fraction).
         *
         * @return Double value of the desired attribute name.
         */
        double convertValueStringToDouble(const std::string& value) const {
          std::string myString = value;
          while (myString.find(" ") != std::string::npos){
            myString.erase(myString.find(" "),1);
          }
          double result = std::numeric_limits<double>::signaling_NaN();
          if ( myString.find("/") == std::string::npos ) {  //no "/" in
            if (myString.find("pi") == std::string::npos) { //no "pi" in
              std::istringstream in( myString);
              if ( !(in >> result) ) {
                result = std::numeric_limits<double>::signaling_NaN();
              }
            }
            else {                                       // "pi" in
              std::string left  = myString.substr( 0, myString.find("pi") );
              std::string right = myString.substr( myString.find("pi")+2, myString.size()-myString.find("pi")-1 );
              double leftValue;
              if (left!="") leftValue = convertValueStringToDouble( left );
              else          leftValue = 1.0;
              double rightValue;
              if (right!="") rightValue = convertValueStringToDouble( right );
              else           rightValue = 1.0;
              result = leftValue * tarch::la::PI * rightValue;
            }
          }
          else {                                           // "/" in
            std::string left  = myString.substr( 0,myString.find("/") );
            std::string right = myString.substr( myString.find("/")+1, myString.size()-myString.find("/")-1 );
            double leftValue  = convertValueStringToDouble( left );
            double rightValue = convertValueStringToDouble( right );
            if (rightValue==0) {
              result = std::numeric_limits<double>::signaling_NaN();
            }
            else if (   rightValue==std::numeric_limits<double>::signaling_NaN()
                || leftValue==std::numeric_limits<double>::signaling_NaN() ) {
              result = std::numeric_limits<double>::signaling_NaN();
            }
            else {
              result = leftValue / rightValue;
            }
          }
          return result;
        }

        /**
         * Evaluates a string to find out if it represents a bool.
         *
         * @author Tobias Weinzierl
         */
        bool convertValueStringToBool( const std::string& value ) const {
          if ( value=="1" ) return true;
          if ( value=="yes" ) return true;
          if ( value=="true" ) return true;
          if ( value=="on" ) return true;

          return false;
        }


        /**
         * Reads a attribute with the given name which is of (double) vector type from
         * the xmlReader and stores it in the vector. The dimension of the vector is
         * determined by the <i>Dim</i> defines of the preprocessor.
         *
         * @param attributeName - name of attribute
         * attributeName = [0.0 1.0 1.5]
         */
        template<int num>
        inline tarch::la::Vector<num,double> getAttributeValueAsDoubleVector(
            const std::string& attributeName
        ) {
          return getAttributeValueAsDoubleVector<num>(attributeName.c_str());
        }

        /**
         * Reads a attribute with the given name which is of (double) vector type from
         * the xmlReader and stores it in the vector. The dimension of the vector is
         * determined by the <i>Dim</i> defines of the preprocessor.
         *
         * @param attributeName - name of attribute
         * attributeName = [0.0 1.0 1.5]
         */
        template<int num>
        tarch::la::Vector<num,double> getAttributeValueAsDoubleVector(
            const char_type* attributeName
        ){
          tarch::la::Vector<num,double> vec(0.0);
          const char_type* attributeValue = getAttributeValue(attributeName);
          if (attributeValue==0) return vec;

          std::string valueString(attributeValue);
          for (int i = 0; i < num; i++){
            std::string currentEntry = valueString;

            size_t found = valueString.find(";");
            if (found != std::string::npos){
              valueString.erase(0,found+1);
              currentEntry.erase(found);
            }
            else if (i != num-1) {
              logError( "getAttributeValueAsDoubleVector(...)", "expected " << num << " ; within value " << std::string(getAttributeValue(attributeName)) );
              return 0.0;
            }

            vec(i) = convertValueStringToDouble(currentEntry);
          }
          return vec;
        }

        /**
         * Reads an attribute with the given name which is of tarch::la::DynamicVector<double>
         * type from the xmlReader and stores it in the vector.
         *
         * @param attributeName - name of attribute
         * attributeName = [0.0 1.0 1.5]
         */
        inline tarch::la::DynamicVector<double> getAttributeValueAsDynamicDoubleVector(
            const std::string& attributeName
        ){
          return getAttributeValueAsDynamicDoubleVector(attributeName.c_str());
        }
        /**
         * Reads an attribute with the given name which is of tarch::la::DynamicVector<double>
         * type from the xmlReader and stores it in the vector.
         *
         * @param attributeName - name of attribute
         * attributeName = [0.0 1.0 1.5]
         */
        tarch::la::DynamicVector<double> getAttributeValueAsDynamicDoubleVector(
            const char_type* attributeName
        ){
          tarch::la::DynamicVector<double> vec;
          const char_type* attributeValue = getAttributeValue(attributeName);
          if (attributeValue==0) return vec;

          std::string valueString(attributeValue);
          bool componentsLeft = true;
          int i =0;
          while (componentsLeft){
            std::string tmp1(valueString);
            // erase entries before i-th entry
            for (int j = 0; j < i; j++){
              if (tmp1.find(";") != std::string::npos){
                tmp1.erase(0,tmp1.find(";")+1);
              }
              else {
                componentsLeft = false;
              }
            }
            // if we are not in the last vector component...
            if (tmp1.find(";") != std::string::npos){
              // ..., erase entries after i-th entry
              tmp1.erase(tmp1.find(";"),tmp1.size());
            }
            if (componentsLeft){
              vec.append(convertValueStringToDouble(tmp1));
            }
            i++;
          }
          return vec;
        }

        /**
         * Reads an attribute with the given name which is of (double) vector type from
         * the xmlReader and stores it in the vector. The dimension of the vector is
         * determined by the <i>Dim</i> defines of the preprocessor.
         */
        template<int num>
        tarch::la::Vector<num,int> getAttributeValueAsIntVector(
            const char_type* attributeName
        ){
          tarch::la::Vector<num,int> vec(0);
          int count = 0;
          bool continueReading = true;
          const char_type* attributeValue = getAttributeValue(attributeName);
          if (attributeValue==0) return vec;
          std::string valueString(attributeValue);
          std::stringstream ss(std::stringstream::in | std::stringstream::out);
          ss.clear();
          ss.str("");
          ss << valueString;
          while ( (!ss.eof()) && continueReading){
            if (count < num){
              std::string tmp1(ss.str());
              std::string tmp2(tmp1);

              // remove ';'-parts
              if (count < num-1){
                if (tmp2.find(";") == std::string::npos){
                  return vec;
                }
                tmp2.erase(0,tmp2.find(";"));

                // determine string for this coordinate
                tmp1 = tmp1.erase(tmp1.find(tmp2), tmp2.size());

                tmp2.erase(0,tmp2.find(";")+1);
              }

              // put value back to stringstream and read it into vector
              ss.clear();
              ss.str(tmp1);
              ss >> vec(count);

              // get remaining string into stringstream
              ss.clear();
              ss.str(tmp2);
            } else {
              continueReading = false;
            }
            count++;
          }
          return vec;
        }

        /**
         * Added by Michael Lieb 08.08.2012
         */
        inline bool getAttributeValueAsBool(
            const std::string& name
            ) const {
          return getAttributeValueAsBool(name.c_str());
        }

        /**
         * Added by Tobias Weinzierl 1.10.2007
         */
        virtual bool getAttributeValueAsBool(const char_type* name) const = 0;

        /**
         * Added by Tobias Weinzierl 1.10.2007
         */
        virtual bool getAttributeValueAsBool(int idx) const = 0;

        //! Returns the name of the current node.
        /** Only non null, if the node type is EXN_ELEMENT.
		\return Name of the current node or 0 if the node has no name. */
        virtual const char_type* getNodeName() const = 0;

        //! Returns data of the current node.
        /** Only non null if the node has some
		data and it is of type EXN_TEXT or EXN_UNKNOWN. */
        virtual const char_type* getNodeData() const = 0;

        //! Returns if an element is an empty element, like "smaller_than" foo />
        virtual bool isEmptyElement() const = 0;

        //! Returns format of the source xml file.
        /** It is not necessary to use
		this method because the parser will convert the input file format
		to the format wanted by the user when creating the parser. This
		method is useful to get/display additional informations. */
        virtual ETEXT_FORMAT getSourceFormat() const = 0;

        //! Returns format of the strings returned by the parser.
        /** This will be UTF8 for example when you created a parser with
		IrrXMLReaderUTF8() and UTF32 when it has been created using
		IrrXMLReaderUTF32. It should not be necessary to call this
		method and only exists for informational purposes. */
        virtual ETEXT_FORMAT getParserFormat() const = 0;
      };

      template<class char_type, class super_class>
      tarch::logging::Log IIrrXMLReader<char_type, super_class>::_log("tarch::logging::Log");


      //! defines the utf-16 type.
      /** Not using wchar_t for this because
	wchar_t has 16 bit on windows and 32 bit on other operating systems. */
      typedef unsigned short char16;

      //! defines the utf-32 type.
      /** Not using wchar_t for this because
	wchar_t has 16 bit on windows and 32 bit on other operating systems. */
      typedef unsigned long char32;

      //! A UTF-8 or ASCII character xml parser.
      /** This means that all character data will be returned in 8 bit ASCII or UTF-8 by this parser.
	The file to read can be in any format, it will be converted to UTF-8 if it is not
	in this format.
	Create an instance of this with createIrrXMLReader();
	See IIrrXMLReader for description on how to use it. */
      typedef IIrrXMLReader<char, IXMLBase> IrrXMLReader;

      //! A UTF-16 xml parser.
      /** This means that all character data will be returned in UTF-16 by this parser.
	The file to read can be in any format, it will be converted to UTF-16 if it is not
	in this format.
	Create an instance of this with createIrrXMLReaderUTF16();
	See IIrrXMLReader for description on how to use it.  */
      //	typedef IIrrXMLReader<char16, IXMLBase> IrrXMLReaderUTF16;

      //! A UTF-32 xml parser.
      /** This means that all character data will be returned in UTF-32 by this parser.
	The file to read can be in any format, it will be converted to UTF-32 if it is not
	in this format.
	Create an instance of this with createIrrXMLReaderUTF32();
	See IIrrXMLReader for description on how to use it. */
      //	typedef IIrrXMLReader<char32, IXMLBase> IrrXMLReaderUTF32;


      //! Creates an instance of an UFT-8 or ASCII character xml parser.
      /** This means that all character data will be returned in 8 bit ASCII or UTF-8.
	The file to read can be in any format, it will be converted to UTF-8 if it is not in this format.
	If you are using the Irrlicht Engine, it is better not to use this function but
	IFileSystem::createXMLReaderUTF8() instead.
	\param filename: Name of file to be opened.
	\return Returns a pointer to the created xml parser. This pointer should be
	deleted using 'delete' after no longer needed. Returns 0 if an error occured
	and the file could not be opened. */
      IrrXMLReader* createIrrXMLReader(const char* filename);
      IrrXMLReader* createIrrXMLReaderFromString(std::string);
      //! Creates an instance of an UFT-8 or ASCII character xml parser.
      /** This means that all character data will be returned in 8 bit ASCII or UTF-8. The file to read can
	be in any format, it will be converted to UTF-8 if it is not in this format.
	If you are using the Irrlicht Engine, it is better not to use this function but
	IFileSystem::createXMLReaderUTF8() instead.
	\param file: Pointer to opened file, must have been opened in binary mode, e.g.
	using fopen("foo.bar", "wb"); The file will not be closed after it has been read.
	\return Returns a pointer to the created xml parser. This pointer should be
	deleted using 'delete' after no longer needed. Returns 0 if an error occured
	and the file could not be opened. */
      IrrXMLReader* createIrrXMLReader(FILE* file);

      //! Creates an instance of an UFT-8 or ASCII character xml parser.
      /** This means that all character data will be returned in 8 bit ASCII or UTF-8. The file to read can
	 be in any format, it will be converted to UTF-8 if it is not in this format.
	 If you are using the Irrlicht Engine, it is better not to use this function but
	 IFileSystem::createXMLReaderUTF8() instead.
	 \param callback: Callback for file read abstraction. Implement your own
	 callback to make the xml parser read in other things than just files. See
	 IFileReadCallBack for more information about this.
	 \return Returns a pointer to the created xml parser. This pointer should be
	 deleted using 'delete' after no longer needed. Returns 0 if an error occured
	 and the file could not be opened. */
      IrrXMLReader* createIrrXMLReader(IFileReadCallBack* callback);

      //! Creates an instance of an UFT-16 xml parser.
      /** This means that
	all character data will be returned in UTF-16. The file to read can
	be in any format, it will be converted to UTF-16 if it is not in this format.
	If you are using the Irrlicht Engine, it is better not to use this function but
	IFileSystem::createXMLReader() instead.
	\param filename: Name of file to be opened.
	\return Returns a pointer to the created xml parser. This pointer should be
	deleted using 'delete' after no longer needed. Returns 0 if an error occured
	and the file could not be opened. */
      //	IrrXMLReaderUTF16* createIrrXMLReaderUTF16(const char* filename);

      //! Creates an instance of an UFT-16 xml parser.
      /** This means that all character data will be returned in UTF-16. The file to read can
	be in any format, it will be converted to UTF-16 if it is not in this format.
	If you are using the Irrlicht Engine, it is better not to use this function but
	IFileSystem::createXMLReader() instead.
	\param file: Pointer to opened file, must have been opened in binary mode, e.g.
	using fopen("foo.bar", "wb"); The file will not be closed after it has been read.
	\return Returns a pointer to the created xml parser. This pointer should be
	deleted using 'delete' after no longer needed. Returns 0 if an error occured
	and the file could not be opened. */
      //	IrrXMLReaderUTF16* createIrrXMLReaderUTF16(FILE* file);

      //! Creates an instance of an UFT-16 xml parser.
      /** This means that all character data will be returned in UTF-16. The file to read can
	be in any format, it will be converted to UTF-16 if it is not in this format.
	If you are using the Irrlicht Engine, it is better not to use this function but
	IFileSystem::createXMLReader() instead.
	\param callback: Callback for file read abstraction. Implement your own
	callback to make the xml parser read in other things than just files. See
	IFileReadCallBack for more information about this.
	\return Returns a pointer to the created xml parser. This pointer should be
	deleted using 'delete' after no longer needed. Returns 0 if an error occured
	and the file could not be opened. */
      //	IrrXMLReaderUTF16* createIrrXMLReaderUTF16(IFileReadCallBack* callback);


      //! Creates an instance of an UFT-32 xml parser.
      /** This means that all character data will be returned in UTF-32. The file to read can
	be in any format, it will be converted to UTF-32 if it is not in this format.
	If you are using the Irrlicht Engine, it is better not to use this function but
	IFileSystem::createXMLReader() instead.
	\param filename: Name of file to be opened.
	\return Returns a pointer to the created xml parser. This pointer should be
	deleted using 'delete' after no longer needed. Returns 0 if an error occured
	and the file could not be opened. */
      //	IrrXMLReaderUTF32* createIrrXMLReaderUTF32(const char* filename);

      //! Creates an instance of an UFT-32 xml parser.
      /** This means that all character data will be returned in UTF-32. The file to read can
	be in any format, it will be converted to UTF-32 if it is not in this format.
	if you are using the Irrlicht Engine, it is better not to use this function but
	IFileSystem::createXMLReader() instead.
	\param file: Pointer to opened file, must have been opened in binary mode, e.g.
	using fopen("foo.bar", "wb"); The file will not be closed after it has been read.
	\return Returns a pointer to the created xml parser. This pointer should be
	deleted using 'delete' after no longer needed. Returns 0 if an error occured
	and the file could not be opened. */
      //	IrrXMLReaderUTF32* createIrrXMLReaderUTF32(FILE* file);

      //! Creates an instance of an UFT-32 xml parser.
      /** This means that
	all character data will be returned in UTF-32. The file to read can
	be in any format, it will be converted to UTF-32 if it is not in this format.
	If you are using the Irrlicht Engine, it is better not to use this function but
	IFileSystem::createXMLReader() instead.
	\param callback: Callback for file read abstraction. Implement your own
	callback to make the xml parser read in other things than just files. See
	IFileReadCallBack for more information about this.
	\return Returns a pointer to the created xml parser. This pointer should be
	deleted using 'delete' after no longer needed. Returns 0 if an error occured
	and the file could not be opened. */
      //	IrrXMLReaderUTF32* createIrrXMLReaderUTF32(IFileReadCallBack* callback);
    } // end namespace io
  } // end namespace irr
}

#endif // __IRR_XML_H_INCLUDED__

