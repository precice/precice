#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
// Copyright (C) 2002-2005 Nikolaus Gebhardt
// This file is part of the "Irrlicht Engine" and the "irrXML" project.
// For conditions of distribution and use, see copyright notice in irrlicht.h and/or irrXML.h

#include "tarch/irr/XML.h"
#include "tarch/irr/String.h"
#include "tarch/irr/Array.h"
#include "tarch/irr/FastStringConversion.h"
#include "tarch/irr/CXMLReaderImpl.h"
#include <string.h>


//
// Added by the Peano project
//
namespace tarch {
namespace irr
{
namespace io
{

//! Implementation of the file read callback for ordinary files
class CStringReadCallBack : public IFileReadCallBack
{
private:
	char *_stringPointer;
	std::string _xmlString;
	int _bytesRead;
	int Size;
	void getFileSize()
	{

		Size = static_cast<int>(_xmlString.size());

	}
public:
	//construct from string
	CStringReadCallBack(std::string xmlString)
		:_xmlString(xmlString),_bytesRead(0)
	{
		getFileSize();
		_stringPointer=(char*)_xmlString.c_str();
	}
	virtual ~CStringReadCallBack()
	{
	}
	//! Reads an amount of bytes from the file.
	virtual int read(void* buffer, int sizeToRead)
	{

		if(_bytesRead>=Size)
			return 0;
		if(_bytesRead+sizeToRead>Size)
			sizeToRead=Size-_bytesRead;
		_bytesRead+=sizeToRead;
		memcpy(buffer,_stringPointer,sizeToRead);
		_stringPointer+=sizeToRead;
		return sizeToRead;
	}
	//! Returns size of string in bytes
	virtual int getSize()
	{
		return Size;
	}

};
class CFileReadCallBack : public IFileReadCallBack
{
public:

	//! construct from filename
	CFileReadCallBack(const char* filename)
		: File(0), Size(0), Close(true)
	{
		// open file
		File = fopen(filename, "rb");

		if (File)
			getFileSize();
	}

	//! construct from FILE pointer
	CFileReadCallBack(FILE* file)
		: File(file), Size(0), Close(false)
	{
		if (File)
			getFileSize();
	}


	//! destructor
	virtual ~CFileReadCallBack()
	{
		if (Close && File)
			fclose(File);
	}

	//! Reads an amount of bytes from the file.
	virtual int read(void* buffer, int sizeToRead)
	{
		if (!File)
			return 0;

		return (int)fread(buffer, 1, sizeToRead, File);
	}

	//! Returns size of file in bytes
	virtual int getSize()
	{
		return Size;
	}

private:

	//! retrieves the file size of the open file
	void getFileSize()
	{
		fseek(File, 0, SEEK_END);
		Size = static_cast<int>(ftell(File));
		fseek(File, 0, SEEK_SET);
	}

	FILE* File;
	int Size;
	bool Close;

}; // end class CFileReadCallBack



// FACTORY FUNCTIONS:


IrrXMLReader* createIrrXMLReader(const char* filename)
{
	return new CXMLReaderImpl<char, IXMLBase>(new CFileReadCallBack(filename));
}

IrrXMLReader* createIrrXMLReaderFromString(std::string xmlString)
{
	return new CXMLReaderImpl<char, IXMLBase>(new CStringReadCallBack(xmlString));
}
IrrXMLReader* createIrrXMLReader(FILE* file)
{
	return new CXMLReaderImpl<char, IXMLBase>(new CFileReadCallBack(file));
}


IrrXMLReader* createIrrXMLReader(IFileReadCallBack* callback)
{
	return new CXMLReaderImpl<char, IXMLBase>(callback, false);
}

} // end namespace io
} // end namespace irr
}
