#include "STEP2Vrml.hpp"

STEP2Vrml:: STEP2Vrml ( void )
{

   mySTEPFileName = "";
}

STEP2Vrml:: ~STEP2Vrml ( void )
{
}

bool STEP2Vrml:: readFile ( const std::string & inputFileName )
{
   mySTEPFileName = inputFileName.c_str ();
   if ( !mySTEPFileName || ( mySTEPFileName[0] == '\0' ) ) {
     return false;
   }
   IFSelect_ReturnStatus stat = myReader.ReadFile ( mySTEPFileName );
   if ( stat != IFSelect_RetDone )  {
      return false;
   }
   return true;
}

void STEP2Vrml:: loadReport ( void )
{
   if ( !mySTEPFileName || ( mySTEPFileName [0] == '\0' ) ) {
      cout << "no input Filename" << '\n';
      return;
   }

   // get list of all entities
   Handle( TColStd_HSequenceOfTransient ) myList =
         myReader.GiveList ( "xst-transferrable-all" );

   myReader.TransferList ( myList );

   myReader.PrintStatsTransfer ( true, 3 );
   myReader.ClearShapes ();
}


bool STEP2Vrml:: transfer ( const std::string & outputFileName )
{
   // get list of all entities
   Handle ( TColStd_HSequenceOfTransient ) myList =
         myReader.GiveList ( "xst-transferrable-all" );

   myReader.TransferList ( myList );

   VrmlAPI_Writer myWriter;
   TopoDS_Shape   myOneShape;

   Standard_CString vrmlOutputFile = outputFileName.c_str ();

   for ( int i =  1; i <= myList->Length (); i++ ) {
      myOneShape =  myReader.Shape ( i );
      if ( myOneShape.ShapeType () == TopAbs_SOLID )
         break;
   }

   myWriter.SetRepresentation ( VrmlAPI_ShadedRepresentation );
   myWriter.Write ( myOneShape,vrmlOutputFile );

   return true;
}
