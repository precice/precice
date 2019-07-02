#include "IGES2Vrml.hpp"

IGES2Vrml:: IGES2Vrml ( void )
{
   myIGESFileName = "";
}

IGES2Vrml:: ~IGES2Vrml ( void )
{
}

bool IGES2Vrml:: readFile (  const std::string & inputFileName )
{
   myIGESFileName = inputFileName.c_str ();
   if ( !myIGESFileName || ( myIGESFileName[0] == '\0' ) ) {
      return false;
   }

   IFSelect_ReturnStatus stat = myReader.ReadFile( myIGESFileName );

   if ( stat != IFSelect_RetDone)  {
      return false;
   }

   return true;

}

void IGES2Vrml:: loadReport ( void )
{
   if ( !myIGESFileName || ( myIGESFileName [0] == '\0' ) ) {
      std::cout << "no input file" << '\n';
      return;
   }

   // get list of all entities
   Handle ( TColStd_HSequenceOfTransient ) myList =
         myReader.GiveList( "xst-transferrable-all" );

   myReader.TransferList ( myList );

   myReader.PrintTransferInfo ( IFSelect_FailAndWarn, IFSelect_Mapping );

   myReader.ClearShapes ();

}

bool IGES2Vrml:: transfer( const std::string & outputFileName )
{
   // get list of all entities
   Handle ( TColStd_HSequenceOfTransient ) myList =
         myReader.GiveList ( "xst-transferrable-all" );

   myReader.TransferList ( myList );

   VrmlAPI_Writer myWriter;
   TopoDS_Shape   myOneShape;

   for ( int i =  1; i < myList-> Length (); i++ ) {

      myOneShape =  myReader.Shape ( i );
      if ( myOneShape.ShapeType () == TopAbs_COMPOUND ) {
         break;
      }
   }

   Standard_CString vrmlOutputFile = outputFileName.c_str ();

   myWriter.SetRepresentation ( VrmlAPI_ShadedRepresentation );
   myWriter.Write ( myOneShape, vrmlOutputFile );
   return true;
}

