// ---- Open-CASCADE Includes ----
// Geometric Data Structure
#include <gp_Pnt.hxx>							// Geometric Point
#include <gp_Lin.hxx>							// Geometric Line
// Topologic Data Structure
#include <TopoDS.hxx>							// Topology 
#include <TopoDS_Vertex.hxx>					// Topological node
#include <TopoDS_Edge.hxx>						// Topological edge
#include <TopoDS_Wire.hxx>						// Topological wire
#include <TopoDS_Face.hxx>						// Topological face
#include <TopoDS_Shell.hxx>						// Topological shell
#include <TopoDS_Solid.hxx>						// Topological solid
#include <TopoDS_CompSolid.hxx>					// Topological compound solid
#include <TopoDS_Compound.hxx>					// Topological compound
#include <TopTools_HSequenceOfShape.hxx>		// Handle sequence of shapes
// Import/Export-Tools
#include <Interface_Static.hxx>					// Parser tool
#include <BRepTools.hxx>						// Parser fuer BREP files
#include <BRep_Builder.hxx>						// Parser fuer BREP files
#include <IGESControl_Controller.hxx>			// Parser fuer Igs files
#include <IGESControl_Reader.hxx>				// Parser fuer Igs files
#include <IGESControl_Writer.hxx>				// Parser fuer Igs files
#include <STEPControl_Reader.hxx>				// Parser fuer Step files
#include <STEPControl_Writer.hxx>				// Parser fuer Step files
#include <FSD_File.hxx>							// Parser fuer CSFDB files
#include <ShapeSchema.hxx>						// Parser fuer CSFDB files
#include <Storage_Data.hxx>						// Parser fuer CSFDB files
#include <Storage_Root.hxx>						// Parser fuer CSFDB files
#include <Storage_HSeqOfRoot.hxx>				// Parser fuer CSFDB files
#include <PTopoDS_HShape.hxx>					// Parser fuer CSFDB files
#include <PTColStd_PersistentTransientMap.hxx>	// Parser fuer CSFDB files
#include <PTColStd_TransientPersistentMap.hxx>	// Parser fuer CSFDB files
#include <MgtBRep.hxx>							// Parser fuer CSFDB files
#include <MgtBRep_TriangleMode.hxx>				// Parser fuer CSFDB files
#include <StlAPI_Writer.hxx>					// Parser fuer STL files
#include <VrmlAPI_Writer.hxx>					// Parser fuer VRML files
#include <TCollection_ExtendedString.hxx>
// Explorer-Tools
#include <BRep_Tool.hxx>
#include <TopExp_Explorer.hxx>					// Explorer fuer TopoDS	








/***************************************************************************
 **** class Geo_Shape                                                   ****
 ***************************************************************************/
class Geo_Shape
{
private:
// ---- Private variables ----
	// Open-CASCADE objects
	TopoDS_Shape						my_TopoDS_Shape;
	Handle(TopTools_HSequenceOfShape)	my_HSequenceOfShape;
	
public:
// ---- Constructors and Destructors ----
	Geo_Shape();

public:
	// Import
	bool import_IgesFile(char* filename);
	bool import_StepFile(char* filename);
	// Export
	bool export_BrepFile(char* filename);
	bool export_IgesFile(char* filename);
	bool export_StepFile(char* filename);
};





/***************************************************************************
 **** class Geo_Shape : Constructor                                     ****
 ***************************************************************************/
Geo_Shape::Geo_Shape()
{
// ---- Initialize Geo_Shape ----
	// Open-CASCADE objects (TopoDS_Face)
	my_TopoDS_Shape.Nullify();
	my_HSequenceOfShape.Nullify();
}


/***************************************************************************
 **** class Geo_Shape : Public member funcion                        ****
 ***************************************************************************/
bool Geo_Shape::import_IgesFile(char* filename)
{
    IGESControl_Reader aReader;
    int status = aReader.ReadFile(filename);
    if ( status != 0 )
    {
        aReader.TransferRoots();
        my_TopoDS_Shape = aReader.OneShape();
		my_HSequenceOfShape = new TopTools_HSequenceOfShape();
		my_HSequenceOfShape->Append( my_TopoDS_Shape );
    }
	return true;
}


/***************************************************************************
 **** class Geo_Shape : Public member funcion                        ****
 ***************************************************************************/
bool Geo_Shape::import_StepFile(char* filename)
{
	STEPControl_Reader aReader;
    IFSelect_ReturnStatus status = aReader.ReadFile( filename );
    if ( status != 0 )
    {
        aReader.TransferRoots();
        my_TopoDS_Shape = aReader.OneShape();
		my_HSequenceOfShape = new TopTools_HSequenceOfShape();
		my_HSequenceOfShape->Append( my_TopoDS_Shape );
    }
	return true;
}


/***************************************************************************
 **** class Geo_Shape : Public member funcion                        ****
 ***************************************************************************/
bool Geo_Shape::export_IgesFile(char* filename)
{
// ---- Check sequence of shapes ----
    if ( my_HSequenceOfShape.IsNull() || my_HSequenceOfShape->IsEmpty() )
        return false;
// ---- Initialize writer ----
	IGESControl_Controller::Init();
	IGESControl_Writer writer( Interface_Static::CVal( "XSTEP.iges.unit" ),
                               Interface_Static::IVal( "XSTEP.iges.writebrep.mode" ) );
 // ---- Write sequence of shapes ----	
	for ( int i = 1; i <= my_HSequenceOfShape->Length(); i++ )
		writer.AddShape ( my_HSequenceOfShape->Value( i ) );
	writer.ComputeModel();
	writer.Write(filename);
	return true;
}


/***************************************************************************
 **** class Geo_Shape : Public member funcion                        ****
 ***************************************************************************/
bool Geo_Shape::export_StepFile(char* filename)
{
// ---- Check sequence of shapes ----
    if ( my_HSequenceOfShape.IsNull() || my_HSequenceOfShape->IsEmpty() )
        return false;
// ---- Set step-control ----
    IFSelect_ReturnStatus status;
	STEPControl_StepModelType type = STEPControl_AsIs;
	//STEPControl_StepModelType type = STEPControl_ManifoldSolidBrep;
	//STEPControl_StepModelType type = STEPControl_FacetedBrep;
	//STEPControl_StepModelType type = STEPControl_ShellBasedSurfaceModel;
	//STEPControl_StepModelType type = STEPControl_GeometricCurveSet;
// ---- Initialize writer ----
	STEPControl_Writer writer;
	for ( int i = 1; i <= my_HSequenceOfShape->Length(); i++ )
    {
		status = writer.Transfer( my_HSequenceOfShape->Value( i ), type );
        if ( status != IFSelect_RetDone )
            return false;
    }
// ---- Falls Sequence of shapes leer
    status = writer.Write(filename);
    return status == IFSelect_RetDone;
}
