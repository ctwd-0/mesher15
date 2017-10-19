#pragma once

#include <opennurbs.h>

#include "bspline_conversion.h"


#include <Precision.hxx>

#include <Geom_BSplineSurface.hxx>

#include <BRep_Builder.hxx>
#include <BRepBuilderAPI.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Solid.hxx>

#include <BRepLib.hxx>


#pragma comment(lib, "TKBRep.lib")
#pragma comment(lib, "TKernel.lib")
#pragma comment(lib, "TKG3d.lib")
#pragma comment(lib, "TKG2d.lib")
#pragma comment(lib, "TKGeomBase.lib")
#pragma comment(lib, "TKMath.lib")
#pragma comment(lib, "TKTopAlgo.lib")


TopoDS_Edge convert_opennurbs_trim_to_occt_edge(ON_BrepTrim* trim, Handle_Geom_BSplineSurface& occt_surface);

TopoDS_Wire convert_opennurbs_loop_to_occt_wire(const ON_BrepLoop* loop, Handle_Geom_BSplineSurface& occt_surface);

TopoDS_Face convert_opennurbs_face_to_occt_face(const ON_BrepFace& face);

TopoDS_Compound convert_opennurbs_solid_to_occt_solid(const ON_Brep* brep);
