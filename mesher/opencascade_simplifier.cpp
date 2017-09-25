#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <time.h>

#include "folder.h"

#define ON_DLL_IMPORTS

#include "opennurbs\opennurbs_dynamic_linking.h"

#include <opennurbs.h>

#include <TColgp_Array2OfPnt.hxx>
#include <TColStd_Array2OfReal.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <BSplCLib.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GeomTools_SurfaceSet.hxx>
#include <Precision.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeShell.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepTools.hxx>
#include <Geom_BSplineCurve.hxx>
#include <GeomTools_CurveSet.hxx>

#include <BRep_Builder.hxx>

////#pragma comment(lib, "TKBin.lib")
////#pragma comment(lib, "TKBinL.lib")
////#pragma comment(lib, "TKBinTObj.lib")
////#pragma comment(lib, "TKBinXCAF.lib")
////#pragma comment(lib, "TKBO.lib")
////#pragma comment(lib, "TKBool.lib")
//#pragma comment(lib, "TKBRep.lib")
////#pragma comment(lib, "TKCAF.lib")
////#pragma comment(lib, "TKCDF.lib")
////#pragma comment(lib, "TKD3DHost.lib")
////#pragma comment(lib, "TKDCAF.lib")
////#pragma comment(lib, "TKDraw.lib")
//#pragma comment(lib, "TKernel.lib")
////#pragma comment(lib, "TKFeat.lib")
////#pragma comment(lib, "TKFillet.lib")
////#pragma comment(lib, "TKG2d.lib")
//#pragma comment(lib, "TKG3d.lib")
////#pragma comment(lib, "TKGeomAlgo.lib")
//#pragma comment(lib, "TKGeomBase.lib")
////#pragma comment(lib, "TKHLR.lib")
////#pragma comment(lib, "TKIGES.lib")
////#pragma comment(lib, "TKIVtk.lib")
////#pragma comment(lib, "TKIVtkDraw.lib")
////#pragma comment(lib, "TKLCAF.lib")
//#pragma comment(lib, "TKMath.lib")
////#pragma comment(lib, "TKMesh.lib")
////#pragma comment(lib, "TKMeshVS.lib")
////#pragma comment(lib, "TKOffset.lib")
////#pragma comment(lib, "TKOpenGl.lib")
////#pragma comment(lib, "TKPrim.lib")
////#pragma comment(lib, "TKQADraw.lib")
////#pragma comment(lib, "TKService.lib")
////#pragma comment(lib, "TKShHealing.lib")
////#pragma comment(lib, "TKStd.lib")
////#pragma comment(lib, "TKStdL.lib")
////#pragma comment(lib, "TKSTEP.lib")
////#pragma comment(lib, "TKSTEP209.lib")
////#pragma comment(lib, "TKSTEPAttr.lib")
////#pragma comment(lib, "TKSTEPBase.lib")
////#pragma comment(lib, "TKSTL.lib")
////#pragma comment(lib, "TKTObj.lib")
////#pragma comment(lib, "TKTObjDRAW.lib")
//#pragma comment(lib, "TKTopAlgo.lib")
////#pragma comment(lib, "TKTopTest.lib")
////#pragma comment(lib, "TKV3d.lib")
////#pragma comment(lib, "TKVCAF.lib")
////#pragma comment(lib, "TKViewerTest.lib")
////#pragma comment(lib, "TKVRML.lib")
////#pragma comment(lib, "TKXCAF.lib")
////#pragma comment(lib, "TKXDEDRAW.lib")
////#pragma comment(lib, "TKXDEIGES.lib")
////#pragma comment(lib, "TKXDESTEP.lib")
////#pragma comment(lib, "TKXMesh.lib")
////#pragma comment(lib, "TKXml.lib")
////#pragma comment(lib, "TKXmlL.lib")
////#pragma comment(lib, "TKXmlTObj.lib")
////#pragma comment(lib, "TKXmlXCAF.lib")
////#pragma comment(lib, "TKXSBase.lib")
////#pragma comment(lib, "TKXSDRAW.lib")

#pragma comment(lib, "TKBRep.lib")
#pragma comment(lib, "TKernel.lib")
#pragma comment(lib, "TKG3d.lib")
#pragma comment(lib, "TKGeomBase.lib")
#pragma comment(lib, "TKMath.lib")
#pragma comment(lib, "TKTopAlgo.lib")



///**
//* @breif Convert OpenNURBS NURBS curve to OpenCASCADE Geom_BSplineCurve.
//* @param [in] theCurve opennurbs nurbs curve;
//* @param [in] theBRepFile the curve is in brep file of opencascade;
//* @note pay attention to the knots of opennurbs nurbs curve/surface.
//*/
//void ConvertCurve(const ON_NurbsCurve& theCurve, const std::string& theBRepFile)
//{
//	TColgp_Array1OfPnt aPoles(1, theCurve.CVCount());
//	TColStd_Array1OfReal aWeights(1, theCurve.CVCount());
//
//	TColStd_Array1OfReal aKnotSequence(1, theCurve.KnotCount() + 2);
//
//	bool IsRational = theCurve.IsRational();
//	bool IsPeriodic = (theCurve.IsPeriodic()) ? true : false;
//
//	// Control point and its weight.
//	for (int i = 0; i < theCurve.CVCount(); ++i)
//	{
//		if (IsRational)
//		{
//			ON_4dPoint aPole;
//
//			theCurve.GetCV(i, aPole);
//
//			aPoles.SetValue(i + 1, gp_Pnt(aPole.x / aPole.w, aPole.y / aPole.w, aPole.z / aPole.w));
//			aWeights.SetValue(i + 1, aPole.w);
//		}
//		else
//		{
//			ON_3dPoint aPole;
//
//			theCurve.GetCV(i, aPole);
//
//			aPoles.SetValue(i + 1, gp_Pnt(aPole.x, aPole.y, aPole.z));
//		}
//	}
//
//	// Knot vector and its multiplicity.
//	for (int i = 0; i < theCurve.KnotCount(); ++i)
//	{
//		aKnotSequence.SetValue(i + 2, theCurve.Knot(i));
//	}
//
//	aKnotSequence.SetValue(aKnotSequence.Lower(), theCurve.Knot(0));
//	aKnotSequence.SetValue(aKnotSequence.Upper(), theCurve.Knot(theCurve.KnotCount() - 1));
//
//	TColStd_Array1OfReal aKnots(1, BSplCLib::KnotsLength(aKnotSequence, IsPeriodic));
//	TColStd_Array1OfInteger aMultiplicities(1, aKnots.Upper());
//
//	BSplCLib::Knots(aKnotSequence, aKnots, aMultiplicities);
//
//	Handle_Geom_BSplineCurve aBSplineCurve = new Geom_BSplineCurve(
//		aPoles, aWeights, aKnots, aMultiplicities,
//		theCurve.Degree(), theCurve.IsPeriodic());
//
//	GeomTools_CurveSet::PrintCurve(aBSplineCurve, std::cout);
//
//	TopoDS_Edge anEdge = BRepBuilderAPI_MakeEdge(aBSplineCurve);
//
//	BRepTools::Write(anEdge, theBRepFile.c_str());
//}
//
///**
//* @breif Convert OpenNURBS NURBS surface to OpenCASCADE Geom_BSplineSurface.
//* @param [in] theSurface opennurbs nurbs surface;
//* @param [in] theBRepFile the surface is in the brep file of opencascade;
//* @note pay attention to the knots of opennurbs nurbs curve/surface.
//*/
//void ConvertSurface(const ON_NurbsSurface& theSurface, const std::string& theBRepFile)
//{
//	TColgp_Array2OfPnt aPoles(1, theSurface.CVCount(0), 1, theSurface.CVCount(1));
//	TColStd_Array2OfReal aWeights(1, theSurface.CVCount(0), 1, theSurface.CVCount(1));
//
//	TColStd_Array1OfReal aUKnotSequence(1, theSurface.KnotCount(0) + 2);
//	TColStd_Array1OfReal aVKnotSequence(1, theSurface.KnotCount(1) + 2);
//
//	bool IsRational = theSurface.IsRational();
//	bool IsUPeriodic = (theSurface.IsPeriodic(0)) ? true : false;
//	bool IsVPeriodic = (theSurface.IsPeriodic(1)) ? true : false;
//
//	// control point and its weight.
//	for (int i = 0; i < theSurface.CVCount(0); ++i)
//	{
//		for (int j = 0; j < theSurface.CVCount(1); ++j)
//		{
//			if (IsRational)
//			{
//				ON_4dPoint aPole;
//
//				theSurface.GetCV(i, j, aPole);
//
//				aPoles.SetValue(i + 1, j + 1, gp_Pnt(aPole.x / aPole.w, aPole.y / aPole.w, aPole.z / aPole.w));
//				aWeights.SetValue(i + 1, j + 1, aPole.w);
//			}
//			else
//			{
//				ON_3dPoint aPole;
//
//				theSurface.GetCV(i, j, aPole);
//
//				aPoles.SetValue(i + 1, j + 1, gp_Pnt(aPole.x, aPole.y, aPole.z));
//			}
//
//		}
//	}
//
//	// Knot vector and its multiplicity.
//	for (int i = 0; i < theSurface.KnotCount(0); ++i)
//	{
//		aUKnotSequence.SetValue(i + 2, theSurface.Knot(0, i));
//	}
//
//	aUKnotSequence.SetValue(aUKnotSequence.Lower(), theSurface.Knot(0, 0));
//	aUKnotSequence.SetValue(aUKnotSequence.Upper(), theSurface.Knot(0, theSurface.KnotCount(0) - 1));
//
//	TColStd_Array1OfReal aUKnots(1, BSplCLib::KnotsLength(aUKnotSequence, IsUPeriodic));
//	TColStd_Array1OfInteger aUMults(1, aUKnots.Upper());
//
//	BSplCLib::Knots(aUKnotSequence, aUKnots, aUMults);
//
//	for (int j = 0; j < theSurface.KnotCount(1); ++j)
//	{
//		aVKnotSequence.SetValue(j + 2, theSurface.Knot(1, j));
//	}
//
//	aVKnotSequence.SetValue(aVKnotSequence.Lower(), theSurface.Knot(1, 0));
//	aVKnotSequence.SetValue(aVKnotSequence.Upper(), theSurface.Knot(1, theSurface.KnotCount(1) - 1));
//
//	TColStd_Array1OfReal aVKnots(1, BSplCLib::KnotsLength(aVKnotSequence, IsVPeriodic));
//	TColStd_Array1OfInteger aVMults(1, aVKnots.Upper());
//
//	BSplCLib::Knots(aVKnotSequence, aVKnots, aVMults);
//
//	Handle_Geom_BSplineSurface aBSplineSurface = new Geom_BSplineSurface(aPoles,
//		aWeights, aUKnots, aVKnots, aUMults, aVMults,
//		theSurface.Degree(0), theSurface.Degree(1),
//		IsUPeriodic, IsVPeriodic);
//
//	GeomTools_SurfaceSet::PrintSurface(aBSplineSurface, std::cout);
//
//	TopoDS_Face aFace = BRepBuilderAPI_MakeFace(aBSplineSurface, Precision::Confusion());
//
//	BRepTools::Write(aFace, theBRepFile.c_str());
//}
//


Handle_Geom_BSplineCurve convert_opennurbs_curve_to_occt_curve(const ON_NurbsCurve* curve) {
	TColgp_Array1OfPnt poles(1, curve->CVCount());
	TColStd_Array1OfReal weights(1, curve->CVCount());

	TColStd_Array1OfReal knot_sequence(1, curve->KnotCount() + 2);

	// Control point and its weight.
	for (int i = 0; i < curve->CVCount(); ++i) {
		if (curve->IsRational()) {
			ON_4dPoint aPole;
			curve->GetCV(i, aPole);
			poles.SetValue(i + 1, gp_Pnt(aPole.x / aPole.w, aPole.y / aPole.w, aPole.z / aPole.w));
			weights.SetValue(i + 1, aPole.w);
		}
		else {
			ON_3dPoint aPole;
			curve->GetCV(i, aPole);
			poles.SetValue(i + 1, gp_Pnt(aPole.x, aPole.y, aPole.z));
			weights.SetValue(i + 1, 1.);
		}
	}

	// Knot vector and its multiplicity.
	for (int i = 0; i < curve->KnotCount(); ++i) {
		knot_sequence.SetValue(i + 2, curve->Knot(i));
	}

	knot_sequence.SetValue(knot_sequence.Lower(), curve->Knot(0));
	knot_sequence.SetValue(knot_sequence.Upper(), curve->Knot(curve->KnotCount() - 1));

	TColStd_Array1OfReal knots(1, BSplCLib::KnotsLength(knot_sequence, curve->IsPeriodic()));
	TColStd_Array1OfInteger muls(1, knots.Upper());

	BSplCLib::Knots(knot_sequence, knots, muls);

	Handle_Geom_BSplineCurve occt_curve = new Geom_BSplineCurve(
		poles, weights, knots, muls,
		curve->Degree(), curve->IsPeriodic());

	return occt_curve;
}

Handle_Geom_BSplineSurface convert_opennurbs_surface_to_occt_surface(const ON_NurbsSurface* surface) {
	TColgp_Array2OfPnt poles(1, surface->CVCount(0), 1, surface->CVCount(1));
	TColStd_Array2OfReal weights(1, surface->CVCount(0), 1, surface->CVCount(1));

	TColStd_Array1OfReal u_knot_sequence(1, surface->KnotCount(0) + 2);
	TColStd_Array1OfReal v_knot_sequence(1, surface->KnotCount(1) + 2);

	bool is_u_periodic = (bool)surface->IsPeriodic(0);
	bool is_v_periodic = (bool)surface->IsPeriodic(1);

	// control point and its weight.
	for (int i = 0; i < surface->CVCount(0); ++i) {
		for (int j = 0; j < surface->CVCount(1); ++j) {
			if (surface->IsRational()) {
				ON_4dPoint aPole;
				surface->GetCV(i, j, aPole);
				poles.SetValue(i + 1, j + 1, gp_Pnt(aPole.x / aPole.w, aPole.y / aPole.w, aPole.z / aPole.w));
				weights.SetValue(i + 1, j + 1, aPole.w);
			}
			else {
				ON_3dPoint aPole;
				surface->GetCV(i, j, aPole);
				poles.SetValue(i + 1, j + 1, gp_Pnt(aPole.x, aPole.y, aPole.z));
				weights.SetValue(i + 1, j + 1, 1);
			}
		}
	}

	// Knot vector and its multiplicity.
	for (int i = 0; i < surface->KnotCount(0); ++i) {
		u_knot_sequence.SetValue(i + 2, surface->Knot(0, i));
	}
	u_knot_sequence.SetValue(u_knot_sequence.Lower(), surface->Knot(0, 0));
	u_knot_sequence.SetValue(u_knot_sequence.Upper(), surface->Knot(0, surface->KnotCount(0) - 1));
	TColStd_Array1OfReal u_knots(1, BSplCLib::KnotsLength(u_knot_sequence, is_u_periodic));
	TColStd_Array1OfInteger u_mults(1, u_knots.Upper());
	BSplCLib::Knots(u_knot_sequence, u_knots, u_mults);

	for (int j = 0; j < surface->KnotCount(1); ++j) {
		v_knot_sequence.SetValue(j + 2, surface->Knot(1, j));
	}
	v_knot_sequence.SetValue(v_knot_sequence.Lower(), surface->Knot(1, 0));
	v_knot_sequence.SetValue(v_knot_sequence.Upper(), surface->Knot(1, surface->KnotCount(1) - 1));
	TColStd_Array1OfReal v_knots(1, BSplCLib::KnotsLength(v_knot_sequence, is_v_periodic));
	TColStd_Array1OfInteger v_mults(1, v_knots.Upper());
	BSplCLib::Knots(v_knot_sequence, v_knots, v_mults);

	Handle_Geom_BSplineSurface occt_surface = new Geom_BSplineSurface(poles,
		weights, u_knots, v_knots, u_mults, v_mults,
		surface->Degree(0), surface->Degree(1),
		is_u_periodic, is_v_periodic);

	return occt_surface;
}

TopoDS_Edge convert_opennurbs_edge_to_occt_edge(ON_BrepEdge* edge) {
	BRep_Builder occt_brep_builder;
	auto occt_curve = convert_opennurbs_curve_to_occt_curve(edge->NurbsCurve());

	TopoDS_Edge occt_edge = BRepBuilderAPI_MakeEdge(occt_curve);

	auto p0 = edge->Vertex(0)->Point();
	auto p1 = edge->Vertex(1)->Point();

	TopoDS_Vertex occt_v0 = BRepBuilderAPI_MakeVertex(gp_Pnt(p0.x, p0.y, p0.z));
	TopoDS_Vertex occt_v1 = BRepBuilderAPI_MakeVertex(gp_Pnt(p1.x, p1.y, p1.z));
	occt_brep_builder.Add(occt_edge, occt_v0);
	occt_brep_builder.Add(occt_edge, occt_v1);

	return occt_edge;
}

TopoDS_Edge convert_opennurbs_trim_to_occt_edge(ON_BrepTrim* trim) {
	auto edge = trim->Edge();
	if (edge != nullptr) {
		return convert_opennurbs_edge_to_occt_edge(edge);
	}
	else {
		cout << "trim type:"<< trim->m_type <<endl;
		return TopoDS_Edge();
	}
}

TopoDS_Wire convert_opennurbs_loop_to_occt_wire(const ON_BrepLoop* loop) {
	BRep_Builder occt_brep_builder;

	TopoDS_Wire occt_wire;
	bool wire_init = false;

	for (int i = 0; i < loop->TrimCount(); i++) {
		if (loop->Trim(i)->Edge() != nullptr) {
			auto occt_edge = convert_opennurbs_trim_to_occt_edge(loop->Trim(i));
			if (wire_init) {
				occt_brep_builder.Add(occt_wire, occt_edge);
			}
			else {
				occt_wire = BRepBuilderAPI_MakeWire(occt_edge);
				wire_init = true;
			}
			//occt_brep_builder.Add(occt_wire, occt_edge);
		}
		else {
			cout << "trim type:" << loop->Trim(i)->m_type << endl;
		}
	}

	return occt_wire;
}

TopoDS_Face convert_opennurbs_face_to_occt_face(const ON_BrepFace& face) {
	BRep_Builder occt_brep_builder;

	//这里遇到个问题。一定要先加外环。
	//如果直接遍历环往上加，会出现面裁剪出问题的情况。
	auto occt_surface = convert_opennurbs_surface_to_occt_surface(face.NurbsSurface());

	if (face.OuterLoop()) {
		auto outer_loop = face.OuterLoop();
		auto occt_outer_wire = convert_opennurbs_loop_to_occt_wire(outer_loop);

		TopoDS_Face occt_face = BRepBuilderAPI_MakeFace(occt_surface, occt_outer_wire);

		for (int i = 0; i < face.LoopCount(); i++) {
			if (face.Loop(i) == outer_loop) {
				continue;
			}
			auto occt_wire = convert_opennurbs_loop_to_occt_wire(face.Loop(i));
			occt_face = BRepBuilderAPI_MakeFace(occt_surface, occt_wire);
		}

		return occt_face;
	}
	else {

		TopoDS_Face occt_face = BRepBuilderAPI_MakeFace(occt_surface, Precision::Confusion());
		for (int i = 0; i < face.LoopCount(); i++) {
			auto occt_wire = convert_opennurbs_loop_to_occt_wire(face.Loop(i));
			occt_face = BRepBuilderAPI_MakeFace(occt_surface, occt_wire);
		}

		return occt_face;
	}
}

TopoDS_Compound convert_opennurbs_solid_to_occt_solid(const ON_Brep* brep) {
	BRep_Builder occt_brep_builder;
	TopoDS_Compound occt_compound;
	occt_brep_builder.MakeCompound(occt_compound);

	auto & faces = brep->m_F;
	for (int i = 0; i < faces.Count(); i++) {
		auto occt_face = convert_opennurbs_face_to_occt_face(faces[i]);
		occt_brep_builder.Add(occt_compound, occt_face);
	}

	return occt_compound;
}

int main() {
	ON::Begin();

	ONX_Model model;
	FILE* archive_fp = ON::OpenFile("full.3dm", "rb");
	ON_BinaryFile archive(ON::read3dm, archive_fp);
	bool rc = model.Read(archive);
	if (rc) {
		cout << "reading 3dm file success" << endl;
	}
	else {
		cout << "failed" << endl;
		return -1;
	}
	ON::CloseFile(archive_fp);

	BRep_Builder occt_brep_builder;
	TopoDS_Compound occt_compound_all;
	occt_brep_builder.MakeCompound(occt_compound_all);

	for (int obj_id = 0; obj_id < model.m_object_table.Count(); obj_id++) {
		ONX_Model_Object object = model.m_object_table[obj_id];

		const ON_Brep* brep = ON_Brep::Cast(object.m_object);
		if (brep != nullptr) {

			//if (brep->IsSolid()) {
				auto occt_compound = convert_opennurbs_solid_to_occt_solid(brep);

			//}
			//else {
			//	cout << "none solid." <<endl;
			//}
				occt_brep_builder.Add(occt_compound_all, occt_compound);
		}
	}

	BRepTools::Write(occt_compound_all,"full.release.brep");
}