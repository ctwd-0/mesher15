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
#include <TColgp_Array2OfPnt2d.hxx>
#include <TColStd_Array2OfReal.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <BSplCLib.hxx>
#include <Geom2d_BSplineCurve.hxx>
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

#include<BRepLib.hxx>

#include <gp_Pnt.hxx> 
#include <BRepMesh.hxx>
#include <BRep_Tool.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepMesh_FastDiscret.hxx>
#include <BRepMesh_FastDiscretFace.hxx>
#include <TopExp_Explorer.hxx>
#include <Poly_Triangulation.hxx>

#include <BRepProj_Projection.hxx>

#include <BRepBndLib.hxx>

#include <ProjLib.hxx>

#include <ProjLib_CompProjectedCurve.hxx>

#include <Geom2dAPI_ProjectPointOnCurve.hxx>


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
#pragma comment(lib, "TKG2d.lib")
#pragma comment(lib, "TKGeomBase.lib")
#pragma comment(lib, "TKMath.lib")
#pragma comment(lib, "TKTopAlgo.lib")

#pragma comment(lib, "TKMesh.lib") 

using namespace std;

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


Handle_Geom2d_BSplineCurve convert_opennurbs_2d_curve_to_occt_2d_curve(const ON_NurbsCurve* curve) {
	TColgp_Array1OfPnt2d poles(1, curve->CVCount());
	TColStd_Array1OfReal weights(1, curve->CVCount());

	TColStd_Array1OfReal knot_sequence(1, curve->KnotCount() + 2);

	// Control point and its weight.
	for (int i = 0; i < curve->CVCount(); ++i) {
		if (curve->IsRational()) {
			ON_4dPoint aPole;
			curve->GetCV(i, aPole);
			poles.SetValue(i + 1, gp_Pnt2d(aPole.x / aPole.w, aPole.y / aPole.w));
			weights.SetValue(i + 1, aPole.w);
		}
		else {
			ON_3dPoint aPole;
			curve->GetCV(i, aPole);
			poles.SetValue(i + 1, gp_Pnt2d(aPole.x, aPole.y));
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

	Handle_Geom2d_BSplineCurve occt_curve = new Geom2d_BSplineCurve(
		poles, weights, knots, muls,
		curve->Degree(), curve->IsPeriodic());

	return occt_curve;
}

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

//TopoDS_Edge convert_opennurbs_edge_to_occt_edge(ON_BrepEdge* edge, bool reverse) {
//	BRep_Builder brep_builder;
//	
//	auto occt_curve = convert_opennurbs_curve_to_occt_curve(edge->NurbsCurve());
//
//	TopoDS_Edge occt_edge = BRepBuilderAPI_MakeEdge(occt_curve);
//
//	auto p0 = edge->Vertex(0)->Point();
//	auto p1 = edge->Vertex(1)->Point();
//
//
//
//	//cout << edge->ProxyCurveIsReversed() << endl;
//
//	//马勒戈壁的。
//	//loop中边的顺序可能会有反的。
//	//而且看起来一定会有。
//	//原因是在trim层，可能把3d curve反转。
//	//问题仍然没有解决。
//
//	TopoDS_Vertex occt_v0 = BRepBuilderAPI_MakeVertex(gp_Pnt(p0.x, p0.y, p0.z));
//	TopoDS_Vertex occt_v1 = BRepBuilderAPI_MakeVertex(gp_Pnt(p1.x, p1.y, p1.z));
//	brep_builder.Add(occt_edge, occt_v0);
//	brep_builder.Add(occt_edge, occt_v1);
//
//	//cout << p0.x+300 << " " << p0.y+2200 << " " << 0 << endl;
//	//cout << p1.x+300 << " " << p1.y+2200 << " " << 0 << endl;
//	//cout<<endl;
//
//
//	if (reverse) {
//		occt_edge.Reverse();
//	}
//
//	
//	return occt_edge;
//}

//static int curve_proxy_c = 0;
//static int on_surface_c = 0;
//static int line_c = 0;
//static int arc_c = 0;
//static int polyline_c = 0;
//static int polycurve_c = 0;
//static int nurbs_c = 0;

//int opennurbs_curve2d_tester(ON_Curve* curve_2d) {
//
//	//auto curve_proxy = dynamic_cast<ON_CurveProxy*>(curve_2d); //checked
//	//auto curve_on_surface = dynamic_cast<ON_CurveOnSurface*>(curve_2d); //checked
//	//auto line = dynamic_cast<ON_LineCurve*>(curve_2d); //checked
//	//auto arc = dynamic_cast<ON_ArcCurve*>(curve_2d); //checked
//	//auto polyline = dynamic_cast<ON_PolylineCurve*>(curve_2d); //checked
//	//auto polycurve = dynamic_cast<ON_PolyCurve*>(curve_2d); //checked
//	//auto nurbs = dynamic_cast<ON_NurbsCurve*>(curve_2d); //checked
//
//	//if (curve_proxy != nullptr) {
//	//	curve_proxy_c++;
//	//}
//	//if (curve_on_surface != nullptr) {
//	//	on_surface_c++;
//	//}
//	//if (line != nullptr) {
//	//	line_c++;
//	//}
//	//if (arc != nullptr) {
//	//	arc_c++;
//	//}
//	//if (polyline != nullptr) {
//	//	polyline_c++;
//	//}
//	//if (polycurve != nullptr) {
//	//	polycurve_c++;
//	//}
//	//if (nurbs != nullptr) {
//	//	nurbs_c++;
//	//}
//
//
//
//	return 0;
//}

TopoDS_Edge convert_opennurbs_trim_to_occt_edge(ON_BrepTrim* trim, Handle_Geom_BSplineSurface& occt_surface) {
	BRep_Builder brep_builder;

	auto brep = trim->Brep();
	auto edge = trim->Edge();



	if (trim->m_c2i != -1) {

		auto interval = trim->Domain();
		auto curve_2d = brep->m_C2[trim->m_c2i];
		auto nurbs = dynamic_cast<ON_NurbsCurve*>(curve_2d); //checked
		if (nurbs == nullptr) {
			nurbs = curve_2d->NurbsCurve();
		}
		
		auto occt_bspline_2d = convert_opennurbs_2d_curve_to_occt_2d_curve(nurbs);
		
		TopoDS_Edge occt_edge;

		occt_edge = BRepBuilderAPI_MakeEdge(occt_bspline_2d, occt_surface, interval.Min(), interval.Max());

		//if (trim->m_bRev3d) {
		//	occt_edge.Reverse();
		//}

		return occt_edge;
	}
	else {
		cout <<"no 2d curve"<<endl;
		return TopoDS_Edge();
	}


	//if (trim->m_c2i != -1) {
	//	auto curve_2d = brep->m_C2[trim->m_c2i];
	//	auto nurbs = dynamic_cast<ON_NurbsCurve*>(curve_2d); //checked
	//	if (nurbs == nullptr) {
	//		nurbs = curve_2d->NurbsCurve();
	//	}

	//	auto occt_bspline_2d = convert_opennurbs_2d_curve_to_occt_2d_curve(nurbs);
	//	
	//	auto p0 = trim->Vertex(0)->Point();
	//	auto p1 = trim->Vertex(1)->Point();

	//	gp_Pnt occt_p0(p0.x, p0.y, p0.z);
	//	gp_Pnt occt_p1(p1.x, p1.y, p1.z);

	//	TopoDS_Vertex occt_v0, occt_v1;

	//	brep_builder.MakeVertex(occt_v0, occt_p0, Precision::Confusion());
	//	brep_builder.MakeVertex(occt_v1, occt_p0, Precision::Confusion());

	//	TopoDS_Edge occt_edge;

	//	if (edge != nullptr) {
	//		auto occt_bspline = convert_opennurbs_curve_to_occt_curve(edge->NurbsCurve());
	//		brep_builder.MakeEdge(occt_edge, occt_bspline, Precision::Confusion());
	//	}
	//	else {
	//		brep_builder.MakeEdge(occt_edge);
	//	}

	//	brep_builder.UpdateEdge(occt_edge, occt_bspline_2d, occt_face, Precision::Confusion());
	//	brep_builder.Add(occt_edge, occt_v0);
	//	brep_builder.Add(occt_edge, occt_v1);
	//	
	//	brep_builder.Range(occt_edge, interval.Min(), interval.Max());

	//	//cout<<endl;
	//	//cout<< trim->m_bRev3d <<endl;
	//	//cout<< interval.Min()  << " "<< interval.Max() <<endl;
	//	
	//	if (trim->m_bRev3d) {
	//		occt_edge.Reverse();
	//		brep_builder.Range(occt_edge, -interval.Max(), -interval.Min(), true);
	//	}

	//	return occt_edge;
	//}
	//else {
	//	cout <<"no 2d curve"<<endl;
	//	return TopoDS_Edge();
	//}
}

TopoDS_Wire convert_opennurbs_loop_to_occt_wire(const ON_BrepLoop* loop, Handle_Geom_BSplineSurface& occt_surface) {
	BRep_Builder occt_brep_builder;

	BRepBuilderAPI_MakeWire occt_wire;
	bool wire_init = false;

	for (int i = 0; i < loop->TrimCount(); i++) {
		if (loop->Trim(i)->Edge() != nullptr) {
			auto occt_edge = convert_opennurbs_trim_to_occt_edge(loop->Trim(i), occt_surface);
			//if (wire_init) {
			//	occt_wire.Add(occt_edge);
			//}
			//else {
			//	occt_wire = BRepBuilderAPI_MakeWire(occt_edge);
			//	wire_init = true;
			//}
			occt_wire.Add(occt_edge);
			//cout << occt_wire.IsDone() << endl;
			//occt_brep_builder.Add(occt_wire, occt_edge);
		}
		else {
			//这里 trim 的类型没有详细的研究。
			//当前忽略singular的trim。
			//如果之后有时间的话，应该详细研究下保证不出问题。
			//cout << "trim type:" << loop->Trim(i)->m_type << endl;
		}
	}

	return occt_wire;
}

TopoDS_Face convert_opennurbs_face_to_occt_face(const ON_BrepFace& face) {
	BRep_Builder brep_builder;

	//这里遇到个问题。一定要先加外环。
	//如果直接遍历环往上加，会出现面裁剪出问题的情况。
	auto occt_surface = convert_opennurbs_surface_to_occt_surface(face.NurbsSurface());
	TopoDS_Face occt_face;// = BRepBuilderAPI_MakeFace(occt_surface, Precision::Confusion());
	brep_builder.MakeFace(occt_face, occt_surface, Precision::Confusion());

	if (face.OuterLoop()) {
		auto outer_loop = face.OuterLoop();
		auto occt_outer_wire = convert_opennurbs_loop_to_occt_wire(outer_loop, occt_surface);

		brep_builder.Add(occt_face, occt_outer_wire);

		for (int i = 0; i < face.LoopCount(); i++) {
			if (face.Loop(i) == outer_loop) {
				continue;
			}
			auto occt_wire = convert_opennurbs_loop_to_occt_wire(face.Loop(i), occt_surface);
			//occt_face = BRepBuilderAPI_MakeFace(occt_surface, occt_wire);
			brep_builder.Add(occt_face, occt_wire);
		}
	}
	else {
		for (int i = 0; i < face.LoopCount(); i++) {
			auto occt_wire = convert_opennurbs_loop_to_occt_wire(face.Loop(i) , occt_surface);
			//occt_face = BRepBuilderAPI_MakeFace(occt_surface, occt_wire);
			brep_builder.Add(occt_face, occt_wire);
		}
	}

	BRepLib::BuildCurves3d(occt_face);
	return occt_face;
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

void write_obj(vector<gp_Pnt>& vertices, vector<vector<int>>& faces, const char* file_name) {
	FILE* f;
	auto error = fopen_s(&f, file_name, "wb");

	char line[1000] = { 0 };

	for (int i = 0; i < vertices.size(); i++) {
		sprintf(line, "v %lf %lf %lf\n\0", vertices[i].XYZ().X(), vertices[i].XYZ().Y(), vertices[i].XYZ().Z());
		fwrite(line, strlen(line), 1, f);
	}

	for (int i = 0; i < faces.size(); i++) {
		sprintf(line, "f %d %d %d\n\0", faces[i][0] + 1, faces[i][1] + 1, faces[i][2] + 1);
		fwrite(line, strlen(line), 1, f);
	}

	fclose(f);
}

void append_mesh(vector<gp_Pnt>& vs, vector<vector<int>>& ts, const opencascade::handle<Poly_Triangulation> triangulation) {
	auto n_triangles = triangulation->NbTriangles();
	auto n_nodes = triangulation->NbNodes();
	auto triangles = triangulation->Triangles();
	auto nodes = triangulation->Nodes();

	Standard_Integer i1, i2, i3;

	int offset = vs.size();
	for (auto i = 1; i <= n_triangles; i++) {
		auto triangle = triangles.Value(i);
		triangle.Get(i1, i2, i3);
		ts.push_back(vector<int>{i1 + offset - 1, i2 + offset - 1, i3 + offset - 1});
	}

	for (auto i = 1; i <= n_nodes; i++) {
		vs.push_back(nodes.Value(i));
	}
}

int main() {
	ON::Begin();

	ONX_Model model;
	FILE* archive_fp = ON::OpenFile("small12.3dm", "rb");
	ON_BinaryFile archive(ON::read3dm, archive_fp);
	bool rc = model.Read(archive);
	if (rc) {
		cout << "reading 3dm file success" << endl;
	}
	else {
		cout << "failed" << endl;
		system("pause");
		return -1;
	}
	ON::CloseFile(archive_fp);

	//ON_TextLog log;
	//model.Dump(log);

	

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

	cout << "convert done" << endl;

	ON::End();


	//cout << "curve_proxy: " << curve_proxy_c << endl;
	//cout << "on_surface: " << on_surface_c << endl;
	//cout << "line: " << line_c << endl;
	//cout << "arc: " << arc_c << endl;
	//cout << "polyline: " << polyline_c << endl;
	//cout << "polycurve: " << polycurve_c << endl;
	//cout << "nurbs: " << nurbs_c << endl;
	
	BRepTools::Write(occt_compound_all, "small12.brep");
	
	//system("pause");
	//return -2;

	Bnd_Box b;

	BRepBndLib::Add(occt_compound_all, b);
	//B.Get(Xmin, Ymin, Zmin, Xmax, Ymax, Zmax);

	BRepMesh_FastDiscret::Parameters p;
	auto occt_discret = BRepMesh_FastDiscret(b, p);
	occt_discret.Perform(occt_compound_all);

	auto occt_mesh = BRepMesh_IncrementalMesh(occt_compound_all, 1);

	if (occt_mesh.IsDone()) {
		cout << "mesh done" << endl;
	}

	vector<gp_Pnt> vs;
	vector<vector<int>> fs;

	for (TopExp_Explorer face_exp(occt_compound_all, TopAbs_FACE); face_exp.More(); face_exp.Next()) {
		TopLoc_Location location;
		auto &face = TopoDS::Face(face_exp.Current());

		auto triangulation = BRep_Tool::Triangulation(face, location);

		if (!triangulation.IsNull()) {
			append_mesh(vs, fs, triangulation);
		}
		else {

		}
	}
	
	write_obj(vs, fs, "small12.obj");

	ON::End();
	system("pause");
}