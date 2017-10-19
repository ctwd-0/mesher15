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

#include <BRepBndLib.hxx>

#pragma comment(lib, "TKBRep.lib")
#pragma comment(lib, "TKernel.lib")
#pragma comment(lib, "TKG3d.lib")
#pragma comment(lib, "TKG2d.lib")
#pragma comment(lib, "TKGeomBase.lib")
#pragma comment(lib, "TKMath.lib")
#pragma comment(lib, "TKTopAlgo.lib")
#pragma comment(lib, "TKMesh.lib")


using namespace std;

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

		return occt_edge;
	}
	else {
		cout <<"no 2d curve"<<endl;
		return TopoDS_Edge();
	}
}

TopoDS_Wire convert_opennurbs_loop_to_occt_wire(const ON_BrepLoop* loop, Handle_Geom_BSplineSurface& occt_surface) {
	BRep_Builder occt_brep_builder;

	BRepBuilderAPI_MakeWire occt_wire;
	bool wire_init = false;

	for (int i = 0; i < loop->TrimCount(); i++) {
		if (loop->Trim(i)->Edge() != nullptr) {
			auto occt_edge = convert_opennurbs_trim_to_occt_edge(loop->Trim(i), occt_surface);
			occt_wire.Add(occt_edge);
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

static bool IsLayerVisible(const ON_ObjectArray<ON_Layer>& layer_table, int layer_index)
{
	bool rc = false;
	// Validate the layer index
	if (layer_index >= 0 && layer_index < layer_table.Count())
	{
		// Get the layer
		const ON_Layer& layer = layer_table[layer_index];
		// Get the layer's visibility
		rc = layer.IsVisible();
		// If the layer is visible, see if the layer has a parent. If so,
		// check to see if the layer's parent is visible. If not, then
		// the layer is also not visible.
		if (rc && ON_UuidIsNotNil(layer.m_parent_layer_id))
		{
			int i, layer_count = layer_table.Count();
			for (i = 0; i < layer_count; i++)
			{
				if (0 == ON_UuidCompare(layer.m_parent_layer_id, layer_table[i].m_layer_id))
					return IsLayerVisible(layer_table, i); // recursive
			}
		}
	}
	return rc;
}

static bool IsModelObjectVisible(const ONX_Model& model, const ONX_Model_Object& model_object)
{
	bool rc = false;
	switch (model_object.m_attributes.Mode())
	{
	case ON::normal_object:
	case ON::idef_object:
	case ON::locked_object:
	{
		// Get the object's layer
		int layer_index = model_object.m_attributes.m_layer_index;
		// See if the layer is visible
		rc = IsLayerVisible(model.m_layer_table, layer_index);
	}
	break;
	}
	return rc;
}

int main() {

	string case_name;

	case_name = "full";

	ON::Begin();

	ONX_Model model;
	FILE* archive_fp = ON::OpenFile((case_name + ".3dm").c_str(), "rb");
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

		if (!IsModelObjectVisible(model, object)) {
			continue;
		}
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
	
	BRepTools::Write(occt_compound_all, (case_name + ".brep").c_str());
	

	Bnd_Box b;

	//BRepBndLib::Add(occt_compound_all, b);
	//B.Get(Xmin, Ymin, Zmin, Xmax, Ymax, Zmax);

	//BRepMesh_FastDiscret::Parameters p;
	//auto occt_discret = BRepMesh_FastDiscret(b, p);
	//occt_discret.Perform(occt_compound_all);

	auto occt_mesh_10 = BRepMesh_IncrementalMesh(occt_compound_all, 10, false, 5, true, false);

	if (occt_mesh_10.IsDone()) {
		cout << "mesh done" << endl;
	}

	vector<gp_Pnt> vs;
	vector<vector<int>> fs;

	int counter = 0;
	for (TopExp_Explorer face_exp(occt_compound_all, TopAbs_FACE); face_exp.More(); face_exp.Next()) {
		TopLoc_Location location;
		auto &face = TopoDS::Face(face_exp.Current());

		auto triangulation = BRep_Tool::Triangulation(face, location);

		counter++;
		if (!triangulation.IsNull()) {
			append_mesh(vs, fs, triangulation);
		}
		else {
			cout <<counter << " is null" << endl;
		}
	}
	
	write_obj(vs, fs, (case_name + ".obj").c_str());

	ON::End();
	system("pause");
}