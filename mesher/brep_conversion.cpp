#include "brep_conversion.h"

TopoDS_Edge convert_opennurbs_trim_to_occt_edge(ON_BrepTrim* trim, Handle_Geom_BSplineSurface& occt_surface) {
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
		cout << "no 2d curve" << endl;
		return TopoDS_Edge();
	}
}

TopoDS_Wire convert_opennurbs_loop_to_occt_wire(const ON_BrepLoop* loop, Handle_Geom_BSplineSurface& occt_surface) {
	BRepBuilderAPI_MakeWire occt_wire;
	for (int i = 0; i < loop->TrimCount(); i++) {
		auto occt_edge = convert_opennurbs_trim_to_occt_edge(loop->Trim(i), occt_surface);
		occt_wire.Add(occt_edge);
	}

	return occt_wire;
}

TopoDS_Face convert_opennurbs_face_to_occt_face(const ON_BrepFace& face) {
	BRep_Builder brep_builder;

	auto occt_surface = convert_opennurbs_surface_to_occt_surface(face.NurbsSurface());
	TopoDS_Face occt_face;
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
			brep_builder.Add(occt_face, occt_wire);
		}
	}
	else {
		for (int i = 0; i < face.LoopCount(); i++) {
			auto occt_wire = convert_opennurbs_loop_to_occt_wire(face.Loop(i), occt_surface);
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
