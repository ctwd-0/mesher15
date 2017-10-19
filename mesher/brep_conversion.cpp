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
