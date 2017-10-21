#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <time.h>

using namespace std;

#define ON_DLL_IMPORTS

#include "opennurbs\opennurbs_dynamic_linking.h"

#include "opennurbs_utils.h"

#include "obj_model_output.h"

#include "bspline_conversion.h"

#include "brep_conversion.h"

#include "folder.h"

#include <Precision.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <BRepTools.hxx>

#include <TopExp_Explorer.hxx>
#include <Poly_Triangulation.hxx>

#include <BRepMesh.hxx>
#include <BRepMesh_IncrementalMesh.hxx>

#pragma comment(lib, "TKTopAlgo.lib")
#pragma comment(lib, "TKMesh.lib")


const int MAX_V = 50 * 10000;
const int MAX_F = 100 * 10000;

int main() {

	string case_name;

	string output_prefix = "D:\\garbage\\";

	case_name = "full";


	mkdir((output_prefix + case_name).c_str());
	mkdir((output_prefix + case_name + "_opennurbs").c_str());
	mkdir((output_prefix + case_name + "_groups").c_str());

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


	OpennurbsGroupInfo group_info(model);

	cout <<"number of groups: "<< group_info.group_ids.size() << endl;

	map<int, TopoDS_Compound> breps;

	map<int, RhinoMesh> obj_mesh;

	for (int obj_id = 0; obj_id < model.m_object_table.Count(); obj_id++) {
		ONX_Model_Object object = model.m_object_table[obj_id];

		const ON_Brep* brep = ON_Brep::Cast(object.m_object);

		if (brep == nullptr) {
			continue;
		}

		ON_MeshParameters mp;
		ON_SimpleArray<const ON_Mesh*> meshes;
		int mc = brep->GetMesh(ON::render_mesh, meshes);
		int fc = 0;

		if (mc) {
			obj_mesh[obj_id] = RhinoMesh();
			for (int i = 0; i < meshes.Count(); i++) {
				auto mesh = meshes[i];
				append_mesh(obj_mesh[obj_id].v, obj_mesh[obj_id].f, mesh);
			}
			obj_mesh[obj_id].merge_vertices();
			//obj_mesh[obj_id].write_obj((case_name + "_opennurbs\\obj_" + to_string(obj_id) + ".obj").c_str());
		}

		if (IsModelObjectVisible(model, object)) {
			auto occt_compound = convert_opennurbs_solid_to_occt_solid(brep);
			breps[obj_id] = occt_compound;
		}
	}

	model.Destroy();
	cout << "convert opennurbs to opencascad finished." << endl;

	ON::End();

	map<int, map<int, RhinoMesh>> grp_meshs_1;
	map<int, map<int, OcctMesh>> grp_meshs;

	for (auto grp_id : group_info.group_ids) {
		auto& obj_ids = group_info.map_group_id_objects[grp_id];
		int vs = 0;
		int fs = 0;

		for (auto obj_id : obj_ids) {
			vs += obj_mesh[obj_id].v.size();
			fs += obj_mesh[obj_id].f.size();
		}

		if (vs <= MAX_V && fs < MAX_F) {
			//use original mesh;
			grp_meshs_1[grp_id] = map<int, RhinoMesh>();
			for (auto obj_id : obj_ids) {
				grp_meshs_1[grp_id][obj_id] = obj_mesh[obj_id];
			}
		}
		else {
			//generate mesh;
			double line_tolerance = 1;
			double angle_tolerance = 0.5;

			double multiplier = fs / MAX_F;
			line_tolerance *= multiplier;
			angle_tolerance *= multiplier;

			for (auto obj_id : obj_ids) {
				auto x = BRepMesh_IncrementalMesh(breps[obj_id],
					line_tolerance, false, angle_tolerance, true, false);
				grp_meshs[grp_id][obj_id] = OcctMesh();
				for (TopExp_Explorer face_exp(breps[obj_id], TopAbs_FACE); face_exp.More(); face_exp.Next()) {
					TopLoc_Location location;

					auto &face = TopoDS::Face(face_exp.Current());
					auto triangulation = BRep_Tool::Triangulation(face, location);

					if (!triangulation.IsNull()) {
						grp_meshs[grp_id][obj_id].append_mesh(triangulation);
					}
				}
			}
		}
	}


	//for (auto &x : breps) {
	//	auto obj_id = x.first;
	//	auto& brep = x.second;

	//	auto x  = BRepMesh_IncrementalMesh(brep, 1);

	//	vector<gp_Pnt> vs;
	//	vector<vector<int>> fs;
	//	for (TopExp_Explorer face_exp(brep, TopAbs_FACE); face_exp.More(); face_exp.Next()) {
	//		TopLoc_Location location;
	//		auto &face = TopoDS::Face(face_exp.Current());

	//		auto &triangulation = BRep_Tool::Triangulation(face, location);

	//		if (!triangulation.IsNull()) {
	//			append_mesh(vs, fs, triangulation);
	//		}
	//		else {
	//		}
	//	}
	//	merge_vertices(vs, fs);
	//	
	//	write_obj(vs, fs, (output_prefix + case_name + "\\obj_" + to_string(obj_id) + ".obj").c_str());
	//}
	
	system("pause");

	return 0;
}
