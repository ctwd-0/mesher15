#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <time.h>

#include <omp.h>

#include <math.h>

using namespace std;

#define ON_DLL_IMPORTS

#include "opennurbs\opennurbs_dynamic_linking.h"

#include "opennurbs_utils.h"

#include "obj_model_output.h"

#include "bspline_conversion.h"

#include "brep_conversion.h"

#include "folder.h"

#include <Precision.hxx>

#pragma comment(lib, "TKTopAlgo.lib")
#pragma comment(lib, "TKMesh.lib")


const int MAX_V = 50 * 10000;
const int MAX_F = 100 * 10000;

const int thread_count = 8;

int main() {

	string case_name;

	string output_prefix = "D:\\garbage\\";

	case_name = "full";

	mkdir((output_prefix + case_name).c_str());
	mkdir((output_prefix + case_name + "_opennurbs").c_str());
	mkdir((output_prefix + case_name + "_groups").c_str());

	ON::Begin();

	ONX_Model *model = new ONX_Model();;
	FILE* archive_fp = ON::OpenFile((case_name + ".3dm").c_str(), "rb");
	ON_BinaryFile archive(ON::read3dm, archive_fp);
	bool rc = model->Read(archive);
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

	vector<TopoDS_Compound> breps(model->m_object_table.Count());

	vector<RhinoMesh> obj_mesh(model->m_object_table.Count());

#pragma omp parallel for  
	for (int obj_id = 0; obj_id < model->m_object_table.Count(); obj_id++) {
		ONX_Model_Object &object = model->m_object_table[obj_id];

		const ON_Brep* brep = ON_Brep::Cast(object.m_object);

		if (brep == nullptr) {
			continue;
		}

		ON_SimpleArray<const ON_Mesh*> meshes;
		int mc = brep->GetMesh(ON::render_mesh, meshes);

		if (mc) {
			for (int i = 0; i < meshes.Count(); i++) {
				auto mesh = meshes[i];
				append_mesh(obj_mesh[obj_id].v, obj_mesh[obj_id].f, mesh);
			}
			obj_mesh[obj_id].merge_vertices();
			//obj_mesh[obj_id].write_obj((case_name + "_opennurbs\\obj_" + to_string(obj_id) + ".obj").c_str());
		}

		if (IsModelObjectVisible(*model, object)) {
			breps[obj_id] = convert_opennurbs_solid_to_occt_solid(brep);
		}
	}

	model->Destroy();
	model->DestroyCache();
	delete model;
	cout << "convert opennurbs to opencascad finished." << endl;

	ON::End();

	map<int, map<int, OcctMesh>> grp_meshs;

	map<int, BRepMesh_FastDiscret::Parameters> grp_para;

	for (auto grp_id : group_info.group_ids) {
		auto& obj_ids = group_info.group_id_objects[grp_id];
		int vs = 0;
		int fs = 0;

		for (auto obj_id : obj_ids) {
			vs += obj_mesh[obj_id].v.size();
			fs += obj_mesh[obj_id].f.size();
		}

		if (vs <= MAX_V && fs <= MAX_F) {
			//use original mesh;
		}
		else {
			BRepMesh_FastDiscret::Parameters p;

			auto multiplier =  ceil((double)fs / MAX_F);
			p.Deflection = 2.0 * multiplier;
			p.Angle = 1.0 * multiplier;

			grp_para[grp_id] = p;
			grp_meshs[grp_id] = map<int, OcctMesh>();
			for (auto obj_id : obj_ids) {
				grp_meshs[grp_id][obj_id] = OcctMesh();
			}
		}
	}


//#pragma omp parallel for
	for (auto i = 0; i < group_info.object_ids.size(); i++) {
		auto obj_id = group_info.object_ids[i];
		cout << obj_id << endl;
		for (auto grp_id : group_info.object_id_groups[obj_id]) {
			if (grp_para.count(grp_id) != 0) {
				grp_meshs[grp_id][obj_id] = generate_occt_mesh(breps[obj_id], grp_para[grp_id]);
			}
		}
	}


	//for (auto grp_id : group_info.group_ids) {
	//	auto& obj_ids = group_info.map_group_id_objects[grp_id];
	//	int vs = 0;
	//	int fs = 0;
	//
	//	for (auto obj_id : obj_ids) {
	//		vs += obj_mesh[obj_id].v.size();
	//		fs += obj_mesh[obj_id].f.size();
	//	}
	//
	//	if (vs <= MAX_V && fs < MAX_F) {
	//		//use original mesh;
	//	}
	//	else {
	//		//generate mesh;
	//		double line_tolerance = 1;
	//		double angle_tolerance = 0.5;
	//
	//		double multiplier = fs / MAX_F;
	//		line_tolerance *= multiplier;
	//		angle_tolerance *= multiplier;
	//		BRepMesh_FastDiscret::Parameters p;
	//
	//		p.Deflection = line_tolerance;
	//		p.Angle = angle_tolerance;
	//		for (auto obj_id : obj_ids) {
	//			grp_meshs[grp_id][obj_id] = generate_occt_mesh(breps[obj_id], p);
	//		}
	//	}
	//}
	

	//这段代码用来测试以不同的精度离散obj
	//for (auto &x : breps) {
	//	auto obj_id = x.first;
	//	auto &brep = x.second;
	//	BRepMesh_FastDiscret::Parameters p;
	//	p.Angle = 0.5;
	//	p.Deflection = 0.1;
	//	auto mesh1 = generate_occt_mesh(brep, p);
	//
	//	p.Angle = 0.5;
	//	p.Deflection = 1;
	//	auto mesh2 = generate_occt_mesh(brep, p);
	//	mesh1.write_obj((output_prefix + case_name + "\\obj_" + to_string(obj_id) + "_0.1.obj").c_str());
	//	mesh2.write_obj((output_prefix + case_name + "\\obj_" + to_string(obj_id) + "_1.0.obj").c_str());
	//}

	system("pause");

	return 0;
}
