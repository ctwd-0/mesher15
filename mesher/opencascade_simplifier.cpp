#include <iostream>
#include <vector>
#include <map>
#include <thread>
#include <mutex>
#include <set>

#include <math.h>
#include <omp.h>
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

#pragma comment(lib, "TKTopAlgo.lib")
#pragma comment(lib, "TKMesh.lib")


const int MAX_V = 50 * 10000;
const int MAX_F = 100 * 10000;

const int thread_count = 8;

class nurbs_conversion_data {
public:
	mutex lock;
	volatile int now_obj_id;
	volatile int finished_thread;
	ONX_Model* model;
	vector<Mesh>* obj_mesh;
	vector<TopoDS_Compound>* breps;
};

int process_nurbs(nurbs_conversion_data* data) {
	int obj_id;

	while (true)
	{
		data->lock.lock();
		if (data->now_obj_id < data->model->m_object_table.Count()) {
			obj_id = data->now_obj_id;
			(data->now_obj_id)++;
			data->lock.unlock();
		}
		else {
			data->finished_thread++;
			data->lock.unlock();
			break;
		}

		ONX_Model_Object &object = data->model->m_object_table[obj_id];

		const ON_Brep* brep = ON_Brep::Cast(object.m_object);

		if (brep == nullptr) {
			continue;
		}

		ON_SimpleArray<const ON_Mesh*> meshes;
		int mc = brep->GetMesh(ON::render_mesh, meshes);

		if (mc) {
			(*(data->obj_mesh))[obj_id].name = "object_" + to_string(obj_id);
			for (int i = 0; i < meshes.Count(); i++) {
				auto mesh = meshes[i];
				(*(data->obj_mesh))[obj_id].append_mesh(mesh);
			}
			//(*(data->obj_mesh))[obj_id].merge_vertices();
			//obj_mesh[obj_id].write_obj((case_name + "_opennurbs\\obj_" + to_string(obj_id) + ".obj").c_str());
		}

		if (IsModelObjectVisible(*(data->model), object)) {
			(*(data->breps))[obj_id] = convert_opennurbs_solid_to_occt_solid(brep);
		}
	}

	return 0;
}

int main() {
	string case_name;

	string output_prefix = "D:\\garbage\\";

	case_name = "full";

	mkdir((output_prefix + case_name).c_str());
	mkdir((output_prefix + case_name + "_opennurbs").c_str());
	mkdir((output_prefix + case_name + "_groups").c_str());

	clock_t start, end;

	start = clock();
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

	end = clock();

	cout << "finished reading file. " << (end - start) / 1000.0 << "s used." << endl;

	OpennurbsGroupInfo group_info(model);

	cout <<"number of groups: "<< group_info.group_ids.size() << endl;

	vector<TopoDS_Compound> breps(model->m_object_table.Count());

	vector<Mesh> obj_mesh(model->m_object_table.Count());

	nurbs_conversion_data data;
	data.model = model;
	data.obj_mesh = &obj_mesh;
	data.breps = &breps;
	
	start = clock();

	for (int i = 0; i < thread_count; i++) {
		thread thd(&process_nurbs, &data);
		thd.detach();
	}

	while (true) {
		Sleep(50);
		data.lock.lock();
		if (data.now_obj_id >= model->m_object_table.Count()
			&& data.finished_thread == thread_count) {
			data.lock.unlock();
			break;
		}
		else {
			data.lock.unlock();
		}
	}
	end = clock();

	cout << "convert opennurbs to opencascad finished. " << (end - start) / 1000.0 << "s used." << endl;

	start = clock();
	
	model->Destroy();
	model->DestroyCache();
	delete model;
	
	end = clock();
	cout << "delete opennurbs model finished. " << (end - start) / 1000.0 << "s used." << endl;
	
	ON::End();

	map<int, map<int, Mesh>> grp_meshs;

	map<int, Mesh> grp_mesh;

	map<int, BRepMesh_FastDiscret::Parameters> grp_para;

	start = clock();
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
			grp_meshs[grp_id] = map<int, Mesh>();
			for (auto obj_id : obj_ids) {
				grp_meshs[grp_id][obj_id] = Mesh();
			}
		}
	}

#pragma omp parallel for
	for (auto i = 0; i < group_info.object_ids.size(); i++) {
		auto obj_id = group_info.object_ids[i];
		//cout << obj_id << endl;
		for (auto grp_id : group_info.object_id_groups[obj_id]) {
			if (grp_para.count(grp_id) != 0) {
				grp_meshs[grp_id][obj_id] = generate_occt_mesh(breps[obj_id], grp_para[grp_id]);
			}
		}
	}
	end = clock();
	cout << "generate mesh finished. " << (end - start) / 1000.0 << "s used." << endl;
	
	start = clock();
	vector<int> grp_ids;
	for (auto grp_id : group_info.group_ids) {
		grp_mesh[grp_id] = Mesh();
		grp_ids.push_back(grp_id);
	}

#pragma omp parallel for
	for (auto i = 0; i < grp_ids.size(); i++) {
		auto grp_id = grp_ids[i];
		grp_mesh[grp_id] = Mesh();
		grp_mesh[grp_id].name = "group_" + to_string(grp_id);
		if (grp_meshs.count(grp_id)) {
			for (auto& x : grp_meshs[grp_id]) {
				auto& obj_mesh = x.second;
				grp_mesh[grp_id].append_mesh(obj_mesh);
			}
		}
		else {
			auto& obj_ids = group_info.group_id_objects[grp_id];
			for (auto obj_id : obj_ids) {
				grp_mesh[grp_id].append_mesh(obj_mesh[obj_id]);
			}
		}

		//cout <<"group "<< grp_id << "gets "<< grp_mesh[grp_id].v.size() <<" vertices."<< endl;
		//grp_mesh[grp_id].merge_vertices();
	}

	end = clock();
	cout << "combine mesh finished. " << (end - start)/1000.0 << "s used." << endl;

	start = clock();
	breps.swap(vector<TopoDS_Compound>());
	obj_mesh.swap(vector<Mesh>());
	grp_meshs.swap(map<int, map<int, Mesh>>());
	grp_mesh.swap(map<int, Mesh>());
	grp_para.swap(map<int, BRepMesh_FastDiscret::Parameters>());
	end = clock();
	cout << "release memory finished. " << (end - start) / 1000.0 << "s used." << endl;

	system("pause");

	return 0;
}
