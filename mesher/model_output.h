#pragma once

#include<vector>
#include<map>
#include <algorithm>
#include <set>

#include <opennurbs.h>

#include <gp_Pnt.hxx>
#include <Poly_Triangulation.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>

#include <TopExp_Explorer.hxx>
#include <Poly_Triangulation.hxx>

#include <BRepTools.hxx>
#include <BRepMesh.hxx>
#include <BRepMesh_IncrementalMesh.hxx>


#pragma comment(lib, "TKMesh.lib")

#include "xbj.h"
#include "opennurbs_utils.h"
#include "md5.h"

using namespace std;

class Mesh;

void append_mesh(vector<gp_Pnt>& vertices, vector<vector<int>>& faces, const opencascade::handle<Poly_Triangulation> triangulation);

void append_mesh(vector<gp_Pnt>& vertices, vector<vector<int>>& faces, const ON_Mesh* mesh);

void append_mesh(vector<gp_Pnt>& vertices, vector<vector<int>>& faces, const Mesh& mesh);

void merge_vertices(vector<gp_Pnt>& vertices, vector<vector<int>>& faces, double eps = 1e-7);

void write_obj(vector<gp_Pnt>& vertices, vector<vector<int>>& faces, const char* file_name);

class Mesh {
public:
	vector<gp_Pnt> v;
	vector<vector<int>> f;
	string name;
	char* xbj = nullptr;
	int len_xbj = 0;
	void append_mesh(const opencascade::handle<Poly_Triangulation> triangulation) {
		::append_mesh(v, f, triangulation);
	}

	void append_mesh(const ON_Mesh* mesh) {
		::append_mesh(v, f, mesh);
	}
	void append_mesh(const Mesh& other) {
		::append_mesh(v, f, other);
	}
	void write_obj(const char* file_name) {
		::write_obj(v, f, file_name);
	}

	void merge_vertices() {
		::merge_vertices(v, f);
	}
	bool empty() {
		return (v.size() == 0) || (f.size() == 0);
	}
	void generate_xbj(int byte_float = 4) {
		if (name.length() == 0 || v.size() == 0 || f.size() == 0) {
			return;
		}
		if (xbj != nullptr) {
			delete[] xbj;
		}
		len_xbj = calc_xbj_len(byte_float);
		xbj = new char[len_xbj + 4];
		int pos = 0;
		int name_length = name.length();
		memcpy_s(xbj + pos, 4, &name_length, 4);
		pos += 4;
		memcpy_s(xbj + pos, name_length, name.c_str(), name_length);
		pos += name_length;
		while (pos % 4 != 0) pos++;
		memcpy_s(xbj + pos, 1, &float_bytes, 1);
		pos++;
		memcpy_s(xbj + pos, 1, &int_bytes, 1);
		pos++;
		pos += 2;
		int number_vertices = v.size();
		memcpy_s(xbj + pos, 4, &number_vertices, 4);
		pos += 4;
		for (int i = 0; i < number_vertices; i++) {
			if (float_bytes == 8) {
				double xyz[] = { v[i].X(), v[i].Y(), v[i].Z() };
				memcpy_s(xbj + pos, float_bytes * 3, xyz, float_bytes * 3);
				pos += float_bytes * 3;
			}
			else {
				float xyz[] = { (float)v[i].X(), (float)v[i].Y(), (float)v[i].Z() };
				memcpy_s(xbj + pos, float_bytes * 3, xyz, float_bytes * 3);
				pos += float_bytes * 3;
			}
		}
		int number_faces = f.size();
		memcpy_s(xbj + pos, 4, &number_faces, 4);
		pos += 4;
		for (int i = 0; i < number_faces; i++) {
			int abc[] = { f[i][0]+1, f[i][1] + 1, f[i][2] + 1 };
			memcpy_s(xbj + pos, int_bytes, abc, int_bytes);
			pos += int_bytes;
			memcpy_s(xbj + pos, int_bytes, abc + 1, int_bytes);
			pos += int_bytes;
			memcpy_s(xbj + pos, int_bytes, abc + 2, int_bytes);
			pos += int_bytes;
		}

		if (pos != len_xbj) {
			cout << "xbj error at " << name << "." << endl;
		}
	}

	int calc_xbj_len(int byte_float = 4) {
		if (name.length() == 0 || v.size() == 0 || f.size() == 0) {
			return 0;
		}

		int max_int = v.size();
		if (max_int <= 255) {
			int_bytes = 1;
		}
		else if (max_int <= 65535) {
			int_bytes = 2;
		}
		else if (max_int < 256 * 256 * 256 - 1) {
			int_bytes = 3;
		}
		if (byte_float == 8) {
			float_bytes = 8;
		}
		int result = 0;
		result += 4;												// 名字的长度
		result += name.length() / 4 * 4;
		if (name.length() % 4 != 0) result += 4;					//名字本身，4字节对齐
		result += 4;												//float 和 int 的大小
		result += 4;												//点数
		result += float_bytes * v.size() * 3;						//点本身
		result += 4;												//面数
		result += int_bytes * f.size() * 3;							//面本身

		return result;
	}

	~Mesh() {
		if (xbj != nullptr) {
			delete[] xbj;
		}
	}
private:
	BYTE float_bytes = 4;
	BYTE int_bytes = 4;
};

Mesh generate_occt_mesh(const TopoDS_Shape& shape, const BRepMesh_FastDiscret::Parameters& p);

class OutputXbj {
public:
	map<string, string> meta_data;
	map<string, pair<string, int>> look_up;
	int groups;
	map<string, string> xbjs;
	map<string, string> files;
	OutputXbj(OpennurbsGroupInfo& info,
		map<int, Mesh>& obj_meshs,
		map<int, map<int, Mesh>>& grp_obj_meshs,
		map<int, map<int, Mesh>>& grp_grp_meshs):
	_info(info),_obj_meshs(obj_meshs),
		_grp_obj_meshs(grp_obj_meshs), _grp_grp_meshs(grp_grp_meshs){

		process_xbjs();
		process_group();
		process_lookup();
		process_meta();
	}

	void process_group() {
		files["group_info"] = string();
		auto& x = files["group_info"];
		x += (to_string(_info.group_ids.size())+ " \n");
		x += (to_string(_info.root_group_ids.size()) + " \n");
		for (auto id : _info.root_group_ids) {
			x += (to_string(id) + " ");
		}
		x += "\n";
		for (auto group_id : _info.group_ids) {
			x += (to_string(group_id) + " \n");
			x += (to_string(_info.group_composed_of_groups[group_id].size()) + " \n");
			for (auto sub_group_id : _info.group_composed_of_groups[group_id]) {
				x += (to_string(sub_group_id) + " ");
			}
			x += "\n";
			x += (to_string(_info.group_composed_of_objects[group_id].size()) + " \n");
			for (auto object_id : _info.group_composed_of_objects[group_id]) {
				x += (to_string(object_id) + " ");
			}
			x += "\n";
		}
		MD5 md5(x);
		meta_data["group_info"] = md5.toStr();
	}

	void process_xbjs() {
		//objects
		process(_obj_meshs, "object");
		//groups and objects
		for (auto &x : _grp_grp_meshs) {
			process(x.second, "group_" + to_string(x.first) + "group");
		}
		for (auto &x : _grp_obj_meshs) {
			process(x.second, "group_" + to_string(x.first) + "object");
		}
	}

	void process_meta() {
		files["metadata"] = "";
		auto& x = files["metadata"];
		for (auto & y : meta_data) {
			x += y.first + " " + y.second + "\n";
		}
	}

	void process_lookup() {
		files["lookup"] = "";
		auto& x = files["lookup"];
		for (auto& y : look_up) {
			x += y.first + " " + y.second.first + " " + to_string(y.second.second) + "\n";
		}
		MD5 md5(x);

		meta_data["lookup"] = md5.toStr();
	}

	void output(const string& prefix) {
		for (auto& x : files) {
			FILE* f;
			fopen_s(&f, (prefix + "\\" + x.first).c_str(), "wb");
			fwrite(x.second.c_str(), x.second.length(), 1, f);
			fclose(f);
		}
	}

private:
	OpennurbsGroupInfo& _info;
	map<int, Mesh>& _obj_meshs;
	map<int, map<int, Mesh>>& _grp_obj_meshs;
	map<int, map<int, Mesh>>& _grp_grp_meshs;

	void process(map<int, Mesh> meshes, string prefix) {
		vector<pair<int, int>> base;
		for (auto &x : _obj_meshs) {
			base.push_back(pair<int, int>(x.second.len_xbj, x.first));
		}
		sort(base.begin(), base.end());
		reverse(base.begin(), base.end());
		int postfix = 1;
		int left = 0; int right = 0;
		while (left < base.size()) {
			right = left + 1;
			int now_size = base[left].first;
			while (right < base.size() && (now_size + base[right].first) < MAX_XBJ_SIZE) {
				now_size += base[right].first;
				right++;
			}
			int file_size = now_size + 12;
			int mesh_count = right - left;
			cout << left << " " << right << endl;
			cout << file_size << endl;
			string file_name = prefix + "_" + to_string(postfix) + ".xbj";
			files[file_name] = string(file_size, ' ');
			auto &file = files[file_name];
			cout << file.size() << endl;
			file.replace(0, 8, "XBJV0001");
			file.replace(8, 4, (char*)&mesh_count);
			int lookup_pos = 12;
			for (int i = left; i < right; i++) {
				cout << lookup_pos<<" " << base[i].first <<endl;
				look_up[_obj_meshs[base[i].second].name] = pair<string, int>(file_name, lookup_pos);
				file.replace(lookup_pos, base[i].first, _obj_meshs[base[i].second].xbj);
				lookup_pos += base[i].first;
			}
			MD5 md5(file);
			meta_data[file_name] = md5.toStr();
			postfix++;
			left = right;
		}
	}
};
