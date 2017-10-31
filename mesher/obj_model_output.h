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
