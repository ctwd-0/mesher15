#pragma once

#include<vector>
#include<map>

#include <gp_Pnt.hxx>
#include <Poly_Triangulation.hxx>

#include <opennurbs.h>

#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>

#include <TopExp_Explorer.hxx>
#include <Poly_Triangulation.hxx>

#include <BRepTools.hxx>
#include <BRepMesh.hxx>
#include <BRepMesh_IncrementalMesh.hxx>

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
};

Mesh generate_occt_mesh(const TopoDS_Shape& shape, const BRepMesh_FastDiscret::Parameters& p);
