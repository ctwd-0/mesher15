#pragma once

#include<vector>
#include<map>

#include <gp_Pnt.hxx>
#include <Poly_Triangulation.hxx>

#include <opennurbs.h>

using namespace std;


void append_mesh(vector<gp_Pnt>& vs, vector<vector<int>>& ts, const opencascade::handle<Poly_Triangulation> triangulation);

void append_mesh(vector<ON_3fPoint>& vertices, vector<vector<int>>& faces, const ON_Mesh* mesh);

void merge_vertices(vector<gp_Pnt>& vertices, vector<vector<int>>& faces, double eps = 1e-7);

void merge_vertices(vector<ON_3fPoint>& vertices, vector<vector<int>>& faces, double eps = 1e-7);

void write_obj(vector<gp_Pnt>& vertices, vector<vector<int>>& faces, const char* file_name);

void write_obj(vector<ON_3fPoint>& vertices, vector<vector<int>>& faces, const char* file_name);