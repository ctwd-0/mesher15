#pragma once

#include<vector>

#include <gp_Pnt.hxx>
#include <Poly_Triangulation.hxx>

using namespace std;


void write_obj(vector<gp_Pnt>& vertices, vector<vector<int>>& faces, const char* file_name);

void append_mesh(vector<gp_Pnt>& vs, vector<vector<int>>& ts, const opencascade::handle<Poly_Triangulation> triangulation);