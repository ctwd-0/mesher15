#include "obj_model_output.h"


void write_obj(vector<gp_Pnt>& vertices, vector<vector<int>>& faces, const char* file_name) {
	FILE* f;
	auto error = fopen_s(&f, file_name, "wb");

	char line[1000] = { 0 };

	for (int i = 0; i < vertices.size(); i++) {
		sprintf(line, "v %lf %lf %lf\n\0", vertices[i].XYZ().X(), vertices[i].XYZ().Y(), vertices[i].XYZ().Z());
		fwrite(line, strlen(line), 1, f);
	}

	for (int i = 0; i < faces.size(); i++) {
		sprintf(line, "f %d %d %d\n\0", faces[i][0] + 1, faces[i][1] + 1, faces[i][2] + 1);
		fwrite(line, strlen(line), 1, f);
	}

	fclose(f);
}

void append_mesh(vector<gp_Pnt>& vs, vector<vector<int>>& ts, const opencascade::handle<Poly_Triangulation> triangulation) {
	auto n_triangles = triangulation->NbTriangles();
	auto n_nodes = triangulation->NbNodes();
	auto triangles = triangulation->Triangles();
	auto nodes = triangulation->Nodes();

	Standard_Integer i1, i2, i3;

	int offset = vs.size();
	for (auto i = 1; i <= n_triangles; i++) {
		auto triangle = triangles.Value(i);
		triangle.Get(i1, i2, i3);
		ts.push_back(vector<int>{i1 + offset - 1, i2 + offset - 1, i3 + offset - 1});
	}

	for (auto i = 1; i <= n_nodes; i++) {
		vs.push_back(nodes.Value(i));
	}
}

void merge_vertices(vector<gp_Pnt>& vertices, vector<vector<int>>& faces, double eps) {
	eps *= eps;
	map<int, int> duplicates;
	vector<int> used_id;
	map<int, int> used_id_new;

	for (int i = 0; i < vertices.size(); i++) {
		gp_Pnt point = vertices[i];
		if (used_id.size() == 0) {
			used_id.push_back(i);
			used_id_new[i] = used_id_new.size();
		}
		else {
			bool is_new = true;
			int dup = -1;
			for (int j = 0; j < used_id.size(); j++) {
				double square_distance = point.SquareDistance(vertices[used_id[j]]);
				if (square_distance < eps) {
					is_new = false;
					dup = used_id[j];
					break;
				}
			}
			if (is_new) {
				used_id.push_back(i);
				used_id_new[i] = used_id_new.size();
			}
			else {
				duplicates[i] = dup;
			}
		}
	}

	for (int i = 0; i < faces.size(); i++) {
		for (int j = 0; j < faces[i].size(); j++) {
			int old_id = faces[i][j];
			if (duplicates.count(old_id)) {
				old_id = duplicates[old_id];
			}
			faces[i][j] = used_id_new[old_id];
		}
	}
	for (int i = 0; i < used_id.size(); i++) {
		vertices[used_id_new[used_id[i]]] = vertices[used_id[i]];
	}
	while (vertices.size() > used_id.size()) {
		vertices.pop_back();
	}
}