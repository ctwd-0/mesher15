#include "obj_model_output.h"




void append_mesh(vector<gp_Pnt>& vertices, vector<vector<int>>& faces,
	const opencascade::handle<Poly_Triangulation> triangulation) {
	auto n_triangles = triangulation->NbTriangles();
	auto n_nodes = triangulation->NbNodes();
	auto triangles = triangulation->Triangles();
	auto nodes = triangulation->Nodes();

	Standard_Integer i1, i2, i3;

	int offset = vertices.size();
	for (auto i = 1; i <= n_triangles; i++) {
		auto triangle = triangles.Value(i);
		triangle.Get(i1, i2, i3);
		faces.push_back(vector<int>{i1 + offset - 1, i2 + offset - 1, i3 + offset - 1});
	}

	for (auto i = 1; i <= n_nodes; i++) {
		vertices.push_back(nodes.Value(i));
	}
}

void append_mesh(vector<gp_Pnt>& vertices, vector<vector<int>>& faces, const Mesh& mesh) {
	int offset = vertices.size();

	for (auto i = 0; i < mesh.f.size(); i++) {
		auto t = mesh.f[i];
		faces.push_back(vector<int>{t[0] + offset, t[1] + offset, t[2] + offset});
	}

	for (auto i = 0; i < mesh.v.size(); i++) {
		vertices.push_back(mesh.v[i]);
	}
}

void append_mesh(vector<gp_Pnt>& vertices, vector<vector<int>>& faces, const ON_Mesh* mesh) {
	map<int, int> m_index;
	vector<int> used_point;
	int offset = vertices.size();
	for (int j = 0; j < mesh->m_F.Count(); j++) {
		auto face = mesh->m_F[j];
		for (int k = 0; k < 4; k++) {
			auto index = face.vi[k];
			if (m_index.count(index) == 0) {
				m_index[index] = offset + m_index.size();
				used_point.push_back(index);
			}
		}
		if (face.vi[2] == face.vi[3]) {
			vector<int> fi{ m_index[face.vi[0]], m_index[face.vi[1]], m_index[face.vi[2]]};
			faces.push_back(fi);
		}
		else {
			vector<int> f1{ m_index[face.vi[0]], m_index[face.vi[1]], m_index[face.vi[2]] };
			faces.push_back(f1);
			vector<int> f2{ m_index[face.vi[0]], m_index[face.vi[2]], m_index[face.vi[3]] };
			faces.push_back(f2);
		}
	}

	for (int j = 0; j < used_point.size(); j++) {
		auto point = mesh->m_V[used_point[j]];
		gp_Pnt pnt(point.x, point.y, point.z);
		vertices.push_back(pnt);
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

Mesh generate_occt_mesh(
	const TopoDS_Shape& shape,
	const BRepMesh_FastDiscret::Parameters& p) {
	Mesh result;

	auto x = BRepMesh_IncrementalMesh(shape, p);
	for (TopExp_Explorer face_exp(shape, TopAbs_FACE); face_exp.More(); face_exp.Next()) {
		TopLoc_Location location;

		auto &face = TopoDS::Face(face_exp.Current());
		auto triangulation = BRep_Tool::Triangulation(face, location);
		if (!triangulation.IsNull()) {
			result.append_mesh(triangulation);
		}
		BRepTools::Clean(face);
	}
	return result;
}