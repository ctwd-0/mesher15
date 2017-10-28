#include "obj_model_output.h"

using namespace std;

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
			vector<int> fi{ m_index[face.vi[0]], m_index[face.vi[1]], m_index[face.vi[2]] };
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

class comp_x {
public:
	comp_x(vector<gp_Pnt>& vertices) :_vertices(vertices) {}
	bool operator()(const int &a, const int &b) {
		if (_vertices[a].X() <= _vertices[b].X())
			return true;
		return false;
	}
private:
	vector<gp_Pnt>& _vertices;
};

class comp_y {
public:
	comp_y(vector<gp_Pnt>& vertices) :_vertices(vertices) {}
	bool operator()(const int &a, const int &b) {
		if (_vertices[a].Y() <= _vertices[b].Y())
			return true;
		return false;
	}
private:
	vector<gp_Pnt>& _vertices;
};

class comp_z {
public:
	comp_z(vector<gp_Pnt>& vertices) :_vertices(vertices) {}
	bool operator()(const int &a, const int &b) {
		if (_vertices[a].Z() <= _vertices[b].Z())
			return true;
		return false;
	}
private:
	vector<gp_Pnt>& _vertices;
};

#define __lb(x, X) auto iter_##x = lower_bound(_v##x.begin(), _v##x.end(), _vertices[index].X(), comp_##x(_vertices));
#define __dc(x, X) while (iter_##x > _v##x.begin() && abs(_vertices[*iter_##x].X() - _vertices[index].X())< eps) iter_##x--;
#define __ad(x, X) \
		while (iter_##x < _v##x.end() && !(_vertices[*iter_##x].X() - _vertices[index].X() > eps)){\
			if (abs(_vertices[*iter_##x].X() - _vertices[index].X()) < eps) {\
				c_##x.insert(*iter_##x);\
			}\
		iter_##x++;\
		}

class sort_for_merge {
public:
	sort_for_merge(vector<gp_Pnt>& vertices,
		vector<int>& vx,
		vector<int>& vy,
		vector<int>& vz)
		:_vertices(vertices), _vx(vx), _vy(vy), _vz(vz) {
		exec();
	}
	vector<int> candidate(int index, double eps) {
		eps = abs(eps);
		set<int> c_x, c_y, c_z;
		//__lb(x, X);
		//__lb(y, Y);
		//__lb(z, Z);
		//__dc(x, X);
		//__dc(y, Y);
		//__dc(z, Z);
		//__ad(x, X);
		//__ad(y, Y);
		//__ad(z, Z);

		auto iter_x = lower_bound(_vx.begin(), _vx.end(), _vertices[index].X(), comp_x(_vertices));;
		auto iter_y = lower_bound(_vy.begin(), _vy.end(), _vertices[index].Y(), comp_y(_vertices));;
		auto iter_z = lower_bound(_vz.begin(), _vz.end(), _vertices[index].Z(), comp_z(_vertices));;
		while (iter_x > _vx.begin() && abs(_vertices[*iter_x].X() - _vertices[index].X()) < eps) {
			iter_x--;;
		}
		while (iter_y > _vy.begin() && abs(_vertices[*iter_y].Y() - _vertices[index].Y()) < eps) {
			iter_y--;;
		}
		while (iter_z > _vz.begin() && abs(_vertices[*iter_z].Z() - _vertices[index].Z()) < eps) {
			iter_z--;;
		}
		while (iter_x < _vx.end() && !(_vertices[*iter_x].X() - _vertices[index].X() > eps)) {
			if (abs(_vertices[*iter_x].X() - _vertices[index].X()) < eps) {
				c_x.insert(*iter_x);
			}
			iter_x++;
		};
		while (iter_y < _vy.end() && !(_vertices[*iter_y].Y() - _vertices[index].Y() > eps)) {
			if (abs(_vertices[*iter_y].Y() - _vertices[index].Y()) < eps) {
				c_y.insert(*iter_y);
			} iter_y++;
		};
		while (iter_z < _vz.end() && !(_vertices[*iter_z].Z() - _vertices[index].Z() > eps)) {
			if (abs(_vertices[*iter_z].Z() - _vertices[index].Z()) < eps) {
				c_z.insert(*iter_z);
			}
			iter_z++;
		};

		vector<int> result;
		eps *= eps;
		for (auto i : c_x) {
			if (i != index && c_y.count(i) != 0 && c_x.count(i) != 0) {
				if (_vertices[i].SquareDistance(_vertices[index]) < eps) {
					result.push_back(i);
				}
			}
		}
		return result;
	}
private:
	void exec() {
		_vx.swap(vector<int>(_vertices.size()));
		_vy.swap(vector<int>(_vertices.size()));
		_vz.swap(vector<int>(_vertices.size()));
		for (int i = 0; i < _vertices.size(); i++) {
			_vx[i] = i;
			_vy[i] = i;
			_vz[i] = i;
		}
		sort(_vx.begin(), _vx.end(), comp_x(_vertices));
		sort(_vy.begin(), _vy.end(), comp_y(_vertices));
		sort(_vz.begin(), _vz.end(), comp_z(_vertices));
	}
	vector<gp_Pnt>& _vertices;
	vector<int>& _vx;
	vector<int>& _vy;
	vector<int>& _vz;
};

void merge_vertices(vector<gp_Pnt>& vertices, vector<vector<int>>& faces, double eps) {
	eps *= eps;
	map<int, int> duplicates;
	vector<int> used_id;
	set<int> used_id_set;
	map<int, int> used_id_new;

	vector<int> vx;
	vector<int> vy;
	vector<int> vz;
	sort_for_merge sfm(vertices, vx, vy, vz);

	for (int i = 0; i < vertices.size(); i++) {
		gp_Pnt point = vertices[i];
		if (used_id.size() == 0) {
			used_id.push_back(i);
			used_id_set.insert(i);
			used_id_new[i] = used_id_new.size();
		}
		else {
			bool is_new = true;
			int dup = -1;
			auto candidates = sfm.candidate(i, eps);
			for (auto candidate : candidates) {
				if (used_id_set.count(candidate) != 0) {
					is_new = false;
					dup = candidate;
					break;
				}
			}
			if (is_new) {
				used_id.push_back(i);
				used_id_set.insert(i);
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