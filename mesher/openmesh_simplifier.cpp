//#include <iostream>
//#include <vector>
//#include <map>
//#include <set>
//#include <time.h>
//
//#include "folder.h"
//
//#define ON_DLL_IMPORTS
//
//#include "opennurbs\opennurbs_dynamic_linking.h"
//#include "openmesh\openmesh_dynamic_linking.h"
//
//#include <OpenMesh\Core\IO\MeshIO.hh>
//#include <OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh>
//#include <OpenMesh\Core\Mesh\TriConnectivity.hh>
//#include <OpenMesh\Core\Mesh\PolyConnectivity.hh>
//#include <OpenMesh\Tools\Decimater\ModQuadricT.hh>
//
//#include <OpenMesh\Tools\Decimater\DecimaterT.hh>
//
//#include <opennurbs.h>
//
//using namespace std;
//using namespace OpenMesh;
//
//const float pi = 3.1415926535897932384626433832795;
//
//const float esp = 1e-6;
//
//typedef OpenMesh::TriMesh_ArrayKernelT<> TriMesh;
//
//typedef OpenMesh::Decimater::DecimaterT<TriMesh> MyDecimater;
//
//typedef OpenMesh::Decimater::ModQuadricT<TriMesh>::Handle HModQuadric;
//
//void write_obj(vector<ON_3fPoint>& vertices, vector<vector<int>>& faces, const char* file_name) {
//	FILE* f;
//	auto error = fopen_s(&f, file_name, "wb");
//
//	char line[100] = { 0 };
//
//	for (int i = 0; i < vertices.size(); i++) {
//		sprintf(line, "v %lf %lf %lf\n\0", vertices[i].x, vertices[i].y, vertices[i].z);
//		fwrite(line, strlen(line), 1, f);
//	}
//
//	for (int i = 0; i < faces.size(); i++) {
//		sprintf(line, "f %d %d %d\n\0", faces[i][0] + 1, faces[i][1] + 1, faces[i][2] + 1);
//		fwrite(line, strlen(line), 1, f);
//	}
//
//	fclose(f);
//}
//
//void protectMeshBoundariesFromDecimation(TriMesh& mesh) {
//	mesh.request_vertex_status();
//	map<int, Vec3f> face_normals;
//
//	//int counter = 0;
//	for (const auto & fh : mesh.faces()) {
//		face_normals[fh.idx()] = mesh.calc_face_normal(fh);
//	}
//
//	for (const auto &vh : mesh.vertices()) {
//		vector<Vec3f> results;
//		for (auto fh = mesh.vf_iter(vh); fh.is_valid(); ++fh) {
//			int fid = fh.handle().idx();
//			auto fn = face_normals[fid];
//			bool new_direction = true;
//			for (auto n : results) {
//				auto temp = n.data()[0] * fn.data()[0]
//					+ n.data()[1] * fn.data()[1]
//					+ n.data()[2] * fn.data()[2];
//				auto angle = acosf(temp) * 180 / pi;
//				//if (vh.idx() == 21) {
//				//	cout << temp << " " << angle << endl;
//				//}
//				if (temp >= 1 || angle < 20 || angle > 160) {
//					new_direction = false;
//					break;
//				}
//			}
//			if (new_direction) {
//				results.push_back(fn);
//			}
//		}
//		if (results.size() >= 3) {
//			mesh.status(vh).set_locked(true);
//			//counter++;
//			//cout << "vertex:" << vh.idx() << endl;
//			//for (auto v : results) {
//			//	cout << v.data()[0] << " " << v.data()[1] << " " << v.data()[2] << " " << endl;
//			//}
//		}
//	}
//
//	//cout << "locked vertices:" << counter << endl;
//}
//
//void append_mesh(vector<ON_3fPoint>& vertices, vector<vector<int>>& faces, const ON_Mesh* mesh) {
//	map<int, int> map_index;
//	vector<int> used_point;
//	int offset = vertices.size();
//	for (int j = 0; j < mesh->m_F.Count(); j++) {
//		auto face = mesh->m_F[j];
//		for (int k = 0; k < 4; k++) {
//			auto index = face.vi[k];
//			if (map_index.count(index) == 0) {
//				map_index[index] = offset + map_index.size();
//				used_point.push_back(index);
//			}
//		}
//		if (face.vi[2] == face.vi[3]) {
//			vector<int> fi;
//			fi.push_back(map_index[face.vi[0]]);
//			fi.push_back(map_index[face.vi[1]]);
//			fi.push_back(map_index[face.vi[2]]);
//			faces.push_back(fi);
//		} else {
//			vector<int> f1;
//			f1.push_back(map_index[face.vi[0]]);
//			f1.push_back(map_index[face.vi[1]]);
//			f1.push_back(map_index[face.vi[2]]);
//			faces.push_back(f1);
//			vector<int> f2;
//			f2.push_back(map_index[face.vi[0]]);
//			f2.push_back(map_index[face.vi[2]]);
//			f2.push_back(map_index[face.vi[3]]);
//			faces.push_back(f2);
//		}
//	}
//
//	for (int j = 0; j < used_point.size(); j++) {
//		vertices.push_back(mesh->m_V[used_point[j]]);
//	}
//}
//
//void append_mesh(vector<ON_3fPoint>& vertices, vector<vector<int>>& faces, TriMesh& mesh) {
//	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
//		vector<int> ids;
//		for (auto fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
//			ids.push_back(fv_it->idx() + vertices.size());
//		}
//
//		int neighbor_count = 0;
//		for (auto ff_it = mesh.ff_iter(f_it); ff_it.is_valid(); ff_it++) {
//			neighbor_count++;
//		}
//		if (neighbor_count >= 2)
//			faces.push_back(ids);
//	}
//
//	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
//		auto d = mesh.point(v_it).data();
//		vertices.push_back(ON_3fPoint(d[0], d[1], d[2]));
//	}
//}
//
//void append_mesh(vector<ON_3fPoint>& dst_vertices, vector<vector<int>>& dst_faces, vector<ON_3fPoint>& src_vertices, vector<vector<int>>& src_faces) {
//	for (int i = 0; i < src_faces.size(); i++) {
//		vector<int> ids = src_faces[i];
//		for (int j = 0; j < ids.size(); j++) {
//			ids[j] += dst_vertices.size();
//		}
//		dst_faces.push_back(ids);
//	}
//
//	for (int i = 0; i < src_vertices.size(); i++) {
//		dst_vertices.push_back(src_vertices[i]);
//	}
//}
//
//float triangle_area(float a, float b, float c) {
//	float m = (a + b + c) / 2;
//	if (a + b > c || a + c > b || b + c > a) {
//		return sqrt(m*(m - a)*(m - b)*(m - c));
//	} else {
//		return 0;
//	}
//}
//
//float triangle_area(TriMesh& mesh, FaceHandle fh) {
//	vector<float> fels;
//	for (auto fe = mesh.fe_begin(fh); fe.is_valid(); fe++) {
//		fels.push_back(mesh.calc_edge_length(fe));
//	}
//	float area = triangle_area(fels[0], fels[1], fels[2]);
//	return area;
//}
//
//Vec3f edge_center(TriMesh& mesh, EdgeHandle &et) {
//	auto to = mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(et, 0)));
//	auto from = mesh.point(mesh.from_vertex_handle(mesh.halfedge_handle(et, 0)));
//
//	return (to + from) / 2;
//}
//
//TriMesh openmesh_simplify_mesh(vector<ON_3fPoint>& vertices, vector<vector<int>>& faces) {
//	//cout << "begin push mesh" << endl;
//	TriMesh mesh;
//	for (int i = 0; i < vertices.size(); i++) {
//		mesh.add_vertex(Vec3f(vertices[i].x, vertices[i].y, vertices[i].z));
//	}
//
//	for (int i = 0; i < faces.size(); i++) {
//		if (faces[i][0] == faces[i][1] || faces[i][1] == faces[i][2] || faces[i][0] == faces[i][2])
//			continue;
//		mesh.add_face(VertexHandle(faces[i][0]), VertexHandle(faces[i][1]), VertexHandle(faces[i][2]));
//	}
//
//	if (mesh.n_faces() <= 1) {
//		return mesh;
//	}
//	//cout << "end push mesh" << endl;
//
//	//for (auto f_it : mesh.faces()) {
//	//	cout << f_it.idx() << endl;
//	//}
//	protectMeshBoundariesFromDecimation(mesh);
//
//	MyDecimater decimater(mesh);
//	HModQuadric hModQuadric;
//	decimater.add(hModQuadric);
//
//	//float max_area = 0;
//	//map<int, float> areas;
//	//for (auto ft : mesh.faces()) {
//	//	float area = triangle_area(mesh, ft);
//	//	if (max_area < area) {
//	//		max_area = area;
//	//	}
//	//	areas[ft.idx()] = area;
//	//}
//
//	Vec3f mins = mesh.calc_face_centroid(mesh.faces_begin());
//	Vec3f maxs = mesh.calc_face_centroid(mesh.faces_begin());
//	map<int, Vec3f> mids;
//	for (auto ft : mesh.faces()) {
//		mids[ft.idx()] = mesh.calc_face_centroid(ft);
//	}
//
//	for (auto vt : mesh.vertices()) {
//		auto p = mesh.point(vt);
//		auto d = p.data();
//		for (int i = 0; i < 3; i++) {
//			if (mins.data()[i] > d[i]) {
//				mins.data()[i] = d[i];
//			}
//			if (maxs.data()[i] < d[i]) {
//				maxs.data()[i] = d[i];
//			}
//		}
//	}
//
//	Vec3f center = (mins + maxs) / 2;
//
//	Vec3f box = (maxs - mins);
//
//	float size = box.data()[0] > box.data()[1] ? box.data()[0] : box.data()[1];
//	size = size > box.data()[2] ? size : box.data()[2];
//
//
//	//cout << box.data()[0] << " " << box.data()[1] << " " << box.data()[2] << endl;
//	bool check = true;
//	for (auto et : mesh.edges()) {
//		auto ec = edge_center(mesh, et);
//		auto dis = (ec - center).length();
//		if (dis < size * 0.2) {
//			check = false;
//			break;
//		}
//	}
//	//cout << "check" << check << endl;
//
//	//cout << "begin simplify mesh" << endl;
//	decimater.initialize();
//	//decimater.decimate();
//	decimater.decimate_to_faces(vertices.size() / 10, faces.size() / 10);
//
//	//for (auto f_it : mesh.faces()) {
//	//	if(!mesh.status(f_it).deleted())
//	//		cout << f_it.idx() << endl;
//	//}
//
//	//Vec3f ds = maxs - mins;
//	//vector<float> dis;
//	//for (int i = 0; i < 3; i++) {
//	//	dis.push_back(ds.data()[i]);
//	//}
//	//sort(dis.begin(), dis.end(), greater<float>());
//	//if (dis[2] < dis[0] * 0.5) {
//	//	dis.pop_back();
//	//}
//	//if (dis[1] < dis[0] * 0.5) {
//	//	dis.pop_back();
//	//}
//	vector<int> removed_faces;
//
//	if (check) {
//		for (auto ft : mesh.faces()) {
//			for (auto fe = mesh.fe_begin(ft); fe.is_valid(); ++fe) {
//				auto mid = edge_center(mesh, fe.handle());
//				auto dis = (mid - center).length();
//				if (dis < size * 0.2) {
//					removed_faces.push_back(ft.idx());
//					break;
//				}
//			}
//		}
//	}
//	//for (auto ft : mesh.faces()) {
//	//	auto center = mesh.calc_face_centroid(ft);
//	//	auto len = (center - mids[ft.idx()]).length();
//	//	if (len > dis.back() * 0.6) {
//	//		removed_faces.push_back(ft.idx());
//	//	}
//	//}
//
//	//for (auto ft : mesh.faces()) {
//	//	vector<int> neightbors;
//	//	for (auto fnt = mesh.ff_begin(ft); fnt.is_valid(); ++fnt) {
//	//		neightbors.push_back(fnt.handle().idx());
//	//	}
//
//	//	if (neightbors.size() == 3) {
//	//		if (neightbors[0] == neightbors[1] || neightbors[0] == neightbors[2] || neightbors[2] == neightbors[1]) {
//	//			removed_faces.push_back(ft.idx());
//	//		}
//	//	}
//	//}
//	//for (auto ft : mesh.faces()) {
//	//	if (mesh.status(ft).deleted()) {
//	//		continue;
//	//	}
//	//	float old_area = areas[ft.idx()];
//	//	float area = triangle_area(mesh, ft);
//	//	float neighbor_area = area;
//	//	for (auto fnt = mesh.ff_begin(ft); fnt.is_valid(); fnt++) {
//	//		double n_area = triangle_area(mesh, fnt);
//	//		if (neighbor_area > n_area) {
//	//			neighbor_area = n_area;
//	//		}
//	//	}
//	//	if (area > old_area * 20 && area > neighbor_area * 20) {
//	//		removed_faces.push_back(ft.idx());
//	//	}
//	//}
//
//	for (auto id : removed_faces) {
//		mesh.delete_face(FaceHandle(id));
//	}
//
//	mesh.garbage_collection();
//
//	//cout << "end simplify mesh" << endl;
//
//	return mesh;
//}
//
//
//#include <CGAL/Orthogonal_k_neighbor_search.h>
//#include <CGAL/Kd_tree.h>
//#include <CGAL/Search_traits.h>
//#include "neighbor_search.hpp"
//void merge_vertices(vector<ON_3fPoint>& vertices, vector<vector<int>>& faces) {
//	typedef CGAL::Dimension_tag<3> D;
//	typedef CGAL::Search_traits<float, Point, const float*, Construct_coord_iterator, D> Traits;
//	typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance> Neighbor_search;
//	typedef Neighbor_search::Tree Tree;
//
//	//cout << "old vertices:" << vertices.size() << endl;
//	map<int, int> duplicates;
//	vector<int> used_id;
//	map<int, int> used_id_new;
//	Tree tree;
//	Distance dist;
//
//	time_t start, end;
//	start = clock();
//	for (int i = 0; i < vertices.size(); i++) {
//		Point point(vertices[i].x, vertices[i].y, vertices[i].z, i);
//		if (used_id.size() == 0) {
//			used_id.push_back(i);
//			used_id_new[i] = used_id_new.size();
//		} else {
//			bool is_new = true;
//			int dup = -1;
//			for (int j = 0; j < used_id.size(); j++) {
//				double dis = dist.transformed_distance(point, Point(vertices[used_id[j]].x, vertices[used_id[j]].y, vertices[used_id[j]].z));
//				if (dist.inverse_of_transformed_distance(dis) < esp) {
//					is_new = false;
//					dup = used_id[j];
//					break;
//				}
//			}
//			if (is_new) {
//				used_id.push_back(i);
//				used_id_new[i] = used_id_new.size();
//			} else {
//				duplicates[i] = dup;
//			}
//		}
//	}
//	//for (int i = 0; i < vertices.size(); i++) {
//	//	Point point(vertices[i].x, vertices[i].y, vertices[i].z, i);
//	//	if (tree.size() == 0) {
//	//		tree.insert(point);
//	//		used_id.push_back(i);
//	//		used_id_new[i] = used_id_new.size();
//	//	} else {
//	//		Neighbor_search search(tree, point, 1);
//	//		auto it = search.begin();
//	//		if (it != search.end() && tr_dist.inverse_of_transformed_distance(it->second) < esp) {
//	//			duplicates[i] = it->first.id;
//	//		} else {
//	//			tree.insert(point);
//	//			used_id.push_back(i);
//	//			used_id_new[i] = used_id_new.size();
//	//		}
//	//	}
//	//}
//	//end = clock();
//	//cout << "search time:" << end - start << "ms"<< endl;
//
//	start = clock();
//	for (int i = 0; i < faces.size(); i++) {
//		for (int j = 0; j < faces[i].size(); j++) {
//			int old_id = faces[i][j];
//			if (duplicates.count(old_id)) {
//				old_id = duplicates[old_id];
//			}
//			faces[i][j] = used_id_new[old_id];
//		}
//	}
//	//end = clock();
//	//cout << "modify time:" << end - start << "ms" << endl;
//
//	start = clock();
//	for (int i = 0; i < used_id.size(); i++) {
//		vertices[used_id_new[used_id[i]]] = vertices[used_id[i]];
//	}
//	while (vertices.size() > used_id.size()) {
//		vertices.pop_back();
//	}
//	end = clock();
//	//cout << "move time:" << end - start << "ms" << endl;
//
//	//cout << "final vertices:" << vertices.size() << endl;
//}
//
//void Wchar_tToString(std::string& szDst, const wchar_t *wchar)
//{
//	const wchar_t * wText = wchar;
//	DWORD dwNum = WideCharToMultiByte(CP_OEMCP, NULL, wText, -1, NULL, 0, NULL, FALSE);// WideCharToMultiByte的运用
//	char *psText;  // psText为char*的临时数组，作为赋值给std::string的中间变量
//	psText = new char[dwNum];
//	WideCharToMultiByte(CP_OEMCP, NULL, wText, -1, psText, dwNum, NULL, FALSE);// WideCharToMultiByte的再次运用
//	szDst = psText;// std::string赋值
//	delete[]psText;// psText的清除
//}
//
//bool contains(set<int> &small, set<int> &big) {
//	if (small.size() > big.size()) {
//		return false;
//	}
//
//	for (auto id : small) {
//		if (big.count(id) == 0) {
//			return false;
//		}
//	}
//
//	return true;
//}
//
//int main() {
//	ON::Begin();
//
//	ONX_Model model;
//	FILE* archive_fp = ON::OpenFile("full.3dm", "rb");
//	ON_BinaryFile archive(ON::read3dm, archive_fp);
//	bool rc = model.Read(archive);
//	if (rc) {
//		cout << "reading 3dm file success" << endl;
//	} else {
//		cout << "failed" << endl;
//		return -1;
//	}
//	ON::CloseFile(archive_fp);
//
//	map<int, vector<ON_3fPoint>> group_vertices;
//	map<int, vector<vector<int>>> group_faces;
//	
//	map<int, vector<ON_3fPoint>> simplified_group_vertices;
//	map<int, vector<vector<int>>> simplified_group_faces;
//
//	vector<ON_3fPoint> fvs;
//	vector<vector<int>> ffs;
//
//	map<int, set<int>> map_group_objects;
//	vector<int> used_group;
//
//	map<int, int> small_groups;
//
//
//	map<int, vector<ON_3fPoint>> object_vertices;
//	map<int, vector<vector<int>>> object_faces;
//
//	map<int, vector<ON_3fPoint>> simplified_object_vertices;
//	map<int, vector<vector<int>>> simplified_object_faces;
//
//	vector<int> used_object;
//	for (int obj_id = 0; obj_id < model.m_object_table.Count(); obj_id++) {
//		ONX_Model_Object object = model.m_object_table[obj_id];
//
//		auto group_count = object.m_attributes.GroupCount();
//		auto group_list = object.m_attributes.GroupList();
//
//		for (int i = 0; i < group_count; i++) {
//			int group = group_list[i];
//			if (map_group_objects.count(group) == 0) {
//				map_group_objects[group] = set<int>();
//				used_group.push_back(group);
//			}
//			map_group_objects[group].insert(obj_id);
//		}
//	}
//
//	for (int i = 0; i < used_group.size(); i++) {
//		for (int j = 0; j < used_group.size(); j++) {
//			int groupi = used_group[i];
//			int groupj = used_group[j];
//			if (groupj != groupi && small_groups.count(groupj) == 0) {
//				if (contains(map_group_objects[groupi], map_group_objects[groupj])) {
//					small_groups[groupi] = groupi;
//					break;
//				}
//			}
//		}
//	}
//
//	for (int i = 0; i < used_group.size(); i++) {
//		int group = used_group[i];
//		if (small_groups.count(group) != 0) {
//			continue;
//		}
//		group_vertices[group] = vector<ON_3fPoint>();
//		group_faces[group] = vector<vector<int>>();
//		simplified_group_vertices[group] = vector<ON_3fPoint>();
//		simplified_group_faces[group] = vector<vector<int>>();
//	}
//
//	//cout << "group count:" << model.m_group_table.Count() << endl;
//	//cout << "object count:" << model.m_object_table.Count() << endl;
//	//cout << "object count:" << model.m_layer_table.Count() << endl;
//
//	//for (int group_id = 0; group_id < model.m_group_table.Count(); group_id++) {
//	//	auto group = model.m_group_table[group_id];
//	//	string name;
//	//	Wchar_tToString(name, group.GroupName());
//	//	cout << group_id << ":" << name  << endl;
//	//}
//
//	for (int obj_id = 0; obj_id < model.m_object_table.Count(); obj_id++) {
//		ONX_Model_Object object = model.m_object_table[obj_id];
//		const ON_Brep* brep = ON_Brep::Cast(object.m_object);
//		
//		if (brep) {
//			ON_MeshParameters mp;
//			ON_SimpleArray<const ON_Mesh*> meshes;
//			int mc = brep->GetMesh(ON::render_mesh, meshes);
//			int fc = 0;
//			if (mc) {
//				vector<ON_3fPoint> vertices;
//				vector<vector<int>> faces;
//
//				clock_t start, end;
//
//				//start = clock();
//				for (int i = 0; i < meshes.Count(); i++) {
//					auto mesh = meshes[i];
//					append_mesh(vertices, faces, mesh);
//				}
//				//end = clock();
//				//cout << "read mesh:" << end - start << "ms" << endl;
//
//				//start = clock();
//				merge_vertices(vertices, faces);
//				//end = clock();
//				//cout << "merge vertices:" << end - start << "ms" << endl;
//
//
//				//char fname[100];
//				//sprintf(fname, "%d.obj\0", obj_id);
//				//write_obj(vertices, faces, fname);
//
//				//start = clock();
//				auto mesh = openmesh_simplify_mesh(vertices, faces);
//				//sprintf(fname, "%d.simplified.obj\0", obj_id);
//				//OpenMesh::IO::write_mesh(mesh, fname);
//				//end = clock();
//				//cout << "simplify mesh:" << end - start << "ms" << endl;
//
//
//				//start = clock();
//				//append_mesh(fvs, ffs, mesh);
//				//end = clock();
//				//cout << "append mesh:" << end - start << "ms" << endl;
//
//				auto group_count = object.m_attributes.GroupCount();
//				auto group_list = object.m_attributes.GroupList();
//
//				for (int idx = 0; idx < group_count; idx++) {
//					auto group = group_list[idx];
//					if (small_groups.count(group) != 0) {
//						continue;
//					}
//
//					append_mesh(group_vertices[group], group_faces[group], vertices, faces);
//					append_mesh(simplified_group_vertices[group], simplified_group_faces[group], mesh);
//				}
//
//				if (group_count == 0 && vertices.size() != 0 && faces.size() != 0) {
//					used_object.push_back(obj_id);
//					object_vertices[obj_id] = vector<ON_3fPoint>();
//					object_faces[obj_id] = vector<vector<int>>();
//					simplified_object_vertices[obj_id] = vector<ON_3fPoint>();
//					simplified_object_faces[obj_id] = vector<vector<int>>();
//					//cout << "object id: " << obj_id << " not added." << endl;
//					append_mesh(object_vertices[obj_id], object_faces[obj_id], vertices, faces);
//					append_mesh(simplified_object_vertices[obj_id], simplified_object_faces[obj_id], mesh);
//				}
//
//			}
//		}
//	}
//
//	model.Destroy();
//	ON::End();
//
//	make_dir("full");
//	make_dir("simplified");
//
//	FILE* full_files;
//	FILE* simplified_files;
//
//	fopen_s(&full_files, "full_files.txt", "wb");
//
//	fopen_s(&simplified_files, "simplified_files.txt", "wb");
//
//
//	for (int i = 0; i < used_group.size(); i++) {
//		auto group = used_group[i];
//		if (small_groups.count(group) != 0) {
//			continue;
//		}
//
//		if (group_vertices[group].size() < 4 || group_faces[group].size() < 3) {
//			continue;
//		}
//		char filename[100];
//		sprintf(filename, "full\\group_%d.obj\0", group);
//		fwrite(filename, strlen(filename), 1, full_files);
//		fwrite("\n", 1, 1, full_files);
//		write_obj(group_vertices[group], group_faces[group], filename);
//
//		sprintf(filename, "simplified\\group_%d.obj\0", group);
//		fwrite(filename, strlen(filename), 1, simplified_files);
//		fwrite("\n", 1, 1, simplified_files);
//		write_obj(simplified_group_vertices[group], simplified_group_faces[group], filename);
//
//		append_mesh(fvs, ffs, simplified_group_vertices[group], simplified_group_faces[group]);
//	}
//
//	for (int i = 0; i < used_object.size(); i++) {
//		auto obj_id = used_object[i];
//		
//		if (object_vertices[obj_id].size() < 4|| object_faces[obj_id].size() < 3) {
//			continue;
//		}
//
//		char filename[100];
//		sprintf(filename, "full\\object_%d.obj\0", obj_id);
//		fwrite(filename, strlen(filename), 1, full_files);
//		fwrite("\n", 1, 1, full_files);
//		write_obj(object_vertices[obj_id], object_faces[obj_id], filename);
//
//		sprintf(filename, "simplified\\object_%d.obj\0", obj_id);
//		fwrite(filename, strlen(filename), 1, simplified_files);
//		fwrite("\n", 1, 1, simplified_files);
//		write_obj(simplified_object_vertices[obj_id], simplified_object_faces[obj_id], filename);
//
//		append_mesh(fvs, ffs, simplified_object_vertices[obj_id], simplified_object_faces[obj_id]);
//	}
//
//	fclose(full_files);
//	fclose(simplified_files);
//
//	//merge_vertices(fvs, ffs);
//	write_obj(fvs, ffs, "full.simplified.20.obj");
//	
//
//
//	system("pause");
//	return 0;
//}