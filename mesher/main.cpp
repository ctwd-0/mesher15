//#include <iostream>
//#include <vector>
//#include <map>
//#include <set>
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
//typedef OpenMesh::TriMesh_ArrayKernelT<> TriMesh;
//
//typedef OpenMesh::Decimater::DecimaterT< TriMesh > MyDecimater;
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
//		sprintf(line, "f %d %d %d\n\0", faces[i][0], faces[i][1], faces[i][2]);
//		fwrite(line, strlen(line), 1, f);
//	}
//
//	fclose(f);
//}
//
//void write_obj(TriMesh& mesh, const char* file_name) {
//	FILE* f;
//	auto error = fopen_s(&f, file_name, "wb");
//
//	map<int, int> map_vid;
//
//	char line[100] = { 0 };
//
//	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
//		auto d = mesh.point(v_it).data();
//		sprintf(line, "v %lf %lf %lf\n\0", d[0], d[1], d[2]);
//		fwrite(line, strlen(line), 1, f);
//		if (map_vid.count(v_it->idx()) == 0) {
//			map_vid[v_it->idx()] = map_vid.size() + 1;
//		}
//	}
//
//	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
//		vector<int> ids;
//
//		for (auto fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
//			ids.push_back(map_vid[fv_it->idx()]);
//		}
//
//		sprintf(line, "f %d %d %d\n\0", ids[0], ids[1], ids[2]);
//		fwrite(line, strlen(line), 1, f);
//	}
//
//	fclose(f);
//}
//
//
//void protectMeshBoundariesFromDecimation(TriMesh& mesh, double size) {
//	mesh.request_vertex_status();
//
//	for (const auto& halfEdgeHandle : mesh.halfedges()) {
//		if (mesh.is_boundary(halfEdgeHandle)) {
//			if (mesh.calc_edge_length(halfEdgeHandle) > size * 0.01) {
//				mesh.status(mesh.to_vertex_handle(halfEdgeHandle)).set_locked(true);
//				mesh.status(mesh.from_vertex_handle(halfEdgeHandle)).set_locked(true);
//			}
//		}
//	}
//}
//
//int hc = 0;
//map<int, vector<vector<ON_3fPoint>>> detect_small_holes(ON_SimpleArray<const ON_Mesh*> &meshes, double size) {
//	map<int, vector<vector<ON_3fPoint>>> result;
//	for (int i = 0; i < meshes.Count(); i++) {
//		auto mesh = meshes[i];
//		auto vc = mesh->m_V.Count();
//		auto fc = mesh->m_F.Count();
//		if (vc >= 10 && vc % 2 == 0 && vc == fc) {
//			//cout << "hole" << endl;
//			hc++;
//			ON_3fPoint p;
//			for (int i = 0; i < mesh->m_V.Count(); i++) {
//				p += mesh->m_V[i];
//			}
//			p.x /= vc;
//			p.y /= vc;
//			p.z /= vc;
//
//			auto dis = p.DistanceTo(mesh->m_V[0]);
//
//			if (size * 0.01 < dis) {
//				//cout << "size" << endl;
//				continue;
//			}
//			bool equals = true;
//			for (int i = 1; i < mesh->m_V.Count(); i++) {
//				auto dd = p.DistanceTo(mesh->m_V[i]);
//				if (abs(dis - dd) > 0.001) {
//					equals = false;
//					break;
//				}
//			}
//			if (!equals){
//				//cout << "equals" << endl; 
//				continue;
//			}
//
//			TriMesh hmesh;
//			for (int i = 0; i < mesh->m_V.Count(); i++) {
//				hmesh.add_vertex(Vec3f(mesh->m_V[i].x, mesh->m_V[i].y, mesh->m_V[i].z));
//			}
//
//			for (int i = 0; i < mesh->m_F.Count(); i++) {
//				if (mesh->m_F[i].vi[2] == mesh->m_F[i].vi[3]) {
//					hmesh.add_face(VertexHandle(mesh->m_F[i].vi[0]), VertexHandle(mesh->m_F[i].vi[1]), VertexHandle(mesh->m_F[i].vi[2]));
//				} else {
//					hmesh.add_face(VertexHandle(mesh->m_F[i].vi[0]), VertexHandle(mesh->m_F[i].vi[1]), VertexHandle(mesh->m_F[i].vi[2]));
//					hmesh.add_face(VertexHandle(mesh->m_F[i].vi[0]), VertexHandle(mesh->m_F[i].vi[2]), VertexHandle(mesh->m_F[i].vi[3]));
//				}
//			}
//
//			bool is_loop_topology = true;
//			for (auto ff = hmesh.faces_begin(); ff != hmesh.faces_end(); ++ff) {
//				int no_boundary = 0;
//				int boundary = 0;
//				for (auto fe = hmesh.fe_iter(ff); fe.is_valid(); ++fe) {
//					if (hmesh.is_boundary(fe)) {
//						boundary++;
//					} else {
//						no_boundary++;
//					}
//				}
//				if (!(boundary == 1 && no_boundary == 2)) {
//					is_loop_topology = false;
//					break;
//				}
//			}
//
//			if (!is_loop_topology) {
//				//cout << "topology" << endl;
//				continue;
//			}
//			
//			set<int> set_a;
//			set<int> set_b;
//
//			set<int> used_face;
//
//			auto ff_id = hmesh.faces_begin().handle();
//
//			used_face.insert(ff_id.idx());
//
//			int next_face = -1;
//
//			for (auto hfe = hmesh.fh_begin(ff_id); hfe.is_valid(); ++hfe) {
//				auto ohe = hmesh.opposite_halfedge_handle(hfe);
//				if (hmesh.is_boundary(ohe)) {
//					set_a.insert(hmesh.opposite_vh(hfe).idx());
//				} else {
//					set_b.insert(hmesh.opposite_vh(hfe).idx());
//					auto id = hmesh.opposite_face_handle(hmesh.opposite_halfedge_handle(hfe)).idx();
//					if (used_face.count(id) == 0) {
//						next_face = id;
//					}
//				}
//			}
//
//			while (next_face != -1 && next_face != hmesh.faces_begin().handle().idx()){
//				set_a.swap(set_b);
//				ff_id = FaceHandle(next_face);
//				used_face.insert(next_face);
//
//				for (auto hfe = hmesh.fh_begin(ff_id); hfe.is_valid(); ++hfe) {
//					auto ohe = hmesh.opposite_halfedge_handle(hfe);
//					if (hmesh.is_boundary(ohe)) {
//						set_a.insert(hmesh.opposite_vh(hfe).idx());
//					} else {
//						set_b.insert(hmesh.opposite_vh(hfe).idx());
//						auto id = hmesh.opposite_face_handle(hmesh.opposite_halfedge_handle(hfe)).idx();
//						if (used_face.count(id) == 0) {
//							next_face = id;
//						}
//					}
//				}
//			}
//
//			result[i] = vector<vector<ON_3fPoint>>();
//			result[i].push_back(vector<ON_3fPoint>());
//			result[i].push_back(vector<ON_3fPoint>());
//			ON_3fPoint mid_a(0, 0, 0);
//			ON_3fPoint mid_b(0, 0, 0);
//			for (auto id : set_a) {
//				result[i][0].push_back(mesh->m_V[id]);
//				mid_a += mesh->m_V[id];
//			}
//			mid_a.x /= set_a.size();
//			mid_a.y /= set_a.size();
//			mid_a.z /= set_a.size();
//			result[i][0].push_back(mid_a);
//
//			for (auto id : set_b) {
//				result[i][1].push_back(mesh->m_V[id]);
//				mid_b += mesh->m_V[id];
//			}
//			mid_b.x /= set_a.size();
//			mid_b.y /= set_a.size();
//			mid_b.z /= set_a.size();
//			result[i][1].push_back(mid_b);
//		}
//	}
//
//	//cout << result.size() <<endl;
//	return result;
//}
//
//int test_remove(map<int, vector<vector<ON_3fPoint>>>& result, ON_3fPoint point) {
//	for (auto & kvp : result) {
//		int i = kvp.first;
//		for (int j = 0; j < 2; j++) {
//			for (int k = 0; k < kvp.second[j].size() - 1; k++) {
//				if (kvp.second[j][k].DistanceTo(point) < 0.0001) {
//					return i * 2 + j;
//				}
//			}
//		}
//	}
//	return -1;
//}
//int main() {
//	ON::Begin();
//
//	ONX_Model model;
//	FILE* archive_fp = ON::OpenFile("small1.3dm", "rb");
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
//	int mesh_size = 0;
//	int mesh_count = 0;
//	int max_fc = 0;
//	int min_fc = INT_MAX;
//	vector<ON_3fPoint> vertices;
//	vector<vector<int>> faces;
//	//for (int i = 0; i < model.m_layer_table.Count(); i++) {
//	//	auto layer_object = model.m_layer_table[i];
//
//	//}
//	auto bbox = model.BoundingBox();
//	auto size = bbox.m_max.x - bbox.m_min.x;
//	if (size > bbox.m_max.y - bbox.m_min.y)
//		size = bbox.m_max.y - bbox.m_min.y;
//	if (size > bbox.m_max.z - bbox.m_min.z)
//		size = bbox.m_max.z - bbox.m_min.z;
//
//	for (int obj_id = 0; obj_id < model.m_object_table.Count(); obj_id++) {
//		ONX_Model_Object model_object = model.m_object_table[obj_id];
//		const ON_Brep* brep = ON_Brep::Cast(model_object.m_object);
//
//		//NurbsSurface_obj obj;
//
//		if (brep) {
//			ON_MeshParameters mp;
//			ON_SimpleArray<const ON_Mesh*> meshes;
//			int mc = brep->GetMesh(ON::render_mesh, meshes);
//			int fc = 0;
//			if (mc) {
//				//cout<<meshes.Count()<<endl;
//				mesh_size += mc;
//				mesh_count++;
//				
//				auto result = detect_small_holes(meshes,size);
//
//				map<int, int> global;
//				for (int i = 0; i < meshes.Count(); i++) {
//					if (result.count(i) == 0) {
//						continue;
//					}
//
//					for (int j = 0; j < 2; j++) {
//						int key = i * 2 + j;
//						int value = vertices.size() + 1;
//						vertices.push_back(result[i][j].back());
//						global[key] = value;
//					}
//				}
//
//				for (int i = 0; i < meshes.Count(); i++) {
//					if (result.count(i)) {
//						continue;
//					}
//
//					auto mesh = meshes[i];
//					map<int, int> map_index;
//					vector<int> used_point;
//					int offset = vertices.size();
//					fc += mesh->m_F.Count();
//					for (int j = 0; j < mesh->m_F.Count(); j++) {
//						auto face = mesh->m_F[j];
//						for (int k = 0; k < 4; k++) {
//							auto index = face.vi[k];
//							if (map_index.count(index) == 0) {
//								auto point = mesh->m_V[index];
//								auto removebyid = test_remove(result, point);
//								if (removebyid != -1) {
//									map_index[index] = removebyid;
//								} else {
//									map_index[index] = offset + map_index.size() + 1;
//									used_point.push_back(index);
//								}
//								
//							}
//						}
//						if (face.vi[2] == face.vi[3]) {
//							vector<int> fi;
//							fi.push_back( map_index[face.vi[0]]);
//							fi.push_back( map_index[face.vi[1]]);
//							fi.push_back( map_index[face.vi[2]]);
//							if(fi[0] != fi[1] && fi[0] != fi[2] && fi[1]!= fi[2])
//								faces.push_back(fi);
//						} else {
//							vector<int> f1;
//							f1.push_back( map_index[face.vi[0]]);
//							f1.push_back( map_index[face.vi[1]]);
//							f1.push_back( map_index[face.vi[2]]);
//							if (f1[0] != f1[1] && f1[0] != f1[2] && f1[1] != f1[2])
//								faces.push_back(f1);
//							vector<int> f2;
//							f2.push_back( map_index[face.vi[0]]);
//							f2.push_back( map_index[face.vi[2]]);
//							f2.push_back( map_index[face.vi[3]]);
//							if (f2[0] != f2[1] && f2[0] != f2[2] && f2[1] != f2[2])
//								faces.push_back(f2);
//						}
//					}
//
//					for (int j = 0; j < used_point.size(); j++) {
//						vertices.push_back(mesh->m_V[used_point[j]]);
//					}
//				}
//				//cout << fc << endl;
//				if (fc < min_fc) {
//					min_fc = fc;
//				}
//				if (fc > max_fc) {
//					max_fc = fc;
//				}
//			}
//		}
//	}
//
//	model.Destroy();
//	ON::End();
//
//	cout << "mesh_count:" << mesh_count << endl;
//	cout << "mesh_size:" << mesh_size << endl;
//
//	cout << "vertices:" << vertices.size() << endl;
//	cout << "faces:" << faces.size() << endl;
//
//	cout << "max_fc:" << max_fc << endl;
//	cout << "min_fc:" << min_fc << endl;
//	//cout << "hc:" << hc << endl;
//
//	//write_obj(vertices, faces, "full.obj");
//
//	cout << "begin push mesh" << endl;
//	TriMesh mesh;
//	for (int i = 0; i < vertices.size(); i++) {
//		mesh.add_vertex(Vec3f(vertices[i].x, vertices[i].y, vertices[i].z));
//	}
//
//	for (int i = 0; i < faces.size(); i++) {
//		mesh.add_face(VertexHandle(faces[i][0] - 1), VertexHandle(faces[i][1] - 1), VertexHandle(faces[i][2] - 1));
//	}
//	
//	write_obj(mesh, "small1.obj");
//	cout << "end push mesh" << endl;
//
//	//protectMeshBoundariesFromDecimation(mesh, size);
//
//	MyDecimater decimater(mesh);
//	HModQuadric hModQuadric;
//
//	decimater.add(hModQuadric);
//
//	cout << "begin simplify mesh" << endl;
//
//	decimater.initialize();
//	decimater.decimate();
//
//	mesh.garbage_collection();
//
//	cout << "end simplify mesh" << endl;
//
//	write_obj(mesh, "simplified_small1.obj");
//
//	//OpenMesh::IO::write_mesh(mesh, "mesIO_simplified_mesh_remove_holes.obj");
//
//	system("pause");
//
//	return 0;
//}