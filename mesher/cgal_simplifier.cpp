//#include <iostream>
//#include <fstream>
//#include <vector>
//#include <map>
//#include <set>
//#define ON_DLL_IMPORTS
//
//const float esp = 1e-6;
//
//#include "opennurbs\opennurbs_dynamic_linking.h"
//
//#include <opennurbs.h>
//
//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Polyhedron_3.h>
//#include <CGAL/Surface_mesh.h>
//#include <CGAL/IO/Polyhedron_iostream.h>
//#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
//#include <CGAL/Orthogonal_k_neighbor_search.h>
//#include <CGAL/Kd_tree.h>
//#include <CGAL/Search_traits.h>
//
//// Simplification function
//#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
//#include <CGAL/Surface_mesh_simplification/>
//
//// Stop-condition policy
//#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
//#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
////Placement wrapper
//#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
//#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
//
//#include "neighbor_search.hpp"
//
//typedef CGAL::Simple_cartesian<float> Kernel;
//typedef Kernel::Point_3 Point_3;
//typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
//typedef CGAL::Dimension_tag<3> D;
//typedef CGAL::Search_traits<float, Point, const float*, Construct_coord_iterator, D> Traits;
//typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance> Neighbor_search;
//typedef Neighbor_search::Tree Tree;
//
//
//typedef Polyhedron::HalfedgeDS HDS;
//
//namespace SMS = CGAL::Surface_mesh_simplification;
//
//using namespace std;
//class Mesh_to_polyhedron : public CGAL::Modifier_base<HDS> {
//public:
//	Mesh_to_polyhedron(vector<ON_3fPoint> & points, vector<vector<int>>& faces) : m_points(points), m_faces(faces) {}
//	void operator()(HDS& hds) {
//		// Postcondition: `hds' is a valid polyhedral surface.
//		CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
//		B.begin_surface(m_points.size(), m_faces.size());
//
//		// vertices
//		typedef typename HDS::Vertex::Point Vertex;
//		typedef typename Vertex Point;
//		for (unsigned i = 0; i < m_points.size(); i++) {
//			HDS::Vertex_handle vh = B.add_vertex(Point(m_points[i].x, m_points[i].y, m_points[i].z));
//			//vh->id = i;
//		}
//
//		// triangles
//		for (unsigned i = 0; i < m_faces.size(); i++) {
//			HDS::Face_handle fh = B.begin_facet();
//			B.add_vertex_to_facet(m_faces[i][0]);
//			B.add_vertex_to_facet(m_faces[i][1]);
//			B.add_vertex_to_facet(m_faces[i][2]);
//			B.end_facet();
//			//fh->id = i;
//		}
//
//		B.end_surface();
//	}
//
//private:
//	vector<ON_3fPoint> & m_points;
//	vector<vector<int>>& m_faces;
//};
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
//void merge_vertices(vector<ON_3fPoint>& vertices, vector<vector<int>>& faces) {
//	map<int, int> duplicates;
//	vector<int> used_id;
//	map<int, int> used_id_new;
//	Tree tree;
//
//	Distance tr_dist;
//
//	for (int i = 0; i < vertices.size(); i++) {
//		Point point(vertices[i].x, vertices[i].y, vertices[i].z, i);
//		if (tree.size() == 0) {
//			tree.insert(point);
//			used_id.push_back(i);
//			used_id_new[i] = used_id_new.size();
//		} else {
//			Neighbor_search search(tree, point, 1);
//			auto it = search.begin();
//			if (it != search.end() && tr_dist.inverse_of_transformed_distance(it->second) < esp) {
//				duplicates[i] = it->first.id;
//			} else {
//				tree.insert(point);
//				used_id.push_back(i);
//				used_id_new[i] = used_id_new.size();
//			}
//		}
//	}
//
//	for (int i = 0; i < faces.size(); i++) {
//		for (int j = 0; j < faces[i].size(); j++) {
//			int old_id = faces[i][j];
//			if (duplicates.count(old_id)) {
//				old_id = duplicates[old_id];
//			}
//			faces[i][j] = used_id_new[old_id];
//		}
//	}
//	for (int i = 0; i < used_id.size(); i++) {
//		vertices[used_id_new[used_id[i]]] = vertices[used_id[i]];
//	}
//
//	while (vertices.size() > used_id.size()) {
//		vertices.pop_back();
//	}
//
//	cout << "vertices:" << vertices.size() << endl;
//	cout << "tree size:" << tree.size() << endl;
//}
//
//struct Border_is_constrained_edge_map {
//	const Polyhedron* sm_ptr;
//	typedef boost::graph_traits<Polyhedron>::edge_descriptor key_type;
//	typedef bool value_type;  typedef value_type reference;
//	typedef boost::readable_property_map_tag category;
//
//	Border_is_constrained_edge_map(const Polyhedron& sm)
//		: sm_ptr(&sm) {}
//	friend bool get(Border_is_constrained_edge_map m, const key_type& edge) {
//		return CGAL::is_border(edge, *m.sm_ptr);
//	}
//};
////// Placement class//
//typedef SMS::Constrained_placement<SMS::Midpoint_placement<Polyhedron>, Border_is_constrained_edge_map > Placement;
//
//int main() {
//	ON::Begin();
//
//	ONX_Model model;
//	FILE* archive_fp = ON::OpenFile("small1.3dm", "rb");
//	ON_BinaryFile archive(ON::read3dm, archive_fp);
//	bool rc = model.Read(archive);
//
//	if (rc) {
//		cout << "reading 3dm file success" << endl;
//	} else {
//		cout << "failed" << endl;
//		return -1;
//	}
//
//	ON::CloseFile(archive_fp);
//
//	cout << "group count:" << model.m_group_table.Count() << endl;
//	cout << "number of objects:" << model.m_object_table.Count() << endl;
//	vector<ON_3fPoint> vertices;
//	vector<vector<int>> faces;
//	int counter = 0;
//	for (int obj_id = 0; obj_id < model.m_object_table.Count(); obj_id++) {
//		ONX_Model_Object model_object = model.m_object_table[obj_id];
//		const ON_Brep* brep = ON_Brep::Cast(model_object.m_object);
//
//		if (brep) {
//			ON_SimpleArray<const ON_Mesh*> meshes;
//			int mc = brep->GetMesh(ON::mesh_type::render_mesh, meshes);
//			if (mc != 0) {
//				for (int i = 0; i < meshes.Count(); i++) {
//					append_mesh(vertices, faces, meshes[i]);
//				}
//				counter++;
//				break;
//			}
//		}
//	}
//	cout << "number of objects:" << counter << endl;
//
//	model.Destroy();
//	ON::End();
//
//	merge_vertices(vertices, faces);
//
//	write_obj(vertices, faces, "small1_render_mesh_merge_vertices.obj");
//
//
//	Polyhedron polyhedron;
//
//	Mesh_to_polyhedron m2p(vertices, faces);
//	polyhedron.delegate(m2p);
//
//
//	cout << "valid:" << polyhedron.is_valid() << endl;
//	cout << "closed:" << polyhedron.is_closed() << endl;
//	cout << "vertices: " << polyhedron.size_of_vertices() << endl;
//	cout << "facets: " << polyhedron.size_of_facets() << endl;
//
//	std::map<Polyhedron::Halfedge_handle, std::pair<Point_3, Point_3> >constrained_edges;
//	std::size_t nb_border_edges = 0;
//	for (Polyhedron::Halfedge_iterator hit = polyhedron.halfedges_begin(), hit_end = polyhedron.halfedges_end(); hit != hit_end; ++hit) {
//		if (hit->is_border()) {
//			constrained_edges[hit] = std::make_pair(hit->opposite()->vertex()->point(), hit->vertex()->point());
//			++nb_border_edges;
//		}
//	}
//
//	Border_is_constrained_edge_map bem(polyhedron);
//
//	SMS::Count_stop_predicate<Polyhedron> stop(300);
//
//	int r = SMS::edge_collapse(
//		polyhedron,
//		stop,
//		CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, polyhedron))
//		.halfedge_index_map(get(CGAL::halfedge_external_index, polyhedron))
//		.edge_is_constrained_map(bem).get_placement(Placement(bem)));
//
//	std::cout << "\nFinished...\n" << r << " edges removed.\n" << (polyhedron.size_of_halfedges() / 2) << " final edges.\n";
//
//	std::ofstream os("small1_simplified_render_mesh.off"); os << polyhedron;
//
//	system("pause");
//
//	return 0;
//}
