#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <time.h>

#define ON_DLL_IMPORTS

#include "opennurbs\opennurbs_dynamic_linking.h"

#include "opennurbs_utils.h"

#include "obj_model_output.h"

#include "bspline_conversion.h"

#include "brep_conversion.h"

#include "folder.h"

#include <Precision.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <BRepTools.hxx>

#include <TopExp_Explorer.hxx>
#include <Poly_Triangulation.hxx>

#include <BRepMesh.hxx>
#include <BRepMesh_IncrementalMesh.hxx>

#pragma comment(lib, "TKTopAlgo.lib")
#pragma comment(lib, "TKMesh.lib")

using namespace std;


int main() {

	string case_name;

	case_name = "full";

	ON::Begin();

	ONX_Model model;
	FILE* archive_fp = ON::OpenFile((case_name + ".3dm").c_str(), "rb");
	ON_BinaryFile archive(ON::read3dm, archive_fp);
	bool rc = model.Read(archive);
	if (rc) {
		cout << "reading 3dm file success" << endl;
	}
	else {
		cout << "failed" << endl;
		system("pause");
		return -1;
	}
	ON::CloseFile(archive_fp);

	BRep_Builder occt_brep_builder;
	TopoDS_Compound occt_compound_all;
	occt_brep_builder.MakeCompound(occt_compound_all);

	for (int obj_id = 0; obj_id < model.m_object_table.Count(); obj_id++) {
		ONX_Model_Object object = model.m_object_table[obj_id];

		if (!IsModelObjectVisible(model, object)) {
			continue;
		}
		const ON_Brep* brep = ON_Brep::Cast(object.m_object);
		if (brep != nullptr) {

			auto occt_compound = convert_opennurbs_solid_to_occt_solid(brep);
			occt_brep_builder.Add(occt_compound_all, occt_compound);
		}
	}

	cout << "convert done" << endl;

	ON::End();
	
	BRepTools::Write(occt_compound_all, (case_name + ".brep").c_str());
	

	Bnd_Box b;

	auto occt_mesh_10 = BRepMesh_IncrementalMesh(occt_compound_all, 10, false, 5, true, false);

	if (occt_mesh_10.IsDone()) {
		cout << "mesh done" << endl;
	}

	vector<gp_Pnt> vs;
	vector<vector<int>> fs;

	int counter = 0;
	for (TopExp_Explorer face_exp(occt_compound_all, TopAbs_FACE); face_exp.More(); face_exp.Next()) {
		TopLoc_Location location;
		auto &face = TopoDS::Face(face_exp.Current());

		auto triangulation = BRep_Tool::Triangulation(face, location);

		counter++;
		if (!triangulation.IsNull()) {
			append_mesh(vs, fs, triangulation);
		}
		else {
			cout <<counter << " is null" << endl;
		}
	}
	
	write_obj(vs, fs, (case_name + ".obj").c_str());

	ON::End();
	system("pause");
}