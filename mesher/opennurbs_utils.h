#pragma once

#define ON_DLL_IMPORTS

#include "opennurbs\opennurbs_dynamic_linking.h"

#include <opennurbs.h>
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <string>
#include <locale>
#include <codecvt>

using namespace std;

//从opennurbs官网教程中复制来的代码
bool IsLayerVisible(const ON_ObjectArray<ON_Layer>& layer_table, int layer_index);

bool IsModelObjectVisible(const ONX_Model& model, const ONX_Model_Object& model_object);

inline bool contains(set<int> &big, set<int> &small) {
	if (small.size() > big.size()) {
		return false;
	}

	for (auto id : small) {
		if (big.count(id) == 0) {
			return false;
		}
	}
		
	return true;
}

inline string wstring_to_string(ON_wString& input) {
	wstring temp;
	for (int i = 0; i < input.Length(); i++) {
		temp.push_back(input[i]);
	}

	wstring_convert<std::codecvt_utf8<wchar_t>> c;

	return c.to_bytes(temp);
}

class OpennurbsGroupInfo {
public:
	OpennurbsGroupInfo(ONX_Model* model) 
		:m_model(model){
		generate();
		compose();
		cout << "root groups:" << endl;
		for (auto group_id : root_group_ids) {
			cout << group_id << " ";
		}
		cout << endl;
	}
	vector<int> object_ids;
	map<int, set<int>> group_id_objects;
	map<int, set<int>> object_id_groups;
	set<int> group_ids;
	map<int, set<int>> group_id_groups;
	set<int> root_group_ids;
	map<int, set<int>> group_composed_of_groups;
	map<int, set<int>> group_composed_of_objects;

	map<int, string> map_group_id_name;
	map<int, string> map_object_id_name;
	map<int, string> map_object_id_uuid;

	void output() {
		for (auto & x: group_ids) {
			cout << "group id: " << x << endl;
			cout << "group ids: ";
			for (auto y : group_composed_of_groups[x]) {
				cout << y <<" ";
			}
			cout << endl;
			cout << "object ids: ";
			for (auto y : group_composed_of_objects[x]) {
				cout << y << " ";
			}
			cout << endl;
			cout << endl;
		}
	}

	void release_mem() {
		object_ids.swap(vector<int>());
		group_id_objects.swap(map<int, set<int>>());
		object_id_groups.swap(map<int, set<int>>());
		group_ids.swap(set<int>());
		group_id_groups.swap(map<int, set<int>>());
		root_group_ids.swap(set<int>());
		group_composed_of_groups.swap(map<int, set<int>>());
		group_composed_of_objects.swap(map<int, set<int>>());
	}

private:
	void generate() {
		//for (int group_id = 0; group_id < m_model->m_group_table.Count(); group_id++) {
		//	ON_Group& group = m_model->m_group_table[group_id];
		//}
		int max_group_size = 0;
		for (int object_id = 0; object_id < m_model->m_object_table.Count(); object_id++) {
			auto object = m_model->m_object_table[object_id];
			if (!IsModelObjectVisible(*m_model, object)) {
				continue;
			}

			const ON_Brep* brep = ON_Brep::Cast(object.m_object);
			if (brep == nullptr) {
				continue;
			}

			object_ids.push_back(object_id);
			object_id_groups[object_id] = set<int>();
			auto name = wstring_to_string(object.m_attributes.m_name);
			char uuid[37];
			ON_UuidToString(object.m_attributes.m_uuid, uuid);
			map_object_id_name[object_id] = name;
			map_object_id_uuid[object_id] = uuid;
			auto group_count = object.m_attributes.GroupCount();
			auto group_list = object.m_attributes.GroupList();

			for (int i = 0; i < group_count; i++) {
				int group_id = group_list[i];
				if (group_id_objects.count(group_id) == 0) {
					group_id_objects[group_id] = set<int>();
					group_ids.insert(group_id);
				}
				group_id_objects[group_id].insert(object_id);
				object_id_groups[object_id].insert(group_id);
			}
		}

		//cout << "object_ids.size() = "<< object_ids .size()<<endl;

		//移除重复的group
		vector<int> groups;
		set<int> removed;
		for (auto group_id : group_ids) {
			groups.push_back(group_id);
		}

		for (int i = 0; i < groups.size(); i++) {
			if (removed.count(i) != 0) {
				continue;
			}
			for (int j = i + 1; j < groups.size(); j++) {
				if (removed.count(j) != 0) {
					continue;
				}
 				if (group_id_objects[i].size() == group_id_objects[j].size()
					&& contains(group_id_objects[i], group_id_objects[j])) {
					removed.insert(j);
				}
			}
		}

		for (auto group_id : removed) {
			group_ids.erase(group_id);
			for (auto object_id : group_id_objects[group_id]) {
				object_id_groups[object_id].erase(object_id);
			}
			group_id_objects.erase(group_id);
		}

		for (auto&x : group_id_objects) {
			if (x.second.size() > max_group_size) {
				max_group_size = x.second.size();
			}
		}
		if (max_group_size < object_ids.size()) {
			group_ids.insert(-1);
			group_id_objects[-1] = set<int>();
			for (auto obj_id : object_ids) {
				group_id_objects[-1].insert(obj_id);
				object_id_groups[obj_id].insert(-1);
			}
		}

		root_group_ids = group_ids;
		for (auto group_id_a : group_ids) {
			for (auto group_id_b : group_ids) {
				if (group_id_a != group_id_b) {
					if (contains(group_id_objects[group_id_a], group_id_objects[group_id_b])) {
						root_group_ids.erase(group_id_b);
						if (group_id_groups.count(group_id_a) == 0) {
							group_id_groups[group_id_a] = set<int>();
						}
						group_id_groups[group_id_a].insert(group_id_b);
					}
				}
			}
		}

	}

	void compose() {
		for (auto group_id : group_ids) {
			group_composed_of_groups[group_id] = set<int>();
			group_composed_of_objects[group_id] = set<int>();
			if (group_id_groups.count(group_id) != 0) {
				auto candidates = group_id_groups[group_id];
				auto obj_ids = group_id_objects[group_id];
				while (candidates.size() != 0 && obj_ids.size() != 0) {
					auto candidate_id = largest_group(candidates);
					candidates.erase(candidate_id);
					if (contains(obj_ids, group_id_objects[candidate_id])) {
						for (auto obj_id : group_id_objects[candidate_id]) {
							obj_ids.erase(obj_id);
						}
						group_composed_of_groups[group_id].insert(candidate_id);
					}
				}
				for (auto obj_id : obj_ids) {
					group_composed_of_objects[group_id].insert(obj_id);
				}
			}
			else {
				for (auto obj_id : group_id_objects[group_id]) {
					group_composed_of_objects[group_id].insert(obj_id);
				}
			}
		}
	}

	int largest_group(set<int> group_ids) {
		if (group_ids.size() == 0) {
			return -1;
		}

		int id = -1;
		int size = 0;
		for (auto group_id : group_ids) {
			if (group_id_objects[group_id].size() > size) {
				size = group_id_objects[group_id].size();
				id = group_id;
			}
		}
		return id;
	}
	ONX_Model* m_model;
};