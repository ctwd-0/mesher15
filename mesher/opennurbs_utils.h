#pragma once

#include <opennurbs.h>
#include <map>
#include <set>
#include <vector>

using namespace std;

//从opennurbs官网教程中复制来的代码
bool IsLayerVisible(const ON_ObjectArray<ON_Layer>& layer_table, int layer_index);

bool IsModelObjectVisible(const ONX_Model& model, const ONX_Model_Object& model_object);

static inline bool contains(set<int> &big, set<int> &small) {
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

class OpennurbsGroupInfo {
public:
	OpennurbsGroupInfo(ONX_Model* model) 
		:m_model(model){
		generate();
		compose();
	}
	map<int, set<int>> group_id_objects;
	set<int> group_ids;
	map<int, set<int>> group_id_groups;
	set<int> root_group_ids;
	map<int, set<int>> group_composed_of_groups;
	map<int, set<int>> group_composed_of_objects;

private:
	void generate() {
		for (int group_id = 0; group_id < m_model->m_group_table.Count(); group_id++) {
			ON_Group& group = m_model->m_group_table[group_id];
		}
		
		for (int object_id = 0; object_id < m_model->m_object_table.Count(); object_id++) {
			ONX_Model_Object object = m_model->m_object_table[object_id];

			auto group_count = object.m_attributes.GroupCount();
			auto group_list = object.m_attributes.GroupList();

			for (int i = 0; i < group_count; i++) {
				int group_id = group_list[i];
				if (group_id_objects.count(group_id) == 0) {
					group_id_objects[group_id] = set<int>();
					group_ids.insert(group_id);
				}
				group_id_objects[group_id].insert(object_id);
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
			if (group_id_groups.count(group_id) != 0) {
				group_composed_of_groups[group_id] = set<int>();
				group_composed_of_objects[group_id] = set<int>();
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