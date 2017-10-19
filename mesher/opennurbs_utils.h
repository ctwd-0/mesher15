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
	OpennurbsGroupInfo(ONX_Model& model) 
		:m_model(model){
		generate();
	}
	map<int, set<int>> map_group_id_objects;
	set<int> group_ids;
	map<int, set<int>> map_group_id_groups;
	set<int> root_group_ids;
private:
	void generate() {
		for (int group_id = 0; group_id < m_model.m_group_table.Count(); group_id++) {
			ON_Group& group = m_model.m_group_table[group_id];
		}
		
		for (int object_id = 0; object_id < m_model.m_object_table.Count(); object_id++) {
			ONX_Model_Object object = m_model.m_object_table[object_id];

			auto group_count = object.m_attributes.GroupCount();
			auto group_list = object.m_attributes.GroupList();

			for (int i = 0; i < group_count; i++) {
				int group_id = group_list[i];
				if (map_group_id_objects.count(group_id) == 0) {
					map_group_id_objects[group_id] = set<int>();
					group_ids.insert(group_id);
				}
				map_group_id_objects[group_id].insert(object_id);
			}
		}

		root_group_ids = group_ids;
		for (auto group_id_a : group_ids) {
			for (auto group_id_b : group_ids) {
				if (group_id_a != group_id_b) {
					if (contains(map_group_id_objects[group_id_a], map_group_id_objects[group_id_b])) {
						root_group_ids.erase(group_id_b);
						if (map_group_id_groups.count(group_id_a) == 0) {
							map_group_id_groups[group_id_a] = set<int>();
						}
						map_group_id_groups[group_id_a].insert(group_id_b);
					}
				}
			}
		}

	}
	ONX_Model& m_model;
};