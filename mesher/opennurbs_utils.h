#pragma once

#include <opennurbs.h>


//��opennurbs�����̳��и������Ĵ���
bool IsLayerVisible(const ON_ObjectArray<ON_Layer>& layer_table, int layer_index);

bool IsModelObjectVisible(const ONX_Model& model, const ONX_Model_Object& model_object);