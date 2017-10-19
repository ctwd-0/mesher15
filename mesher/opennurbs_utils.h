#pragma once

#include <opennurbs.h>


//从opennurbs官网教程中复制来的代码
bool IsLayerVisible(const ON_ObjectArray<ON_Layer>& layer_table, int layer_index);

bool IsModelObjectVisible(const ONX_Model& model, const ONX_Model_Object& model_object);