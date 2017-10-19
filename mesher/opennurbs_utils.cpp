#include "opennurbs_utils.h"

//从opennurbs官网教程中复制来的代码
 bool IsLayerVisible(const ON_ObjectArray<ON_Layer>& layer_table, int layer_index)
 {
	bool rc = false;
	// Validate the layer index
	if (layer_index >= 0 && layer_index < layer_table.Count())
	{
		// Get the layer
		const ON_Layer& layer = layer_table[layer_index];
		// Get the layer's visibility
		rc = layer.IsVisible();
		// If the layer is visible, see if the layer has a parent. If so,
		// check to see if the layer's parent is visible. If not, then
		// the layer is also not visible.
		if (rc && ON_UuidIsNotNil(layer.m_parent_layer_id))
		{
			int i, layer_count = layer_table.Count();
			for (i = 0; i < layer_count; i++)
			{
				if (0 == ON_UuidCompare(layer.m_parent_layer_id, layer_table[i].m_layer_id))
					return IsLayerVisible(layer_table, i); // recursive
			}
		}
	}
	return rc;
}

bool IsModelObjectVisible(const ONX_Model& model, const ONX_Model_Object& model_object)
{
	bool rc = false;
	switch (model_object.m_attributes.Mode())
	{
	case ON::normal_object:
	case ON::idef_object:
	case ON::locked_object:
	{
		// Get the object's layer
		int layer_index = model_object.m_attributes.m_layer_index;
		// See if the layer is visible
		rc = IsLayerVisible(model.m_layer_table, layer_index);
	}
	break;
	}
	return rc;
}