#pragma once

#include <opennurbs.h>

#include <Precision.hxx>
#include <BSplCLib.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColgp_Array2OfPnt2d.hxx>
#include <TColStd_Array2OfReal.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>

#pragma comment(lib, "TKernel.lib")
#pragma comment(lib, "TKG3d.lib")
#pragma comment(lib, "TKG2d.lib")
#pragma comment(lib, "TKGeomBase.lib")
#pragma comment(lib, "TKMath.lib")



Handle_Geom2d_BSplineCurve convert_opennurbs_2d_curve_to_occt_2d_curve(const ON_NurbsCurve* curve);

Handle_Geom_BSplineCurve convert_opennurbs_curve_to_occt_curve(const ON_NurbsCurve* curve);

Handle_Geom_BSplineSurface convert_opennurbs_surface_to_occt_surface(const ON_NurbsSurface* surface);