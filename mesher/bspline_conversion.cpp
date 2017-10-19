#include "bspline_conversion.h"

Handle_Geom2d_BSplineCurve convert_opennurbs_2d_curve_to_occt_2d_curve(const ON_NurbsCurve* curve) {
	TColgp_Array1OfPnt2d poles(1, curve->CVCount());
	TColStd_Array1OfReal weights(1, curve->CVCount());

	TColStd_Array1OfReal knot_sequence(1, curve->KnotCount() + 2);

	// Control point and its weight.
	for (int i = 0; i < curve->CVCount(); ++i) {
		if (curve->IsRational()) {
			ON_4dPoint aPole;
			curve->GetCV(i, aPole);
			poles.SetValue(i + 1, gp_Pnt2d(aPole.x / aPole.w, aPole.y / aPole.w));
			weights.SetValue(i + 1, aPole.w);
		}
		else {
			ON_3dPoint aPole;
			curve->GetCV(i, aPole);
			poles.SetValue(i + 1, gp_Pnt2d(aPole.x, aPole.y));
			weights.SetValue(i + 1, 1.);
		}
	}

	// Knot vector and its multiplicity.
	for (int i = 0; i < curve->KnotCount(); ++i) {
		knot_sequence.SetValue(i + 2, curve->Knot(i));
	}

	knot_sequence.SetValue(knot_sequence.Lower(), curve->Knot(0));
	knot_sequence.SetValue(knot_sequence.Upper(), curve->Knot(curve->KnotCount() - 1));

	TColStd_Array1OfReal knots(1, BSplCLib::KnotsLength(knot_sequence, curve->IsPeriodic()));
	TColStd_Array1OfInteger muls(1, knots.Upper());

	BSplCLib::Knots(knot_sequence, knots, muls);

	Handle_Geom2d_BSplineCurve occt_curve = new Geom2d_BSplineCurve(
		poles, weights, knots, muls,
		curve->Degree(), curve->IsPeriodic());

	return occt_curve;
}

Handle_Geom_BSplineCurve convert_opennurbs_curve_to_occt_curve(const ON_NurbsCurve* curve) {
	TColgp_Array1OfPnt poles(1, curve->CVCount());
	TColStd_Array1OfReal weights(1, curve->CVCount());

	TColStd_Array1OfReal knot_sequence(1, curve->KnotCount() + 2);

	// Control point and its weight.
	for (int i = 0; i < curve->CVCount(); ++i) {
		if (curve->IsRational()) {
			ON_4dPoint aPole;
			curve->GetCV(i, aPole);
			poles.SetValue(i + 1, gp_Pnt(aPole.x / aPole.w, aPole.y / aPole.w, aPole.z / aPole.w));
			weights.SetValue(i + 1, aPole.w);
		}
		else {
			ON_3dPoint aPole;
			curve->GetCV(i, aPole);
			poles.SetValue(i + 1, gp_Pnt(aPole.x, aPole.y, aPole.z));
			weights.SetValue(i + 1, 1.);
		}
	}

	// Knot vector and its multiplicity.
	for (int i = 0; i < curve->KnotCount(); ++i) {
		knot_sequence.SetValue(i + 2, curve->Knot(i));
	}

	knot_sequence.SetValue(knot_sequence.Lower(), curve->Knot(0));
	knot_sequence.SetValue(knot_sequence.Upper(), curve->Knot(curve->KnotCount() - 1));

	TColStd_Array1OfReal knots(1, BSplCLib::KnotsLength(knot_sequence, curve->IsPeriodic()));
	TColStd_Array1OfInteger muls(1, knots.Upper());

	BSplCLib::Knots(knot_sequence, knots, muls);

	Handle_Geom_BSplineCurve occt_curve = new Geom_BSplineCurve(
		poles, weights, knots, muls,
		curve->Degree(), curve->IsPeriodic());

	return occt_curve;
}

Handle_Geom_BSplineSurface convert_opennurbs_surface_to_occt_surface(const ON_NurbsSurface* surface) {
	TColgp_Array2OfPnt poles(1, surface->CVCount(0), 1, surface->CVCount(1));
	TColStd_Array2OfReal weights(1, surface->CVCount(0), 1, surface->CVCount(1));

	TColStd_Array1OfReal u_knot_sequence(1, surface->KnotCount(0) + 2);
	TColStd_Array1OfReal v_knot_sequence(1, surface->KnotCount(1) + 2);

	bool is_u_periodic = (bool)surface->IsPeriodic(0);
	bool is_v_periodic = (bool)surface->IsPeriodic(1);

	// control point and its weight.
	for (int i = 0; i < surface->CVCount(0); ++i) {
		for (int j = 0; j < surface->CVCount(1); ++j) {
			if (surface->IsRational()) {
				ON_4dPoint aPole;
				surface->GetCV(i, j, aPole);
				poles.SetValue(i + 1, j + 1, gp_Pnt(aPole.x / aPole.w, aPole.y / aPole.w, aPole.z / aPole.w));
				weights.SetValue(i + 1, j + 1, aPole.w);
			}
			else {
				ON_3dPoint aPole;
				surface->GetCV(i, j, aPole);
				poles.SetValue(i + 1, j + 1, gp_Pnt(aPole.x, aPole.y, aPole.z));
				weights.SetValue(i + 1, j + 1, 1);
			}
		}
	}

	// Knot vector and its multiplicity.
	for (int i = 0; i < surface->KnotCount(0); ++i) {
		u_knot_sequence.SetValue(i + 2, surface->Knot(0, i));
	}
	u_knot_sequence.SetValue(u_knot_sequence.Lower(), surface->Knot(0, 0));
	u_knot_sequence.SetValue(u_knot_sequence.Upper(), surface->Knot(0, surface->KnotCount(0) - 1));
	TColStd_Array1OfReal u_knots(1, BSplCLib::KnotsLength(u_knot_sequence, is_u_periodic));
	TColStd_Array1OfInteger u_mults(1, u_knots.Upper());
	BSplCLib::Knots(u_knot_sequence, u_knots, u_mults);

	for (int j = 0; j < surface->KnotCount(1); ++j) {
		v_knot_sequence.SetValue(j + 2, surface->Knot(1, j));
	}
	v_knot_sequence.SetValue(v_knot_sequence.Lower(), surface->Knot(1, 0));
	v_knot_sequence.SetValue(v_knot_sequence.Upper(), surface->Knot(1, surface->KnotCount(1) - 1));
	TColStd_Array1OfReal v_knots(1, BSplCLib::KnotsLength(v_knot_sequence, is_v_periodic));
	TColStd_Array1OfInteger v_mults(1, v_knots.Upper());
	BSplCLib::Knots(v_knot_sequence, v_knots, v_mults);

	Handle_Geom_BSplineSurface occt_surface = new Geom_BSplineSurface(poles,
		weights, u_knots, v_knots, u_mults, v_mults,
		surface->Degree(0), surface->Degree(1),
		is_u_periodic, is_v_periodic);

	return occt_surface;
}
