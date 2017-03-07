#pragma once
struct Point {
	int id;
	float vec[3];
	Point() { vec[0] = vec[1] = vec[2] = 0; id = -1;}
	Point(float x, float y, float z, int _id = -1) { vec[0] = x; vec[1] = y; vec[2] = z; id = _id; }
	float x() const { return vec[0]; }
	float y() const { return vec[1]; }
	float z() const { return vec[2]; }
	float& x() { return vec[0]; }
	float& y() { return vec[1]; }
	float& z() { return vec[2]; }
	bool operator==(const Point& p) const
	{
		return (x() == p.x()) && (y() == p.y()) && (z() == p.z());
	}
	bool  operator!=(const Point& p) const { return !(*this == p); }
}; //end of class

namespace CGAL {
	template <>
	struct Kernel_traits<Point> {
		struct Kernel {
			typedef float FT;
			typedef float RT;
		};
	};
}

struct Construct_coord_iterator {
	typedef  const float* result_type;
	const float* operator()(const Point& p) const
	{
		return static_cast<const float*>(p.vec);
	}
	const float* operator()(const Point& p, int)  const
	{
		return static_cast<const float*>(p.vec + 3);
	}
};

struct Distance {
	typedef Point Query_item;
	typedef float FT;
	typedef CGAL::Dimension_tag<3> D;
	float transformed_distance(const Point& p1, const Point& p2) const {
		float distx = p1.x() - p2.x();
		float disty = p1.y() - p2.y();
		float distz = p1.z() - p2.z();
		return distx*distx + disty*disty + distz*distz;
	}
	float min_distance_to_rectangle(const Point& p,
		const CGAL::Kd_tree_rectangle<FT, D>& b) const {
		float distance(0.0), h = p.x();
		if (h < b.min_coord(0)) distance += (b.min_coord(0) - h)*(b.min_coord(0) - h);
		if (h > b.max_coord(0)) distance += (h - b.max_coord(0))*(h - b.max_coord(0));
		h = p.y();
		if (h < b.min_coord(1)) distance += (b.min_coord(1) - h)*(b.min_coord(1) - h);
		if (h > b.max_coord(1)) distance += (h - b.max_coord(1))*(h - b.min_coord(1));
		h = p.z();
		if (h < b.min_coord(2)) distance += (b.min_coord(2) - h)*(b.min_coord(2) - h);
		if (h > b.max_coord(2)) distance += (h - b.max_coord(2))*(h - b.max_coord(2));
		return distance;
	}
	float min_distance_to_rectangle(const Point& p,
		const CGAL::Kd_tree_rectangle<FT, D>& b, std::vector<float>& dists) {
		float distance(0.0), h = p.x();
		if (h < b.min_coord(0)) {
			dists[0] = (b.min_coord(0) - h);
			distance += dists[0] * dists[0];
		}
		if (h > b.max_coord(0)) {
			dists[0] = (h - b.max_coord(0));
			distance += dists[0] * dists[0];
		}
		h = p.y();
		if (h < b.min_coord(1)) {
			dists[1] = (b.min_coord(1) - h);
			distance += dists[1] * dists[1];
		}
		if (h > b.max_coord(1)) {
			dists[1] = (h - b.max_coord(1));
			distance += dists[1] * dists[1];
		}
		h = p.z();
		if (h < b.min_coord(2)) {
			dists[2] = (b.min_coord(2) - h);
			distance += dists[2] * dists[2];
		}
		if (h > b.max_coord(2)) {
			dists[2] = (h - b.max_coord(2));
			distance += dists[2] * dists[2];
		}
		return distance;
	}
	float max_distance_to_rectangle(const Point& p,
		const CGAL::Kd_tree_rectangle<FT, D>& b) const {
		float h = p.x();
		float d0 = (h >= (b.min_coord(0) + b.max_coord(0)) / 2.0) ?
			(h - b.min_coord(0))*(h - b.min_coord(0)) : (b.max_coord(0) - h)*(b.max_coord(0) - h);
		h = p.y();
		float d1 = (h >= (b.min_coord(1) + b.max_coord(1)) / 2.0) ?
			(h - b.min_coord(1))*(h - b.min_coord(1)) : (b.max_coord(1) - h)*(b.max_coord(1) - h);
		h = p.z();
		float d2 = (h >= (b.min_coord(2) + b.max_coord(2)) / 2.0) ?
			(h - b.min_coord(2))*(h - b.min_coord(2)) : (b.max_coord(2) - h)*(b.max_coord(2) - h);
		return d0 + d1 + d2;
	}
	float max_distance_to_rectangle(const Point& p,
		const CGAL::Kd_tree_rectangle<FT, D>& b, std::vector<float>& dists) {
		float h = p.x();
		dists[0] = (h >= (b.min_coord(0) + b.max_coord(0)) / 2.0) ?
			(h - b.min_coord(0)) : (b.max_coord(0) - h);

		h = p.y();
		dists[1] = (h >= (b.min_coord(1) + b.max_coord(1)) / 2.0) ?
			(h - b.min_coord(1)) : (b.max_coord(1) - h);
		h = p.z();
		dists[2] = (h >= (b.min_coord(2) + b.max_coord(2)) / 2.0) ?
			(h - b.min_coord(2)) : (b.max_coord(2) - h);
		return dists[0] * dists[0] + dists[1] * dists[1] + dists[2] * dists[2];
	}
	float new_distance(float& dist, float old_off, float new_off,
		int /* cutting_dimension */)  const {
		return dist + new_off*new_off - old_off*old_off;
	}
	float transformed_distance(float d) const { return d*d; }
	float inverse_of_transformed_distance(float d) { return std::sqrt(d); }
}; // end of struct Distance