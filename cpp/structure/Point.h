#ifndef POINT_H

#define POINT_H// POINT_H
template <class T>
std::vector<T> operator*(const std::vector<T>& vec, const T& scalar) {
	//assumes starts at 1
	return{0, vec[1] * scalar, vec[2] * scalar, vec[3] * scalar };
}
template <class T>
std::vector<T> operator+(const std::vector<T>& vec1, const std::vector<T>& vec2) { 
	//assumes that the vectors both start at 1
	return{0, vec1[1] + vec2[1], vec1[2] + vec2[2], vec1[3] + vec2[3] };
}



class Point {
public:
	double x;
	double y;
	double z;
	Point() {
	};
	Point(const double& x, const double& y, const double& z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	std::vector<double> toVector() {
		//converts to start at 1
		/*std::vector<double> vec;
		vec.push_back(0); vec.push_back(x); vec.push_back(y); vec.push_back(z);*/
		return {0, x, y, z};
	}
};

inline Point operator*(const Point& p1, const double scalar) {
	return Point(p1.x*scalar, p1.y*scalar, p1.z*scalar);
}
inline Point operator/(const Point& p1, const double scalar) {
	return Point(p1.x / scalar, p1.y / scalar, p1.z / scalar);
}
inline Point operator+(Point& p1, const Point p2) {
	return Point(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
}
#endif