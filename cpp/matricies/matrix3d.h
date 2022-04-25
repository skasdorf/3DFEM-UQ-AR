
#ifndef CPP_MATRIX3D_H
#define CPP_MATRIX3D_H

#include <vector>


template <class T>
class matrix3d {

public:
    std::vector<std::vector<std::vector<T>>> matrix;
    matrix3d<T>(const int& a, const int& b, const int& c){
        matrix = std::vector<std::vector<std::vector<T>>>
                (a+1, std::vector<std::vector<T>>(b+1, std::vector<T>(c+1)));
    }

    T& operator()(const int& a, const int& b, const int& c){
        return matrix3d<T>::matrix[a][b][c];
    }
	const T& operator()(const int& a, const int& b, const int& c) const {
		return matrix3d<T>::matrix[a][b][c];
	}
	matrix3d<T>() {
		matrix = std::vector<std::vector<std::vector<T>>>
			(0, std::vector<std::vector<T>>(0, std::vector<T>(0)));
	}

};

#endif //CPP_MATRIX3D_H
