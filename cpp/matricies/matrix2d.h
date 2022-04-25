//
// Created by Blake Troksa on 5/30/18.
//

#ifndef CPP_MATRIX2D_H
#define CPP_MATRIX2D_H

#include <vector>

template <class T>
class matrix2d {

public:
    std::vector<std::vector<T>> matrix;
    matrix2d<T>(int a, int b){
        matrix = std::vector<std::vector<T>>(a+1, std::vector<T>(b+1));
    }

    T& operator()(const int& a, const int& b){
        return matrix2d<T>::matrix[a][b];
    }
	const T& operator()(const int& a, const int& b) const {
		return matrix2d<T>::matrix[a][b];
	}
	std::vector<T>& operator()(int& a) {
		return matrix2d<T>::matrix[a];
	}
	std::vector<T>& operator[](int a) { //for use when index accessing is known at compile time
		return matrix2d<T>::matrix[a];
	}
	matrix2d<T>() {
		matrix = std::vector<std::vector<T>>(0, std::vector<T>(0));
	}
	int size(int dim) {
		if (dim == 1) {
			return matrix.size();
		}
		else {
			return matrix[0].size();
		}
	}
};

#endif //CPP_MATRIX2D_H