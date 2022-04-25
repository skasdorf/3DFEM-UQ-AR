
#ifndef CPP_FOURDMATRIX_H
#define CPP_FOURDMATRIX_H

#include <vector>

template <class T>
class matrix4d {

public:

    std::vector<std::vector<std::vector<std::vector<T>>>> matrix;
    matrix4d<T>(const int& a, const int& b, const int& c, const int& d){
        matrix = std::vector<std::vector<std::vector<std::vector<T>>>>
                (a+1, std::vector<std::vector<std::vector<T>>>(b+1, std::vector<std::vector<T>>(c+1, std::vector<T>(d+1))));
    }

    T& operator()(const int& a, const int& b, const int& c, const int& d){
        return matrix4d<T>::matrix[a][b][c][d];
    }
	const T& operator()(const int& a, const int& b, const int& c, const int& d) const {
		return matrix4d<T>::matrix[a][b][c][d];
	}

    std::vector<T>& operator()(const int& a, const int& b, const int& c){
        return matrix4d<T>::matrix[a][b][c];
    }
	matrix4d<T>() {
		matrix = std::vector<std::vector<std::vector<std::vector<T>>>>
			(0, std::vector<std::vector<std::vector<T>>>(0, std::vector<std::vector<T>>(0, std::vector<T>(0))));
	}

};


#endif //CPP_FOURDMATRIX_H
