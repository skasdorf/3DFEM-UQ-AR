#pragma once
class adBasis {
public:
	int index; //element index - starts at 1  (e->index without adjustment)
	int iuvw; //index of the direction (1,2,3)
	int i, j, k; //indices of the order of the basis
	adBasis(int index, int iuvw, int i, int j, int k) {
		this->index = index;
		this->iuvw = iuvw;
		this->i = i;
		this->j = j;
		this->k = k;
	}
};