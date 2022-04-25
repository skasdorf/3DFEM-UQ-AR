#include <string>
#include <sstream>
#include "Facet.h"
#include <vector>
#include <fstream>

Facet::Facet() {

}
std::vector<std::vector<int>> Facet::MAKEEPT() {
	//result starts at 1
	std::vector<std::vector<int>> ept(6, std::vector<int>(4));
	ept[1][1] = this->designations1.first;
	ept[1][2] = this->designations1.second;
	ept[2][1] = this->designations2.first;
	ept[2][2] = this->designations2.second;
	ept[1][3] = this->local_facet_indices.first;
	ept[2][3] = this->local_facet_indices.second;
	for (int i = 1; i < 3; ++i) {
		switch (abs(ept[1][i])) {
		case 1:
			switch (abs(ept[2][i])) {
			case 1:
				ept[3][i] = 0;
				break;
			case 2:
				if (ept[2][3] == 5 || ept[2][3] == 6) ept[3][i] = 1;
				else ept[3][i] = 0;
				break;
			case 3:
				ept[3][i] = 1;
				break;
			}
			ept[4][i] = 1;
			ept[5][i] = 2;
			break;
		case 2:
			switch (abs(ept[2][i])) {
			case 1:
				if (ept[1][3] == 5 || ept[1][3] == 6) {
					ept[3][i] = 1;
					ept[4][i] = 2;
					ept[5][i] = 1;
					
				}
				else {
					ept[3][i] = 0;
					ept[4][i] = 1;
					ept[5][i] = 2;
				}
				break;
			case 2:
				if (ept[1][3] == 5 || ept[1][3] == 6) {
					ept[4][i] = 2;
					ept[5][i] = 1;
					if (ept[2][3] == 1 || ept[2][3] == 2) ept[3][i]=1;
					else ept[3][i] = 0;
				}
				else {
					ept[4][i] = 1;
					ept[5][i] = 2;
					if (ept[2][3] == 5 || ept[2][3] == 6) ept[3][i] = 1;
					else ept[3][i] = 0;
				}
				break;
			case 3:
				if (ept[1][3] == 1 || ept[1][3] == 2) {
					ept[3][i] = 1;
					ept[4][i] = 1;
					ept[5][i] = 2;
				}
				else {
					ept[3][i] = 0;
					ept[4][i] = 2;
					ept[5][i] = 1;
				}
				break;
			}
			break;
		case 3:
			switch (abs(ept[2][i])) {
			case 1:
				ept[3][i] = 1;
				break;
			case 2:
				if (ept[2][3] == 1 || ept[2][3] == 2) ept[3][i] = 1;
				else ept[3][i] = 0;
				break;
			case 3:
				ept[3][i] = 0;
				break;
			}
			ept[4][i] = 2;
			ept[5][i] = 1;
			break;
		}
	}
	return ept;
}



