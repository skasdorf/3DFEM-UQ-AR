#pragma once
//for saving cBrPAK, cArPAK
std::ofstream file;
file.open("../exampleFiles/debug_data/ " + mesh_name + " /results/cPAKs.txt")
for (auto e = this->elements.begin(); e != this->elements.end(); ++e) {
	//for every element
	for (int i = 0; i < e->cBrPAK.size(); ++i) {
		for (int j = 0; j < e->cBrPAK.size(); ++j) {
			file << i << " " << j << " " << std::setprecision(15) << e->cArPAK(i, j) << " " << e->cBrPAK(i, j) << " ";
		}
	}
	file << std::endl;
}