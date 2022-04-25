#ifndef TABLE_H
#define TABLE_H
//tables
//includes vertex and face tables


//the possibilities for v1 and v2 may be + 1, -1, +2, -2, +3, -3, 
//which mean that the unit vectors correspond to the + u, -u, +v, -v, +w, -w directions, respectively
#include <vector>
class Tables {
public:
	int vertex_table[8][4];
	
	std::vector<std::vector<std::vector<int>>> face_table= std::vector<std::vector<std::vector<int>>>(6);
	std::vector < std::vector<int>> face_table_rel_ind = std::vector<std::vector<int>>(8);
	Tables() {
		//vertex_table initialization
		vertex_table[0][0] = 1;
		vertex_table[0][1] = 1;
		vertex_table[0][2] = 1;
		vertex_table[0][3] = 1;

		vertex_table[1][0] = 2;
		vertex_table[1][1] = 3;
		vertex_table[1][2] = 4;
		vertex_table[1][3] = 5;

		vertex_table[2][0] = 3;
		vertex_table[2][1] = 7;
		vertex_table[2][2] = 13;
		vertex_table[2][3] = 21;

		vertex_table[3][0] = 4;
		vertex_table[3][1] = 9;
		vertex_table[3][2] = 16;
		vertex_table[3][3] = 25;

		vertex_table[4][0] = 5;
		vertex_table[4][1] = 19;
		vertex_table[4][2] = 49;
		vertex_table[4][3] = 101;

		vertex_table[5][0] = 6;
		vertex_table[5][1] = 21;
		vertex_table[5][2] = 52;
		vertex_table[5][3] = 105;

		vertex_table[6][0] = 7;
		vertex_table[6][1] = 25;
		vertex_table[6][2] = 61;
		vertex_table[6][3] = 121;

		vertex_table[7][0] = 8;
		vertex_table[7][1] = 27;
		vertex_table[7][2] = 64;
		vertex_table[7][3] = 125;
		face_table[0] = { { 1, 3, 5, 7, 2, 3 },{ 1, 5, 3, 7, 3, 2 },{ 3, 1, 7, 5,-2, 3 },{ 3, 7, 1, 5, 3,-2 },
		{ 5, 1, 7, 3,-3, 2 },{ 7, 3, 5, 1,-3,-2 },{ 7, 5, 3, 1,-2,-3 } };
		face_table[1] = { { 2, 4, 6, 8, 2, 3 }, { 2, 6, 4, 8, 3, 2 }, { 4, 2, 8, 6,-2, 3 }, 
		{ 4, 8, 2, 6, 3,-2 }, { 6, 2, 8, 4,-3, 2 }, { 6, 8, 2, 4, 2,-3 }, { 8, 4, 6, 2,-3,-2 }, { 8, 6, 4, 2,-2,-3 } };
		face_table[2] = { { 1, 2, 5, 6, 1, 3 }, { 1, 5, 2, 6, 3, 1 }, { 2, 1, 6, 5,-1, 3 },
		{ 2, 6, 1, 5, 3,-1 },{ 5, 1, 6, 2,-3, 1 },{ 5, 6, 1, 2, 1,-3 },{ 6, 2, 5, 1,-3,-1 },{ 6, 5, 2, 1,-1,-3 } };
		face_table[3] = { { 3, 4, 7, 8, 1, 3 } ,{ 3, 7, 4, 8, 3, 1 },{ 4, 3, 8, 7,-1, 3 },
		{ 4, 8, 3, 7, 3,-1 },{ 7, 3, 8, 4,-3, 1 },{ 7, 8, 3, 4, 1,-3 },{ 8, 4, 6, 2,-3,-2 },{ 8, 7, 4, 3,-1,-3 } };
		face_table[4] = { { 1, 2, 3, 4, 1, 2 } ,{ 1, 3, 2, 4, 2, 1 },{ 2, 1, 4, 3,-1, 2 },
		{ 2, 4, 1, 3, 2,-1 },{ 3, 1, 4, 2,-2, 1 },{ 3, 4, 1, 2, 1,-2 },{ 4, 2, 3, 1,-2,-1 },{ 4, 3, 2, 1,-1,-2 } };
		face_table[5] = { { 5, 6, 7, 8, 1, 2 } ,{ 5, 7, 6, 8, 2, 1 } ,{ 6, 5, 8, 7,-1, 2 },
		{ 6, 8, 5, 7, 2,-1 } ,{ 7, 5, 8, 6,-2, 1 } ,{ 7, 8, 5, 6, 1,-2 } ,{ 8, 6, 7, 5,-2,-1 },{ 8, 7, 6, 5,-1,-2 } };

		face_table_rel_ind = { {0,1,2,3,},{0,2,1,3},{1,0,3,2},{1,3,0,2},{2,0,3,1},{2,3,0,1},{3,1,2,0},{3,2,1,0} };


		

	}

};
#endif