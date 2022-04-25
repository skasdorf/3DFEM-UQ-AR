#pragma once
class adjacent {
public:
	int u_dir = 0; //0 no refined neighbors, 1 refined neighbor on -u face, 2 refined neighbor on +u face, 3 refined neighbors on both -u and +u
	int v_dir = 0; //0,1,2,3
	int w_dir = 0; //0,1,2,3

	//additional variables for multiconnections
};