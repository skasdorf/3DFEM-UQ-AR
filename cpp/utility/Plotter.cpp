#include "Plotter.h"
void PlotController::plotterIO(std::string plot_file, std::string coeff_file, bool use_file) {
	std::ifstream file;
	std::string line;
	
	
	//int size = dom->elements.size()*6;

	//this->element_num_for_plot.resize(size);
	//this->u1.resize(size);
	//this->u2.resize(size);
	//this->v1.resize(size);
	//this->v2.resize(size);
	//this->w1.resize(size);
	//this->w2.resize(size);
	//this->Nu.resize(size);
	//this->Nv.resize(size);
	//this->Nw.resize(size);
	this->element_num_for_plot.push_back(0);
	this->u1.push_back(0);
	this->u2.push_back(0);
	this->Nu.push_back(0);
	this->v1.push_back(0);
	this->v2.push_back(0);
	this->Nv.push_back(0);
	this->w1.push_back(0);
	this->w2.push_back(0);
	this->Nw.push_back(0);
	if (use_file) {
		file.open(plot_file);
		getline(file, line); //contains the number of plots to perform
		this->numberOfPlots = std::stoi(functions::split(line, ' ')[0]);
		while (getline(file, line)) {
			auto input_line = functions::split(line, ' ');
			this->element_num_for_plot.push_back(std::stoi(input_line[1]));
			this->u1.push_back(std::stoi(input_line[2]));
			this->u2.push_back(std::stoi(input_line[3]));
			this->Nu.push_back(std::stoi(input_line[4]));
			this->v1.push_back(std::stoi(input_line[5]));
			this->v2.push_back(std::stoi(input_line[6]));
			this->Nv.push_back(std::stoi(input_line[7]));
			this->w1.push_back(std::stoi(input_line[8]));
			this->w2.push_back(std::stoi(input_line[9]));
			this->Nw.push_back(std::stoi(input_line[10]));
		}
		file.close();
	}
	else {
		this->numberOfPlots = dom->elements.size()*6;
		for (int e = 0; e < dom->elements.size(); ++e) {
			for (int i = 1; i <= 6; ++i) {
				this->element_num_for_plot.push_back(e + 1);
				switch (i) {
				case 1:
					this->u1.push_back(-1);
					this->u2.push_back(1);
					this->Nu.push_back(10);
					this->v1.push_back(-1);
					this->v2.push_back(1);
					this->Nv.push_back(10);
					this->w1.push_back(-1);
					this->w2.push_back(-1);
					this->Nw.push_back(1);
					break;

				case 2:
					this->u1.push_back(-1);
					this->u2.push_back(1);
					this->Nu.push_back(10);
					this->v1.push_back(-1);
					this->v2.push_back(1);
					this->Nv.push_back(10);
					this->w1.push_back(1);
					this->w2.push_back(1);
					this->Nw.push_back(1);
					break;
				case 3:
					this->u1.push_back(-1);
					this->u2.push_back(-1);
					this->Nu.push_back(1);
					this->v1.push_back(-1);
					this->v2.push_back(1);
					this->Nv.push_back(10);
					this->w1.push_back(-1);
					this->w2.push_back(1);
					this->Nw.push_back(10);
					break;
				case 4:
					this->u1.push_back(1);
					this->u2.push_back(1);
					this->Nu.push_back(1);
					this->v1.push_back(-1);
					this->v2.push_back(1);
					this->Nv.push_back(10);
					this->w1.push_back(-1);
					this->w2.push_back(1);
					this->Nw.push_back(10);
					break;
				case 5:
					this->u1.push_back(-1);
					this->u2.push_back(1);
					this->Nu.push_back(10);
					this->v1.push_back(-1);
					this->v2.push_back(-1);
					this->Nv.push_back(1);
					this->w1.push_back(-1);
					this->w2.push_back(1);
					this->Nw.push_back(10);
					break;
				case 6:
					this->u1.push_back(-1);
					this->u2.push_back(1);
					this->Nu.push_back(10);
					this->v1.push_back(1);
					this->v2.push_back(1);
					this->Nv.push_back(1);
					this->w1.push_back(-1);
					this->w2.push_back(1);
					this->Nw.push_back(10);
					break;
				}
			}
		}
	}
	file.open(coeff_file);
	std::vector<std::complex<double>> cAlpha_list;
	while (getline(file, line)) {
		auto input_line = functions::split(line, ' ');
		cAlpha_list.push_back(std::complex<double>(std::stod(input_line[0]), std::stod(input_line[1])));
	}
	file.close();
	dom->cAlpha = cAlpha_list;
}
void PlotController::FUODUVWPLOT(const int& m, const int& n, const int& l, const int& i, const int& j,
	const int& k, const matrix2d<double>& uPowers, const matrix2d<double>& fvPowers, const matrix2d<double>& fwPowers,
	const int& nPu, const int& nPv, const int& nPw, double& fu) {
	fu = uPowers(m, i)*fvPowers(n, j)*fwPowers(l, k);
}
void PlotController::FVODUVWPLOT(const int& m, const int& n, const int& l, const int& i, const int& j,
	const int& k, const matrix2d<double>& fuPowers, const matrix2d<double>& vPowers, const matrix2d<double>& fwPowers,
	const int& nPu, const int& nPv, const int& nPw, double& fv) {
	fv = fuPowers(m, i)*vPowers(n, j)*fwPowers(l, k);
}
void PlotController::FWODUVWPLOT(const int& m, const int& n, const int& l, const int& i, const int& j,
	const int& k, const matrix2d<double>& fuPowers, const matrix2d<double>& fvPowers, const matrix2d<double>& wPowers,
	const int& nPu, const int& nPv, const int& nPw, double& fw) {
	fw = fuPowers(m, i)*fvPowers(n, j)*wPowers(l, k);
}
void PlotController::find_ausec(std::vector<double>& av, std::vector<double>& aw, std::vector<double>& ausec) {
	products::cross(av, aw, ausec);
}
void PlotController::find_avsec(std::vector<double>& aw, std::vector<double>& au, std::vector<double>& avsec) {
	products::cross(aw, au, avsec);
}
void PlotController::find_awsec(std::vector<double>& au, std::vector<double>& av, std::vector<double>& awsec) {
	products::cross(au, av, awsec);
}
void PlotController::findPowers(const int& nu, const int& nv, const int& nw, const int& nglu, const int& nglv, const int& nglw, const std::vector<double>& xglu, const std::vector<double>& xglv, const std::vector<double>& xglw,
	matrix2d<double>& uPowers, matrix2d<double>& vPowers, matrix2d<double>& wPowers, matrix2d<double>& fuPowers, matrix2d<double>& fvPowers, matrix2d<double>& fwPowers,
	matrix2d<double>& fpuPowers, matrix2d<double>& fpvPowers, matrix2d<double>& fpwPowers) {
	//upowers goes from (0 to nglu+1, 0 to nu)
	//vpower (0 to nglv + 1, 0 to nv)
	//wPowers (0 to nglw + 1, 0 to nw)
	//fuPowers, fpuPowers same as uPowers
	//fvPowers, fpvPowers same as vPowers
	double x;
	for (int i = 0; i <= nglu + 1; ++i) {
		x = xglu[i];
		uPowers(i, 0) = 1.0;
		uPowers(i, 1) = x;
		fuPowers(i, 0) = 1.0 - x;
		fuPowers(i, 1) = 1.0 + x;
		fpuPowers(i, 0) = -1.0;
		fpuPowers(i, 1) = 1.0;
		for (int j = 2; j <= nu; ++j) {
			uPowers(i, j) = uPowers(i, j - 1)*x;
			if (j % 2 == 0) {
				fuPowers(i, j) = uPowers(i, j) - 1.0;
				fpuPowers(i, j) = j*uPowers(i, j - 1);
			}
			else {
				fuPowers(i, j) = uPowers(i, j) - x;
				fpuPowers(i, j) = j*  uPowers(i, j - 1) - 1.0;
			}
		}
	}
	for (int i = 0; i <= nglv + 1; ++i) {
		x = xglv[i];
		vPowers(i, 0) = 1.0;
		vPowers(i, 1) = x;
		fvPowers(i, 0) = 1.0 - x;
		fvPowers(i, 1) = 1.0 + x;
		fpvPowers(i, 0) = -1.0;
		fpvPowers(i, 1) = 1.0;
		for (int j = 2; j <= nv; ++j) {
			vPowers(i, j) = vPowers(i, j - 1)*x;
			if (j % 2 == 0) {
				fvPowers(i, j) = vPowers(i, j) - 1.0;
				fpvPowers(i, j) = j* vPowers(i, j - 1);
			}
			else {
				fvPowers(i, j) = vPowers(i, j) - x;
				fpvPowers(i, j) = j*vPowers(i, j - 1) - 1.0;
			}
		}
	}
	for (int i = 0; i <= nglw + 1; ++i) {
		x = xglw[i];
		wPowers(i, 0) = 1.0;
		wPowers(i, 1) = x;
		fwPowers(i, 0) = 1.0 - x;
		fwPowers(i, 1) = 1.0 + x;
		fpwPowers(i, 0) = -1.0;
		fpwPowers(i, 1) = 1.0;
		for (int j = 2; j <= nw; ++j) {
			wPowers(i, j) = wPowers(i, j - 1)*x;
			if (j % 2 == 0) {
				fwPowers(i, j) = wPowers(i, j) - 1.0;
				fpwPowers(i, j) = j*wPowers(i, j - 1);
			}
			else {
				fwPowers(i, j) = wPowers(i, j) - x;
				fpwPowers(i, j) = j*wPowers(i, j - 1) - 1.0;
			}
		}
	}

}
void PlotController::unitaryvectorsplot(int& nPu, int& nPv, int& nPw, std::vector<double>& xglu, std::vector<double>& xglv, std::vector<double>& xglw, int& elementOrders,
	matrix2d<double>& rs, matrix4d<double>& rMatrix, matrix4d<double>& auMatrix, matrix4d<double>& avMatrix, matrix4d<double>& awMatrix, matrix3d<double>& jacobianMatrix, const int& nRs) {
	
	std::vector<double> r(4), au(4), av(4), aw(4);
	for (int m = 0; m <= nPu + 1; ++m) {
		for (int n = 0; n <= nPv + 1; ++n) {
			for (int l = 0; l <= nPw + 1; ++l) {
				r = unitVectorsM::findR(elementOrders, xglu[m], xglv[n], xglw[l], rs, nRs);
				rMatrix(m, n, l, 1) = r[1];
				rMatrix(m, n, l, 2) = r[2];
				rMatrix(m, n, l, 3) = r[3];
			}
		}
	}
	for (int m = 0; m <= nPu + 1; ++m) {
		for (int n = 0; n <= nPv + 1; ++n) {
			for (int l = 0; l <= nPw + 1; ++l) {
				au = unitVectorsM::findAU(elementOrders, xglu[m], xglv[n], xglw[l], rs, nRs);
				auMatrix(m, n, l, 1) = au[1];
				auMatrix(m, n, l, 2) = au[2];
				auMatrix(m, n, l, 3) = au[3];
			}
		}
	}
	for (int m = 0; m <= nPu + 1; ++m) {
		for (int n = 0; n <= nPv + 1; ++n) {
			for (int l = 0; l <= nPw + 1; ++l) {
				av = unitVectorsM::findAV(elementOrders, xglu[m], xglv[n], xglw[l], rs, nRs);
				avMatrix(m, n, l, 1) = av[1];
				avMatrix(m, n, l, 2) = av[2];
				avMatrix(m, n, l, 3) = av[3];
			}
		}
	}
	for (int m = 0; m <= nPu + 1; ++m) {
		for (int n = 0; n <= nPv + 1; ++n) {
			for (int l = 0; l <= nPw + 1; ++l) {
				aw = unitVectorsM::findAW(elementOrders, xglu[m], xglv[n], xglw[l], rs, nRs);
				awMatrix(m, n, l, 1) = aw[1];
				awMatrix(m, n, l, 2) = aw[2];
				awMatrix(m, n, l, 3) = aw[3];
			}
		}
	}
	for (int m = 0; m <= nPu + 1; ++m) {
		for (int n = 0; n <= nPv + 1; ++n) {
			for (int l = 0; l <= nPw + 1; ++l) {
				for (int coord = 1; coord <= 3; ++coord) {
					au[coord] = auMatrix(m, n, l, coord);
					av[coord] = avMatrix(m, n, l, coord);
					aw[coord] = awMatrix(m, n, l, coord);
				}
				unitVectorsM::find_jacobian(au, av, aw, jacobianMatrix(m, n, l));
				//jacobianMatrix(m,n,l)
			}
		}
	}

}
void edpointsinterval(double& xbeg, double& xend, int& n, std::vector<double>& x) {
	double h;
	if (n == 1) x[1] = xbeg;
	else if (n > 1) {
		h = (xend - xbeg) / (n - 1);
		for (int i = 0; i < n; ++i) {
			x[i+1] = xbeg + i*h;
		}
	}
}

void PlotController::plotField(std::string mesh_name, bool plot_basis) {
	plot_basis = true;
	BasisEval eval;
	if (eval.basisType == 1) { //legendre basis
		eval.setup_legendre();
	}
	std::ofstream file;
	file.open("../exampleFiles/debug_data/" + mesh_name + "/results/Fieldplot.dat");
	double hu, hv, hw, xu, xv, xw;
	double fu, fv, fw;
	std::vector<double> au(4), av(4), aw(4), ausec(4), avsec(4), awsec(4), fVec(4);
	
	for (int iPlotCount = 1; iPlotCount <= numberOfPlots; ++iPlotCount) {
		if (element_num_for_plot[iPlotCount] == 0) continue;
		matrix4d<double> rMatrix(Nu[iPlotCount] + 1, Nv[iPlotCount] + 1, Nw[iPlotCount] + 2, 3),
			auMatrix(Nu[iPlotCount] + 1, Nv[iPlotCount] + 1, Nw[iPlotCount] + 1, 3), avMatrix(Nu[iPlotCount] + 1,
				Nv[iPlotCount] + 1, Nw[iPlotCount] + 1, 3), awMatrix(Nu[iPlotCount] + 1, Nv[iPlotCount] + 1, Nw[iPlotCount] + 1, 3);
		matrix3d<double> jacobianMatrix(Nu[iPlotCount] + 1, Nv[iPlotCount] + 1, Nw[iPlotCount] + 1);

		matrix4d<std::complex<double>> cField(Nu[iPlotCount] + 1, Nv[iPlotCount] + 1, Nw[iPlotCount] + 1, 3);
		int ehOld = 0;
		
		for (auto e = dom->elements.begin(); e != dom->elements.end(); ++e) {
			if (e->index != element_num_for_plot[iPlotCount]) continue;
			
			int nu = dom->elements[e->index - 1].expansion[0];
			int nv = dom->elements[e->index - 1].expansion[1];
			int nw = dom->elements[e->index - 1].expansion[2];
			auto rs = dom->elements[e->index - 1].rs;
			std::vector<double> xglu(Nu[iPlotCount] + 2), xglv(Nv[iPlotCount] + 2), xglw(Nw[iPlotCount] + 2);
			edpointsinterval(u1[iPlotCount], u2[iPlotCount], Nu[iPlotCount], xglu);
			xglu[0] = -1.0;
			xglu[Nu[iPlotCount] + 1] = 1;
			edpointsinterval(v1[iPlotCount], v2[iPlotCount], Nv[iPlotCount], xglv);
			xglv[0] = -1.0;
			xglv[Nv[iPlotCount] + 1] = 1;
			edpointsinterval(w1[iPlotCount], w2[iPlotCount], Nw[iPlotCount], xglw);
			xglw[0] = -1.0;
			xglw[Nw[iPlotCount] + 1] = 1;
			unitaryvectorsplot(Nu[iPlotCount], Nv[iPlotCount], Nw[iPlotCount], xglu, xglv, xglw,
				dom->elements[e->index - 1].geom_order, rs, rMatrix, auMatrix, avMatrix, awMatrix, jacobianMatrix, dom->elements[e->index - 1].nRs);
			matrix2d<double> uPowers(Nu[iPlotCount] + 1, nu), fuPowers(Nu[iPlotCount] + 1, nu), fpuPowers(Nu[iPlotCount] + 1, nu);
			matrix2d<double> vPowers(Nv[iPlotCount] + 1, nv), fvPowers(Nv[iPlotCount] + 1, nv), fpvPowers(Nv[iPlotCount] + 1, nv);
			matrix2d<double> wPowers(Nw[iPlotCount] + 1, nw), fwPowers(Nw[iPlotCount] + 1, nw), fpwPowers(Nw[iPlotCount] + 1, nw);

			int indicatorBasisType = 1;
			if (indicatorBasisType == 1) {
				eval.finduPowersLeg(nu, Nu[iPlotCount], xglu, uPowers, fuPowers, fpuPowers);
				eval.findvPowersLeg(nv, Nv[iPlotCount], xglv, vPowers, fvPowers, fpvPowers);
				eval.findwPowersLeg(nw, Nw[iPlotCount], xglw, wPowers, fwPowers, fpwPowers);
				/*findPowers(nu, nv, nw, Nu[iPlotCount], Nv[iPlotCount], Nw[iPlotCount], xglu, xglv, xglw,
					uPowers, vPowers, wPowers, fuPowers, fvPowers, fwPowers, fpuPowers, fpvPowers, fpwPowers);*/
				/*uPowers = eval.get_u_samples(xglu.size(), nu)[0];
				vPowers = eval.get_v_samples(xglu.size(), nv)[0];
				wPowers = eval.get_w_samples(xglu.size(), nw)[0];
				fuPowers = eval.get_u_samples(xglu.size(), nu)[1];
				fvPowers = eval.get_v_samples(xglu.size(), nv)[1];
				fwPowers = eval.get_w_samples(xglu.size(), nw)[1];
				fpuPowers = eval.get_u_samples(xglu.size(), nu)[2];
				fpvPowers = eval.get_v_samples(xglu.size(), nv)[2];
				fpwPowers = eval.get_w_samples(xglu.size(), nw)[2];*/
				
			}
			for (int kMat = dom->elements[e->index - 1].unknownsStart; kMat <= dom->elements[e->index - 1].unknownsEnd; ++kMat) {
				int ihCon; 
				plot_basis ? ihCon = kMat : ihCon = abs(dom->vectorD[kMat]);
				if (ihCon == 0) continue;
				int eh = dom->eiUVWijk[kMat][1];
				int iuvwh = dom->eiUVWijk[kMat][2];
				int ih = dom->eiUVWijk[kMat][3];
				int jh = dom->eiUVWijk[kMat][4];
				int kh = dom->eiUVWijk[kMat][5];
				int ehPlot = e->index;
				fVec[1] = 0; fVec[2] = 0; fVec[3] = 0;
				for (int m = 0; m <= Nu[iPlotCount] + 1; ++m) {
					for (int n = 0; n <= Nv[iPlotCount] + 1; ++n) {
						for (int l = 0; l <= Nw[iPlotCount] + 1; ++l) {
							au = auMatrix(m, n, l);
							av = avMatrix(m, n, l);
							aw = awMatrix(m, n, l);
							double jacobian = jacobianMatrix(m, n, l);
							switch (iuvwh) {
							case 1:
								FUODUVWPLOT(m, n, l, ih, jh, kh, uPowers, fvPowers, fwPowers, Nu[iPlotCount], Nv[iPlotCount], Nw[iPlotCount], fu);
								find_ausec(av, aw, ausec);
								fVec[1] = ausec[1] * fu / jacobian;
								fVec[2] = ausec[2] * fu / jacobian;
								fVec[3] = ausec[3] * fu / jacobian;
								break;
							case 2:
								FVODUVWPLOT(m, n, l, ih, jh, kh, fuPowers, vPowers, fwPowers, Nu[iPlotCount], Nv[iPlotCount], Nw[iPlotCount], fv);
								find_avsec(aw, au, avsec);
								fVec[1] = avsec[1] * fv / jacobian;
								fVec[2] = avsec[2] * fv / jacobian;
								fVec[3] = avsec[3] * fv / jacobian;
								break;
							case 3:
								FWODUVWPLOT(m, n, l, ih, jh, kh, fuPowers, fvPowers, wPowers, Nu[iPlotCount], Nv[iPlotCount], Nw[iPlotCount], fw);
								find_awsec(au, av, awsec);
								fVec[1] = awsec[1] * fw / jacobian;
								fVec[2] = awsec[2] * fw / jacobian;
								fVec[3] = awsec[3] * fw / jacobian;
								break;
							}
							if (dom->vectorD[kMat] > 0 || plot_basis) {
								cField(m, n, l, 1) += dom->cAlpha[ihCon] * fVec[1];
								cField(m, n, l, 2) += dom->cAlpha[ihCon] * fVec[2];
								cField(m, n, l, 3) += dom->cAlpha[ihCon] * fVec[3];
							}
							else {
								cField(m, n, l, 1) += -dom->cAlpha[ihCon] * fVec[1];
								cField(m, n, l, 2) += -dom->cAlpha[ihCon] * fVec[2];
								cField(m, n, l, 3) += -dom->cAlpha[ihCon] * fVec[3];
							}
						}
					}
				}
			}
		}
		for (int m = 1; m <= Nu[iPlotCount]; ++m) {
			if (Nu[iPlotCount] == 1) {
				xu = u1[iPlotCount];
			}
			else if (Nu[iPlotCount] > 1) {
				hu = (u2[iPlotCount] - u1[iPlotCount]) / (Nu[iPlotCount] - 1);
				xu = u1[iPlotCount] + (m - 1)*hu;
			}
			for (int n = 1; n <= Nv[iPlotCount]; ++n) {
				if (Nv[iPlotCount] == 1) xv = v1[iPlotCount];
				else if (Nv[iPlotCount] > 1) {
					hv = (v2[iPlotCount] - v1[iPlotCount]) / (Nv[iPlotCount] - 1);
					xv = v1[iPlotCount] + (n - 1)*hv;
				}
				for (int l = 1; l <= Nw[iPlotCount]; ++l) {
					if (Nw[iPlotCount] == 1) xw = w1[iPlotCount];
					else if (Nw[iPlotCount] > 1) {
						hw = (w2[iPlotCount] - w1[iPlotCount]) / (Nw[iPlotCount]);
						xw = w1[iPlotCount] + (l - 1)*hw;
					}
					//write cField to file
					file << rMatrix(m, n, l, 1) << " " << rMatrix(m, n, l, 2) << " " << rMatrix(m, n, l, 3) << " "
						<< cField(m, n, l, 1).real() << " " << cField(m, n, l, 1).imag() << " "
						<< cField(m, n, l, 2).real() << " " << cField(m, n, l, 2).imag() << " "
						<< cField(m, n, l, 3).real() << " " << cField(m, n, l, 3).imag() << " "
						<< abs(cField(m, n, l, 1)) << " " << abs(cField(m, n, l, 2)) << " " << abs(cField(m, n, l, 3)) << std::endl;
				}
			}
		}
	}
}