#include "BasisEval.h"


evaluation & BasisEval::get_u_samples(int nglu, int nu)
{
	std::string key = std::to_string(nglu) + "|" + std::to_string(nu); //create string key

	auto it = this->uMap.find(key); //try to find the evaluation in the map for the given key

	if (it == this->uMap.end()) { //the evaluation for the given key does not exist
		matrix2d<double> uPowers = matrix2d<double>(nglu + 1, nu);
		matrix2d<double> fuPowers = matrix2d<double>(nglu + 1, nu);
		matrix2d<double> fpuPowers = matrix2d<double>(nglu + 1, nu);
		std::vector<double> xglu, wglu;
		functions::gaussk(nglu, xglu, wglu);
		if (basisType == 1) {
			finduPowersLeg(nu, nglu, xglu, uPowers, fuPowers, fpuPowers);
		}
		else {
			finduPowers(nu, nglu, xglu, uPowers, fuPowers, fpuPowers);
		}
		
		uMap.insert({ key, { uPowers,fuPowers,fpuPowers } });
		return uMap.at(key);
	}
	else //the evaluation has been done before
		return it->second; //return its reference
}



//below this is placeholder for now (code just in place to compile)
evaluation & BasisEval::get_v_samples(int nglv, int nv)
{
	std::string key = std::to_string(nglv) + "|" + std::to_string(nv); //create string key

	auto it = this->vMap.find(key); //try to find the evaluation in the map for the given key

	if (it == this->vMap.end()) { //the evaluation for the given key does not exist
		matrix2d<double> vPowers = matrix2d<double>(nglv + 1, nv);
		matrix2d<double> fvPowers = matrix2d<double>(nglv + 1, nv);
		matrix2d<double> fpvPowers = matrix2d<double>(nglv + 1, nv);
		std::vector<double> xglv, wglv;
		functions::gaussk(nglv, xglv, wglv);
		if (basisType == 1) {
			findvPowersLeg(nv, nglv, xglv, vPowers, fvPowers, fpvPowers);
		}
		else {
			findvPowers(nv, nglv, xglv, vPowers, fvPowers, fpvPowers);
		}
		vMap.insert({ key,{ vPowers,fvPowers,fpvPowers } });
		return vMap.at(key);
	}
	else //the evaluation has been done before
		return it->second; //return its reference
}

evaluation & BasisEval::get_w_samples(int nglw, int nw)
{
	std::string key = std::to_string(nglw) + "|" + std::to_string(nw); //create string key

	auto it = this->wMap.find(key); //try to find the evaluation in the map for the given key

	if (it == this->wMap.end()) { //the evaluation for the given key does not exist
		matrix2d<double> wPowers = matrix2d<double>(nglw + 1, nw);
		matrix2d<double> fwPowers = matrix2d<double>(nglw + 1, nw);
		matrix2d<double> fpwPowers = matrix2d<double>(nglw + 1, nw);
		std::vector<double> xglw, wglw;
		functions::gaussk(nglw, xglw, wglw);
		if (basisType == 1) {
			findwPowersLeg(nw, nglw, xglw, wPowers, fwPowers, fpwPowers);
		}
		else {
			findwPowers(nw, nglw, xglw, wPowers, fwPowers, fpwPowers);
		}
		wMap.insert({ key,{ wPowers,fwPowers,fpwPowers } });
		return wMap.at(key);
	}
	else //the evaluation has been done before
		return it->second; //return its reference
}

void BasisEval::matrix_modifiy(matrix2d<double>& mat, int start_x, int end_x, int start_y, int end_y, double & scalar)
{
	for (int i = start_x; i <= end_x; ++i) {
		for (int j = start_y; j <= end_y; ++j) {
			mat(i, j) = mat(i, j)*scalar;
		}
	}
}

void BasisEval::finduPowers(const int& nu, const int& nglu, const std::vector<double>& xglu, matrix2d<double>& uPowers,
	matrix2d<double>& fuPowers, matrix2d<double>& fpuPowers) {
	//upowers goes from (0 to nglu+1, 0 to nu)
	//fuPowers, fpuPowers same as uPowers
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
}

void BasisEval::findvPowers(const int& nv, const int& nglv, const std::vector<double>& xglv, matrix2d<double>& vPowers,
	matrix2d<double>& fvPowers, matrix2d<double>& fpvPowers) {
	//vpower (0 to nglv + 1, 0 to nv)
	//fvPowers, fpvPowers same as vPowers
	double x;

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
}

void BasisEval::findwPowers(const int& nw, const int& nglw, const std::vector<double>& xglw, matrix2d<double>& wPowers,
	matrix2d<double>& fwPowers, matrix2d<double>& fpwPowers) {
	//wPowers (0 to nglw + 1, 0 to nw)
	//fwPowers, fpwPowers same as wPowers
	double x;
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

void BasisEval::finduPowersLeg(const int & nu, const int & nglu, const std::vector<double>& xglu, matrix2d<double>& uPowers, matrix2d<double>& fuPowers, matrix2d<double>& fpuPowers)
{
	double x;
	int kStart = 0;
	matrix2d<double> uPowersL = matrix2d<double>(nglu + 1, nu);
	for (int i = 0; i <= nglu + 1; ++i) {
		x = xglu[i];
		uPowers(i, 0) = 1.0;
		uPowers(i, 1) = x;
		uPowersL(i, 0) = 1.0;
		uPowersL(i, 1) = x;
		fuPowers(i, 0) = 1.0 - x;
		fuPowers(i, 1) = 1.0 + x;
		fpuPowers(i, 0) = -1.0;
		fpuPowers(i, 1) = 1.0;
		for (int j = 2; j <= nu; ++j) {
			uPowers(i, j) = uPowers(i, j - 1)*x;
			uPowersL(i, j) = 0.0;
			fuPowers(i, j) = 0.0;
			fpuPowers(i, j) = 0.0;
			kStart = 1;
			if (j % 2 == 0) {
				kStart = 0;
			}
			for (int k = kStart; k <= j; k += 2) {
				fuPowers(i, j) = fuPowers(i, j) + this->leg2(j, k)*uPowers(i, k);
				uPowersL(i, j) = uPowersL(i, j) + this->leg(j, k)*uPowers(i, k);
				if (k != 0) {
					fpuPowers(i, j) = fpuPowers(i, j) + k * this->leg2(j, k)*uPowers(i, k - 1);
				}
			}
		}
	}
	uPowers = uPowersL;
	double cmhat = sqrt(3.0) / 4.0;
	matrix_modifiy(fuPowers, 0, nglu + 1, 0, 1, cmhat);
	matrix_modifiy(fpuPowers, 0, nglu + 1, 0, 1, cmhat);
	for (int j = 0; j <= nu; ++j) {
		double cm = sqrt(j + 0.5);
		matrix_modifiy(uPowers, 0, nglu + 1, j, j, cm);
		if (j >= 2) {
			cmhat = 0.5*sqrt((2 * j - 3.0) * (2 * j + 1.0) / (2*j - 1.0));
			matrix_modifiy(fuPowers, 0, nglu + 1, j, j, cmhat);
			matrix_modifiy(fpuPowers, 0, nglu + 1, j, j, cmhat);
		}
	}
}

void BasisEval::findvPowersLeg(const int & nv, const int & nglv, const std::vector<double>& xglv, matrix2d<double>& vPowers, matrix2d<double>& fvPowers, matrix2d<double>& fpvPowers)
{
	double x;
	int kStart = 0;
	matrix2d<double> vPowersL = matrix2d<double>(nglv + 1, nv);
	for (int i = 0; i <= nglv + 1; ++i) {
		x = xglv[i];
		vPowers(i, 0) = 1.0;
		vPowers(i, 1) = x;
		vPowersL(i, 0) = 1.0;
		vPowersL(i, 1) = x;
		fvPowers(i, 0) = 1.0 - x;
		fvPowers(i, 1) = 1.0 + x;
		fpvPowers(i, 0) = -1.0;
		fpvPowers(i, 1) = 1.0;
		for (int j = 2; j <= nv; ++j) {
			vPowers(i, j) = vPowers(i, j - 1)*x;
			vPowersL(i, j) = 0.0;
			fvPowers(i, j) = 0.0;
			fpvPowers(i, j) = 0.0;
			kStart = 1;
			if (j % 2 == 0) {
				kStart = 0;
			}
			for (int k = kStart; k <= j; k += 2) {
				fvPowers(i, j) = fvPowers(i, j) + this->leg2(j, k)*vPowers(i, k);
				vPowersL(i, j) = vPowersL(i, j) + this->leg(j, k)*vPowers(i, k);
				if (k != 0) {
					fpvPowers(i, j) = fpvPowers(i, j) + k * this->leg2(j, k)*vPowers(i, k - 1);
				}
			}
		}
	}
	vPowers = vPowersL;
	double cmhat = sqrt(3.0) / 4.0;
	matrix_modifiy(fvPowers, 0, nglv + 1, 0, 1, cmhat);
	matrix_modifiy(fpvPowers, 0, nglv + 1, 0, 1, cmhat);
	for (int j = 0; j <= nv; ++j) {
		double cm = sqrt(j + 0.5);
		matrix_modifiy(vPowers, 0, nglv + 1, j, j, cm);
		if (j >= 2) {
			cmhat = 0.5*sqrt((2 * j - 3.0) * (2 * j + 1.0) / (2 * j - 1.0));
			matrix_modifiy(fvPowers, 0, nglv + 1, j, j, cmhat);
			matrix_modifiy(fpvPowers, 0, nglv + 1, j, j, cmhat);
		}
	}
}

void BasisEval::findwPowersLeg(const int & nw, const int & nglw, const std::vector<double>& xglw, matrix2d<double>& wPowers, matrix2d<double>& fwPowers, matrix2d<double>& fpwPowers)
{
	double x;
	int kStart = 0;
	matrix2d<double> wPowersL = matrix2d<double>(nglw + 1, nw);
	for (int i = 0; i <= nglw + 1; ++i) {
		x = xglw[i];
		wPowers(i, 0) = 1.0;
		wPowers(i, 1) = x;
		wPowersL(i, 0) = 1.0;
		wPowersL(i, 1) = x;
		fwPowers(i, 0) = 1.0 - x;
		fwPowers(i, 1) = 1.0 + x;
		fpwPowers(i, 0) = -1.0;
		fpwPowers(i, 1) = 1.0;
		for (int j = 2; j <= nw; ++j) {
			wPowers(i, j) = wPowers(i, j - 1)*x;
			wPowersL(i, j) = 0.0;
			fwPowers(i, j) = 0.0;
			fpwPowers(i, j) = 0.0;
			kStart = 1;
			if (j % 2 == 0) {
				kStart = 0;
			}
			for (int k = kStart; k <= j; k += 2) {
				fwPowers(i, j) = fwPowers(i, j) + this->leg2(j, k)*wPowers(i, k);
				wPowersL(i, j) = wPowersL(i, j) + this->leg(j, k)*wPowers(i, k);
				if (k != 0) {
					fpwPowers(i, j) = fpwPowers(i, j) + k * this->leg2(j, k)*wPowers(i, k - 1);
				}
			}
		}
	}
	wPowers = wPowersL;
	double cmhat = sqrt(3.0) / 4.0;
	matrix_modifiy(fwPowers, 0, nglw + 1, 0, 1, cmhat);
	matrix_modifiy(fpwPowers, 0, nglw + 1, 0, 1, cmhat);
	for (int j = 0; j <= nw; ++j) {
		double cm = sqrt(j + 0.5);
		matrix_modifiy(wPowers, 0, nglw + 1, j, j, cm);
		if (j >= 2) {
			cmhat = 0.5*sqrt((2 * j - 3.0) * (2 * j + 1.0) / (2 * j - 1.0));
			matrix_modifiy(fwPowers, 0, nglw + 1, j, j, cmhat);
			matrix_modifiy(fpwPowers, 0, nglw + 1, j, j, cmhat);
		}
	}
}

void BasisEval::setup_legendre()
{
	this->leg = matrix2d<double>(11, 11);
	this->leg2 = matrix2d<double>(11, 11);

	leg(0, 0) = 1.0; leg(1, 1) = 1.0; leg(2, 0) = -0.5; leg(2, 2) = 1.5;
	leg(3, 1) = -1.5; leg(3, 3) = 2.5;
	leg(4, 0) = 0.375; leg(4, 2) = -3.75; leg(4, 4) = 4.375;
	leg(5, 1) = 1.875; leg(5, 3) = -8.75; leg(5, 5) = 7.875;
	leg(6, 0) = -0.3125; leg(6, 2) = 6.5625; leg(6, 4) = -19.6875; leg(6, 6) = 14.4375;
	leg(7, 1) = -2.1875; leg(7, 3) = 19.6875; leg(7, 5) = -43.3125; leg(7, 7) = 26.8125;
	leg(8, 0) = 0.2734375; leg(8, 2) = -9.84375; leg(8, 4) = 54.140625; leg(8, 6) = -93.84375; leg(8,8) = 50.2734375;
	leg(9, 1) = 2.4609375; leg(9, 3) = -36.09375; leg(9, 5) = 140.765625; leg(9, 7) = -201.09375; leg(9,9) = 94.9609375;
	leg(10, 0) = -0.24609375; leg(10, 2) = 13.53515625; leg(10, 4) = -117.3046875; leg(10, 6) = 351.9140625; leg(10, 8) = -427.32421875; leg(10,10) = 180.42578125;
	leg(11, 1) = -2.70703125; leg(11, 3) = 58.65234375; leg(11, 5) = -351.9140625; leg(11, 7) = 854.6484375; leg(11, 9) = -902.12890625; leg(11,11) = 344.44921875;

	leg2(0, 0) = 1.0; leg2(0, 1) = -1.0; leg2(1, 0) = 1.0; leg2(1, 1) = 1.0;
	leg2(2, 0) = -1.5; leg2(2, 2) = 1.5; leg2(3, 1) = -2.5; leg2(3, 3) = 2.5;
	leg2(4, 0) = 0.875; leg2(4, 2) = -5.25; leg2(4, 4) = 4.375;
	leg2(5, 1) = 3.375; leg2(5, 3) = -11.25; leg2(5, 5) = 7.875;
	leg2(6, 0) = -0.6875; leg2(6, 2) = 10.3125; leg2(6, 4) = -24.0625; leg2(6,6) = 14.4375;
	leg2(7, 1) = -4.0625; leg2(7, 3) = 28.4375; leg2(7, 5) = -51.1875; leg2(7,7) = 26.8125;
	leg2(8, 0) = 0.5859375; leg2(8, 2) = -16.40625; leg2(8, 4) = 73.828125; leg2(8, 6) = -108.28125; leg2(8,8)= 50.2734375;
	leg2(9, 1) = 4.6484375; leg2(9, 3) = -55.78125; leg2(9, 5) = 184.078125; leg2(9, 7) = -227.90625; leg2(9,9) = 94.9609375;
	leg2(10, 0) = -0.51953125; leg2(10, 2) = 23.37890625; leg2(10, 4) = -171.4453125; leg2(10, 6) = 445.7578125; leg2(10, 8) = -477.59765625; leg2(10,10) = 180.42578125;
	leg2(11, 1) = -5.16796875; leg2(11, 3) = 94.74609375; leg2(11, 5) = -492.6796875; leg2(11, 7) = 1055.7421875; leg2(11, 9) = -997.08984375; leg2(11,11) = 344.44921875;

}

