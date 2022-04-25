#include "qoi_helper.h"
bool qoi_error::is_higher(int& u_order, int& v_order, int& w_order,
	int& iuvwh, int& ih, int& jh, int& kh) {
	if (ih > u_order || jh > v_order || kh > w_order) return true;
	else switch (iuvwh) {
	case 1:
		if (ih > u_order - 1) {
			return true;
		}
		break;
	case 2:
		if (jh > v_order - 1) {

			return true;

		}
		break;
	case 3:
		if (kh > w_order - 1) {

			return true;

		}
		break;
	}
	return false;
}