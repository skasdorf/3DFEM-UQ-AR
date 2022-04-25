#include "refinement_tests.h"

Point tests::midpoint(Point & p1, Point & p2)
{
	Point p3((p1.x+p2.x)/2.0, (p1.y+p2.y)/2.0, (p1.z + p2.z)/2.0);
	return p3;
}

void tests::test(std::vector<int> refine_elements, Domain & dom, double factor)
{
	std::vector<Point> points;
	points = { Point(-1.0, -1.0, -1.0), Point(0.0, -1.0, -1.0),
		Point(1.0, -1.0, -1.0), Point(-1.0, 0.0, -1.0),Point(0.0, 0.0, -1.0),
		Point(1.0, 0.0, -1.0) ,Point(-1.0, 1.0, -1.0) ,Point(0.0, 1.0, -1.0),
		Point(1.0, 1.0, -1.0) };
	for (int i = 0; i < 9; ++i) {
		points.push_back(points[i] + Point(0.0, 0.0, 1.0));
	}
	for (int i = 0; i < 9; ++i) {
		points.push_back(points[i] + Point(0.0, 0.0, 2.0));
	}
	for (auto it = refine_elements.begin(); it != refine_elements.end(); ++it) {
		Element e1 = dom.elements[*it - 1];
		//temp condition b/c material isn't debugged yet!!!!!!!!!!!!!!
		for (int p = 0; p < points.size(); ++p) {
			auto r = unitVectorsM::findR(e1.geom_order, points[p].x, points[p].y, points[p].z, e1.rs, e1.nRs);
			auto r1 = Point(r[1], r[2], r[3]);
			auto r2 = dom.nodes[e1.all_indices[p] - 1];
		}
	}
}

std::vector<Point> tests::get_base_coords(int geom_order)
{
	if (geom_order == 2) {
		std::vector<Point> points = { Point(-1.0, -1.0, -1.0), Point(0.0, -1.0, -1.0),
			Point(1.0, -1.0, -1.0), Point(-1.0, 0.0, -1.0),Point(0.0, 0.0, -1.0),
			Point(1.0, 0.0, -1.0) ,Point(-1.0, 1.0, -1.0) ,Point(0.0, 1.0, -1.0),
			Point(1.0, 1.0, -1.0) };
		for (int i = 0; i < 9; ++i) {
			points.push_back(points[i] + Point(0.0, 0.0, 1.0));
		}
		for (int i = 0; i < 9; ++i) {
			points.push_back(points[i] + Point(0.0, 0.0, 2.0));
		}
		return points;
	}
	else if (geom_order == 1) {

	}
	return std::vector<Point>();
}

std::vector<Point> tests::get_inner_coords(int geom_order, std::vector<Point>& base_coords, double factor)
{
	std::vector<Point> points;
	if (geom_order == 2) {
		for (int i = 0; i < base_coords.size(); ++i) {
			if (i != 13) { //at index 13, the parametric point is zero
				points.push_back(base_coords[i]*factor);
			}
		}
	}
	return points;
}

std::vector<Point> tests::get_inbetween_coords(int geom_order, std::vector<Point>& base_coords, std::vector<Point>& inner_coords)
{
	std::vector<Point> points;
	if (geom_order == 2) {
		for (int i = 0; i < base_coords.size(); ++i) {
			if (i < 13) {
				points.push_back(midpoint(base_coords[i], inner_coords[i]));
			}if (i > 13) {
				points.push_back(midpoint(base_coords[i], inner_coords[i-1])); //-1 shift b/c inner coords has one less (it shares point at index 13)
			}
		}
	}
	return points;
}

void tests::split_elements(Domain & dom, std::vector<int>& refine_elements, double factor)
{
	int geom_order = dom.elements[0].geom_order;
	//get base, inner, and in between coordinates in parametric domain
	auto base_coords = get_base_coords(geom_order);
	auto inner_coords = get_inner_coords(geom_order, base_coords, factor);
	auto inbetween_coords = get_inbetween_coords(geom_order, base_coords, inner_coords);
	for (int i = 0; i < refine_elements.size(); ++i) { //refine elements starts indexing at 0
		Element e1 = dom.elements[refine_elements[i]];
		std::vector<int> coords;
		//add the the elements original vertices
		coords.insert(coords.begin(), e1.all_indices.begin(), e1.all_indices.end());
		//add the new inner coordinates (after finding their real values)
		for (int j = 0; j < inner_coords.size(); ++j) {
			auto r = unitVectorsM::findR(e1.geom_order, inner_coords[j].x, inner_coords[j].y, inner_coords[j].z, e1.rs, e1.nRs);
			auto r1 = Point(r[1], r[2], r[3]);
			int index = dom.nodes.size() + 1;
			dom.nodes.push_back(r1);
			coords.push_back(index);
		}
		//add the in between coordinates
		for (int j = 0; j < inbetween_coords.size(); ++j) {
			auto r = unitVectorsM::findR(e1.geom_order, inbetween_coords[j].x, inbetween_coords[j].y, inbetween_coords[j].z, e1.rs, e1.nRs);
			auto r1 = Point(r[1], r[2], r[3]);
			int index = dom.nodes.size() + 1;
			dom.nodes.push_back(r1);
			coords.push_back(index);
		}
		coord_table index_list;
		auto indices = index_list.indices3;
		for (int j = 0; j < 7; ++j) { //for each piece the element is being split into
			Element e_test;
			//grab the list of indices for this segment
			auto element_indices = indices[j];
			//grab indices from domain corresponding to these indices
			std::vector<int> order_pt_indices;
			for (int k = 0; k < element_indices.size(); ++k) {
				order_pt_indices.push_back(coords[element_indices[k]-1]);
			}
			//get material stuff (right now only works for homogoenous
			std::vector<std::vector<std::complex<double>>> new_epsr;
			std::vector<std::vector<std::complex<double>>> new_mur;
			bool do_extra_mats = false;
			if (e1.materials.hcode == 1 && e1.materials.icode == 1) {
				new_epsr = e1.materials.epsr_list;
				new_mur = e1.materials.mur_list;
			}
			else {
				continue;
			}
			e_test.all_indices = order_pt_indices;
			e_test.expansion = e1.expansion;
			e_test.quadrature = e1.quadrature;
			e_test.materials.mur_list = new_mur;
			e_test.materials.epsr_list = new_epsr;
			e_test.index = dom.elements.size() + 1;
			e_test.materials.hcode = e1.materials.hcode;
			e_test.materials.icode = e1.materials.icode;
			e_test.materials.pmlcode = e1.materials.pmlcode;
			e_test.materials.Kuvw = e1.materials.Kuvw;
			e_test.materials.sym = e1.materials.sym;
			e_test.geom_order = e1.geom_order;
			if (j< 6) dom.elements.push_back(e_test);
			else {
				e_test.index = e1.index;
				dom.elements[e1.index - 1] = e_test;
			}
		}
	}
	//print out all necessary information to then redo the mesh
	std::ofstream file_elem, file_geom, file_bc, file_basic, file_material;
	file_geom.open("../Geometry_refined.dat");
	file_geom << "!Refined geometry!" << std::endl;
	file_elem.open("../Element_refined.dat");
	file_elem << "!Refined elements!" << std::endl;
	file_bc.open("../Boundary_refined.dat");
	file_bc << "!Refined boundary!" << std::endl;
	file_basic.open("../Basic_refined.dat");
	file_basic << "!Basic info file" << std::endl;
	file_material.open("../exampleFiles/" + dom.mesh_name + "_refine/materials_refined.txt");
	file_material << "!Material file" << std::endl;


	for (auto dom_points = dom.nodes.begin(); dom_points != dom.nodes.end(); ++dom_points) {
		file_geom << dom_points->x << " " << dom_points->y << " " << dom_points->z << std::endl;
	}
	file_geom.close();
	for (auto elem_p = dom.elements.begin(); elem_p != dom.elements.end(); ++elem_p) {
		if (elem_p->index == 447) {
			int pause = 0;
		}
		file_elem << "Index: " << elem_p->index << " " << elem_p->materials.pmlcode << " " << elem_p->geom_order;
		file_elem << std::endl << "Basis: " << elem_p->expansion[0] << " " << elem_p->expansion[1] << " " << elem_p->expansion[2];
		file_elem << std::endl << "Quadrature: " << elem_p->quadrature[0] << " " << elem_p->quadrature[1] << " " << elem_p->quadrature[2] << std::endl;
		file_elem << "Nodes: ";

		for (auto elm_nodes = elem_p->all_indices.begin(); elm_nodes != elem_p->all_indices.end(); ++elm_nodes) {

			file_elem << *elm_nodes << " ";
		}


		file_elem << std::endl;
		file_material << elem_p->index << " " << elem_p->materials.hcode << " " << elem_p->materials.icode << " " << elem_p->materials.pmlcode << std::endl;
		file_material << elem_p->materials.Kuvw << " " << elem_p->materials.sym << std::endl;
		file_material << "eps:";// << std::endl;

		for (auto mat_ent = elem_p->materials.epsr_list.begin(); mat_ent != elem_p->materials.epsr_list.end(); ++mat_ent) {
			file_material << std::endl;
			//iterates through nodes
			for (auto val_mat = mat_ent->begin(); val_mat != mat_ent->end(); ++val_mat) {
				//gets values
				file_material << " " << val_mat->real() << "," << val_mat->imag();
			}
		}
		file_material << std::endl << "mu:";
		for (auto mat_ent = elem_p->materials.mur_list.begin(); mat_ent != elem_p->materials.mur_list.end(); ++mat_ent) {
			file_material << std::endl;
			//iterates through nodes
			for (auto val_mat = mat_ent->begin(); val_mat != mat_ent->end(); ++val_mat) {
				//gets values
				file_material << " " << val_mat->real() << "," << val_mat->imag();
			}
		}
		file_material << std::endl;
	}
	file_elem.close();
	file_material.close();
	//save boundary condition info
	int bc_index = 1;
	for (auto fac_bc = dom.facets.begin(); fac_bc != dom.facets.end(); ++fac_bc) {
		if (fac_bc->boundary_condition == 0) continue;
		file_bc << bc_index++ << " " << fac_bc->vertices[0] << " " << fac_bc->vertices[1] << " " << fac_bc->vertices[2] << " " << fac_bc->vertices[3] << " " << fac_bc->boundary_condition << std::endl;
	}
	file_bc.close();

	//MODE - NODES - NOEL - NBC - NGEN - NWAVES - NPORTS - FSTART - FSTOP - NFR
	file_basic << dom.sc.mode_of_operation << " " << dom.nodes.size() << " " << dom.elements.size() << " " << dom.sc.nbc << " " << dom.sc.ngen << " " << dom.sc.numberOfWaves << " " << dom.sc.nports << " " << dom.sc.fstart << " " << dom.sc.fstop << " " << dom.sc.nfr << std::endl;
	file_basic.close();
}

