#include "Mesh.h"

double sq(double a)
{
	return a*a;
}


void Split(string str, string* (&parts), int n, int& i)
{
	//const int N = 256;      //Максимальная длина строки
	//char word[N] = {};          //Буфер для считывания строки
	//stringstream x; //Создание потоковой переменной
	//x << str;                //Перенос строки в поток

	//i = 0;
	//while (x >> word) {
	//	//cout << "word " << word << endl;

	//	//if (i < n) parts[i] = word;
	//	parts[i] = word;
	//	//cout << "parts[i] " << parts[i] << endl;

	//	i++;
	//}
}

Mesh::Mesh(string file)
{
	ReadMesh(file);
	FaceGeom();
	CellGeom();
	FaceNormals();
	CellCoeffs();
}

void Mesh::ReadMesh(string file)
{
	file = "Mesh/" + file;
	ifstream reading;
	reading.open(file);

	if (reading)
	{
		// Поиск слова "(2"
		string word = "(2";
		string tmp, cdim;
		while (reading >> tmp)
		{
			if (word == tmp) {
				reading >> cdim;
				break;
			}
		}
		cout << " cdim = " << cdim << endl;

		// Поиск слова "(10" - для определения числа узлов
		word = "(10";
		while (reading >> tmp)
		{
			if (word == tmp) {
				break;
			}
		}
		reading >> tmp >> tmp;
		reading >> tmp;

		nNodes = stoi(tmp, 0, 16);
		cout << " total_nodes = " << nNodes << endl;

		// Выделение памяти под узлы
		nodes = new Node[nNodes];

		// Cells
		// Поиск слова "(12" - для определения числа ячеек
		word = "(12";
		string w2 = "(12(0";
		while (reading >> tmp)
		{
			if (word == tmp || w2 == tmp) {
				if (word == tmp) reading >> tmp;
				break;
			}
		}
		string s1, s2;
		reading >> s1 >> s2;

		int i1 = stoi(s1, 0, 16);
		int i2 = stoi(s2, 0, 16);

		//cout << "i1= " << i1 << ", i2= " << i2 << endl;
		nCells = i2 - i1 + 1;

		cout << " nCells = " << nCells << endl;

		// Выделение памяти под ячейки
		cells = new Cell[nCells];
		//exit(44);

		// Возвращение к началу файла
		reading.clear();
		reading.seekg(0L, std::ios_base::beg);

		cout << "                         Считывание узлов " << endl;
		// Считывание узлов
		// Перемещение на одно из слов:
		word = "(";
		w2 = "2)(";
		string w3 = "3)(";

		while (reading >> tmp)
		{
			if (word == tmp || w2 == tmp || w3 == tmp) {
				break;
			}
		}
		//cout << "tmp = " << tmp << endl;

		int node_counter = 0;
		while (node_counter < nNodes) {
			reading >> tmp;
			if (tmp == "))") {

			}
			else if (tmp == "(10") {
				for (int j = 0; j < 5; j++) {
					reading >> tmp;

					cout << "(10 tmp = " << tmp << endl;
				}
			}
			else {
				nodes[node_counter].x = atof(tmp.c_str());
				for (int j = 1; j < 2; j++) {
					reading >> nodes[node_counter].y;

				}

				node_counter++;
			}

		}
		cout << "node_counter= " << node_counter << endl;
		//cout << nodes[nNodes - 2].x << " " << nodes[nNodes - 2].y << endl;
		//exit(55);
		//****************************************************************
		// Считывание зон
		// ***************************************************************
		cout << "                         Считывание зон " << endl;
		// Перемещение на одно из слов:
		word = "(39";
		w2 = "(45";
		while (reading >> tmp)
		{
			if (word == tmp || w2 == tmp) {
				break;
			}
		}

		//reading >> tmp;
		//cout << " zones " << endl;
		
		while(!reading.eof())
		{	
			reading >> tmp;
			if (word != tmp && w2 != tmp && !reading.eof())
			{
				Zone z;
				
				tmp = tmp.substr(1);
				z.zone_id = stoi(tmp);
				reading >> z.zone_type >> tmp;
				size_t start = tmp.find(")");
				z.zone_name = tmp.substr(0, start);
				cout << "ZONES: id, type, name: " << z.zone_id << " " << z.zone_type << " " << z.zone_name << "\n";
				zones.push_back(z);
			}
		}

		nZones = static_cast<int>(zones.size());
		cout << "nZones= " << nZones  << endl;
		//exit(44);

		// Возвращение к началу файла
		reading.clear();
		reading.seekg(0L, std::ios_base::beg);


		// Считывание граней + zones
		cout << "                         Считывание граней " << endl;

		// Перемещение на слово:
		word = "(2";
		while (reading >> tmp)
		{
			if (word == tmp) {
				reading >> cdim;
				break;
			}
		}

		// Перемещение на одно из слов:
		word = "(13";
		w2 = "(13(0";
		while (reading >> tmp)
		{
			if (word == tmp || w2 == tmp) {
				if (word == tmp) reading >> tmp;
				break;
			}
		}
		
		// Определение числа граней 
		reading >> s1 >> s2;
		//int i1, i2;
		i1 = stoi(s1, 0, 16);
		i2 = stoi(s2, 0, 16);

		//cout << "s1= " << s1 << ", s2= " << s2 << endl;
		//cout << "i1= " << i1 << ", i2= " << i2 << endl;

		nFaces = i2 - i1 + 1; // общее число граней

		cout << " total_faces = " << nFaces << endl;
		
		// Выделение памяти под грани
		faces = new Face[nFaces];
		
		word = "(13";
		w2 = "(13(0";
		// Цикл по числу зон -1  (одна из зон не используется под грани)
		for (int k = 0; k < nZones-1; k++)
		{
			bool find = false;
			// Поиск одного одного из слов: "(13" или "(13(0"
			while (reading >> tmp)
			{
				if (word == tmp || w2 == tmp) {
					if (word == tmp) reading >> tmp;
					find = true;
					break;
				}
			}

			// если нашли:
			if (find)
			{
				tmp = tmp.substr(1);
				int id = stoi(tmp, 0, 16); // Идентификатор зон
				//cout << " tmp = " << tmp << " id= " << id << endl;
				
				// Ищем номер зоны, у которой идентификатор равен id
				int z_id;
				for (int m = 0; m < nZones; m++)
				{
					if (zones[m].zone_id == id) z_id = m;
				}
				
				//cout << "z_id= " << z_id << ", zones[z_id].zone_name= " << zones[z_id].zone_name
				//	<< ", zones[z_id].zone_type= " << zones[z_id].zone_type << endl;

				// Определение начального и конечного номеров граней, относящихся к этой зоне
				reading >> s1 >> s2 >> tmp >> tmp;
				i1 = stoi(s1, 0, 16);
				i2 = stoi(s2, 0, 16);
				zones[z_id].if1 = i1 - 1;
				zones[z_id].if2 = i2 - 1;

				//cout << "zones[z_id].if1= " << zones[z_id].if1 << ", zones[z_id].if2= "
				//	<< zones[z_id].if2 << endl;

				// Считывание для этих граней номеров узлов и номеров ячеек
				for (int i = i1 - 1; i < i2; i++)
				{
					string sn1, sn2, scr, scl;
					reading >> sn1 >> sn2 >> scr >> scl;
					
					//cout << sn1 << " " << sn2 << " " << scr << " " << scl << endl;
					int n1 = stoi(sn1, 0, 16) - 1;
					int n2 = stoi(sn2, 0, 16) - 1;
					int cr = stoi(scr, 0, 16) - 1;
					int cl = stoi(scl, 0, 16) - 1;
					//cout << n1 << " " << n2 << " " << cr << " " << cl << endl;

					faces[i].v[0] = n1;
					faces[i].v[1] = n2;
					faces[i].il = cl;
					faces[i].ir = cr;
					faces[i].zone_id = zones[z_id].zone_id;
					//cout << i << "\t" << faces[i].zone_id  << "\t" << zones[z_id].zone_id << endl;



					cells[cr].fid.push_back(i);
					if (cl>-1) cells[cl].fid.push_back(i);
				}
			}
			//if (k==1) exit(78);
		}

		//cout << "cells[*].fid.size()= " << cells[10].fid.size() << endl;
		//cout << "fid[0]= " << cells[10].fid[0] << ", fid[1]= " << cells[10].fid[1] 
		//	<< ", fid[2]= " << cells[10].fid[2] << ", fid[3]= " << cells[10].fid[3] << endl;
		//exit(77);

		reading.close();
	}
	else
	{
		cout << "File \"" << file << "\" could not be opened" << endl;
		exit(1);
	}
}

void Mesh::CellGeom()
{
	// Площадь и ц.т. ячеек
	cout << "                  Определение ц.т. и площади ячеек" << endl;
	for (int i = 0; i < nCells; i++) {
		cells[i].num_f = static_cast<int>(cells[i].fid.size());

		int n = cells[i].num_f;

		string d = "cw";
		int* point_ind = new int[n];

		SortNodes(i, d, point_ind);

		cells[i].nodes = new int[n];
		for (int j = 0; j < n; j++) {
			cells[i].nodes[j] = point_ind[j];
		}

		//cells[i].nodes = new int[n];
		int* nds = new int[n];

		for (int j = 0; j < n; j++) {
			nds[j] = point_ind[j];
		}

		double* xn = new double[n + 1];
		double* yn = new double[n + 1];

		for (int j = 0; j < n; j++) {
			xn[j] = nodes[point_ind[j]].x;
			yn[j] = nodes[point_ind[j]].y;
		}
		xn[n] = xn[0];
		yn[n] = yn[0];

		double A = 0, Cx = 0, Cy = 0;
		for (int j = 0; j < n; j++) {
			A += xn[j] * yn[j + 1] - xn[j + 1] * yn[j];
			Cx += (xn[j] + xn[j + 1]) * (xn[j] * yn[j + 1] - xn[j + 1] * yn[j]);
			Cy += (yn[j] + yn[j + 1]) * (xn[j] * yn[j + 1] - xn[j + 1] * yn[j]);
		}
		A = 0.5 * A;
		cells[i].x = Cx / (6. * A);
		cells[i].y = Cy / (6. * A);
		cells[i].S = abs(A);


		//cout << cells[i].x << ", " << cells[i].y << ", " << cells[i].S << endl;

		//exit(9);
		// Минимальный размер ячейки
		cells[i].dl = 1.e10;
		for (int j = 0; j < n; j++)
		{
			int m = cells[i].fid[j];
			double a = faces[m].s;
			if (cells[i].dl > a) cells[i].dl = a;
		}


	}

}

void Mesh::SortNodes(int k, string d, int* (&point_ind))
{
	int n = static_cast<int>(cells[k].fid.size()); // число граней, окружающих k-ую ячейку

	int* nf = new int[n];
	int* cfaces = new int[n];

	// выбираем начальную грань
	int i = 0;
	int f = cells[k].fid[i]; // номер начальной грани

	cfaces[0] = f;
	// узлы этой грани
	int v1 = faces[f].v[0];
	int v2 = faces[f].v[1];

	if (d == "cw") {
		if (faces[f].ir == k) {
			point_ind[0] = v2; // v1;
			point_ind[1] = v1; // v2;
		}
		else {
			point_ind[0] = v1;
			point_ind[1] = v2;
		}
	}
	else if (d == "ccw") {
		if (faces[f].ir == k) {
			point_ind[0] = v1;
			point_ind[1] = v2;
		}
		else {
			point_ind[0] = v2;
			point_ind[1] = v1;
		}
	}

	// узел v1 выбираем в качестве начальной точки обхода

	for (int j = 1; j < (n - 1); j++) {

		// ищем грань, которая содержит point_ind[2], но не относится к вышеперечисленным
		for (int i = 0; i < n; i++) {
			int f = cells[k].fid[i];
			bool b = false;
			for (int m = 0; m < j; m++) {
				if (f == cfaces[m]) {
					b = true;
				}
			}
			if (b) continue;

			int vi = faces[f].v[0];
			int vi1 = faces[f].v[1];

			if (vi == point_ind[j]) {
				cfaces[j] = f;
				point_ind[j + 1] = vi1;
				break;
			}

			if (vi1 == point_ind[j]) {
				cfaces[j] = f;
				point_ind[j + 1] = vi;
				break;
			}

		}

	}

}

void Mesh::FaceGeom()
{
	// длина и центр граней
	for (int i = 0; i < nFaces; i++) {

		//faces[i].centr = new double[2];

		double x0 = nodes[faces[i].v[0]].x;
		double y0 = nodes[faces[i].v[0]].y;

		double x1 = nodes[faces[i].v[1]].x;
		double y1 = nodes[faces[i].v[1]].y;

		faces[i].x = 0.5 * (x0 + x1);
		faces[i].y = 0.5 * (y0 + y1);

		faces[i].s = sqrt(sq(x1 - x0) + sq(y1 - y0));


		//cout << i << endl;
		//cout << faces[i].x << ", " << faces[i].y << ", " << faces[i].s << endl;

		//if (i==10) exit(10);

	}

}

void Mesh::ReadBounds()
{
	string file = "Input/bounds.dat";
	ifstream reading;
	reading.open(file);

	for (int m = 0; m < nZones; m++)
	{
		zones[m].b_id = -1;
		zones[m].tid = 0;
	}

	if (reading)
	{
		reading >> nBounds;
		string com;
		bounds = new Boundary[nBounds];
		for (int i = 0; i < nBounds; i++)
		{
			reading >> bounds[i].zone_id >> bounds[i].tp_id;

			cout << "zone_id= " << bounds[i].zone_id << endl;

			reading >> bounds[i].stp_id >> bounds[i].nVals;
			if (bounds[i].nVals > 0)
			{
				bounds[i].values = new double[bounds[i].nVals];

				for (int j = 0; j < bounds[i].nVals; j++)
					reading >> bounds[i].values[j];

			}
			reading >> com;
			//ищем зону с идентификатором bounds[i].zone_id
			for (int m = 0; m < nZones; m++)
			{
				//cout << "m= " << m << " i= " << i << endl;
				//cout << "zones[m].zone_id " << zones[m].zone_id << " bounds[i].zone_id " << bounds[i].zone_id << endl;
				if (zones[m].zone_id == bounds[i].zone_id)
				{
					zones[m].b_id = i;
					zones[m].tid = bounds[i].tp_id;
				}

			}
		}

		exit(10000);

		for (int m = 0; m < nZones; m++)
			cout << m << " zones[m].name= " << zones[m].zone_name
			<< ", zones[m].b_id= " << zones[m].b_id << endl;

		reading.close();

		//exit(44);
	}
	else
	{
		cout << " Boundaries file could not be opened" << endl;
		exit(2);
	}
}

void Mesh::write_OF()
{
	string work_dir = "test";

	string f = work_dir + "//foam.foam";
	ofstream record(f, ios::out);
	record.close();

	string meshfold = work_dir + "//constant//polyMesh";

	cout << "Writing faces" << endl;
	string f2 = meshfold + "//boundary";
	write_OF_boundary(f2);

	cout << "Writing points" << endl;
	string f3 = meshfold + "//points";
	write_OF_points(f3);

	cout << "Writing faces" << endl;
	string f4 = meshfold + "//faces";
	write_OF_faces(f4);

	cout << "Writing owner/neighbour" << endl;
	string f5 = meshfold + "//owner";
	write_OF_owner(f5);

	string f6 = meshfold + "//neighbour";
	write_OF_neighbour(f6);
}

void Mesh::write_OF_boundary(string f)
{
	ofstream record(f, ios::out);
	if (record) {
		int num_of_boundaries;
		if (mesh_dim == 3) {
			num_of_boundaries = 0;
			for (int m = 0; m < nZones; m++) {
				//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			}
		}
		else {
			num_of_boundaries = 0;
			for (int m = 0; m < nZones; m++) {
				if (zones[m].if1 != 0) num_of_boundaries++;
			}
		}
		cout << "num_of_boundaries = " << num_of_boundaries << endl;

		string S(9, ' ');

		record << "FoamFile" << endl;
		record << "{" << endl;
		record << "version" << S << "2.0;" << endl;
		record << "format" << S << "ascii;" << endl;
		record << "class" << S << "polyBoundaryMesh;" << endl;
		record << "location" << S << "constant/polyMesh;" << endl;
		record << "object" << S << "boundary;" << endl;
		record << "}" << endl;
		record << "//***************************************************" << endl;
		record << setw(3) << num_of_boundaries << endl;
		record << "(" << endl;

		int buff;
		for (int i = 0; i < nZones; i++) {

			if (zones[i].tid > 0) {

				int nFaces = zones[i].if2 - zones[i].if1 + 1;
				int startFace = zones[i].if1 - 1;
				buff = startFace + nFaces;
				record << zones[i].zone_name << endl;
				record << "{" << endl;
				record << "type" << S << "patch;" << endl;

				record << " nFaces" << S << setw(12) << nFaces << ";" << endl;
				record << " startFace" << S << setw(12) << startFace << ";" << endl;
				record << "}" << endl;
				record << " " << endl;
			}

		}

		if (mesh_dim == 2) { // Make artificial symmetry planes for 2d mesh (2d to 3d conversion)
			int nFaces = nCells;
			int startFace = buff;
			record << "Interior_artificial_1" << endl;
			record << "{" << endl;
			record << "type" << S << "empty;" << endl;
			record << "inGroups" << S << "1(empty);" << endl;
			record << " nFaces" << setw(12) << S << nFaces << ";" << endl;
			record << " startFace" << setw(12) << S << startFace << ";" << endl;
			record << "}" << endl;
			record << " " << endl;

			startFace = startFace + nFaces;
			record << "Interior_artificial_2" << endl;
			record << "{" << endl;
			record << "type" << S << "empty;" << endl;
			record << "inGroups" << S << "1(empty);" << endl;
			record << " nFaces" << setw(12) << S << nFaces << ";" << endl;
			record << " startFace" << setw(12) << S << startFace << ";" << endl;
			record << "}" << endl;
			record << " " << endl;
		}
		record << ")" << endl;
		record << "//***************************************************" << endl;
	}
	record.close();
}

void Mesh::write_OF_points(string f)
{
	ofstream record(f, ios::out);
	if (record) {
		int nPoints;
		if (mesh_dim == 3) {
			nPoints = nNodes;
		}
		else {
			nPoints = 2 * nNodes;
		}
		string S(9, ' ');
		record << "FoamFile" << endl;
		record << "{" << endl;
		record << "version" << S << "2.0;" << endl;
		record << "format" << S << "ascii;" << endl;
		record << "class" << S << "vectorField;" << endl;
		record << "location" << S << "constant/polyMesh;" << endl;
		record << "object" << S << "points;" << endl;
		record << "}" << endl;
		record << "//***************************************************" << endl;
		record << setw(10) << nPoints << endl;
		record << "(" << endl;

		if (mesh_dim == 3) {
			for (int i = 0; i < nNodes; i++) {
				record << "(";
				record << scientific << setw(20) << nodes[i].x << " ";
				record << scientific << setw(20) << nodes[i].y << " ";
				record << scientific << setw(20) << 0.0 << ")" << endl;
			}
		}
		else {
			for (int i = 0; i < nNodes; i++) {
				record << "(";
				record << scientific << setw(20) << nodes[i].x << " ";
				record << scientific << setw(20) << nodes[i].y << " ";
				record << scientific << setw(20) << 0.0 << ")" << endl;
			}
			for (int i = 0; i < nNodes; i++) {
				record << "(";
				record << scientific << setw(20) << nodes[i].x << " ";
				record << scientific << setw(20) << nodes[i].y << " ";
				record << scientific << setw(20) << 1.e-6 << ")" << endl;
			}

		}

		record << ")" << endl;
		record << "//***************************************************" << endl;



	}
	record.close();
}

void Mesh::write_OF_owner(string f)
{
	ofstream record(f, ios::out);
	if (record) {
		int nFcs;
		if (mesh_dim == 3) {
			nFcs = nFaces;
		}
		else {
			nFcs = nFaces + 2 * nCells;
		}

		string S(9, ' ');

		record << "FoamFile" << endl;
		record << "{" << endl;
		record << "version" << S << "2.0;" << endl;
		record << "format" << S << "ascii;" << endl;
		record << "class" << S << "labelList;" << endl;
		record << "location" << S << "constant/polyMesh;" << endl;
		record << "object" << S << "owner;" << endl;
		record << "}" << endl;
		record << "//***************************************************" << endl;
		record << setw(10) << nFcs << endl;
		record << "(" << endl;

		if (mesh_dim == 3) {
			for (int i = 0; i < nFcs; i++)
				record << faces[i].ir << endl;
		}
		else {
			for (int i = 0; i < nFaces; i++)
				record << faces[i].ir << endl;
			for (int i = 0; i < nCells; i++)
				record << i << endl;
			for (int i = 0; i < nCells; i++)
				record << i << endl;
		}

		record << ")" << endl;
		record << "//***************************************************" << endl;

	}
	record.close();
}

void Mesh::write_OF_neighbour(string f)
{
	ofstream record(f, ios::out);
	if (record) {
		int nFcs;
		if (mesh_dim == 3) {
			nFcs = nFaces;
		}
		else {
			nFcs = nFaces + 2 * nCells;
		}

		string S(9, ' ');

		record << "FoamFile" << endl;
		record << "{" << endl;
		record << "version" << S << "2.0;" << endl;
		record << "format" << S << "ascii;" << endl;
		record << "class" << S << "labelList;" << endl;
		record << "location" << S << "constant/polyMesh;" << endl;
		record << "object" << S << "neighbour;" << endl;
		record << "}" << endl;
		record << "//***************************************************" << endl;
		record << setw(10) << nFcs << endl;
		record << "(" << endl;

		if (mesh_dim == 3) {
			for (int i = 0; i < nFcs; i++)
				record << setw(12) << faces[i].il << endl;
		}
		else {
			for (int i = 0; i < nFaces; i++)
				record << setw(12) << faces[i].il << endl;

			for (int i = 0; i < nCells; i++)
				record << setw(12) << -1 << endl;

			for (int i = 0; i < nCells; i++)
				record << setw(12) << -1 << endl;
		}

		record << ")" << endl;
		record << "//***************************************************" << endl;

	}
	record.close();
}

void Mesh::FaceNormals()
{
	for (int i = 0; i < nFaces; i++) {

		int n1 = faces[i].v[0];
		int n2 = faces[i].v[1];

		double xNd1 = nodes[n1].x;
		double xNd2 = nodes[n2].x;
		double yNd1 = nodes[n1].y;
		double yNd2 = nodes[n2].y;

		double dl = sqrt(sq(xNd2 - xNd1) + sq(yNd2 - yNd1));

		faces[i].nx = -(yNd2 - yNd1) / dl;
		faces[i].ny = (xNd2 - xNd1) / dl;
	}
}

void Mesh::grad_coeffs_least_sq(Pnt* A7, int q)
{
	// A7 - массив координат соседних точек [nn]
	// q - номер ячейки

	// X7 - координаты центра ячейки
	Pnt X_C;
	//double X7[2];
	X_C.x[0] = cells[q].x;
	X_C.x[1] = cells[q].y;

	// nn - число соседей
	int nn = cells[q].num_f;

	cells[q].coeffs = new double* [nn];
	double** delta_X = new double* [nn];
	for (int i = 0; i < nn; i++) {
		cells[q].coeffs[i] = new double[mesh_dim];
		delta_X[i] = new double[mesh_dim];
	}

	// weight coefs and distance between cells
	double* w = new double[nn];
	double* dist = new double[nn];

	double** A, ** M;	//matrix of the system and its inverse
	A = new double* [mesh_dim];
	M = new double* [mesh_dim];
	for (int i = 0; i < mesh_dim; i++) {
		A[i] = new double[mesh_dim];
		M[i] = new double[mesh_dim];
	}

	//!Initialization!
	for (int i = 0; i < mesh_dim; i++) {
		for (int j = 0; j < mesh_dim; j++) {
			A[i][j] = 0.;
			M[i][j] = 0.;
		}
	}

	for (int i = 0; i < nn; i++) {
		for (int j = 0; j < mesh_dim; j++) {
			cells[q].coeffs[i][j] = 0.;
		}
	}

	// difference between coords of the given cell and its neighbors
	for (int i = 0; i < nn; i++) {
		for (int j = 0; j < mesh_dim; j++) {
			delta_X[i][j] = A7[i].x[j] - X_C.x[j];
		}
	}

	// compute weighting coeffs
	for (int i = 0; i < nn; i++) {
		dist[i] = 0.;
		for (int j = 0; j < mesh_dim; j++)
			dist[i] += delta_X[i][j] * delta_X[i][j];

		dist[i] = sqrt(dist[i]);

		w[i] = 1. / dist[i];

		//cout << "dist[i]" << dist[i] << ", w[i] " << w[i] << endl;
	}

	for (int i = 0; i < mesh_dim; i++) {
		for (int j = 0; j < mesh_dim; j++) {
			A[i][j] = 0.;
			for (int k = 0; k < nn; k++) {
				A[i][j] += w[k] * delta_X[k][i] * delta_X[k][j];
			}
		}
	}

	//cout << "1.A" << endl;
	//cout << A[0][0] << " " << A[0][1] << endl;
	//cout << A[1][0] << " " << A[1][1] << endl;

	inversion(A, mesh_dim);

	//cout << "2.A" << endl;
	//cout << A[0][0] << " " << A[0][1] << endl;
	//cout << A[1][0] << " " << A[1][1] << endl;

	for (int i = 0; i < nn; i++) {

		double* g = new double[mesh_dim];

		for (int j = 0; j < mesh_dim; j++) {
			g[j] = 0.;
			for (int k = 0; k < mesh_dim; k++) {
				g[j] += A[j][k] * delta_X[i][k];
			}

			cells[q].coeffs[i][j] = w[i] * g[j];
		}
		/*		cout << " q= " << q << " i= " << i << endl;
				cout << cells[q].coeffs[i][0] << " " << cells[q].coeffs[i][1] << endl*/;
	}
	//if(q=
}

void Mesh::CellCoeffs()
{
	for (int q = 0; q < nCells; q++) {
		int nn = cells[q].num_f;

		//cells[q].Yw = 1.e3;

		Pnt X_C;
		X_C.x[0] = cells[q].x;
		X_C.x[1] = cells[q].y;

		//cout << "q= " << q << endl;
		//cout << X_C.x[0] << " " << X_C.x[1] << endl;

		Pnt* A7 = new Pnt[nn];
		for (int i = 0; i < nn; i++) 
		{
			//cout << "   faces= " << cells[q].fid[i] << endl;
			int f = cells[q].fid[i];

			int ir = faces[f].ir;
			int il = faces[f].il;

			//cout << "   ir= " << ir << "   il= " << il << endl;

			if (il == q)
			{
				A7[i].x[0] = cells[ir].x;
				A7[i].x[1] = cells[ir].y;
			}
			else if (ir == q)
			{
				if (il >= 0)
				{
					A7[i].x[0] = cells[il].x;
					A7[i].x[1] = cells[il].y;
				}
				else
				{
					A7[i].x[0] = faces[f].x;
					A7[i].x[1] = faces[f].y;
				}
			}
			else
			{
				cout << "Something wrong!" << endl;
				exit(13);
			}

			//cout << A7[i].x[0] << " " << A7[i].x[1] << endl;

		}	//for (int i = 0; i <nn; i++) {

		grad_coeffs_least_sq(A7, q);

		//cout << cells[q].coeffs[1][0] << " " << cells[q].coeffs[1][1] << endl;

		//if (q == 2) exit(34);

	}		// for (int q = 0; q < total_cells; q++) {

}

void Mesh::write_OF_faces(string f)
{
	int* f_vert;
	int total_faces;
	string format_line(20, ' ');

	ofstream record(f, ios::out);
	if (record) {

		cout << "total_cells = " << nCells << endl;

		int nFcs;
		if (mesh_dim == 3)
		{
			nFcs = nFaces;
			total_faces = nFaces;
		}
		else
		{
			nFcs = nFaces + 2 * nCells;
			total_faces = nFaces + 2 * nCells;
		}

		string S(9, ' ');

		record << "FoamFile" << endl;
		record << "{" << endl;
		record << "version" << S << "2.0;" << endl;

		record << "format" << S << "ascii;" << endl;
		record << "class" << S << "faceList;" << endl;
		record << "location" << S << "constant/polyMesh;" << endl;
		record << "object" << S << "faces;" << endl;
		record << "}" << endl;
		record << "//***************************************************" << endl;
		record << setw(10) << total_faces << endl;
		record << "(" << endl;

		if (mesh_dim == 3) {
			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		}
		else {

			for (int i = 0; i < nFaces; i++) {
				int f_vert[4];
				f_vert[0] = faces[i].v[0]; // -1;
				f_vert[1] = faces[i].v[1]; // -1;
				f_vert[2] = faces[i].v[1] + nNodes; // -1;
				f_vert[3] = faces[i].v[0] + nNodes; // -1;
				record << "4(";
				for (int m = 0; m < 4; m++)
					record << setw(10) << f_vert[m];
				record << ")" << endl;
			}

			for (int i = 0; i < nCells; i++) {
				int num_f = cells[i].num_f;

				//cout << i << ",  num_f= " << num_f << endl;
				////exit(65);

				int* f_vert = new int[num_f];
				for (int m = 0; m < num_f; m++) {
					f_vert[m] = cells[i].nodes[m];
				}
				record << num_f << "(";
				for (int m = 0; m < num_f; m++) {
					int m1 = num_f - 1 - m;
					record << setw(10) << f_vert[m1];
				}
				record << ")" << endl;
				delete[] f_vert;
			}

			for (int i = 0; i < nCells; i++) {
				int num_f = cells[i].num_f;

				//cout << i << ",  num_f= " << num_f << endl;
				////exit(65);

				int* f_vert = new int[num_f];
				for (int m = 0; m < num_f; m++) {
					f_vert[m] = cells[i].nodes[m] + nNodes;
				}
				record << num_f << "(";
				for (int m = 0; m < num_f; m++)
					record << setw(10) << f_vert[m];

				record << ")" << endl;
				delete[] f_vert;
			}

		}
		record << ")" << endl;
		record << "//***************************************************" << endl;
	}


	record.close();
}
