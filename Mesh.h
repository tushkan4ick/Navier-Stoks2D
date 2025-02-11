#pragma once

#include <iostream>
#include "fstream"
#include "string"
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

extern const double R; //

struct Cell
{
	double x, y;
	double S;
	vector <int> fid;
	int num_f;
	double dl;		// характерный размер €чейки
	int* nodes;
	//int* faces;

	double** coeffs;	//  коэффициенты дл€ градиента :   
	// элемент 1 - номер точки = 1, ..., n; скорее всего : (num_f)
	// элемент 2 - номер координаты (2 элемента)

	double du[2], dv[2], dh[2];
};

struct Node {		// n
	double x, y;
	//int i, j;
};

struct Face {       // n
	int v[2];		// Face vertices(nodes)
	//int ir0, jr0;		// Right cell indexes
	//int il0, jl0;		// Left cell indexes
	int ir, il;
	int zone_id;		//id зоны, к которой принадлежит грань
	int b_id;		//id  границы, к которой принадлежит грань

	double s;		// face length
	double x, y;		// face center

	double nx, ny;
};

struct Zone {   // n
	string zone_name;
	string zone_type;
	int tid;
	int zone_id;
	int if1;
	int if2;
	int b_id;
};

struct Boundary {
	int tp_id;
	int stp_id;
	int nVals;
	double* values;
	int zone_id;
};

struct Gas
{
	double ro, u, v, E, u_mag;
	double p, T, e, h, p0;
	double Cp, gam, a, alfa;

	double jp, jm;

	double U[4], U1[4], dU[4];

	double gr_u[2], gr_v[2], gr_h[2];

	double mu, Pr;
};

struct Pnt {
	double x[2];
};

double sq(double a);

void Split(string str, string* (&parts), int n, int& i);

class Mesh
{
private:

protected:

public:

	int nNodes;	// number of nodes in a mesh
	Node* nodes;
	int nFaces;	// number of faces in a mesh
	/*int nFacesInternal, nFacesLeft, nFacesBottom;*/
	Face* faces;
	int nCells;	// number of cells in a mesh
	Cell* cells;

	int nZones;
	vector <Zone> zones;
	//Zone* zones;

	int mesh_dim = 2;

	int nBounds;
	Boundary* bounds;

	Mesh();
	Mesh(string);

	void ReadMesh(string);
	void CellGeom();
	void SortNodes(int k, string d, int* (&point_ind));
	void FaceGeom();

	void ReadBounds();

	void write_OF();
	void write_OF_faces(string filename);
	void write_OF_boundary(string filename);
	void write_OF_points(string filename);

	void write_OF_owner(string filename);
	void write_OF_neighbour(string filename);

	void FaceNormals();

	void grad_coeffs_least_sq(Pnt* A7, int q);
	void CellCoeffs();

};

void WriteResults(Mesh mesh, double* f);
void WriteResults(Mesh mesh, Gas* g);
void WriteResults(Mesh mesh, double* f, string par);
void rec_scalar(string par, double* f, string units, string work_dir, int nCells, Mesh mesh);

void CreateGas(Mesh mesh, Gas* (&gasb), Gas* (&g));

void Matrix_Matrix(double G[4][4], double S[4][4], int Nm, double R[4][4]);
void Matrix_Diag(double G[4][4], double L[4], int Nm, double R[4][4]);
void Matrix_Vector(double G[4][4], double V[4], int Nm, double R[4]);

void PrintMatrix(double G[4][4]);

void MatrixA(Gas g, double nx, double ny, double Ap[4][4], double Am[4][4]);

void inversion(double** A, int N);

void Gradients(Mesh mesh, Gas* (&g), Gas* (&gasb));

