//#include "Functions.h"

#include "Mesh.h"
#include <iomanip>


void WriteResults(Mesh mesh, double* f)
{
	string units, par;

	string work_dir = "test//1";

	int nCells = mesh.nCells;

	// Temperature
	units = "0 0 0 1 0 0 0 ";
	par = "T";

	rec_scalar(par, f, units, work_dir, nCells, mesh);
}

void WriteResults(Mesh mesh, Gas* g)
{
	string units, par;

	string work_dir = "test//1";

	int nCells = mesh.nCells;
	double* f = new double[nCells];

	// Temperature
	units = "0 0 0 1 0 0 0 ";
	par = "T";
	for (int c = 0; c < nCells; c++)
		f[c] = g[c].T;

	rec_scalar(par, f, units, work_dir, nCells, mesh);

	// u
	units = "0 0 0 1 0 0 0 ";
	par = "u";
	for (int c = 0; c < nCells; c++)
		f[c] = g[c].u;

	rec_scalar(par, f, units, work_dir, nCells, mesh);

	// v
	units = "0 0 0 1 0 0 0 ";
	par = "v";
	for (int c = 0; c < nCells; c++)
		f[c] = g[c].v;

	rec_scalar(par, f, units, work_dir, nCells, mesh);

	// ro
	units = "0 0 0 1 0 0 0 ";
	par = "ro";
	for (int c = 0; c < nCells; c++)
		f[c] = g[c].ro;

	rec_scalar(par, f, units, work_dir, nCells, mesh);

	// p
	units = "0 0 0 1 0 0 0 ";
	par = "p";
	for (int c = 0; c < nCells; c++)
		f[c] = g[c].p;

	rec_scalar(par, f, units, work_dir, nCells, mesh);

	// du_dy
	units = "0 0 0 1 0 0 0 ";
	par = "du_dy";
	for (int c = 0; c < nCells; c++)
		f[c] = g[c].gr_u[1];

	rec_scalar(par, f, units, work_dir, nCells, mesh);

	// du_dx
	units = "0 0 0 1 0 0 0 ";
	par = "du_dx";
	for (int c = 0; c < nCells; c++)
		f[c] = g[c].gr_u[0];

	rec_scalar(par, f, units, work_dir, nCells, mesh);

	// dv_dx
	units = "0 0 0 1 0 0 0 ";
	par = "dv_dx";
	for (int c = 0; c < nCells; c++)
		f[c] = g[c].gr_v[0];

	rec_scalar(par, f, units, work_dir, nCells, mesh);

	// dv_dy
	units = "0 0 0 1 0 0 0 ";
	par = "dv_dy";
	for (int c = 0; c < nCells; c++)
		f[c] = g[c].gr_v[1];

	rec_scalar(par, f, units, work_dir, nCells, mesh);

	// dh_dx
	units = "0 0 0 1 0 0 0 ";
	par = "dh_dx";
	for (int c = 0; c < nCells; c++)
		f[c] = g[c].gr_h[0];

	rec_scalar(par, f, units, work_dir, nCells, mesh);

	// dh_dy
	units = "0 0 0 1 0 0 0 ";
	par = "dh_dy";
	for (int c = 0; c < nCells; c++)
		f[c] = g[c].gr_h[1];

	rec_scalar(par, f, units, work_dir, nCells, mesh);


}

void WriteResults(Mesh mesh, double* f, string par)
{
	string units;

	string work_dir = "test//1";

	int nCells = mesh.nCells;

	units = "0 0 0 1 0 0 0 ";

	rec_scalar(par, f, units, work_dir, nCells, mesh);
}

void rec_scalar(string par, double* f, string units, string work_dir, int nCells, Mesh mesh)
{
	string filename = work_dir + "//" + par;

	ofstream record(filename, ios::out);
	if (record) 
	{

		string S(9, ' ');

		record << "FoamFile" << endl;
		record << "{" << endl;
		record << "version" << S << "2.0;" << endl;
		record << "format" << S << "ascii;" << endl;
		record << "class" << S << "volScalarField;" << endl;
		record << "location" << S << "\"1\";" << endl;
		record << "object" << S << par << ";" << endl;
		record << "}" << endl;
		record << "//***************************************************" << endl;


		record << "dimensions" << S << setw(10) << "[" << units << "];" << endl;
		record << " internalField  nonuniform List<scalar>" << endl;
		record << setw(12) << nCells << endl;
		record << "(" << endl;


		for (int i = 0; i < nCells; i++) {
			record << f[i] << endl;
		}

		record << ")" << endl;
		record << ";" << endl;

//****************************************************
		int nZones = mesh.nZones;
		//int num_of_boundaries = 0;
		//for (int i = 0; i < nZones; i++)
		//{
		//	if (mesh.zones[i].if1 > 0) num_of_boundaries++;
		//}

		int nBounds = mesh.nBounds;
		/*cout << "nBounds= " << nBounds << "num_of_boundaries= " << num_of_boundaries << endl;*/

		record << "boundaryField" << endl;
		record << "{" << endl;

		for (int i = 0; i < nZones; i++)
		{
			if (mesh.zones[i].b_id >= 0)
			{
				record << mesh.zones[i].zone_name << endl;
				record << "{" << endl;
				record << "type" << S << "zeroGradient;" << endl;
				record << "}" << endl;
				record << " " << endl;
			}
		}

		record << "Interior_artificial_1" << endl;
		record << "{" << endl;
		record << "type" << S << "empty;" << endl;
		record << "}" << endl;
		record << " " << endl;

		record << "Interior_artificial_2" << endl;
		record << "{" << endl;
		record << "type" << S << "empty;" << endl;
		record << "}" << endl;
		record << " " << endl;

		record << S << "defaultFaces" << endl;
		record << S << "{" << endl;
		record << S << S << "type" << S << "empty;" << endl;
		record << S << "}" << endl;

		record << "}" << endl;
		record << "//***************************************************" << endl;

		record.close();
	}
	
}

void CreateGas(Mesh mesh, Gas* (&gasb), Gas* (&g))
{

	cout << endl << "           Create FLUID" << endl;

	clock_t start_time = clock();


	// Параметры на границах
	auto nBounds = mesh.nBounds;
	/*Gas* gasb = new Gas[nBounds];*/

	int bc = 0;

	for (int b = 0; b < nBounds; b++)
	{


		//cout << "create bound fluid" << endl;

		if (mesh.bounds[b].tp_id == 2)	// Supersonic Inlet
		{
			if (mesh.bounds[b].zone_id == 12) bc = b;
			//cout << "bc = " << bc << endl;
			gasb[b].u = mesh.bounds[b].values[0];
			gasb[b].v = mesh.bounds[b].values[1];
			gasb[b].p = mesh.bounds[b].values[2];
			gasb[b].T = mesh.bounds[b].values[3];
			cout << "gasb[b].u= " << gasb[b].u << "; gasb[b].p= " << gasb[b].p << "; gasb[b].T= " << gasb[b].T << endl;

			gasb[b].Cp = 1010.;
			gasb[b].gam = 1.4;
			
			gasb[b].mu = 4.e-5;
			gasb[b].Pr = 0.7;

			gasb[b].a = sqrt(gasb[b].gam * gasb[b].T * R);
			gasb[b].ro = gasb[b].p / (R * gasb[b].T);
			gasb[b].e = gasb[b].Cp / gasb[b].gam * gasb[b].T;

			gasb[b].h = gasb[b].Cp * gasb[b].T;

			gasb[b].E = gasb[b].e + 0.5 * (sq(gasb[b].u) + sq(gasb[b].v));

			//cout << b << "; gasb[b].u = " << gasb[b].u << "; gasb[b].p = " << gasb[b].p << "; gasb[b].T = " << gasb[b].T
			//	<< "; gasb[b].ro = " << gasb[b].ro << "; gasb[b].e = " << gasb[b].e
			//	<< "; gasb[b].E = " << gasb[b].E << endl;

			gasb[b].U[0] = gasb[b].ro;
			//cout << "gasb[b].U[0] = " << gasb[b].U[0] << endl;
			gasb[b].U[1] = gasb[b].ro * gasb[b].u;
			//cout << "gasb[b].U[1] = " << gasb[b].U[1] << endl;
			gasb[b].U[2] = gasb[b].ro * gasb[b].v;
			gasb[b].U[3] = gasb[b].ro * gasb[b].E;

			for (int i = 0; i < 4; i++)
			{
				gasb[b].U1[i] = gasb[b].U[i];
				gasb[b].dU[0] = 0;
				//cout << "gasb[b].U1[i] = " << gasb[b].U1[i] << endl;
			}
		}
		if (mesh.bounds[b].tp_id == 5)		//subsonic inlet
		{
			if (mesh.bounds[b].zone_id == 12) bc = b;	//zone_id for initialization

			//cout << "subsonic inlet" << endl;

			if (mesh.bounds[b].stp_id == 1)	//velocity and temperature
			{
				//cout << "bc = " << bc << endl;
				gasb[b].u = mesh.bounds[b].values[0];
				gasb[b].v = mesh.bounds[b].values[1];
				gasb[b].p = mesh.bounds[b].values[2];		//только для начальных условий
				gasb[b].T = mesh.bounds[b].values[3];
				//cout << "gasb[b].u= " << gasb[b].u << "; gasb[b].p= " << gasb[b].p << "; gasb[b].T= " << gasb[b].T << endl;


				gasb[b].Cp = 1010.;
				gasb[b].gam = 1.4;

				gasb[b].mu = 4.e-5;
				gasb[b].Pr = 0.7;

				gasb[b].u_mag = sqrt(sq(gasb[b].u) + sq(gasb[b].v));

				gasb[b].ro = gasb[b].p / (R * gasb[b].T);
				gasb[b].e = gasb[b].Cp / gasb[b].gam * gasb[b].T;

				gasb[b].h = gasb[b].Cp * gasb[b].T;

				gasb[b].E = gasb[b].e + 0.5 * (sq(gasb[b].u) + sq(gasb[b].v));

				//cout << b << "; gasb[b].u = " << gasb[b].u << "; gasb[b].p = " << gasb[b].p << "; gasb[b].T = " << gasb[b].T
				//	<< "; gasb[b].ro = " << gasb[b].ro << "; gasb[b].e = " << gasb[b].e
				//	<< "; gasb[b].E = " << gasb[b].E << endl;

				gasb[b].U[0] = gasb[b].ro;
				//cout << "gasb[b].U[0] = " << gasb[b].U[0] << endl;
				gasb[b].U[1] = gasb[b].ro * gasb[b].u;
				//cout << "gasb[b].U[1] = " << gasb[b].U[1] << endl;
				gasb[b].U[2] = gasb[b].ro * gasb[b].v;
				gasb[b].U[3] = gasb[b].ro * gasb[b].E;

				for (int i = 0; i < 4; i++)
				{
					gasb[b].U1[i] = gasb[b].U[i];
					gasb[b].dU[0] = 0;
					//cout << "gasb[b].U1[i] = " << gasb[b].U1[i] << endl;
				}
			}
			if (mesh.bounds[b].stp_id == 2)		//total pressure, pressure and temperature
			{
				//cout << "bc = " << bc << endl;
				gasb[b].p0 = mesh.bounds[b].values[0];
				gasb[b].alfa = mesh.bounds[b].values[1];
				gasb[b].p = mesh.bounds[b].values[2];		//только для начальных условий
				gasb[b].T = mesh.bounds[b].values[3];
				//cout << "gasb[b].u= " << gasb[b].u << "; gasb[b].p= " << gasb[b].p << "; gasb[b].T= " << gasb[b].T << endl;


				gasb[b].Cp = 1010.;
				gasb[b].gam = 1.4;

				gasb[b].mu = 4.e-5;
				gasb[b].Pr = 0.7;

				double g = gasb[b].gam;
				double p0 = gasb[b].p0;
				double p = gasb[b].p;
				double T = gasb[b].T;

				gasb[b].u = sqrt(2 * g * R * T / (g - 1.) *
								(pow(p0 / p, (g - 1.) / g) - 1.)) * cos(gasb[b].alfa);
				gasb[b].v = sqrt(2 * g * R * T / (g - 1.) *
								(pow(p0 / p, (g - 1.) / g) - 1.)) * sin(gasb[b].alfa);

				gasb[b].ro = gasb[b].p / (R * gasb[b].T);
				gasb[b].e = gasb[b].Cp / gasb[b].gam * gasb[b].T;

				gasb[b].h = gasb[b].Cp * gasb[b].T;

				gasb[b].E = gasb[b].e + 0.5 * (sq(gasb[b].u) + sq(gasb[b].v));

				//cout << b << "; gasb[b].u = " << gasb[b].u << "; gasb[b].p = " << gasb[b].p << "; gasb[b].T = " << gasb[b].T
				//	<< "; gasb[b].ro = " << gasb[b].ro << "; gasb[b].e = " << gasb[b].e
				//	<< "; gasb[b].E = " << gasb[b].E << endl;

				gasb[b].U[0] = gasb[b].ro;
				//cout << "gasb[b].U[0] = " << gasb[b].U[0] << endl;
				gasb[b].U[1] = gasb[b].ro * gasb[b].u;
				//cout << "gasb[b].U[1] = " << gasb[b].U[1] << endl;
				gasb[b].U[2] = gasb[b].ro * gasb[b].v;
				gasb[b].U[3] = gasb[b].ro * gasb[b].E;

				for (int i = 0; i < 4; i++)
				{
					gasb[b].U1[i] = gasb[b].U[i];
					gasb[b].dU[0] = 0;
					//cout << "gasb[b].U1[i] = " << gasb[b].U1[i] << endl;
				}
			}

		}

	}

	// Начальные параметры в ячейках
	auto nCells = mesh.nCells;
	//Граница, по которой задаются начальные условия
	//int b = 0;   // "rect3.msh"; "Head41.msh";
	//int b = 0;
	//int b = bc;

	cout << "b = " << bc << endl;

	for (int c = 0; c < nCells; c++)
	{
		for (int i = 0; i < 4; i++)
		{
			g[c].U[i] = gasb[bc].U1[i];
			//cout << "g[c].U[i] = " << g[c].U[i] << endl;
		}

		g[c] = gasb[bc];
		//cout << "g[c] " << g[c].u << endl;
	}

	clock_t end_time = clock();
	double seconds = (double)(end_time - start_time) / CLOCKS_PER_SEC;

	cout << "           FLUID created at " << seconds << " sec" << endl << endl;

}

void Matrix_Matrix(double G[4][4], double S[4][4], int Nm, double R[4][4])
{
	for (int k = 0; k < Nm; k++) {
		for (int l = 0; l < Nm; l++) {
			R[k][l] = 0.;
			for (int m = 0; m < Nm; m++)
				R[k][l] += G[k][m] * S[m][l];
		}
	}
}

void Matrix_Diag(double G[4][4], double L[4], int Nm, double R[4][4])
{
	for (int k = 0; k < Nm; k++)
		for (int m = 0; m < Nm; m++)
			R[k][m] = G[k][m] * L[m];
}

void Matrix_Vector(double G[4][4], double V[4], int Nm, double R[4])
{
	for (int k = 0; k < Nm; k++) {
		R[k] = 0;
		for (int l = 0; l < Nm; l++)
			R[k] += G[k][l] * V[l];
	}
}

void PrintMatrix(double G[4][4])
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << setw(12) << G[i][j] << " ";

		cout << endl;
	}
}

void MatrixA(Gas g, double nx, double ny, double Ap[4][4], double Am[4][4])
{

	double p = g.p;
	double ro = g.ro;
	double Et = g.E;
	double Ht = Et + p / ro;

	double gam = g.gam;

	double beta = gam - 1;

	double u = g.u;
	double v = g.v;

	double alfa = 0.5 * (sq(u) + sq(v));

	double h = Ht - alfa;
	//касательный вектор
	double ly = nx;
	double lx = -ny;

	double un = u * nx + v * ny;

	//double A[4][4];

	//A[0][0] = 0.;
	//A[0][1] = nx;
	//A[0][2] = ny;                   //v - кoмпонента скорости по оси y
	//A[0][3] = 0.;

	//A[1][0] = alfa * beta * nx - un * u;
	//A[1][1] = u * nx * (1. - beta) + un;
	//A[1][2] = -beta * v * nx + u * ny;
	//A[1][3] = beta * nx;

	//A[2][0] = alfa * beta * ny - un * v;
	//A[2][1] = -beta * u * ny + v * nx;
	//A[2][2] = v * ny * (1. - beta) + un;
	//A[2][3] = beta * ny;

	//A[3][0] = un * (alfa * beta - Ht);
	//A[3][1] = -beta * u * un + Ht * nx;
	//A[3][2] = -beta * v * un + Ht * ny;
	//A[3][3] = gam * un;

	//cout << "Matrix A" << endl;
	//PrintMatrix(A);

	// S^-1
	double S_[4][4];
	double a2 = gam * p / ro;
	double a = sqrt(a2);
	double V = u * lx + v * ly;

	S_[0][0] = 1. / a2;
	S_[0][1] = 0;
	S_[0][2] = 1. / (2. * a2);
	S_[0][3] = 1. / (2. * a2);

	S_[1][0] = u / a2;
	S_[1][1] = lx;
	S_[1][2] = (u + a * nx) / (2. * a2);
	S_[1][3] = (u - a * nx) / (2. * a2);

	S_[2][0] = v / a2;
	S_[2][1] = ly;
	S_[2][2] = (v + a * ny) / (2. * a2);
	S_[2][3] = (v - a * ny) / (2. * a2);

	S_[3][0] = alfa / a2;
	S_[3][1] = V;
	S_[3][2] = (Ht + a * un) / (2. * a2);
	S_[3][3] = (Ht - a * un) / (2. * a2);

	//cout << "Matrix S^-1" << endl;
	//PrintMatrix(S_);

	// S
	double S[4][4];
	S[0][0] = a2 - alfa * beta;           //a, alfa, beta
	S[0][1] = beta * u;                   //u - кoмпонента скорости по оси x
	S[0][2] = beta * v;                  //v - кoмпонента скорости по оси y
	S[0][3] = -beta;

	S[1][0] = -V;
	S[1][1] = lx;
	S[1][2] = ly;
	S[1][3] = 0;

	S[2][0] = alfa * beta - un * a;
	S[2][1] = a * nx - beta * u;
	S[2][2] = a * ny - beta * v;
	S[2][3] = beta;

	S[3][0] = alfa * beta + un * a;
	S[3][1] = -a * nx - beta * u;
	S[3][2] = -a * ny - beta * v;
	S[3][3] = beta;

	double L[4];
	L[0] = un;
	L[1] = un;
	L[2] = un + a;
	L[3] = un - a;

	//double T[4][4], A1[4][4];
	//Matrix_Diag(S_, L, 4, T);
	//Matrix_Matrix(T, S, 4, A1);
	//cout << "Matrix A1" << endl;
	//PrintMatrix(A1);

	double z_static = 0.; // 0.25;
	double epsilon = z_static * (abs(un) + a);

	double Lp[4], Lm[4];
	for (int i = 0; i < 4; i++)
	{
		//Lp[i] = max(L[i], 0.);
		//Lm[i] = min(L[i], 0.);

		Lp[i] = 0.5 * (L[i] + sqrt(sq(L[i]) + sq(epsilon)));
		Lm[i] = 0.5 * (L[i] - sqrt(sq(L[i]) + sq(epsilon)));
	}

	double T[4][4];
	Matrix_Diag(S_, Lp, 4, T);
	Matrix_Matrix(T, S, 4, Ap);
	//cout << "Matrix Ap" << endl;
	//PrintMatrix(Ap);

	Matrix_Diag(S_, Lm, 4, T);
	Matrix_Matrix(T, S, 4, Am);
	//cout << "Matrix Am" << endl;
	//PrintMatrix(Am);

}

void inversion(double** A, int N)
{
	double temp;

	double** E = new double* [N];

	for (int i = 0; i < N; i++)
		E[i] = new double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			E[i][j] = 0.0;

			if (i == j)
				E[i][j] = 1.0;
		}

	for (int k = 0; k < N; k++)
	{
		temp = A[k][k];

		for (int j = 0; j < N; j++)
		{
			A[k][j] /= temp;
			E[k][j] /= temp;
		}

		for (int i = k + 1; i < N; i++)
		{
			temp = A[i][k];

			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int k = N - 1; k > 0; k--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			temp = A[i][k];

			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			A[i][j] = E[i][j];

	for (int i = 0; i < N; i++)
		delete[] E[i];

	delete[] E;
}

void Gradients(Mesh mesh, Gas* (&g), Gas* (&gasb))
{
	int nCells = mesh.nCells;
	int Neq = 4;
	int nZones = mesh.nZones;

	for (int i = 0; i < nCells; i++)
	{
		// Газодинамические параметры в данной ячейке
		double u = g[i].u;
		double v = g[i].v;
		double h = g[i].h;
		
		//cout << u << " " << v << " " << h << endl;

		// nn - число соседей
		int nn = mesh.cells[i].num_f;

		double gr_u[2], gr_v[2], gr_h[2];

		for (int m = 0; m < 2; m++)
		{
			gr_u[m] = 0.;
			gr_v[m] = 0.;
			gr_h[m] = 0.;
		}

		for (int k = 0; k < nn; k++)
		{
			double coefs[2];
			coefs[0] = mesh.cells[i].coeffs[k][0];
			coefs[1] = mesh.cells[i].coeffs[k][1];

			//cout << coefs[0] << " " << coefs[1] << endl;
			//exit(10);

			int f = mesh.cells[i].fid[k];


			int ir = mesh.faces[f].ir;
			int il = mesh.faces[f].il;


			double du, dv, dh;
			// Внутренние грани
			if (il >= 0)
			{
				int ii;
				if (ir == i)
				{
					ii = il;
				}
				else
				{
					ii = ir;
				}
				du = g[ii].u - u;
				dv = g[ii].v - v;
				dh = g[ii].h - h;
			}
			else     //пограничье
			{
				//cout << "Boundaries" << endl;
				//int tp_id = 0;
				//int stp_id = 0;
				//int b = 0;
				//for (b; b < mesh.nBounds; b++)
				//{
				//	//cout << mesh.faces[f].zone_id << endl;
				//	if (mesh.faces[f].zone_id == mesh.bounds[b].zone_id)
				//	{
				//		tp_id = mesh.bounds[b].tp_id;
				//		stp_id = mesh.bounds[b].stp_id;
				//		break;
				//	}
				//}
				// упрощенный вариант
				du = 0;
				dv = 0;
				dh = 0;

				//if (tp_id == 1)
				//{
				//	double xcc = mesh.cells[i].x;
				//	double ycc = mesh.cells[i].y;

				//	double xfc = mesh.faces[f].x;
				//	double yfc = mesh.faces[f].y;

				//	double dr = sqrt(sq(xcc - xfc) + sq(ycc - yfc));

				//	if (stp_id == 1)
				//	{
				//		cout << "gasb[b].Cp = " << gasb[b].Cp << "\tgasb[b].T = " << gasb[b].T <<
				//			    "\th = " << h << endl;
				//		//cout << "b = " << b << endl;
				//		du = 0;
				//		dv = 0;
				//		dh = (gasb[b].Cp * gasb[b].T - h) * dr / mesh.faces[b].s;		// Где брать параметры Cp и Т? Нужно привести dr к безрамерному виду
				//		//cout << fzid << endl;
				//		//cout << "dh= " << dh << endl;
				//	}

				//	if (stp_id == 2)
				//	{

				//	}

				//	if (stp_id == 3)
				//	{

				//	}
				//}
			}
			//cout << "du= " << du << endl;
			for (int m = 0; m < 2; m++)
			{
				gr_u[m] += coefs[m] * du;
				gr_v[m] += coefs[m] * dv;
				gr_h[m] += coefs[m] * dh;
			}
			//cout << " gr_u[0]= " << gr_u[0] << " gr_u[1]= " << gr_u[1] << endl;
		}
		//cout << "2 gr_u[0]= " << gr_u[0] << " gr_u[1]= " << gr_u[1] << endl;
		//exit(66);

		for (int m = 0; m < 2; m++)
		{
			g[i].gr_u[m] = gr_u[m];
			g[i].gr_v[m] = gr_v[m];
			g[i].gr_h[m] = gr_h[m];
		}
		//cout << " gr_u[0]= " << g[i].gr_u[0] << " gr_u[1]= " << g[i].gr_u[1] << endl;

		//exit(66);
	}

}

