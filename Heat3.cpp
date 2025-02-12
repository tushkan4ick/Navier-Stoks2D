#include "Mesh.h"
//#include "Functions.h"

const double R = 287.;

int main()
{
    clock_t start_time = clock();

    setlocale(LC_CTYPE, "rus");
    //cout << "Привет!\n";
    
    string file;
    //file = "ICM_test.msh";
    //file = "FFF.1.msh";
    //file = "nozzleRect1.msh";
    file = "sphere.msh";
    //file = "rect3.msh";
    //file = "Head41.msh";
    //file = "SYS-1.1.msh";
    //file = "Andr.msh";
    //file = "exhaust.msh";

    cout << "MeshFile:\t" << file << endl;
    
    Mesh mesh(file);


    //exit(33);

    mesh.ReadBounds();

    auto nFaces = mesh.nFaces;
    auto nNodes = mesh.nNodes;
    auto nCells = mesh.nCells;
    auto nZones = mesh.nZones;

    auto nBounds = mesh.nBounds;

    Gas* gasb = new Gas[nBounds];
    Gas* g = new Gas[nCells];

    CreateGas(mesh, gasb, g);

    int c = 10;
    cout << c << " " << g[c].u << " " << g[c].p << " " << g[c].T
        << " " << g[c].ro << " " << g[c].e
        << " " << g[c].E << endl << endl;

    //exit(15);


    double resmax = 5.e-6;
    double res = 1.;
    int it = 0;
    while (res > resmax)
    {
        it++;
        // переопределение параметров на след. слое
        for (int c = 0; c < nCells; c++)
        {
            for (int i = 0; i < 4; i++)
            {
                g[c].U[i] = g[c].U1[i];
                g[c].dU[i] = 0.;
                //cout << "g[c].U[i]= " << g[c].U[i] << endl;
            }
        }

        // Определение шага dt по минимал. размеру ячейки и CFL
        double CFL = 0.1;
        double dt = 1.e10;
        for (int c = 0; c < nCells; c++)
        {
            double u = abs(g[c].u);
            double dti = CFL * mesh.cells[c].dl / u;

            if (dt > dti) dt = dti;
        }

        Gradients(mesh, g, gasb);
        //exit(1100);

        double tauXX, tauXY, tauYY, qx, qy;
        double du_dx, du_dy, dv_dx, dv_dy, dh_dx, dh_dy;


        for (int k = 0; k < nFaces; k++)
        {
            // Принадлежность к зоне
            int mz = -1;
            for (int m = 0; m < nZones; m++)
            {
                if (mesh.zones[m].if1 <= k && k <= mesh.zones[m].if2)
                {
                    mz = m;
                    break;
                }
            }




            // идентификатор и имя зоны
            string zone_name = mesh.zones[mz].zone_name;
            int zone_id = mesh.zones[mz].zone_id;
            int b_id = mesh.zones[mz].b_id;

            if (b_id == -1) // internal
            {
                //cout << "internal" << endl;
                int ir = mesh.faces[k].ir;
                int il = mesh.faces[k].il;

                Gas gm;
                // var1
                gm.u = 0.5 * (g[ir].u + g[il].u);
                gm.v = 0.5 * (g[ir].v + g[il].v);
                gm.ro = 0.5 * (g[ir].ro + g[il].ro);
                gm.p = 0.5 * (g[ir].p + g[il].p);
                gm.E = 0.5 * (g[ir].E + g[il].E);

                gm.gam = 0.5 * (g[ir].gam + g[il].gam);

                //cout << "gm.u " << gm.u <<  "; gm.ro " << gm.ro << "; gm.p " << gm.p << "; gm.E " << gm.E << endl;

                double nx = mesh.faces[k].nx;
                double ny = mesh.faces[k].ny;
                double Ap[4][4];
                double Am[4][4];

                MatrixA(gm, nx, ny, Ap, Am);

                double F[4], F1[4], F2[4], V[4];
                for (int i = 0; i < 4; i++)
                    V[i] = g[il].U[i];

                Matrix_Vector(Ap, V, 4, F1);
              
                for (int i = 0; i < 4; i++)
                    V[i] = g[ir].U[i];
                Matrix_Vector(Am, V, 4, F2);

                for (int i = 0; i < 4; i++)
                    F[i] = F1[i] + F2[i];

                //cout << "F var1" << endl;
                //for (int i = 0; i < 4; i++)
                //    cout << setw(12) << F[i] << " ";
                //cout << endl;

                ////exit(6);
                //double H = gm.E + gm.p / gm.ro;

                //double un = gm.u * mesh.faces[k].nx + gm.v * mesh.faces[k].ny;

                ////double F[4];
                //F[0] = gm.ro * un;
                //F[1] = gm.ro * un * gm.u + gm.p * mesh.faces[k].nx;
                //F[2] = gm.ro * un * gm.v + gm.p * mesh.faces[k].ny;
                //F[3] = gm.ro * un * H;

                //cout << "F var2" << endl;
                //for (int i = 0; i < 4; i++)
                //    cout << setw(12) << F[i] << " ";
                //cout << endl;

                //if(k==5) exit(6);
                //cout << un << " u= " << gm.u << " nx= " << mesh.faces[k].nx 
                //    << " ny= " << mesh.faces[k].ny << endl;
                //exit(76);
                //////double un_plus = max(un, 0.);
                //////double un_minus = min(un, 0.);

                ////double F = un_plus * U[il] + un_minus * U[ir];
                for (int i = 0; i < 4; i++)
                {
                    g[ir].dU[i] += dt * mesh.faces[k].s / mesh.cells[ir].S * F[i];
                    g[il].dU[i] -= dt * mesh.faces[k].s / mesh.cells[il].S * F[i];

                    //cout << "i = " << i << "\tg[ir].dU[i] = " << g[ir].dU[i] <<
                    //        "\tg[il].dU[{i}] = " << g[il].dU[i] << endl << endl;
                }


            }
            else
            {
                int tp_id = mesh.bounds[b_id].tp_id;
                int stp_id = mesh.bounds[b_id].stp_id;
                int nVals = mesh.bounds[b_id].nVals;
                if (tp_id == 1)
                {
                    if (stp_id == 1)
                    {
                        //int ir = mesh.faces[k].ir;
                        //double dx = mesh.cells[ir].x - mesh.faces[k].x;
                        //double dy = mesh.cells[ir].y - mesh.faces[k].y;
                        //double dl = sqrt(sq(dx) + sq(dy));
                        //double q = -(lam * (T[ir] - TT) / dl);

                        //double Q = q * mesh.faces[k].s;

                        //dT[ir] += dt / (ro * c * mesh.cells[ir].S) * Q;
                        int ir = mesh.faces[k].ir;

                        double pb = g[ir].p;

                        double F[4];
                        F[0] = 0.;
                        F[1] = pb * mesh.faces[k].nx;
                        F[2] = pb * mesh.faces[k].ny;
                        F[3] = 0.;

                        for (int i = 0; i < 4; i++)
                        {
                            g[ir].dU[i] += dt * mesh.faces[k].s / mesh.cells[ir].S * F[i];
                        }

                    }

                    if (stp_id == 2)
                    {
                        //cout << "wall_Q" << endl;
                        int ir = mesh.faces[k].ir;

                        double pb = g[ir].p;

                        double F[4];
                        F[0] = 0.;
                        F[1] = pb * mesh.faces[k].nx;
                        F[2] = pb * mesh.faces[k].ny;
                        F[3] = 0.;

                        for (int i = 0; i < 4; i++)
                        {
                            g[ir].dU[i] += dt * mesh.faces[k].s / mesh.cells[ir].S * F[i];
                        }
                    }

                    if (stp_id == 3)
                    {
                        //int ir = mesh.faces[k].ir;
                        //double dx = mesh.cells[ir].x - mesh.faces[k].x;
                        //double dy = mesh.cells[ir].y - mesh.faces[k].y;
                        //double dl = sqrt(sq(dx) + sq(dy));

                        //double Tf = mesh.bounds[b_id].values[0];
                        //double alf = mesh.bounds[b_id].values[1];
                        //double z = alf * dl / lam;
                        //double Tw = (T[ir] + z * Tf)/(1+z);
                        //double q = -(lam * (T[ir] - Tw) / dl);

                        //double Q = q * mesh.faces[k].s;
                        //
                        //dT[ir] += dt / (ro * c * mesh.cells[ir].S) * Q;
                    }
                }

                if (tp_id == 2)
                {
                    //cout << mz << " " << mesh.zones[mz].zone_name << endl;
                    //cout << b_id << " " << mesh.bounds[b_id].values[0] << endl;
                    //cout << b_id << " " << mesh.bounds[b_id].values[1] << endl;
                    //cout << "inlet" << endl;
                    //exit(10);

                    int ir = mesh.faces[k].ir;
                    double rob = gasb[b_id].ro;
                    double ub = gasb[b_id].u;
                    double vb = gasb[b_id].v;
                    double pb = gasb[b_id].p;
                    double Eb = gasb[b_id].E;

                    double Hb = Eb + pb / rob;

                    double un = ub * mesh.faces[k].nx + vb * mesh.faces[k].ny;

                    double F[4];
                    F[0] = rob * un;
                    F[1] = rob * un * ub + pb * mesh.faces[k].nx;
                    F[2] = rob * un * vb + pb * mesh.faces[k].ny;
                    F[3] = rob * un * Hb;

                    for (int i = 0; i < 4; i++)
                    {
                        g[ir].dU[i] += dt * mesh.faces[k].s / mesh.cells[ir].S * F[i];
                    }

                }

                if (tp_id == 4)
                {
                    //cout << "outlet" << endl;
                    //exit(11);

                    int ir = mesh.faces[k].ir;
                    double rob = g[ir].ro;
                    double ub = g[ir].u;
                    double vb = g[ir].v;
                    double pb = g[ir].p;
                    double Eb = g[ir].E;

                    double Hb = Eb + pb / rob;

                    double un = ub * mesh.faces[k].nx + vb * mesh.faces[k].ny;

                    double F[4];
                    F[0] = rob * un;
                    F[1] = rob * un * ub + pb * mesh.faces[k].nx;
                    F[2] = rob * un * vb + pb * mesh.faces[k].ny;
                    F[3] = rob * un * Hb;

                    for (int i = 0; i < 4; i++)
                    {
                        g[ir].dU[i] += dt * mesh.faces[k].s / mesh.cells[ir].S * F[i];
                    }
                }

                if (tp_id == 3)
                {
                    //cout << "axis" << endl;
                    //exit(12);
                    int ir = mesh.faces[k].ir;

                    double pb = g[ir].p;

                    double F[4];
                    F[0] = 0.;
                    F[1] = pb * mesh.faces[k].nx;
                    F[2] = pb * mesh.faces[k].ny;
                    F[3] = 0.;

                    for (int i = 0; i < 4; i++)
                    {
                        g[ir].dU[i] += dt * mesh.faces[k].s / mesh.cells[ir].S * F[i];
                    }
                }

                if (tp_id == 5)
                {
                    if (stp_id == 1)    //velocity and temperature fixed
                    {

                        int ir = mesh.faces[k].ir;

                        double rob = gasb[b_id].ro;
                        double ub = gasb[b_id].u;
                        double vb = gasb[b_id].v;
                        double pb = gasb[b_id].p;
                        double Eb = gasb[b_id].E;
                        double Tb = mesh.bounds[b_id].values[3];
                        double ub_mag = gasb[b_id].u_mag;
                        double ab = gasb[b_id].a;

                        double roc = g[ir].ro;
                        double pc = g[ir].p;
                        double uc_mag = g[ir].u_mag;
                        double ac = g[ir].a;


                        double aro = 0.5 * (rob * ab + roc * ac);

                        gasb[b_id].p = pc + aro * (ub_mag - uc_mag);
                        gasb[b_id].ro = gasb[b_id].p / R / Tb;

                        pb = gasb[b_id].p;
                        rob = gasb[b_id].ro;

                        cout << "uc_mag= " << uc_mag << endl;

                        double Hb = Eb + pb / rob;

                        double un = ub * mesh.faces[k].nx + vb * mesh.faces[k].ny;

                        double F[4];

                        F[0] = rob * un;
                        F[1] = rob * un * ub + pb * mesh.faces[k].nx;
                        F[2] = rob * un * vb + pb * mesh.faces[k].ny;
                        F[3] = rob * un * Hb;

                        for (int i = 0; i < 4; i++)
                        {
                            g[ir].dU[i] += dt * mesh.faces[k].s / mesh.cells[ir].S * F[i];
                        }

                    }

                    if (stp_id == 2)    //static pressure and temperature fixed
                    {


                    }


                }

                if (tp_id == 6)     //subsonic inlet
                {

                }

            }

        }
        cout << endl;



        int iViscous = 1;
        if (iViscous == 1)
        {
            // Приращения за счет вязких потоков
            for (int k = 0; k < nFaces; k++)
            {
                // Принадлежность к зоне
                int mz = -1;
                for (int m = 0; m < nZones; m++)
                {
                    if (mesh.zones[m].if1 <= k && k <= mesh.zones[m].if2)
                    {
                        mz = m;
                        break;
                    }
                }

                //cout << "mz= " << mz << endl; 
                //exit(6);
                // идентификатор и имя зоны
                string zone_name = mesh.zones[mz].zone_name;
                int zone_id = mesh.zones[mz].zone_id;
                int b_id = mesh.zones[mz].b_id;

                if (b_id == -1) // internal
                {
                    int ir = mesh.faces[k].ir;
                    int il = mesh.faces[k].il;

                    // Определение градиентов на грани
                    // 1) Расстояние от центров ячеек
                    double dr, dl;
                    dr = sqrt(sq(mesh.faces[k].x - mesh.cells[ir].x) +
                        sq(mesh.faces[k].y - mesh.cells[ir].y));
                    dl = sqrt(sq(mesh.faces[k].x - mesh.cells[il].x) +
                        sq(mesh.faces[k].y - mesh.cells[il].y));

                    //cout << dr << " " << dl << endl;
                    //exit(22);
                    // 2) Интерполяция градиентов
                    double zr = dl / (dr + dl);
                    double zl = dr / (dr + dl);
                    du_dx = g[ir].gr_u[0] * zr + g[il].gr_u[0] * zl;
                    du_dy = g[ir].gr_u[1] * zr + g[il].gr_u[1] * zl;
                    dv_dx = g[ir].gr_v[0] * zr + g[il].gr_v[0] * zl;
                    dv_dy = g[ir].gr_v[1] * zr + g[il].gr_v[1] * zl;

                    dh_dx = g[ir].gr_h[0] * zr + g[il].gr_h[0] * zl;
                    dh_dy = g[ir].gr_h[1] * zr + g[il].gr_h[1] * zl;

                    // 3) средние значения переносных свойств и скорости
                    double mu = 0.5 * (g[ir].mu + g[il].mu);
                    double Pr = 0.5 * (g[ir].Pr + g[il].Pr);
                    double u = 0.5 * (g[ir].u + g[il].u);
                    double v = 0.5 * (g[ir].v + g[il].v);

                    // 4) tauXX, tauXY, tauYY, qx, qy;
                    tauXX = 2. / 3. * mu * (2. * du_dx - dv_dy);
                    tauXY = mu * (du_dy + dv_dx);
                    tauYY = 2. / 3. * mu * (2. * dv_dy - du_dx);

                    qx = -mu / Pr * dh_dx;
                    qy = -mu / Pr * dh_dy;

                    // вязкие потоки на грани
                    double F[4];
                    F[0] = 0.;
                    F[1] = tauXX * mesh.faces[k].nx + tauXY * mesh.faces[k].ny;
                    F[2] = tauXY * mesh.faces[k].nx + tauYY * mesh.faces[k].ny;
                    F[3] = (u * tauXX + v * tauXY - qx) * mesh.faces[k].nx +
                           (u * tauXY + v * tauYY - qy) * mesh.faces[k].ny;

                    for (int i = 0; i < 4; i++)
                    {
                        g[ir].dU[i] -= dt * mesh.faces[k].s / mesh.cells[ir].S * F[i];
                        g[il].dU[i] += dt * mesh.faces[k].s / mesh.cells[il].S * F[i];
                    }

                }
                else
                {
                    int tp_id = mesh.bounds[b_id].tp_id;
                    int stp_id = mesh.bounds[b_id].stp_id;
                    int nVals = mesh.bounds[b_id].nVals;
                    if (tp_id == 1)
                    {

                        if (stp_id == 1)        // Проверить
                        {

                            int ir = mesh.faces[k].ir;
                            //cout << "tauYY = " << tauYY << endl;
                            //exit(101);

                            double pb = g[ir].p;

                            double lam = g[ir].Cp * g[ir].mu / g[ir].Pr;

                            double F[4];
                            F[0] = 0.;
                            F[1] = 0.;
                            F[2] = 2. / 3. * g[ir].mu * (2. * g[ir].gr_v[1] - g[ir].gr_u[0]);
                            F[3] = - lam * (g[ir].T - gasb[b_id].T) / (mesh.faces[k].nx + mesh.faces[k].ny);
                            for (int i = 0; i < 4; i++)
                            {
                                g[ir].dU[i] -= dt * mesh.faces[k].s / mesh.cells[ir].S * F[i];
                            }
                        }

                        if (stp_id == 2)        // Проверить
                        {
                            int ir = mesh.faces[k].ir;

                            //double pb = g[ir].p;

                            double F[4];
                            F[0] = 0.;
                            F[1] = 0;
                            F[2] = 2. / 3. * g[ir].mu * (2. * g[ir].gr_v[1] - g[ir].gr_u[0]);
                            F[3] = - 2 * mesh.bounds[b_id].values[0];

                            ////cout << "F[3] = " << F[3] << endl;

                            //for (int i = 0; i < 4; i++)
                            //{
                            //    g[ir].dU[i] -= dt * mesh.faces[k].s / mesh.cells[ir].S * F[i];
                            //}
                        }

                        if (stp_id == 3)
                        {

                        }
                    }

                    if (tp_id == 2)
                    {

                    }

                    if (tp_id == 4)
                    {

                    }

                    if (tp_id == 3)
                    {
                        int ir = mesh.faces[k].ir;

                        //double pb = g[ir].p;

                        double F[4];
                        F[0] = 0.;
                        F[1] = 0;
                        F[2] = 2./3. * g[ir].mu * (2.* g[ir].gr_v[1]- g[ir].gr_u[0]);
                        F[3] = 0.;

                        for (int i = 0; i < 4; i++)
                        {
                            g[ir].dU[i] -= dt * mesh.faces[k].s / mesh.cells[ir].S * F[i];
                        }
                    }

                }

            }
        }

        for (int c = 0; c < nCells; c++)
        {
            for (int i = 0; i < 4; i++)
                g[c].U1[i] = g[c].U[i] + g[c].dU[i];

            //cout << "dU[0]= " << g[c].dU[0] << "\tdU[1]= " << g[c].dU[1] << "\tdU[2]= " << g[c].dU[2] << "\tdU[3]= " << g[c].dU[3] << endl;

            g[c].ro = g[c].U1[0];
            g[c].u = g[c].U1[1] / g[c].ro;
            g[c].v = g[c].U1[2] / g[c].ro;
            g[c].E = g[c].U1[3] / g[c].ro;
            g[c].u_mag = sqrt(sq(g[c].u) + sq(g[c].v));

            g[c].e = g[c].E - 0.5 * (sq(g[c].u) + sq(g[c].v));

            g[c].T = g[c].e * g[c].gam / g[c].Cp;

            if (g[c].T > 6000.)
            {
                g[c].T = 6000.;
                //cout << "in cell N " << c << " T = 6000" << endl;
            }

            g[c].h = g[c].T * g[c].Cp;

            g[c].p = g[c].ro * R * g[c].T;
            g[c].a = sqrt(g[c].gam * R * g[c].T);

            //nasa format
            if (g[c].T < 1000)
            {
                g[c].Cp = 2898903. * pow(g[c].T, -2) - 56496.26 * pow(g[c].T, -1) + 1437.799 - 
                          1.653609 * g[c].T + 0.0030632254 * pow(g[c].T, 2) - 2.279138e-6 * pow(g[c].T, 3) + 
                          6.272365e-10 * pow(g[c].T, 4);
            }
            else
            {
                g[c].Cp = 6.932494e7 * pow(g[c].T, -2) - 361053.2 * pow(g[c].T, -1) + 1476.665 -
                          0.06138349 * g[c].T + 2.027963e-5 * pow(g[c].T, 2) - 3.075525e-9 * pow(g[c].T, 3) + 
                          1.888054e-13 * pow(g[c].T, 4);
            }
            //cout << "Cp = " << g[c].Cp << endl;
            //exit(11111);

        }
        

        double r_density = 0.;
        double r_velo_X = 0.;
        double r_velo_Y = 0.;
        double r_E = 0.;
        res = 0.;

        for (int c = 0; c < nCells; c++)
        {
            double den = abs(g[c].dU[0] / g[c].U1[0]);
            if (r_density < den) r_density = den;
            double velo_X = abs(g[c].dU[1] / g[c].U1[1]);
            if (r_velo_X < velo_X) r_velo_X = velo_X;
            double velo_Y = abs(g[c].dU[2] / g[c].U1[2]);
            if (r_velo_Y < velo_Y) r_velo_Y = velo_Y;
            double E = abs(g[c].dU[3] / g[c].U1[3]);
            if (r_E < E) r_E = E;
            double r = max({ r_density, r_velo_X, r_velo_Y, r_E });
            if (res < r) res = r;
            //double r = abs(g[c].dU[3] / g[c].U1[3]);
            //if (res < r) res = r;
        }


        if (it % 10 == 0)
        {
            cout << "iteration = " << it << endl;
            cout << "r_density = " << r_density << "; r_velocityX = " << r_velo_X << 
                    "; r_velocityY = " << r_velo_Y << "; r_E = " << r_E << endl;
            //cout << "P= " << max({ g->p }) << endl;
            //cout << "ro = " << max({ g->ro }) << ";\tu = " << max({ g->u }) << ";\tv = " << max({ g->v }) << 
            //    ";\tE = " << max({ g->E }) << endl;
        }
        
        if (it > 50000) break;
    }

    cout << "it= " << it << " res= " << res << endl;

    clock_t end_time = clock();
    double seconds = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    cout << "time= " << seconds << " seconds" << endl;

    //WriteResults(mesh, T);

    WriteResults(mesh, g);

    //WriteResults(mesh, U1, "U");

    mesh.write_OF();

}
