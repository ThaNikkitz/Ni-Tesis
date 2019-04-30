/*  
  Copyright (C) 2013 by Christopher Cooper, Lorena Barba

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/ 

#include <cmath>
#include <stdio.h>
#include <iostream>
#include <time.h>

#define REAL double

REAL norm(REAL *x)
{
    return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

void cross(REAL *x, REAL *y, REAL *z) // z is the resulting array
{
    z[0] = x[1]*y[2] - x[2]*y[1];
    z[1] = x[2]*y[0] - x[0]*y[2];
    z[2] = x[0]*y[1] - x[1]*y[0];
}

void MV(REAL *M, REAL *V, REAL *res) // 3x3 mat-vec
{
    REAL V2[3] = {V[0], V[1], V[2]};
    for (int i=0; i<3; i++)
    {
        REAL sum = 0.;
        for (int j=0; j<3; j++)
        {
            sum += M[3*i+j]*V2[j]; 
        }
        res[i] = sum;
    }
}

REAL dot_prod(REAL *x, REAL *y) // len(3) vector dot product
{
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

void axpy(REAL *x, REAL *y, REAL *z, REAL alpha, int sign, int N)
{
    for(int i=0; i<N; i++)
    {
        z[i] = sign*alpha*x[i] + y[i];
    }
}

void ax(REAL *x, REAL *y, REAL alpha, int N)
{
    for(int i=0; i<N; i++)
    {
        y[i] = alpha*x[i];
    }

}

void lineInt(REAL *PHI, REAL z, REAL x, REAL v1, REAL v2, REAL E, REAL *xk, REAL *wk, REAL *vert_i, REAL *vert_f, REAL *center, int K, int n, REAL Area)
{
    REAL theta1 = atan2(v1,x);
    REAL theta2 = atan2(v2,x);
    REAL dtheta = theta2 - theta1;
    REAL thetam = (theta2 + theta1)/2;

    REAL epsilon_w = 80.;
    REAL a = 100.;
    REAL fact = (E - epsilon_w)/(epsilon_w + E);

    REAL Z_unit[3], Z_u_norm[3], true_center[3];
    REAL absZ = fabs(z), signZ;
    if (absZ<1e-10) signZ = 0;
    else            signZ = z/absZ;

    cross(vert_i, vert_f, Z_unit); //Revisar el signo del producto cruz, puede estar invertido.
    ax(Z_unit, Z_u_norm, norm(Z_unit), 3);
    for (int j = 0; j < 3; j++){true_center[j] = center[j] - Z_u_norm[j] * absZ;}

    // Loop over gauss points

        for (int nn = -(n - 1)/2; nn < (n + 1)/2; nn++)
        {   
            REAL R;
            if (nn==0){
                for (int i=0; i<K; i++){            
                    REAL thetak, Rtheta;
                    thetak = dtheta/2 * xk[i] + thetam;
                    Rtheta = x/cos(thetak);
                    R      = sqrt(Rtheta * Rtheta + z*z);
                    PHI[0]+= wk[i] * (R-absZ) * dtheta/2;
                    PHI[1]+= wk[i] * (z/R - signZ) * dtheta/2;
                }
            }
            else
            {
                REAL Q_i, z_i;
                Q_i = pow(fact, fabs(nn));
                z_i = nn * a + (pow(-1.0, nn)) * true_center[2];
                R = sqrt((z_i - center[2]) * (z_i - center[2]) + (Z_u_norm[0] * absZ)*(Z_u_norm[0] * absZ) + (Z_u_norm[1] * absZ)*(Z_u_norm[1]));
                PHI[0] += Area * (1./(4.*3.1415962*E)) * Q_i/R;
                PHI[1] += Area * (1./(4.*3.1415962*E)) * Q_i * (center[2] - z_i)*Z_u_norm[2]/(R * R * R);
            }
        }
        for (int i=0; i<K; i++){
            REAL thetak, Rtheta, R;
            thetak = dtheta/2 * xk[i] + thetam;
            Rtheta = x/cos(thetak);
            R = sqrt(Rtheta * Rtheta + z * z);
            PHI[2]+= wk[i] * (R-absZ) * dtheta/2;
            PHI[3]+= wk[i] * (z/R - signZ) * dtheta/2;
        }
}

void intSide(REAL *PHI, REAL *v1, REAL *v2, REAL p, REAL E, REAL *xk, REAL *wk, REAL *vert_i, REAL *vert_f, REAL *center, int K, int n, REAL Area)
{
    REAL v21[3];
    for (int i=0; i<3; i++)
    {
        v21[i] = v2[i] - v1[i];
    }

    REAL L21 = norm(v21);
    REAL v21u[3];
    ax(v21, v21u, 1/L21, 3);

    REAL unit[3] = {0.,0.,1.};
    REAL orthog[3];
    cross(unit, v21u, orthog);

    REAL alpha = dot_prod(v21,v1)/(L21*L21);

    REAL rOrthog[3];
    axpy(v21, v1, rOrthog, alpha, -1, 3); //No entiendo que cresta pasa aca

    REAL d_toEdge = norm(rOrthog); //Segun yo esto sobra
    REAL v1_neg[3];
    ax(v1, v1_neg, -1, 3);
    
    REAL side_vec[3];
    cross(v21, v1_neg, side_vec); //Esto tampoco lo usa nunca xD

    REAL rotateToVertLine[9];

    for(int i=0; i<3; i++)
    {
        rotateToVertLine[3*i] = orthog[i];
        rotateToVertLine[3*i+1] = v21u[i];
        rotateToVertLine[3*i+2] = unit[i];
    }

    REAL v1new[3];
    MV(rotateToVertLine,v1,v1new);

    if (v1new[0]<0)
    {
        ax(v21u, v21u, -1, 3);
        ax(orthog, orthog, -1, 3);
        ax(rotateToVertLine, rotateToVertLine, -1, 9);
        rotateToVertLine[8] = 1.;
        MV(rotateToVertLine,v1,v1new);
    }

    REAL v2new[3], rOrthognew[3];
    MV(rotateToVertLine,v2,v2new);
    MV(rotateToVertLine,rOrthog,rOrthognew);
    REAL x = v1new[0];

    if ((v1new[1]>0 && v2new[1]<0) || (v1new[1]<0 && v2new[1]>0))
    {
        REAL PHI1[4] = {0.,0.,0.,0.} , PHI2[4] = {0.,0.,0.,0.};
        lineInt(PHI1, p, x, 0, v1new[1], E, xk, wk, vert_i, vert_f, center, K, n, Area);
        lineInt(PHI2, p, x, v2new[1], 0, E, xk, wk, vert_i, vert_f, center, K, n, Area);

        for(int i=0; i<4; i++)
            PHI[i] += PHI1[i] + PHI2[i];
    }   
    else
    {
        REAL PHI_aux[4] = {0.,0.,0.,0.};
        lineInt(PHI_aux, p, x, v1new[1], v2new[1], E, xk, wk, vert_i, vert_f, center, K, n, Area);

        for(int i=0; i<4; i++)
            PHI[i] -= PHI_aux[i];
    }
}


void SA_wrap(REAL *PHI, REAL *y, REAL *x, REAL E, 
            int same, REAL *xk, int xkSize, REAL *wk, int n, REAL Area)
{
    // Put first panel at origin
    REAL y0_panel[3], y1_panel[3], y2_panel[3], x_panel[3];
    REAL vert0[3], vert1[3], vert2[3];
    REAL X[3], Y[3], Z[3];
    for (int i=0; i<3;i++)
    {
        x_panel[i] = x[i] - y[i];
        y0_panel[i] = 0.; 
        y1_panel[i] = y[3+i] - y[i]; 
        y2_panel[i] = y[6+i] - y[i]; 
        X[i] = y1_panel[i];

        vert0[i] = y[i];
        vert1[i] = y[i+3];
        vert2[i] = y[i+6];
    }

    // Find panel coordinate system X: 0->1
    cross(y1_panel, y2_panel, Z);
    REAL Xnorm = norm(X); 
    REAL Znorm = norm(Z); 
    for (int i=0; i<3; i++)
    {
        X[i] /= Xnorm;
        Z[i] /= Znorm;
    }

    cross(Z,X,Y);

    // Rotate the coordinate system to match panel plane
    REAL rot_matrix[9];
    for (int i=0; i<3; i++)
    {
        rot_matrix[i] = X[i];
        rot_matrix[i+3] = Y[i];
        rot_matrix[i+6] = Z[i];
    }
    
    REAL panel0_plane[3], panel1_plane[3], panel2_plane[3], x_plane[3];
    MV(rot_matrix, y0_panel, panel0_plane);
    MV(rot_matrix, y1_panel, panel1_plane);
    MV(rot_matrix, y2_panel, panel2_plane);
    MV(rot_matrix, x_panel, x_plane);

    // Shift origin so it matches collocation point
    REAL panel0_final[3], panel1_final[3], panel2_final[3];
    for (int i=0; i<3; i++)
    {
        if (i<2)
        {
            panel0_final[i] = panel0_plane[i] - x_plane[i]; 
            panel1_final[i] = panel1_plane[i] - x_plane[i]; 
            panel2_final[i] = panel2_plane[i] - x_plane[i];
        }
        else
        {
            panel0_final[i] = panel0_plane[i]; 
            panel1_final[i] = panel1_plane[i]; 
            panel2_final[i] = panel2_plane[i];
        }
    }

    // Loop over sides
    intSide(PHI, panel0_final, panel1_final, x_plane[2], E, xk, wk, vert0, vert1, x, xkSize, n, Area); // Side 0
    intSide(PHI, panel1_final, panel2_final, x_plane[2], E, xk, wk, vert1, vert2, x, xkSize, n, Area); // Side 1
    intSide(PHI, panel2_final, panel0_final, x_plane[2], E, xk, wk, vert2, vert0, x, xkSize, n, Area); // Side 2

    if (same==1)
    {
//        PHI[1] = 0.;
        PHI[3] = 0.;
    }

//    printf("PHI: %f, %f, %f, %f\n",PHI[0],PHI[1],PHI[2], PHI[3]);
}

void SA_wrap_arr(REAL *y, int ySize, REAL *x, int xSize, 
            REAL *phi_Y, int pYSize, REAL *dphi_Y, int dpYSize, 
            REAL *phi_L, int pLSize, REAL *dphi_L, int dpLSize,
            REAL E, int *same, int sameSize, REAL *xk, int xkSize, REAL *wk, int wkSize, int n, REAL Area)
{
    int N = pYSize;

    for(int i=0; i<N; i++)
    {
        REAL PHI_1[4] = {0.,0.,0.,0.};
        REAL x_1[3] = {x[3*i], x[3*i+1], x[3*i+2]};
        SA_wrap(PHI_1, y, x_1, E, same[i], xk, xkSize, wk, n, Area);
        phi_Y[i]  = PHI_1[0];
        dphi_Y[i] = PHI_1[1];
        phi_L[i]  = PHI_1[2];
        dphi_L[i] = PHI_1[3];

    }
}


void P2P_c(REAL *MY, int MYSize, REAL *dMY, int dMYSize, REAL *ML, int MLSize, REAL *dML, int dMLSize,    
        REAL *triangle, int triangleSize, int *tri, int triSize, int *k, int kSize, 
        REAL *xi, int xiSize, REAL *yi, int yiSize, REAL *zi, int ziSize,
        REAL *s_xj, int s_xjSize, REAL *s_yj, int s_yjSize, REAL *s_zj, int s_zjSize,
        REAL *xt, int xtSize, REAL *yt, int ytSize, REAL *zt, int ztSize,
        REAL *m, int mSize, REAL *mx, int mxSize, REAL *my, int mySize, REAL *mz, int mzSize, REAL *mclean, int mcleanSize, int *target, int targetSize,
        REAL *Area, int AreaSize, REAL *xk, int xkSize, 
        REAL *wk, int wkSize, REAL kappa, REAL threshold, REAL eps, REAL w0, REAL *aux, int auxSize, int n, REAL Area_i)
{
    time_t start, stop;
    int N_target = targetSize;
    int N_source = s_xjSize;
    REAL dx, dy, dz, dx_tri, dy_tri, dz_tri, R, R2, R3, R_tri, expKr;
    bool L_d, same, condition_an, condition_gq;

    for(int i_aux=0; i_aux<N_target; i_aux++)
    {  
        int i = target[i_aux];
        for(int j=0; j<N_source; j++)
        {   
            // Check if panels are far enough for Gauss quadrature
            dx_tri = xt[i_aux] - xi[tri[j]];
            dy_tri = yt[i_aux] - yi[tri[j]];
            dz_tri = zt[i_aux] - zi[tri[j]];
            R_tri  = sqrt(dx_tri*dx_tri + dy_tri*dy_tri + dz_tri*dz_tri);
            
            L_d  = (sqrt(2*Area[tri[j]])/(R_tri+eps)>=threshold);
            same = (i==tri[j]);
            condition_an = (L_d && (k[j]==0));
            condition_gq = (!L_d);

            if(condition_gq)
            {
                dx = xt[i_aux] - s_xj[j];
                dy = yt[i_aux] - s_yj[j];
                dz = zt[i_aux] - s_zj[j];
                R  = sqrt(dx*dx + dy*dy + dz*dz + eps);
                R2 = R*R;
                R3 = R2*R;
                expKr = exp(-kappa*R);
                MY[i_aux]  += m[j]*expKr/R;
                dMY[i_aux] += -expKr/R2*(kappa+1/R) * (dx*mx[j] + dy*my[j] + dz*mz[j]);
                ML[i_aux]  += m[j]/R;
                dML[i_aux] += -1/R3*(dx*mx[j] + dy*my[j] + dz*mz[j]);
                // The minus sign from the chain rule is being considered in the multiplying vector
                // since we are deriving respect to r' in 1/|r-r'|
            }
            
            if(condition_an)
            {
                aux[0] += 1;
                REAL center[3] = {xt[i_aux], yt[i_aux], zt[i_aux]};
                REAL panel[9]  = {triangle[9*tri[j]], triangle[9*tri[j]+1], triangle[9*tri[j]+2],
                                triangle[9*tri[j]+3], triangle[9*tri[j]+4], triangle[9*tri[j]+5],
                                triangle[9*tri[j]+6], triangle[9*tri[j]+7], triangle[9*tri[j]+8]};
                REAL PHI_1[4] = {0., 0., 0., 0.};
                
                start = clock();
                SA_wrap(PHI_1, panel, center, kappa, same, xk, xkSize, wk, n, Area_i);
                stop = clock();

                MY[i_aux]  += PHI_1[0] * m[j]/(w0*Area[tri[j]]);
                dMY[i_aux] += PHI_1[1] * mclean[j];
                ML[i_aux]  += PHI_1[2] * m[j]/(w0*Area[tri[j]]);
                dML[i_aux] += PHI_1[3] * mclean[j];

                aux[1] += ((REAL)(stop-start))/CLOCKS_PER_SEC;
            }
        }
    }

}

