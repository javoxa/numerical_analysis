#include <stdio.h>
#include <stdlib.h>
/*para usar las definiciones*/
#define _USE_MATH_DEFINES
#include <math.h>

#define Cte -4.*M_PI*M_PI
/*ecuacion para el halley*/

double fx(double t, double vx, double vy, double vz);
double fy(double t, double vx, double vy, double vz);
double fz(double t, double vx, double vy, double vz);
double gx(double t, double x, double y, double z);
double gy(double t, double x, double y, double z);
double gz(double t, double x, double y, double z);

int main(int argc, char *argv[])
{
    double x,y,z,vx,vy,vz,t,tau,dt;

    double k1[6],k2[6],k3[6],k4[6];
    
    /*control para */
    double r,r0,err;
    
    err = 1.E-2;
    
    FILE *archivo,*plot;

    /*condiciones iniciales*/
    /*pos*/
    x=0.325514;
    y=-0.459460;
    z=0.166229;
    r0 = sqrt(x*x+y*y+z*z);
    /*vel*/
    vx=-9.096111;
    vy=-6.916686;
    vz=-1.305721;
    
    /*tiempo*/
    t=0.;
    dt = 1.E-6;
    tau = 250.;
    
    archivo=fopen("salida.txt","w");
    plot=fopen("plot.txt","w");

    fprintf(archivo,"%.3g\t%g\t%g\t%g\t%g\t%g\t%g\n",t,x,y,z,vx,vy,vz);
    fprintf(plot,"%.3g\t%g\t%g\n",t, sqrt(x*x+y*y+z*z), sqrt(vx*vx+vy*vy+vz*vz) );
while (t < tau)
	{
    /*Notese que primero calculo los k1 para todas las ecuaciones*/
    k1[0] = fx(t,vx,vy,vz);
    k1[1] = fy(t,vx,vy,vz);
    k1[2] = fz(t,vx,vy,vz);
    k1[3] = gx(t,x,y,z);
    k1[4] = gy(t,x,y,z);
    k1[5] = gz(t,x,y,z);
    
    /*Posteriormente utilizo k1 para calcular k2*/
    k2[0] = fx(t+0.5*dt, vx+0.5*dt*k1[0], vy+0.5*dt*k1[1], vz+0.5*dt*k1[2]);
    k2[1] = fy(t+0.5*dt, vx+0.5*dt*k1[0], vy+0.5*dt*k1[1], vz+0.5*dt*k1[2]);
    k2[2] = fz(t+0.5*dt, vx+0.5*dt*k1[0], vy+0.5*dt*k1[1], vz+0.5*dt*k1[2]);
    k2[3] = gx(t+0.5*dt, x+0.5*dt*k1[3], y+0.5*dt*k1[4], z+0.5*dt*k1[5]);
    k2[4] = gy(t+0.5*dt, x+0.5*dt*k1[3], y+0.5*dt*k1[4], z+0.5*dt*k1[5]);
    k2[5] = gz(t+0.5*dt, x+0.5*dt*k1[3], y+0.5*dt*k1[4], z+0.5*dt*k1[5]);
    
    
    /*ahora uso los k2 para calcular k3*/
    k3[0] = fx(t+0.5*dt, vx+0.5*dt*k2[0], vy+0.5*dt*k2[1], vz+0.5*dt*k2[2]);
    k3[1] = fy(t+0.5*dt, vx+0.5*dt*k2[0], vy+0.5*dt*k2[1], vz+0.5*dt*k2[2]);
    k3[2] = fz(t+0.5*dt, vx+0.5*dt*k2[0], vy+0.5*dt*k2[1], vz+0.5*dt*k2[2]);
    k3[3] = gx(t+0.5*dt, x+0.5*dt*k2[3], y+0.5*dt*k2[4], z+0.5*dt*k2[5]);
    k3[4] = gy(t+0.5*dt, x+0.5*dt*k2[3], y+0.5*dt*k2[4], z+0.5*dt*k2[5]);
    k3[5] = gz(t+0.5*dt, x+0.5*dt*k2[3], y+0.5*dt*k2[4], z+0.5*dt*k2[5]);
    
    /*finalmente calculos los k4*/
    k4[0] = fx(t+dt, vx+dt*k3[0], vy+dt*k3[1], vz+dt*k3[2]);
    k4[1] = fy(t+dt, vx+dt*k3[0], vy+dt*k3[1], vz+dt*k3[2]);
    k4[2] = fz(t+dt, vx+dt*k3[0], vy+dt*k3[1], vz+dt*k3[2]);
    k4[3] = gx(t+dt, x+dt*k3[3], y+dt*k3[4], z+dt*k3[5]);
    k4[4] = gy(t+dt, x+dt*k3[3], y+dt*k3[4], z+dt*k3[5]);
    k4[5] = gz(t+dt, x+dt*k3[3], y+dt*k3[4], z+dt*k3[5]);
    
    /*finalmente calculo los valors x e y*/
    x = x + (dt/6.)*(k1[0] + 2.*k2[0] + 2.*k3[0] + k4[0]);
    y = y + (dt/6.)*(k1[1] + 2.*k2[1] + 2.*k3[1] + k4[1]);
    z = z + (dt/6.)*(k1[2] + 2.*k2[2] + 2.*k3[2] + k4[2]);
    vx = vx + (dt/6.)*(k1[3] + 2.*k2[3] + 2.*k3[3] + k4[3]);
    vy = vy + (dt/6.)*(k1[4] + 2.*k2[4] + 2.*k3[4] + k4[4]);
    vz = vz + (dt/6.)*(k1[5] + 2.*k2[5] + 2.*k3[5] + k4[5]);
	
    t = t + dt;
	r = sqrt(x*x+y*y+z*z);
        if (fabs(r-r0)>err){
        fprintf(archivo,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n",t,x,y,z,vx,vy,vz);
        fprintf(plot,"%g\t%g\t%g\n",t, sqrt(x*x+y*y+z*z), sqrt(vx*vx+vy*vy+vz*vz) );
        r0=r;
        }
    }
fclose(archivo);
fclose(plot);
//system("echo 'source(archivor)'|r -persist");

    //system("echo \"source(\'grafico.txt\')\"| r --vanilla");

exit(0);

}
double fx(double t, double vx, double vy, double vz){
    return vx;
    }
double fy(double t, double vx, double vy, double vz){
    return vy;
    }
double fz(double t, double vx, double vy, double vz){
    return vz;
    }

double gx(double t, double x, double y, double z){
    float r;
    r = sqrt(x*x + y*y + z*z);
    return Cte*x/(r*r*r);
    }
double gy(double t, double x, double y, double z){
    float r;
    r = sqrt(x*x + y*y + z*z);
    return Cte*y/(r*r*r);
    }
double gz(double t, double x, double y, double z){
    float r;
    r = sqrt(x*x + y*y + z*z);
    return Cte*z/(r*r*r);
    }
