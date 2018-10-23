#include <stdio.h>
#include <stdlib.h>
/*para usar las definiciones*/
#define _USE_MATH_DEFINES
#include <math.h>

float funcion(float V,float R,float h);
int main(int argc, char *argv[])
{
	float V,R,h,hi,hf,dh;
    int n=0;
	
    printf("Mediante biseccion se calcula la\n");
	printf("profundidad h de un tanque esferico de agua \n");
	printf("Ingrese el valor de V: \n");scanf("%f",&V);
	printf("Ingrese el valor de R: \n");scanf("%f",&R);
    
    
    if (V > (4.*M_PI*R*R*R/3.)){
    usage: fprintf(stderr,"El Volumen de agua deseado es mayor al que se puede contener dentro del tanque de radio R \n");
    return 1;
    }
    
    if(R<=0.){ goto usage;}
    if(V<=0.){ goto usage;}
    
    int result = 1;
    double err = 1.E-6;
	
    /*valores de alturas*/
    dh = 1.E-3;
	hi = dh; //no puede ser cero la profundidad
    hf = 2.*R; // la profundidad no puede exceder a 2R = D
	
	/*la mitad del intervalo de busqueda es el radio del tanque*/
	h=R;
	
	while (fabs(hf-h)>err)
	{
        n++;
		if (funcion(V,R,h)*funcion(V,R,hf)<0.) {hi=h;}
		else {hf=h;}
		
		h=(hi+hf)/2.;
	}
    printf("La altura es %.3f metros ; con %d interacciones \n",h,n);
    result = 0 ;
    
	return result;
}
float funcion(float V,float R,float h)
{
    // V = \pi * h^2 *(3R -h) / 3
	return V - M_PI*h*h*( 3.*R - h)/3.;
}




