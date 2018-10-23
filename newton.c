#include <stdio.h>
#include <stdlib.h>
/*para usar las definiciones*/
#define _USE_MATH_DEFINES
#include <math.h>

/*defino funciones a usar*/
float funcion(float V,float R,float h){
    // f(V,R;h) = V - \pi h^2 (3R - h)/3.
	return V - M_PI*h*h*(3.*R - h)/3.;
}
float dfuncion(float R,float h){
    // -\pi h (2R - h)
	return -M_PI*h*(2.*R - h);
}

/*Defino el nucleo central de mi programa*/
int main(int argc, char *argv[])
{
	float V,R,hi,h=0.;
    int n=0;
	
    printf("Mediante Newton se calcula la\n");
	printf("profundidad h de un tanque esferico de agua \n");
	printf("Ingrese el valor de V: \n");scanf("%f",&V);
	printf("Ingrese el valor de R: \n");scanf("%f",&R);
    /*valor de semilla inicial*/
    hi = 1.5*R;
    if (V > (4.*M_PI*R*R*R/3.)){
    usage: fprintf(stderr,"El Volumen de agua deseado es mayor al que se puede contener dentro del tanque de radio R \n");
    return 1;
    }
    
    if(R<=0.){ goto usage;}
    if(V<=0.){ goto usage;}
    
    int result = 1;
    double err = 1.E-6;
    
    /*ciclo para el calculo de la raiz*/
    
    while (fabs(hi-h)>err)
    {
        n++;
        h = hi;
        if ( fabs(dfuncion(R,h)) > err*err ){
            hi = h - funcion(V,R,h)/dfuncion(R,h);
            //fprintf(stderr,"raiz: %.3f\n",hi);
            }
        else { fprintf(stderr,"hay un cero en la semilla, elija otro valor\n"); break;}
        if (hi>(2.*R + 0.1)){break;}
    }
    
    if (hi>(2.*R)){fprintf(stderr,"error? altura: %f ; interacciones : %d \n",hi,n);}
    else {fprintf(stderr,"La altura es %.3f metros; interacciones : %d \n",hi,n);}
    result = 0;
    
	return result;
}
