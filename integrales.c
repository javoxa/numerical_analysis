#include <stdio.h>
#include <stdlib.h>
/*para usar las definiciones*/
#define _USE_MATH_DEFINES
#include <math.h>

/*defino la variable de bessel como global*/
double x_bessel;

double funcion(double x);
/*trapecio*/
double trapecio(double a, double b, double dx,double (*funcion)(double));
/*metodo de s 1/3*/
double simpson(double a, double b, double dt,double (*funcion)(double));
/*metodo de simpson 3/8*/
double simpson_3_8(double a, double b, double dt,double (*funcion)(double));

int main (int argc, char *argv[])
{

	/*---------------------------------*/
	
	double a,b,I0,I,dt;
    double pres;
    
    /*defino un contador*/
	int i=0;
    /*trapecio */
    /*doy valores al intervalo [a,b]*/
	a=0.;
	b=M_PI;
    /*asigno mi paso de integracion*/
    dt =(b-a)/2.;
    printf("Integraciones simples \n");
    
    /*integro*/
    /*un trapecio pedorro*/
    x_bessel = 1.;
    I0 = trapecio(a,b,dt,funcion);
	printf ("1/2: dt = %.3g \tvalor: %.16g\n",dt,I0);
    /*---------------------------------*/
    /*---------------------------------*/
    /*---------------------------------*/
    
    /*simple 1/3 */
    I0 = simpson(a,b,dt,funcion);
	printf ("1/3: dt = %.3g \tvalor: %.16g\n",dt,I0);
    /*---------------------------------*/
    /*---------------------------------*/
    /*---------------------------------*/
    /*simple 3/8*/
    I0 = simpson_3_8(a,b,dt,funcion);
    printf ("3/8: dt = %.3g \tvalor: %.16g\n",dt,I0);
    
    
    
    /*---------------------------------*/
    /*---------------------------------*/
    /*---------------------------------*/
    /* Romberg para traprecio simple, 1/3 y 3/8*/
    /*---------------------------------*/
    /*---------------------------------*/
    /*---------------------------------*/
    printf("Romberg con error = 1E-9\n");
    pres=1.E-9;
    i=0;

    /*asigno mi paso de integracion*/
    dt =(b-a)/2.;
    
    /*asigno valores (distintos) a I e I0*/
    I0 = 100.;
    /*integro*/
    while(fabs(I-I0)>pres){
    /*asigno el valor de I0 a I para comparar*/
    I = I0;
    //hago un homero simpson
    I0 = trapecio(a,b,dt,funcion);
    /*hago el intervalo mas chico por si hay una nueva interacción*/
    dt = dt/2.;
    //contador para saber cuantas veces hizo la interaccion
    i++;
    }
	
	printf ("1/2: iteraciones %d \tvalor: %.16g\n",i,I0);
    /*---------------------------------*/
    /*---------------------------------*/
    /*---------------------------------*/
    

    i=0;

    /*asigno mi paso de integracion*/
    dt =(b-a)/2.;
    
    /*asigno valores (distintos) a I e I0*/
    I0 = 100.;
    /*integro*/
    while(fabs(I-I0)>pres){
    /*asigno el valor de I0 a I para comparar*/
    I = I0;
    //hago un homero simpson
    I0 = simpson(a,b,dt,funcion);
    /*hago el intervalo mas chico por si hay una nueva interacción*/
    dt = dt/2.;
    //contador para saber cuantas veces hizo la interaccion
    i++;
    }
	
	printf ("1/3: iteraciones %d \tvalor: %.16g\n",i,I0);
    /*---------------------------------*/
    /*---------------------------------*/
    /*---------------------------------*/
    /*defino un contador*/
	i=0;
    /*asigno mi paso de integracion*/
    dt =(b-a)/2.;
    
    /*asigno valores (distintos) a I e I0*/
    I0 = 100.;
    /*integro*/
    while(fabs(I-I0)>pres){
    /*asigno el valor de I0 a I para comparar*/
    I = I0;
    //hago un homero simpson
    I0 = simpson_3_8(a,b,dt,funcion);
    /*hago el intervalo mas chico por si hay una nueva interacción*/
    dt = dt/2.;
    //contador para saber cuantas veces hizo la interaccion
    i++;
    }
	
    printf ("3/8: iteraciones %d \tvalor: %.16g\n",i,I0);

	/*---------------------------------*/
	
	return 0;
}
/*aqui pueden cambiar como quieran la funcion
    yo aqui uso una cuadratica*/
double funcion(double x){
    return cos(x_bessel*sin(x))/M_PI;
}
double trapecio(double a, double b, double dx,double (*funcion)(double)){
    /*defino mis variables a usar aqui*/
    double x,salida;
    /*el primer valor de x tiene que ser a */
    x=a;
    /*hago salida cero*/
    salida = 0.;
    
	do{
        /*integro en el intervalo [x,x+dx]
           pero ahora sumo el punto intermedio x+0.5*dx, multiplicado
           por cuatro*/
		salida = salida + dx*( funcion(x) + funcion(x+dx) )/2.;
        /*actualizo el intervalo a integrar*/
        x = x + dx;
	}
	while(x<b);
    return salida;
}

double simpson(double a, double b, double dx,double (*funcion)(double)){
    /*defino mis variables a usar aqui*/
    double x,salida;
    /*el primer valor de x tiene que ser a */
    x=a;
    /*hago salida cero*/
    salida = 0.;
    
	do{
        /*integro en el intervalo [x,x+dx]
           pero ahora sumo el punto intermedio x+0.5*dx, multiplicado
           por cuatro*/
		salida = salida + dx*( funcion(x)+ 4.*funcion(x+0.5*dx) + funcion(x+dx) )/6.;
        /*actualizo el intervalo a integrar*/
        x = x + dx;
	}
	while(x<b);
    return salida;
}

double simpson_3_8(double a, double b, double dx,double (*funcion)(double)){
    /*defino mis variables a usar aqui*/
    double x,salida;
    /*el primer valor de x tiene que ser a */
    x=a;
    /*hago salida cero*/
    salida = 0.;
    
	do{
        /*integro en el intervalo [x,x+dx]
           pero ahora sumo el punto intermedio x+0.5*dx, multiplicado
           por cuatro*/
		salida = salida + dx*( funcion(x)+ 3.*funcion(x+dx/3.) + 3.*funcion(x+2.*dx/3.) + funcion(x+dx) )/8.;
        /*actualizo el intervalo a integrar*/
        x = x + dx;
	}
	while(x<b);
    return salida;
}
