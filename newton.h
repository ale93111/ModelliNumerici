#include <iostream>
#include <math.h>

//#ifndef NEWTON_H
//#define NEWTON_H

//#include "ensemble_henon.h"

//creare una classe con la funzione e i suoi coefficienti e poi inizializzare la funzione
//con un puntatore al membro della classe


//Using  a  combination  of  Newton-Raphson  and  bisection,  find  the  root  of  a  function  bracketed
//between x1 and x2.  The  root,  returned  as  the  function  value rtsafe,  will  be  refined  until
//its  accuracy  is  known  within ±xacc. funcd is  a  user-supplied  routine  that  returns  both  the
//function  value  and  the  first  derivative  of  the  function.
double rtsafe(void funcd(double, double *, double *, double), double x1, double x2, double xacc, double offset = 0.0)
{
	//Maximum allowed number of iterations.
	const int MAXIT = 100;
	int j;
	double df,dx,dxold,f,fh,fl;
	double temp,xh,xl,rts;
	
	funcd(x1,&fl,&df,offset);
	funcd(x2,&fh,&df,offset);
	
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
		std::cout << "Root must be bracketed in rtsafe" << std::endl;
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		//Orient the search so that f(xl)<0. 
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	rts=0.5*(x1+x2);
	//Initialize the  guess  for root,
	dxold=fabs(x2-x1);
	//the “stepsize  before  last,”
	dx=dxold;
	//and  the  last  step.
	funcd(rts,&f,&df,offset);
	for (j=1;j<=MAXIT;j++) {
		//Loop  over  allowed  iterations.
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
		//Bisect if Newton out of range,
		|| (fabs(2.0*f) > fabs(dxold*df))) {
		//or not decreasing fast enough.
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts;
			//Change  in root  is negligible.
		} else {
			//Newton  step acceptable.  Take it.
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		//Convergence  criterion.
		funcd(rts,&f,&df,offset);
		//The  one  new  function  evaluation  per  iteration.
		if (f < 0.0)
		//Maintain the  bracket  on the root.
			xl=rts;
		else
			xh=rts;
	}
	std::cout << "Maximum number of iterations exceeded in rtsafe" << std::endl;
	return 0.0;
	//Never  get  here.
}
//#endif