/*
double* tridiagonal_solver(int N, double** matrix, double* d)
{
	
	matrix[0][1] = matrix[0][1]/matrix[0][0]; //i=0;
	d[0] = d[0]/matrix[0][0];//i=0;
	for(int i=1; i<N-1; i++)
	{
		matrix[i][i+1] = matrix[i][i+1]/(matrix[i][i] - matrix[i][i-1]*matrix[i-1][i] );	
		d[i] = (d[i] - matrix[i][i-1]*d[i-1])/(matrix[i][i] - matrix[i][i-1]*matrix[i-1][i]);
	}
	
	d[N-1] = (d[N-1] - matrix[N-1][N-2]*d[N-2])/(matrix[N-1][N-1] - matrix[N-1][N-2]*matrix[N-2][N-1]);
	
	double* x;
	x = new double[N];
	x[N-1] = d[N-1];
	
	for(int i=N-2; i>=0; i--)
		x[i] = d[i] - matrix[i][i+1]*x[i+1];
	
	return x;
}

double* cranknicolson(double* D, double* ro, const int N, const int Npassi, const double h, const double dt)
{
	double h_h_dt = 2.0*h*h/dt;
	
	double **A,**B;
	double *r ; //r = B*ro, right side of equation 
	//A,B, r allocation
	A = new double*[N];
	B = new double*[N];
	for(int i=0; i<N; i++)
	{
		A[i] = new double[N];
		B[i] = new double[N];
	}
	r  = new double[N];	
	
	//A,B,r initialization
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
		{
			A[i][j] = 0;
			B[i][j] = 0;
		}
		r[i] = 0;
	}
	
	
	for(int i=1; i<N; i++)
	{
		A[i][i-1] = - 0.25*D[i+2] - 0.75*D[i];
		A[i][i]   = h_h_dt + D[i+2] + D[i];
		A[i][i+1] = - 0.75*D[i+2] - 0.25*D[i];
		
		B[i][i-1] =   0.25*D[i+2] + 0.75*D[i];
		B[i][i]   = h_h_dt - D[i+2] - D[i];
		B[i][i+1] =   0.75*D[i+2] + 0.25*D[i];
	}
	
	//boundary conditions ro(0)=0; ro(L)=0;
	A[0][0] = h_h_dt + D[2] + D[0];
	A[0][1] = - 0.75*D[2] - 0.25*D[0];
	B[0][0] = h_h_dt - D[2] - D[0];
	A[N-1][N-1] = h_h_dt + D[N+1] + D[N-1];
	A[N-1][N-2] = - 0.25*D[N+1] - 0.75*D[N-1];
	B[N-1][N-1] = h_h_dt - D[N+1] - D[N-1];
	
	
	for(double t=0; t<Npassi*dt; t+=dt)
	{
		for(int i=0; i<N; i++) r[i] = 0.0;
		//find the right-hand side for the solver
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++)
				r[i] += B[i][j]*ro[j];
		
		//boundary conditions
		r[0] = 0;
		ro[0] = 0;
		r[N-1] = 0;
		ro[N-1] = 0;
		
		ro = tridiagonal_solver(N,A,r);
	}
	
	return ro;
}
*/
double* tridiagonal_solver(int N, const double* a, const double *b, const double *c, const double *d)
{
	double *c1 = new double[N];
	double *d1 = new double[N];
	
	c1[0] = c[0]/b[0]; //i=0;
	d1[0] = d[0]/b[0]; //i=0;
	
	for(int i=1; i<N; i++)
	{
		double m = 1.0/(b[i] - a[i]*c1[i-1] );
		c1[i] = m * c[i];	
		d1[i] = m *(d[i] - a[i]*d1[i-1]);
	}
	
	double *x = new double[N];
	x[N-1] = d1[N-1];
	for(int i=N-2; i>=0; i--)
		x[i] = d1[i] - c1[i]*x[i+1];
	
	return x;
}

double* cranknicolson(double* D, double* ro, const int N, const int Npassi, const double h, const double dt)
{
	double h_h_dt = 2.0*h*h/dt;
	
	double *a = new double[N];
	double *b = new double[N];
	double *c = new double[N];
	double *d = new double[N];
	
	double **B;

	B = new double*[N];
	for(int i=0; i<N; i++) B[i] = new double[N];
	
	//B initialization
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			B[i][j] = 0;
	
	for(int i=0; i<N; i++)
	{
		a[i] = - 0.25*D[i+2] - 0.75*D[i];
		b[i]   = h_h_dt + D[i+2] + D[i];
		c[i] = - 0.75*D[i+2] - 0.25*D[i];
	}	
	//put coefficients outside of the matrix = 0.0
	a[0] = 0.0;
	c[N-1] = 0.0;
	
	for(int i=1; i<N-1; i++)
	{
		B[i][i-1] = - a[i];
		B[i][i]   = h_h_dt - D[i+2] - D[i];
		B[i][i+1] = - c[i];
	}	
	
	B[0][0] = h_h_dt - D[2] - D[0];
	B[0][1] = - c[0];
	B[N-1][N-1] = h_h_dt - D[N+1] - D[N-1];
	B[N-1][N-2] = - a[N-1];
	
	for(double t=0; t<Npassi*dt; t+=dt)
	{
		for(int i=0; i<N; i++) d[i] = 0.0;
		//find the right-hand side for the solver
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++)
				d[i] += B[i][j]*ro[j];
		
		//boundary conditions
		d[0] = 0;
		//ro[0] = 0;
		d[N-1] = 0;
		//ro[N-1] = 0;
		
		ro = tridiagonal_solver(N,a,b,c,d);
	}
	
	return ro;
}

double* cranknicolson2(double* D, double* ro, const int N, const int Npassi, const double h, const double dt)
{
	double h_h_dt = 2.0*h*h/dt;
	
	double *a = new double[N];
	double *b = new double[N];
	double *c = new double[N];
	double *d = new double[N];
	
	double **B;

	B = new double*[N];
	for(int i=0; i<N; i++) B[i] = new double[N];
	
	//B initialization
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			B[i][j] = 0;
	
	for(int i=0; i<N; i++)
	{
		a[i] =   0.25*D[i+2] - D[i+1] - 0.25*D[i];
		b[i] = h_h_dt + 2.0*D[i+1];
		c[i] = - 0.25*D[i+2] - D[i+1] + 0.25*D[i];
	}	
	//put coefficients outside of the matrix = 0.0
	a[0] = 0.0;
	c[N-1] = 0.0;
	
	for(int i=1; i<N-1; i++)
	{
		B[i][i-1] = - a[i];
		B[i][i]   = h_h_dt - 2.0*D[i+1];
		B[i][i+1] = - c[i];
	}	
	
	B[0][0] = h_h_dt - 2.0*D[1];
	B[0][1] = - c[0];
	B[N-1][N-1] = h_h_dt - 2.0*D[N];
	B[N-1][N-2] = - a[N-1];
	
	for(double t=0; t<Npassi*dt; t+=dt)
	{
		for(int i=0; i<N; i++) d[i] = 0.0;
		//find the right-hand side for the solver
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++)
				d[i] += B[i][j]*ro[j];
		
		//boundary conditions
		d[0] = 0;
		//ro[0] = 0;
		d[N-1] = 0;
		//ro[N-1] = 0;
		
		ro = tridiagonal_solver(N,a,b,c,d);
	}
	
	return ro;
}

double* cranknicolson_fokkerplanck(double * mu, double* D, double* ro, const int N, const int Npassi, const double h, const double dt)
{
	double h_h_dt = 2.0*h*h/dt;
	
	double *a = new double[N];
	double *b = new double[N];
	double *c = new double[N];
	double *d = new double[N];
	
	double **B;

	B = new double*[N];
	for(int i=0; i<N; i++) B[i] = new double[N];
	
	//B initialization
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			B[i][j] = 0;
	
	for(int i=0; i<N; i++)
	{
		a[i] = - D[i]   - 0.25*mu[i];
		b[i] = h_h_dt + 2.0*D[i+1];
		c[i] = - D[i+2] + 0.25*mu[i];
	}	
	//put coefficients outside of the matrix = 0.0
	a[0] = 0.0;
	c[N-1] = 0.0;
	
	for(int i=1; i<N-1; i++)
	{
		B[i][i-1] = - a[i];
		B[i][i]   = h_h_dt - 2.0*D[i+1];
		B[i][i+1] = - c[i];
	}	
	
	B[0][0] = h_h_dt - 2.0*D[1];
	B[0][1] = - c[0];
	B[N-1][N-1] = h_h_dt - 2.0*D[N];
	B[N-1][N-2] = - a[N-1];
	
	for(double t=0; t<Npassi*dt; t+=dt)
	{
		for(int i=0; i<N; i++) d[i] = 0.0;
		//find the right-hand side for the solver
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++)
				d[i] += B[i][j]*ro[j];
		
		//boundary conditions
		d[0] = 0;
		//ro[0] = 0;
		d[N-1] = 0;
		//ro[N-1] = 0;
		
		ro = tridiagonal_solver(N,a,b,c,d);
	}
	
	return ro;
}