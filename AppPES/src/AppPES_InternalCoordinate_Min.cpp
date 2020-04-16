# include "Rtlib/src/Rtlib.h"
# include <eigen-3.3.7/Eigen/Dense>

using namespace Eigen;

typedef struct BaseInfo_{
	int* exponents; // vector of exponent of each elem
	int argc; // number of functions consists base function ( vector size of exponents )
	int basetype; // the type of basis function ( 0 : exponential of distance, 1 : inverse of distance, 2 : internal coordinaite )
} BaseInfo;


double BaseFunc( BaseInfo* b, double* argv ) {
	double v = 1;
	if( b->basetype == 0 ) { for(int i = 0;i < b->argc;i++) { v *= pow( 1 - exp( - 0.5 * argv[i] ), b->exponents[i] ); } }
	else if( b->basetype == 1 ) {}
	else if( b->basetype == 2 ) {
		for(int i = 0;i < b->argc;i++) {
			if( i %3 == 0 || i == 1 ) { v *= pow( 1 - exp( - 0.5 * log( argv[i] ) ), b->exponents[i] ); } // distance
			else if( i %3 == 1 || i == 2 ) { v *= pow( log( argv[i] ), b->exponents[i] ); } // angle
			else if( i %3 == 2 ) {} // dihedral
		}
	}
	return v;
}


void BaseFunc_Print( BaseInfo* b ) {}


typedef struct AppPESInfo_{
	int natom; // npair = natom * ( natom - 1) / 2, 
	int order; // nbase = Combination( npair + order, order )
	double *coeff;
	FILE *fp;
	BaseInfo **b;
	int basetype;
} AppPESInfo;


typedef struct AppPES_MinInfo {
	AppPESInfo *app;
	int maxoptitr;
	double stepsize;
} AppPES_MinInfo;


double AppPES( AppPESInfo *a, Atom* m ) {
	double e = 0, *argv;
	const int dfree = 3 * a->natom - 6; 
	const int nbase = Combination( dfree + a->order, a->order );
	// const int npair = a->natom * ( a->natom - 1 ) / 2;
	// const int nbase = Combination( npair + a->order, a->order );

	argv = ( double* )malloc( sizeof( double ) * dfree );
	
	if( a->basetype == 0 ) { ConvertMOLtoDIST( a->natom, m, argv );}
	else if( a->basetype == 2 ) { CartesianToInternal( a->natom, m, argv ); }

	/*
	for(int i = 0;i < dfree;i++) { fprintf( stdout, "argv[%d] %lf\t", i, argv[i] ); }
	fprintf( stdout, "\n" );
	*/

	for(int i = 0;i < nbase;i++) {
		e += a->coeff[i] * BaseFunc( a->b[i], argv ) ;

		/*
		fprintf( stdout, "coeff[%d] : %17.12lf, BaseFunc : %17.12lf", i, a->coeff[i], BaseFunc( a->b[i], argv ) );
		for(int j = 0;j < dfree;j++) { fprintf( stdout, "\t%d(%17.12lf)", a->b[i]->exponents[j], argv[j] ); }
		fprintf( stdout, "\n" );
		*/
	}

	free( argv );
	return e;
} // get energy	// ok


void AppPES_Malloc( AppPESInfo* a ) {
	const int dfree = 3 * a->natom - 6; // degree of freedom
	// const int npair = a->natom * ( a->natom - 1 ) / 2;
	const int nbase = Combination( dfree + a->order, a->order );

	a->b = ( BaseInfo** )malloc( sizeof( BaseInfo* ) * nbase );
	for(int i = 0;i < nbase;i++) {
		a->b[i] = ( BaseInfo* ) malloc( sizeof( BaseInfo ) );
		a->b[i]->argc = dfree;
		a->b[i]->basetype = a->basetype;
	}

	a->coeff = ( double* ) malloc( sizeof( double ) * nbase );
}


void AppPES_Free( AppPESInfo* a ) {
	const int dfree = 3 * a->natom - 6;
	// npair = a->natom * ( a->natom - 1 ) / 2;
	const int nbase = Combination( dfree + a->order, a->order );

	free( a->coeff );
	for(int i = 0;i < nbase;i++) {
		free( a->b[i]->exponents );
		free( a->b[i] );
	}
	free( a->b );
}


void AppPES_Fitting_LS( AppPESInfo *app, Atom **mols, double *enes, int nref ) {
	int i, j, k;
	double *argv, **mat_dd;
	const int dfree = 3 * app->natom - 6;
	const int nbase = Combination( dfree + app->order, app->order );
	// const int npair = ( app->natom ) * ( app->natom - 1 ) / 2;
	// const int nbase = Combination( npair + app->order, app->order );

	InternalCoordinate *icrd;
	VectorXf a( nbase ), F( nref );
	MatrixXf g( nref, nbase );

	mat_dd = ( double** )malloc( sizeof( double* ) * nref ); // matrix of ( double, double ) <-- value of basis sets
	for(i = 0;i < nref;i++) { mat_dd[i] = ( double* )malloc( sizeof( double ) * nbase ); }

	icrd = ( InternalCoordinate* )malloc( sizeof( InternalCoordinate ) * app->natom );

	// Make matrix for fit start
	for(i = 0;i < nref;i++) {
		argv = ( double* )malloc( sizeof( double ) * dfree );
		if( app->basetype == 0 ) { ConvertMOLtoDIST( app->natom, mols[i], argv ); }
		else if( app->basetype == 2 ) { CartesianToInternal( app->natom, mols[i], argv ); }
		for(j = 0;j < nbase;j++) { mat_dd[i][j] = BaseFunc( app->b[j], argv ); }
		free( argv );
	}

	// Fitting start
	for(i = 0;i < nref;i++) { F( i ) = enes[i]; } // energy
	for(i = 0;i < nref;i++) { for(j = 0;j < nbase;j++) { g( i, j ) = mat_dd[i][j]; } } // 
	a = g.colPivHouseholderQr().solve( F );
	for(i = 0;i < nbase;i++) { app->coeff[i] = a( i ); }

	fprintf( stdout, "Coeff.\n" );
	for(i = 0;i < nbase;i++) { fprintf( stdout, "%17.12lf\n", app->coeff[i] ); }

	for(i = 0;i < nref;i++) { free( mat_dd[i] ); }
	free( mat_dd );

} // Least square fitting


void AppPES_Print( AppPESInfo* app ) {} // print fitting func


void AppPES_SetBaseFunc( AppPESInfo* app ) {
	int i, j, **mat_ii;
	const int dfree = 3 * app->natom - 6;
	const int nbase = Combination( dfree + app->order, app->order );
	// const int npair = ( app->natom ) * ( app->natom - 1 ) / 2;
	// const int nbase = Combination( npair + app->order, app->order );

	mat_ii = ( int** )malloc( sizeof( int* ) * nbase );
	for(i = 0;i < nbase;i++) { mat_ii[i] = ( int* )malloc( sizeof( int ) * dfree ); }
	
	MakeCombination( dfree, app->order, mat_ii );

	for(i = 0;i < nbase;i++) {
		app->b[i]->exponents = mat_ii[i];
		app->b[i]->argc = dfree;
	}

	free( mat_ii );
}


void AppPES_Minimization( AppPES_MinInfo *min, Atom *mol ) {
	int i, j, k, l, itr, convergechk = 0, dfree = 3 * min->app->natom - 6; // degree of freedom
	double *grad, **ene_diff, ***icrd_diff, threshold = 0.000001, alpha = 0.0001, icrd[dfree]; // internal coordinate
	Atom ***mol_diff;

	grad = ( double* )malloc( sizeof( double ) * dfree );
	ene_diff = ( double** )malloc( sizeof( double* ) * dfree );
	icrd_diff = ( double*** )malloc( sizeof( double** ) * dfree );
	mol_diff = ( Atom*** )malloc( sizeof( Atom** ) * dfree );
	for(i = 0;i < dfree;i++) {
		ene_diff[i] = ( double* )malloc( sizeof( double ) * 2 );
		icrd_diff[i] = ( double** )malloc( sizeof( double* ) * 2 );
		mol_diff[i] = ( Atom** )malloc( sizeof( Atom* ) * 2 );
		for(j = 0;j < 2;j++) {
			icrd_diff[i][j] = ( double* )malloc( sizeof( double) * dfree );
			mol_diff[i][j] = new Atom [min->app->natom];
			for(k = 0;k < min->app->natom;k++) { mol_diff[i][j][k] = mol[k]; }
		}
	}

	CartesianToInternal( min->app->natom, mol, icrd ); // bohr

	for(itr = 0;itr < min->maxoptitr && convergechk == 0;itr++) {
		convergechk = 1;
		InternalToCartesian( min->app->natom, mol, icrd ); // bohr

		for(i = 0;i < dfree;i++) {
			for(j = 0;j < 2;j++) {
				for(k = 0;k < dfree;k++) { icrd_diff[i][j][k] = icrd[k]; } // copy
				if( j == 0 ) { icrd_diff[i][j][i] += min->stepsize; }
				else { icrd_diff[i][j][i] -= min->stepsize; }
				InternalToCartesian( min->app->natom, mol_diff[i][j], icrd_diff[i][j] );
				ene_diff[i][j] = AppPES( min->app, mol_diff[i][j] );
			}
			grad[i] = ( ene_diff[i][0] - ene_diff[i][1] ) / ( 2 * min->stepsize );
			icrd[i] -= alpha * grad[i];
			if( grad[i] > threshold ) { convergechk *= 0; } // if all grad is lower than threshold, convergechk = 1
		}

		for(i = 0;i < min->app->natom;i++) { mol[i].bohrtoang(); } // bohr --> ang

		fprintf( stdout, "# ITR. %d\n", itr );
		for(i = 0;i < min->app->natom;i++) {
			fprintf( stdout, "%s", mol[i].GetElm().c_str() );
			for(j = 0;j < 3;j++) { fprintf( stdout, "\t%17.12lf", mol[i].GetCrd( j ) ); }
			fprintf( stdout, "\n" );
		}	
		fprintf( stdout, "Item    Value\nENERGY   %17.12lf\n", AppPES( min->app, mol ) );
		for(i = 0;i < dfree;i++) {
			fprintf( stdout, "gradient[%d]\t%17.12lf\n", i, grad[i] );
			for(j = 0;j < 2;j++) {
				for(k = 0;k < min->app->natom;k++) {
					mol_diff[i][j][k].bohrtoang();
					fprintf( stdout, "%s", mol_diff[i][j][k].GetElm().c_str() );
					for(l = 0;l < 3;l++) { fprintf( stdout, "\t%17.12lf", mol_diff[i][j][k].GetCrd( l ) ); }
					fprintf( stdout, "\n" );
				}
				fprintf( stdout, "energy   %17.12lf\n", ene_diff[i][j] );
			}
		}
		fprintf( stdout,"\n" );

	} // optcycle

	for(i = 0;i < dfree;i++) { for(j = 0;j < 2;j++) { free( icrd_diff[i][j] ); delete [] mol_diff[i][j]; }
	free( ene_diff[i] ); free( icrd_diff[i] ); free( mol_diff[i] ); }
	free( grad ); free( ene_diff ); free( icrd_diff ); free( mol_diff );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {
	int natom, nref, i, j, k, cnt = 0, npoint, maxnum = 10000, fitchk = 0;
	double enes[maxnum];
	FILE *fp_ref, *fp_input;
	char *pt, line[256], ref[256], input[256], infile[256];
	Atom *mols[maxnum], mol[100];
	AppPESInfo *app = new AppPESInfo;
	AppPES_MinInfo *min = new AppPES_MinInfo;

	fp_input = fopen( argv[1], "r" );
	if( !fp_input ) { return -1; }
	fgets( line, 256, fp_input ); // line 0
	if( strstr( line, "%" ) ) { sprintf( infile, "%s", line ); fgets( line, 256, fp_input ); }
	fgets( line, 256, fp_input ); // line 1 (comment)
	fgets( line, 256, fp_input ); // line 2 (charge and multiplicity)
	for(natom = 0;fgets( line, 256, fp_input );natom++) {
		if( strstr( line, "Op") || strstr( line, "op") ) { break; }
		mol[natom].SetfromString( line );		
	}
	while( fgets( line, 256, fp_input ) ) {
		pt = strstr( line ,"=" );
		if( strstr( line, "fitorder") ) { sscanf( pt + 1, "%d", &app->order ); }
		else if( strstr( line, "reference") ) { sscanf( pt + 1, "%s", ref ); }
		else if( strstr( line, "opt=internal") ) { app->basetype = 2;}
		else if( strstr( line, "maxoptitr") ) { sscanf( pt + 1, "%d", &min->maxoptitr ); }
		else if( strstr( line, "stepsize") ) { sscanf( pt + 1, "%lf", &min->stepsize ); }
		else if( strstr( line, "fit=leastsquare") ) { fitchk = 0 ; }
	}
	fclose( fp_input );

	min->app = app;
	app->natom = natom;
	AppPES_Malloc( app );
	AppPES_SetBaseFunc( app );

	// Get reference point
	for(i = 0;i < maxnum;i++) { mols[i] = new Atom [natom]; }
	fp_ref = fopen( ref, "r" ); nref = GetPointfromXYZFILE( fp_ref, mols, enes ); fclose( fp_ref ); // get ref. dat
	for(i = 0;i < nref;i++) { for(j = 0;j < natom;j++) { mols[i][j].angtobohr(); } } // convert Ang --> Bohr

	// Fitting
	if( fitchk == 0 ) { AppPES_Fitting_LS( app, mols, enes, nref ); } // Main routine ( Fitting )

	// Minimization
	for(i = 0;i < natom;i++) { mol[i].angtobohr(); }
	AppPES_Minimization( min, mol );

	// Free
	AppPES_Free( app );
	for(i = 0;i < maxnum;i++) { delete [] mols[i]; }
	delete app, min;

	return 0;
}
