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


void ConvertMOLtoDIST( int natom, Atom* m, double* vec ) {
	int index = 0;
	for(int i = 0;i < natom;i++) {
		for(int j = i + 1;j < natom;j++, index++) {
			vec[index] = Dist( m[i], m[j] );

			/* chk */
			/*
			m[i].Print();
			m[j].Print();
			fprintf( stdout, "dist[%d, %d] %lf\n", i, j, vec[index] );
			*/
			/* chk end */
		}
	}

	/* chk */
	/*
	for(i = 0;i < natom;i++) {
		fprintf( stdout, "m : %p,", &m[i] );
		m[i].Print();
	} // not copied
	*/
	/* chk ok */

} // Convert MOL ( Atom* ) to DIST ( double* ) // ok


double AppPES( AppPESInfo *a, Atom* m ) {
	double e = 0, *argv;
	const int dfree = 3 * a->natom - 6; 
	const int nbase = Combination( dfree + a->order, a->order );
	// const int npair = a->natom * ( a->natom - 1 ) / 2;
	// const int nbase = Combination( npair + a->order, a->order );

	argv = ( double* )malloc( sizeof( double ) * dfree );
	
	if( a->basetype == 0 ) { ConvertMOLtoDIST( a->natom, m, argv );}
	else if( a->basetype == 2 ) { CartesianToInternal( a->natom, m, argv ); }

	for(int i = 0;i < dfree;i++) { fprintf( stdout, "argv[%d] %lf\t", i, argv[i] ); }
	fprintf( stdout, "\n" );

	for(int i = 0;i < nbase;i++) {
		e += a->coeff[i] * BaseFunc( a->b[i], argv ) ;
		fprintf( stdout, "coeff[%d] : %17.12lf, BaseFunc : %17.12lf", i, a->coeff[i], BaseFunc( a->b[i], argv ) );
		for(int j = 0;j < dfree;j++) { fprintf( stdout, "\t%d(%17.12lf)", a->b[i]->exponents[j], argv[j] ); }
		fprintf( stdout, "\n" );
	}	// energy calculated here	// error

	/* FILE OUT START */
	/*
	fprintf( a->fp, "Energy : %17.12lf\t", e );	// energy chk error	// here
	for(i = 0;i < npair;i++) { fprintf( a->fp, "dist[%d] : %17.12lf\t", i, vec[i] ); }	// data chk ok
	fprintf( a->fp, "\n" );
	*/
	/* FILE OUT END */

	free( argv );
	return e;
} // get energy	// ok


void LoadfromXYZ( Atom* m, AppPESInfo* a, FILE *fp ) {
	int n, i;
	char line[256];
	Atom atom;

	fgets( line, 256, fp );
	sscanf( line, "%d", &a->natom );
	fgets( line, 256, fp );

	for(i = 0;i < a->natom;i++) {
		fgets( line, 256, fp );
		atom.SetfromString( line );
		m[i] = atom;
	}

	/* chk */
	/*
	for(i = 0;i < a->natom;i++) {
		fprintf( stdout, "m[%d] : %p,", i, &m[i] );
		m[i].Print();
	}
	*/
	/* chk ok */
}


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

		for(j = 0;j < dfree;j++) { fprintf( stdout, "argv[%d] %lf\t", j, argv[j] ); }
		fprintf( stdout, "\n" );

		for(j = 0;j < nbase;j++) { mat_dd[i][j] = BaseFunc( app->b[j], argv ); }
		free( argv );
	}
	// Make matrix for fit end



	/*
	for(i = 0;i < nref;i++) {
		argv = ( double* )malloc( sizeof( double ) * dfree );
		CartesianToInternal( app->natom, mols[i], argv );
		for(j = 0;j < dfree;j++) { fprintf( stdout, "%lf\n", argv[j] ); }
		for(j = 0;j < nbase;j++) { mat_dd[i][j] = BaseFunc( app->b[j], argv ); } // calc. value

		for(j = 0;j < dfree;j++) {
			if( j == 0 ) { argv[j] = Dist( mols[i][1], mols[i][0] ); }
			else if( j == 1 ) { argv[j] = Dist( mols[i][2], mols[i][0] ); }
			else if( j == 2 ) { argv[j] = Angle( mols[i][2], mols[i][0], mols[i][1] ); }
			else if( j %3 == 0 ) { argv[j] = Dist( mols[i][j/3 + 2], mols[i][0] ); }
			else if( j %3 == 1 ) { argv[j] = Angle( mols[i][j/3 + 2], mols[i][0], mols[i][1] ); }
			else if( j %3 == 2 ) { argv[j] = Dihedral( mols[i][j/3 + 2], mols[i][0], mols[i][1], mols[i][2] ); }
		} // Cartesian --> Internal Coordinate
		// for(j = 0;j < nbase;j++) { mat_dd[i][j] = BaseFunc( app->b[j], argv ); } // calc. value
		// for(j = 0;j < dfree;j++) { fprintf( stdout, "%lf\n", argv[j] ); }


		free( argv );
	}
	*/

	// Fitting start
	for(i = 0;i < nref;i++) { F( i ) = enes[i]; } // energy
	for(i = 0;i < nref;i++) { for(j = 0;j < nbase;j++) { g( i, j ) = mat_dd[i][j]; } } // 
	a = g.colPivHouseholderQr().solve( F );
	for(i = 0;i < nbase;i++) { app->coeff[i] = a( i ); }
	// Fitting end

	fprintf( stdout, "Coeff.\n" );
	for(i = 0;i < nbase;i++) { fprintf( stdout, "%17.12lf\n", app->coeff[i] ); }

	for(i = 0;i < nref;i++) { free( mat_dd[i] ); }
	free( mat_dd );

} // Least square fitting


void AppPES_Print( AppPESInfo* app ) {
	int i;
	const int npair = ( app->natom ) * ( app->natom - 1 ) / 2;
	const int nbase = Combination( npair + app->order, app->order );

	fprintf( stdout, "app() = " );
	for(i = 0;i < nbase;i++) {
		fprintf( stdout, "( %17.12lf ) * ", app->coeff[i] );
		BaseFunc_Print( app->b[i] );
		if(i < nbase - 1) { fprintf( stdout, " + " ); }
	}
	fprintf( stdout, "\n" );

} // print fitting func


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


int GetPointfromXYZFILE( FILE* fp, Atom** mols_ref, double* enes_ref ) {
	int natom, nref, i, j;
	char *pt, line[256];
	Atom a;

	nref = 0;
	while( fgets( line, 256, fp ) ) {
		sscanf( line, "%d", &natom );
		fgets( line, 256, fp );
		pt = strstr( line, "/" );
		sscanf( pt + 1, "%17lf", &enes_ref[nref] );

		for(i = 0;i < natom;i++) {
			fgets( line, 256, fp );
			a.SetfromString( line );
			mols_ref[nref][i] = a;
		}

		nref++;
	}

	return nref;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main( int argc, char* argv[] ) {
	AppPESInfo app;
	Atom **mols;
	double *enes;
	int nref, i, j, k, cnt = 0, npoint, maxnum = 10000;
	FILE *fp_xyz, *fp_infile, *fp_in;
	char line[256], xyz[256], *pt, infile[256];

	fp_infile = fopen( argv[1], "r");
	if( !fp_infile ) { return -1; }

	app.basetype = 2;
	while( fgets( line, 256, fp_infile ) ) {
		pt = strstr( line ,"=" );
		if( strstr( line, "ref file") ) { sscanf( pt + 1, "%s", xyz ); cnt++;  }
		else if( strstr( line, "inp file") ) { sscanf( pt + 1, "%s", infile ); cnt++;  }
		else if( strstr( line, "atom number") ) { sscanf( pt + 1, "%d", &app.natom ); cnt++; }
		else if( strstr( line, "ref number") ) { sscanf( pt + 1, "%d", &nref ); cnt++; }
		else if( strstr( line, "order") ) { sscanf( pt + 1, "%d", &app.order ); cnt++; }
	}

	fclose( fp_infile );
	if( cnt < 4 ) { return -1; }

	mols = new Atom* [maxnum];
	enes = new double [maxnum];
	for(i = 0;i < maxnum;i++) { mols[i] = new Atom [app.natom]; }

	fp_xyz = fopen( xyz, "r" );
	nref = GetPointfromXYZFILE( fp_xyz, mols, enes ); // ref.
	for(i = 0;i < nref;i++) {
		for(j = 0;j < app.natom;j++) {
			for(k = 0;k < 3;k++) { mols[i][j].SetCrd( k, ang_to_bohr( mols[i][j].GetCrd(k) ) ); }
		}
	} // convert Ang --> Bohr

	AppPES_Malloc( &app );
	AppPES_SetBaseFunc( &app );
	if( "LS" ) { AppPES_Fitting_LS( &app, mols, enes, nref ); } // Main routine ( Fitting )

	fp_in = fopen( infile, "r" );
	npoint = GetPointfromXYZFILE( fp_in, mols, enes );
	for(i = 0;i < npoint;i++) {
		for(j = 0;j < app.natom;j++) {
			for(k = 0;k < 3;k++) { mols[i][j].SetCrd( k, ang_to_bohr( mols[i][j].GetCrd(k) ) ); }
		}
	} // convert Ang --> Bohr
	for(i = 0;i < npoint;i++) { fprintf( stdout, "Point %d : %17.12lf\n", i, AppPES( &app, mols[i] ) ); } // chk ok

	// AppPES_Print( &app );
	AppPES_Free( &app );

	fclose( fp_xyz );
	fclose( fp_in );
	for(i = 0;i < maxnum;i++) { delete [] mols[i]; }
	delete [] mols;
	delete [] enes;

	return 0;
}
