# include "Rtlib/src/Rtlib.h"
# include <eigen-3.3.7/Eigen/Dense>

using namespace Eigen;

typedef struct BaseInfo_{
	int* exp;
	int npair; // number of pair ( vector size )
} BaseInfo;


double BaseFunc( BaseInfo* b, double* vec ) {
	int i;
	double v = 1;

	for(i = 0;i < b->npair;i++) { v *= pow( 1 - exp( - 0.5 * vec[i] ), b->exp[i] ); }

	/* chk */
	/*
	fprintf( stdout, "BaseFunc\n" );
	for(i = 0;i < b->npair;i++) { fprintf( stdout, "vec[%d] : %17.12lf, b->exp[%d] : %d\n", i, vec[i], i, b->exp[i] ); }
	*/
	// chk end //

	return v;
}


void BaseFunc_Print( BaseInfo* b ) {
	int i;
	for(i = 0;i < b->npair;i++) {
		fprintf( stdout, "( (1 - exp(-0.5 * x%d) ) ** %d )", i, b->exp[i] );
		if( i < b->npair - 1 ) fprintf( stdout, " * ");
	}
}


typedef struct AppPESInfo_{
	int natom; // npair = natom * ( natom - 1) / 2, 
	int order; // nbase = Combination( npair + order, order )
	double *coeff;
	FILE *fp;
	BaseInfo **b;
} AppPESInfo;


void ConvertMOLtoDIST( int natom, Atom* m, double* vec ) {
	int i ,j, index;
	
	index = 0;
	for(i = 0;i < natom;i++) {
		for(j = i + 1;j < natom;j++, index++) {
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
	int i, nbase, npair;
	double e = 0, *vec;

	npair = a->natom * ( a->natom - 1 ) / 2;
	nbase = Combination( npair + a->order, a->order );
	vec = ( double* )malloc( sizeof( double ) * npair );
	ConvertMOLtoDIST( a->natom, m, vec );

	for(i = 0;i < nbase;i++) {
		e += a->coeff[i] * BaseFunc( a->b[i], vec ) ;
		fprintf( stdout, "coeff[%d] : %17.12lf, BaseFunc : %17.12lf\n", i, a->coeff[i], BaseFunc( a->b[i], vec ) );
	}	// energy calculated here	// error

	/* FILE OUT START */
	/*
	fprintf( a->fp, "Energy : %17.12lf\t", e );	// energy chk error	// here
	for(i = 0;i < npair;i++) { fprintf( a->fp, "dist[%d] : %17.12lf\t", i, vec[i] ); }	// data chk ok
	fprintf( a->fp, "\n" );
	*/
	/* FILE OUT END */

	free( vec );

	fprintf( stdout, "AppPES() :%lf\n", e );
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
	int nbase, **mat, npair, i, j;
	npair = a->natom * ( a->natom - 1 ) / 2;
	nbase = Combination( npair + a->order, a->order );

//	mat = ( int** )malloc( sizeof( int* ) * nbase );
//	for(i = 0;i < nbase;i++) { mat[i] = ( int* )malloc( sizeof( int ) * npair ); }
//	MakeCombination( npair, a->order, mat ); // here

	a->b = ( BaseInfo** )malloc( sizeof( BaseInfo* ) * nbase );
	for(i = 0;i < nbase;i++) {
		a->b[i] = ( BaseInfo* ) malloc( sizeof( BaseInfo ) );
//		a->b[i]->exp = mat[i];
		a->b[i]->npair = npair;
	}

	a->coeff = ( double* ) malloc( sizeof( double ) * nbase );

//	for(i = 0;i < nbase;i++) { free( mat[i] ); }
//	free( mat );
}


void AppPES_Free( AppPESInfo* a ) {
	int nbase, npair, i;
	npair = a->natom * ( a->natom - 1 ) / 2;
	nbase = Combination( npair + a->order, a->order );

	free( a->coeff );
	for(i = 0;i < nbase;i++) {
		free( a->b[i]->exp );
		free( a->b[i] );
	}
	free( a->b );
}


void AppPES_Fitting_LS( AppPESInfo *app, Atom **mols, double *enes, int nref ) {
	int i, j, k;
	double *dist, **mat_dd;
	const int npair = ( app->natom ) * ( app->natom - 1 ) / 2;
	const int nbase = Combination( npair + app->order, app->order );
	VectorXf a( nbase ), F( nref );
	MatrixXf g( nref, nbase );

	mat_dd = ( double** )malloc( sizeof( double* ) * nref ); // basis sets
	for(i = 0;i < nref;i++) { mat_dd[i] = ( double* )malloc( sizeof( double ) * nbase ); }

	// Make matrix for fit start
	for(i = 0;i < nref;i++) {
		dist = ( double* )malloc( sizeof( double ) * npair );
		ConvertMOLtoDIST( app->natom, mols[i], dist );
		for(j = 0;j < nbase;j++) { mat_dd[i][j] = BaseFunc( app->b[j], dist ); }
		free( dist );
	}
	// Make matrix for fit end

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
	const int npair = ( app->natom ) * ( app->natom - 1 ) / 2;
	const int nbase = Combination( npair + app->order, app->order );

	mat_ii = ( int** )malloc( sizeof( int* ) * nbase );
	for(i = 0;i < nbase;i++) { mat_ii[i] = ( int* )malloc( sizeof( int ) * npair ); }
	
	MakeCombination( npair, app->order, mat_ii );

	for(i = 0;i < nbase;i++) {
		app->b[i]->exp = mat_ii[i];
		app->b[i]->npair = npair;
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
