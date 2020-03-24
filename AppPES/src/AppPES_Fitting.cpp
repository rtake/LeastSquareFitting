# include "Rtlib/Rtlib.h"
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
		// fprintf( stdout, "coeff[%d] : %17.12lf, BaseFunc : %17.12lf\n", i, a->coeff[i], BaseFunc( a->b[i], vec ) );
	}	// energy calculated here	// error

	/* FILE OUT START */
	/*
	fprintf( a->fp, "Energy : %17.12lf\t", e );	// energy chk error	// here
	for(i = 0;i < npair;i++) { fprintf( a->fp, "dist[%d] : %17.12lf\t", i, vec[i] ); }	// data chk ok
	fprintf( a->fp, "\n" );
	*/
	/* FILE OUT END */

	free( vec );

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

	mat = ( int** )malloc( sizeof( int* ) * nbase );
	for(i = 0;i < nbase;i++) { mat[i] = ( int* )malloc( sizeof( int ) * npair ); }
	MakeCombination( npair, a->order, mat ); // here

	a->b = ( BaseInfo** )malloc( sizeof( BaseInfo* ) * nbase );
	for(i = 0;i < nbase;i++) {
		a->b[i] = ( BaseInfo* ) malloc( sizeof( BaseInfo ) );
		a->b[i]->exp = mat[i];
		a->b[i]->npair = npair;
	}

	a->coeff = ( double* ) malloc( sizeof( double ) * nbase );

//	for(i = 0;i < nbase;i++) { free( mat[i] ); }
	free( mat );
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


void AppPES_Fitting_LS( AppPESInfo* app, char* infilename ) {	
	int i, j, **mat_ii;
	double *vec_d, **mat_dd;
	const int npair = ( app->natom ) * ( app->natom - 1 ) / 2, order = 8;
	const int nbase = Combination( npair + order, order );
	VectorXf a( nbase ), F( npoint );
	MatrixXf g( npoint, nbase );

	vec_d = ( double* )malloc( sizeof( double ) * npoint );
	mat_dd = ( double** )malloc( sizeof( double* ) * npoint );
	mat_ii = ( int** )malloc( sizeof( int* ) * nbase );
	for(i = 0;i < npoint;i++) { mat_dd[i] = ( double* )malloc( sizeof( double ) * nbase ); }
	for(i = 0;i < nbase;i++) { mat_ii[i] = ( int* )malloc( sizeof( int ) * npair ); }

	// Make combination for fitting func start
	MakeCombination( npair, order, mat_ii );
	// Make combination for fitting func end


	// Set BaseFunc for matrix start
	for(i = 0;i < npoint;i++) {
		for(j = 0;j < nabse;j++) {
			// here
		}
	}
	// Set BaseFunc for matrix end


	// Make matrix for fit start
	for(i = 0;i < npoint;i++) {
		for(j = 0;j < nbase;j++) {
			mat[i][j] = BaseFunc(); // here
		}
	}
	// Make matrix for fit end


	// Fitting start
	for(i = 0;i < npoint;i++) { F( i ) = vec[i]; }
	for(i = 0;i < npoint;i++) { for(j = 0;j < nbase;j++) { g( i, j ) = mat[i][j]; } }
	a = g.colPivHouseholderQr().solve( F );
	for(i = 0;i < nbase;i++) { app.coeff[i] = a(i); }
	// Fitting end


	for(i = 0;i < npoint;i++) { free( mat_dd[i] ); }
	for(i = 0;i < nbase;i++) { free( mat_ii[i] ); }
	free( vec );
	free( mat_dd );
	free( mat_ii );

} // Least square fitting


void AppPES_Print( AppPESInfo* app ) {
} // print fitting func


void AppPES_Fitting( AppPESInfo* app, char* infilename ) {
	char ftype[256];
	sprintf( ftype, "LS" ); // default : Least Square fitting
	if( strstr( ftype, "LS" ) ) { AppPES_Fitting_LS( app, argv ); }
	AppPES_Print( app );
} // Main routine for fitting


int main( int argc, char* argv[] ) {
	AppPESInfo *app;
	AppPES_Malloc( app );
	AppPES_Fitting( app, argv[1] );
	AppPES_Free( app );
	return 0;
}
