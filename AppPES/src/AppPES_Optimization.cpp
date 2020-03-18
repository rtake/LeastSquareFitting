# include "Rtlib/Rtlib.h"

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


// void AppPES_LeastSquareFitting() {}
// void AppPES_SetMatrixXf() { }
// void AppPES_SetBaseFunc() {}

void AppPES_SetCoeff( AppPESInfo *a, FILE* fp ) {
	int i, nbase, npair;
	char line[256];

	npair = a->natom * ( a->natom - 1 ) / 2;
	nbase = Combination( npair + a->order, a->order );

	while( fgets( line, 256, fp ) ) { if( strstr( line, "Coeff." ) ) { break; } }
	for(i = 0;i < nbase;i++) {
		fgets( line, 256, fp );
		sscanf( line, "%17lf", &a->coeff[i] );
	}

	/* chk */
	/*
	for(i = 0;i < nbase;i++) { fprintf( a->fp, "%d\t%17.12lf\n", i, a->coeff[i] ); }
	*/
	/* chk end */
}


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


int main( int argc, char* argv[] ) {
	FILE *fp, *fp_in;
	char line[256], infiles[256], fit[256], infile[256], log[256];
	int i, j, k, num, cycle, maxcycle = 10000;
	double ene_diff, ene0, ene1, ssize = 0.002, threshold = 0.0000001, alpha = 0.0000001, *diff;

	AppPESInfo a;
	Atom *m, *m_diff;

	fp = fopen( argv[1], "r" );
	while( fgets( line, 256, fp ) ) {
		const char *pt = strstr( line, ":" );
		if( strstr( line, "order" ) ) { sscanf( pt + 2, "%d", &a.order ); }
		if( strstr( line, "nfile" ) ) { sscanf( pt + 2, "%d", &num ); }
		if( strstr( line, "infiles" ) ) { sscanf( pt + 2, "%s", infiles ); }
		if( strstr( line, "fit" ) ) { sscanf( pt + 2, "%s", fit ); }
		if( strstr( line, "atom" ) ) { sscanf( pt + 2, "%d", &a.natom ); }
		if( strstr( line, "log" ) ) { sscanf( pt + 2, "%s", log ); }
	}
	fclose( fp );

	fp = fopen( infiles, "r" );
	a.fp = fopen( log, "w" );
	for(i = 0;i < num;i++) {
		m = new Atom[a.natom];

		fprintf( a.fp, "%d start\n", i );

		fgets( line, 256, fp );
		sscanf( line, "%s", infile );
		fp_in = fopen( infile, "r" );
		LoadfromXYZ( m, &a, fp_in );
		fclose( fp_in );

		fprintf( a.fp, "AppPES_Malloc() start\n");
		AppPES_Malloc( &a );

		fprintf( a.fp, "AppPES_SetCoeff() start\n" );
		fp_in = fopen( fit, "r" );
		AppPES_SetCoeff( &a, fp_in );
		fclose( fp_in );

		/* Optimization start */

		for(cycle = 0;cycle < maxcycle;cycle++) {
			diff = new double[a.natom * 3]; // energy diff
			ene0 = AppPES( &a, m );

			for(j = 0;j < a.natom * 3;j++) {
				m_diff = new Atom[a.natom];

				for(k = 0;k < a.natom;k++) { m_diff[k] = m[k]; }
				m_diff[j / 3].SetCrd( j % 3, m_diff[j / 3].GetCrd( j % 3 ) + ssize );
				diff[j] = ( AppPES( &a, m_diff ) - ene0 ) / ssize;

				delete [] m_diff;
			} // calc. diff

			for(j = 0;j < a.natom * 3;j++) { m[j / 3].SetCrd( j % 3, m[j / 3].GetCrd( j % 3 ) - alpha * diff[j] ); } // update

			ene1 = AppPES( &a, m );
			ene_diff = ene1 - ene0;

			if( abs( ene_diff ) < threshold ) { break; }

			/* FILE OUT START */
			///*
			fprintf( stdout, "%d\n", a.natom );
			fprintf( stdout, "# Cycle. %d/%17.12lf\n", cycle, ene1 );
			for(j = 0;j < a.natom;j++) {
				fprintf( stdout, "%s", m[j].GetElm().c_str() );
				for(k = 0;k < 3;k++) { fprintf( stdout, "\t%17.12lf", m[j].GetCrd(k) ); }
				fprintf( stdout, "\n" );
			}
			//*/
			/* FILE OUT END */

			delete [] diff;
		} // while not converged 

		/* Optimization end */

		// fprintf( a.fp, "AppPES() start\n" );
		// AppPES( &a, m );

		fprintf( a.fp, "AppPES_Free() start\n" );
		AppPES_Free( &a );
		delete [] m;		
	}
	fclose( fp );
	fclose( a.fp );
	
	return 0;
}
