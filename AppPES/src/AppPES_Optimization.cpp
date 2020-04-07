# include "Rtlib/src/Rtlib.h"

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

	// fprintf( stdout, "AppPES() :%lf\n", e );
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
	FILE *fp, *fp_in, *fp_out = stdout;
	char line[256], input[256], fit[256];
	int i, j, k, l, m, num, itr, maxitr = 10000, chk = 0;
	double ene0, ene1, ssize = 0.002, threshold = 0.0000001, alpha = 0.002, *grad, *ene_fwd_diff, *ene_bck_diff;

	AppPESInfo a;
	Atom *mol, **mols_fwd_diff, **mols_bck_diff;

	fp_in = fopen( argv[1], "r" );
	while( fgets( line, 256, fp_in ) ) {
		const char *pt = strstr( line, "=" );
		if( strstr( line, "input" ) ) { sscanf( pt + 1, "%s", input ); }
		else if( strstr( line, "fit_order" ) ) { sscanf( pt + 1, "%d", &a.order ); }
		else if( strstr( line, "fit_param" ) ) { sscanf( pt + 1, "%s", fit ); }
		else if( strstr( line, "maxoptitr" ) ) { sscanf( pt + 1, "%d", &maxitr ); }
		else if( strstr( line, "stepsize" ) ) { sscanf( pt + 1, "%d", &ssize ); }
	}
	fclose( fp_in );

	mol = new Atom[1000]; // memory alloc for molcule

	fp_in = fopen( input, "r" );
	LoadfromXYZ( mol, &a, fp_in ); // load input
	for(j = 0;j < a.natom;j++) { for(k = 0;k < 3;k++) { mol[j].SetCrd( k, ang_to_bohr( mol[j].GetCrd(k) ) ); } } // convert Ang --> Bohr
	fclose( fp_in );

	AppPES_Malloc( &a ); // memory alloc

	fp_in = fopen( fit, "r" );
	AppPES_SetCoeff( &a, fp_in ); // get param for AppPES()
	fclose( fp_in );

	chk = 1;
	for(itr = 0;itr < maxitr && chk > 0;itr++) {
		grad = new double[a.natom * 3]; // gradient
		ene_fwd_diff = new double [a.natom * 3]; // energy for diff. str.
		ene_bck_diff = new double [a.natom * 3];
		mols_fwd_diff = new Atom* [a.natom * 3];
		mols_bck_diff = new Atom* [a.natom * 3];
		for(j = 0;j < a.natom;j++) {
			for(k = 0;k < 3;k++) {
				mols_fwd_diff[3*j + k] = new Atom[a.natom];
				mols_bck_diff[3*j + k] = new Atom[a.natom];
			}
		} 

		ene0 = AppPES( &a, mol ); // energy for current str.
		
		for(j = 0;j < a.natom;j++) {
			for(k = 0;k < 3;k++) {
				for(l = 0;l < a.natom;l++) {
					mols_fwd_diff[3*j + k][l] = mol[l];
					mols_bck_diff[3*j + k][l] = mol[l];
				} // copy
				mols_fwd_diff[3*j + k][j].SetCrd( k, mol[j].GetCrd( k ) + ang_to_bohr( ssize ) );
				mols_bck_diff[3*j + k][j].SetCrd( k, mol[j].GetCrd( k ) - ang_to_bohr( ssize ) );
				ene_fwd_diff[3*j + k] = AppPES( &a, mols_fwd_diff[3*j + k] );
				ene_bck_diff[3*j + k] = AppPES( &a, mols_bck_diff[3*j + k] );	
				grad[3*j + k] = ( ene_fwd_diff[3*j + k] - ene_bck_diff[3*j + k] ) / ang_to_bohr( ssize * 2 );
				if( abs( grad[3*j + k] ) < threshold ) { chk = -1; }
				if( k > 0 ) { grad[3*j + k] = 0; }
			}
		} // calc. grad.


		// FILE OUT START
		fprintf( fp_out, "# ITR. %d\n", itr);
		for(j = 0;j < a.natom;j++) {
			fprintf( fp_out, "%s", mol[j].GetElm().c_str() );
			for(k = 0;k < 3;k++) { fprintf( fp_out, "\t%17.12lf", bohr_to_ang( mol[j].GetCrd(k) ) ); }
			fprintf( fp_out, "\n" );
		} // geometry
		fprintf( fp_out, "Item          Value\nENERGY          %17.12lf\n\n", ene0 ); 
		for(j = 0;j < a.natom;j++) {
			for(k = 0;k < 3;k++) {
				fprintf( fp_out, "FWD DIFF Str. %d\n", 3*j + k );
				for(l = 0;l < a.natom;l++) {
					fprintf( fp_out, "%s", mols_fwd_diff[3*j + k][l].GetElm().c_str() );
					for(m = 0;m < 3;m++) { fprintf( fp_out, "\t%17.12lf", bohr_to_ang( mols_fwd_diff[3*j + k][l].GetCrd( m ) ) ); }
					fprintf( fp_out, "\n" );
				}
				fprintf( fp_out, "ENERGY          %17.12lf\n", ene_fwd_diff[3*j + k] );

				fprintf( fp_out, "BCK DIFF Str. %d\n", 3*j + k );
				for(l = 0;l < a.natom;l++) {
                                        fprintf( fp_out, "%s", mols_bck_diff[3*j + k][l].GetElm().c_str() );
                                        for(m = 0;m < 3;m++) { fprintf( fp_out, "\t%17.12lf", bohr_to_ang( mols_bck_diff[3*j + k][l].GetCrd( m ) ) ); }
                                        fprintf( fp_out, "\n" );
                                }
				fprintf( fp_out, "ENERGY          %17.12lf\n", ene_bck_diff[3*j + k] );
				fprintf( fp_out, "GRADIENT[%d]\t%17.12lf(%17.12lf - %17.12lf / %17.12lf)\n\n", 3*j + k, grad[3*j + k], ene_fwd_diff[3*j + k], ene_bck_diff[3*j + k], ang_to_bohr( ssize * 2 ) );
			}
		} // gradient
		fprintf( stdout, "\n" );
		// FILE OUT END


		for(j = 0;j < a.natom;j++) { for(k = 0;k < 3;k++) { mol[j].SetCrd( k, mol[j].GetCrd( k ) - ang_to_bohr( alpha ) * grad[3*j + k] ); } } // update

		for(j = 0;j < a.natom;j++) { for(k = 0;k < 3;k++) { delete [] mols_fwd_diff[3*j + k], mols_bck_diff[3*j + k]; } }
		delete [] ene_fwd_diff, ene_bck_diff, mols_fwd_diff, mols_bck_diff, grad;
	} // optimize


	AppPES_Free( &a );
	delete [] mol;
	
	return 0;
}
