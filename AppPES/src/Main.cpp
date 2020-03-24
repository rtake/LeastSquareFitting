
/* Ideal Main routine */

int main( int argc, char* argv[] ) {
	FILE *infilept;
	char jtype[256], line[256];

	AppPESInfo app;

	infilept = fopen( argv[1], "r" );
	while( fgets( line, 256, infilept ) ) {
		if( strstr( line, "#" ) ) {
			if( strstr( line, "FIT") || strstr( line, "Fit") || strstr( line, "fit") ) { sprintf( jtype, "FIT"); }
			else if( strstr( line, "MIN") || strstr( line, "Min") || strstr( line, "min") ) { sprintf( jtype, "MIN"); }
			break;
		}
	}
	fclose( infilept );

	if( strstr( jtype, "FIT" ) ) { AppPES_Fitting( app, argv[1] ); }
	else if( strstr( jtype, "MIN" ) ) { AppPES_Minimization( app, argv[1] ); }
	else return -1; // no job type

	return 0;
}
