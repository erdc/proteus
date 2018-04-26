// This is C:/JF/Software/Cpp/Include/jf.h
#define	ANSI
#include	<math.h>
#include	<stdio.h>
#include	<conio.h>
#include	<stdlib.h>
#include	<string.h>

#define D(x,y)			printf(" "#x":%"#y, x)
#define Diag(x,y)		if(diagon=='y') {fprintf(diag,"\n" #x"\t%"#y"\t", x); fflush(diag);}
#define Diag_Maple(x,y)	fprintf(diag,"\n\# "#x" = %"#y";", x);
#define For(x,y,z)		for((x)=(y);(x)<=(z);(x)++)

#define Write(x,y,z)	fprintf(out,"\n" #x" %"#z" ",y)
#define print(y,z)		Write( ,y,z)

#define nl(x)			fprintf((x),"\n")
#define Heading(x,y)	fprintf((x),"\n\n" #y "\n")

#define Asterisks(x)	For(i,0,72)fprintf(x,"*")
#define c(i)			putchar(i)
#define w(x)			fprintf(results,#x)
#define pf(x)			printf("%8.5g",x)
#define pd(x)			printf("%ld",x)

#define Readtext(stream,x)	fgets(x,400,stream); x[strlen(x)-1] = '\0'
#define Readtab(x)		fscanf(in,"%[^\t]%*c",x);

#define Skip(stream)	fgets(dummy,400,stream)
#define Read(stream,x,y) {fscanf(stream,"%"#y, &x);Skip(stream);}

#define endfile(x)		if(feof(x)) break

#define pause(x)		{printf("\nPaused: " x); fgetchar();}

#define IN				fscanf(in,
#define OUT				fprintf(out,
#define DIAG			fprintf(diag,
#define SCREEN			fprintf(stdout,
#define RESULTS			fprintf(results,

#define END				);

#define sign(x,y)		(fabs(x) * fabs(y)/y)
#define Sign(x)			(x > 0. ? 1. : -1.)

#define iff(x,y)		if(strcmp(x,#y)==0)
#define Screenfull(x,y)	if((x) % ((y)+1) == (y)) pause

#define pi				3.14159265358979324
#define twopi			6.2831853071795864769

int *ivector(long, long), **imatrix(long,long,long,long);
char *cvector(long, long), **cmatrix(long, long, long, long);

//complex<double> *complex_vector(long, long);

float		*vector(long , long);
float **matrix(long , long , long , long);
unsigned long *lvector(long , long );

double *dvector(long, long);
char *cvector(long, long);
double **dmatrix(long , long , long , long );
void free_ivector(int *, long, long);
void free_dvector(double *, long, long);
void free_lvector(unsigned long *, long, long);
void free_dmatrix(double **, long , long , long , long );

#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))

