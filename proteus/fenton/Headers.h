#define Skip		fgets(dummy,100,Input1)
#define pi			3.14159265358979324
#define iff(x,y)	if(strcmp(x,#y)==0)
#define HI			"\t%10.6f"
#define LO			"\t%7.4f"
#define Readtext(x)		fgets(x,400,Input1); x[strlen(x)-1] = '\0'
#define Read(x,y)		{fscanf(Input1,"%"#y, &x);Skip;}

#ifdef Main
FILE
	*monitor, *Input1, *Input2, *Elevation, *Flowfield, *Solution;
char
	Title[100], dummy[100], Case[20], Currentname[10], Current1[10]="Euler", Current2[10]="Stokes",
	Method[100], Convergence_file[50], Points_file[50], Theory[10], Diagname[30];
#endif

#ifdef Subprograms
extern FILE
	*monitor, *Input1, *Input2, *Elevation, *Flowfield, *Solution;
extern char
	Case[], dummy[], Title[], Currentname[], Current1[], Current2[], Method[],
	Convergence_file[], Points_file[], Theory[], Diagname[];
#endif

Int

Current_criterion,
n,
Nprofiles,
ns,
nstep,
num,
number,
Points,
Surface_points,
wave;

Double

**sol,
*B,
*coeff,
*cosa,
*rhs1,
*rhs2,
*sina,
*Tanh,
*Y,
*z;

Double

Bernoulli_check,
c,
ce,
crit,
criter,
cs,
Current,
dphidt,
dudt,
dvdt,
Eta,
f,
H,
Height,
Highest,
Hoverd,
height,
kd,
ke,
L,
MaxH,
Mean_half_keta_squared,
pe,
phi,
Pressure,
psi,
pulse,
Q,
q,
R,
r,
s,
sum,
sxx,
T,
U,
u,
ub2,
ubar,
ut,
ux,
uy,
v,
vt,
vx,
vy;
