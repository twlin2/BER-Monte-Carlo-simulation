//This is a MD program for periodic system.

#include <iostream>    //header for basic io
#include <cmath>       //header for math functions
#include <fstream>     //header for file io
#include <cstdlib>     //header for srand
 
using namespace std;

//This is the main driver of a MD code
class MODEL
{
   public:
   MODEL();
   ~MODEL();
   //model description
   int npart;               //number of particles
   double (*x)[3];      //position of particles
   double (*v)[3];      //velocity
   double (*f)[3];       //force
   double (*m);         //mass
   double cell[3];      //unit cell length
   bool   period;       //periodicity 0:nonperiodic 1:periodic

   //MD parameters
   double dt;             //time step
   double t;               //current time
   double Tset;         //temperature used in velocity initialization
   int nstep;              //total steps to run  
   int istep;               //current step 
   double Rcut;         //cutoff distance for vdw interaction
   double rc2;            //Rcut*Rcut;
   bool   nvt;             //flag for nvt dynamics
   
   //MC parameters
   int icycl; //currnet MC step
   int ncycl; //total MC step
   double delx; //step size in a MC move
   double beta; //1/kT
   int accpt; //numder of accepted move

   //force field parameters
   double eps;  // well depth of LJ potential
   double Ro;  // sigma 
   double p0;  //p0= 12*eps*pow(Ro,12)
   double p1;  //p1= 12*eps*pow(Ro,6)
   
   //model properties
   double T;                 //temperature
   int df;                       //degree of freedom
   double Ek;               //kinetic energy
   double Ep;               //potential energy
   double Etot;             //Ek+Ep
   double P;                 //pressure
   double V;                 //volume      (in sigma^3)
   double Ptail;             //tail correction for pressure
   double Eptail;           //tail correction for P.E.
   double strs[6]; //stress tensor


   //class functions
   int init();                    //initialization of variables
   int sample();             //calculation of properties
   double myrand01(); //generate a random number between 0 and 1 [0,1]
   int dump_trj();    //output trajectory data to a file
   int mcmov();			//MC move
   int ener(int, double &); //calculate change in potential energy due to particle displacement
   int cal_eng();		// 
   ofstream outf;     //file stream of the trajectory
   ofstream logf;     //file stream for the log
};


int main()
{
    MODEL model;                                 //declare a system
      
    model.init();                                       //initialization
    model.sample();
	        
    while(model.icycl<model.ncycl) {  //MC loop
        model.mcmov();
        //model.integrate();                         //integrate equation of motion
        if (model.icycl>3000000){
			if(model.icycl%1000==0) model.sample();                            //determine system property
    	}
	}
    
    cout << "All done\n";
       
    return 0;
}

int MODEL::init()
{   
	double rho = 0.78;		// density 
    istep=0;                                  //current step
    npart=500;                              //number of particles
    Tset=0.9;                                //target temperature
    dt=0.01;                                 //time step
    nvt=1;                           //true for NVT simulations
    
    //MC parameter init
	icycl = 0;
	ncycl = 5000000;
	delx = 0.2; //
	beta = 1.0/(Tset);
	accpt = 0;

    //allocation of memerory    
    x=new double [npart][3];        //position in reduced unit
    v=new double [npart][3];         //velocity in reduced unit
    f=new double [npart][3];          //force in reduced unit
    m=new double [npart];             //mass in reduced unit

    //give force field parameters
    eps= 1.0;			// epsilon: energy unit
    Ro= 1.0;			// bead radius (sigma): length unit
    p0= 12*eps*pow(Ro,12);
    p1= 12*eps*pow(Ro,6);

	
	double box_len = pow(double(npart)/rho, 1.0/3.0);
    cell[0]=cell[1]=cell[2]=box_len;		//length of unit cell in Ro
    period=1;                                //flag for periodicity
    Rcut=3.0;                          //cutoff distance for vdw interaction
    rc2=Rcut*Rcut;                       //Rcut2

    //assign mass of particles
    int i,j,k;
    for(i=0;i<npart;i++) m[i]=1.0; 
    
    //One simple way to place particles in space
    int nt=0;

    if(period) { //periodic system
      //A rough method to place particles in the box
      double delta[3];
      int pps;
      pps= int(pow(npart,1.0/3.0))+1;            //particles per side
      for(k=0;k<3;k++) delta[k]=cell[k]/pps; //spacing of particles
      for(i=0;i<pps;i++) {
          for(j=0;j<pps;j++) {
              for(k=0;k<pps;k++) {
                  x[nt][0]=i*delta[0];
                  x[nt][1]=j*delta[1];
                  x[nt][2]=k*delta[2];
                  nt++;                                       //number of particles placed
                  if(nt==npart) i=j=k=pps;        //nt has reached npart, exit loops
              }
          }
      }    
    } 
    else { // previous code for non-periodic system
      double sep=0.8;  //separation       
	  for(i=0;i<npart;i++) {
          for(j=0;j<=i;j++) {
             for(k=0;k<=j;k++) {
                x[nt][0]=i*sep;
                x[nt][1]=j*sep;
                x[nt][2]=k*sep;
                nt++;
                if(nt==npart) i=j=k=npart;
             }
         }
      }
    }

    outf.open("mytrj_T0.9_rho0.78.txt",ios::out); // log file: print out pressure and energy 
    logf.open("mytrj_T0.9_rho0.78.log",ios::out); // trajectory file: atom position
    
   //calculate volume
    V=cell[0]*cell[1]*cell[2];   //in A^3
    
    //Tail correction for Ep
    double Ro_Rc3=pow(Ro/Rcut,3.0);
    double Ro_Rc6=Ro_Rc3*Ro_Rc3;
    double Ro_Rc9=Ro_Rc6*Ro_Rc3;
    double Ro3=Ro*Ro*Ro;
    Eptail = (8.0/3.0)*3.1415926*npart*(npart/V)*eps*Ro3*(Ro_Rc9/3.0-Ro_Rc3);

    //Tail correction for pressure
    double unitc=1.0; // pressure unit
	Ptail = (32.0)*3.1415926*npart*npart/V/V*eps*Ro3*(Ro_Rc9/9.0-Ro_Rc3/6.0)*unitc; 
    return 0;
}

int MODEL::cal_eng(){
	    
    int i,j,k;
    
    for(i=0;i<6;i++) strs[i]=0;       //set stress to zero 
    double r2,r2i,r6i,ff,xr[3],redu;
    Ep=0;
    for(i=0;i<npart;i++) {
        for(j=i+1;j<npart;j++) {
            for(k=0;k<3;k++) xr[k]= x[i][k] - x[j][k];  //distance vector

            if(period==1) { //periodic system, find distance within one cell  
                for(k=0;k<3;k++) {
                     redu= (xr[k]/cell[k]);        //reduced coordinates
                     redu= redu - round (redu);   //between -0.5 and 0.5
                     xr[k] = redu*cell[k];        //real coordinates 
                 } 
            }                              

            r2= xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];   //distance squared
            
            if(r2<rc2 || period==0) {  //within cutoff distance

               r2i= 1/r2;
               r6i= r2i*r2i*r2i;
               ff = r6i*(4.0*p0*r6i-2.0*p1)*r2i;
			   
			   //the stress tensor 
               strs[0]+= ff*xr[0]*xr[0];  //xx in unit of kcal/mol
               strs[1]+= ff*xr[1]*xr[1];  //yy
               strs[2]+= ff*xr[2]*xr[2];  //zz
               strs[3]+= ff*xr[0]*xr[1];  //xy
               strs[4]+= ff*xr[0]*xr[2];  //xz
               strs[5]+= ff*xr[1]*xr[2];  //yz          
               //Ep += (p0*r6i - p1*2)*r6i/12.0;
               Ep += 4*eps*(r6i*r6i-r6i) ;
            }
       }
   }       
   return 0;
}

int MODEL::mcmov()
{ //attempts to displace a particle
	int k,o,n;
	double eno,enn,xo[3];
	o= int(rand()*double(npart)/double(RAND_MAX+1.0)); //a random # btwn 0 and npart-1
	ener(o,eno); //calculates interaction energy of o with other particles
	
	for(k=0;k<3;k++) {
		xo[k]=x[o][k]; //store the original coordinates
		x[o][k] = x[o][k] + (myrand01()-0.5)*delx; //give particle random movement
	}
	ener(o,enn); //interaction energy of particle o with other particles
	
	if( myrand01() > exp(-beta*(enn-eno)) ) { //reject the move
		for(k=0;k<3;k++) x[o][k] = xo[k];
	}
	else accpt++; //accpt records the number of successful attempts
	
	icycl++; //icycl records the total number of attempts
	t = icycl*dt;
	return 0;
}

int MODEL::ener(int i,double &energy)
{ //this is almost the same as the energy calculation in force()
	int j,k,nbox;
	double r2,r2i,r6i,xr[3], redu;
	energy=0;
	for(j=0;j<npart;j++) {
		if(j==i) continue;
		for(k=0;k<3;k++) xr[k]= x[i][k] - x[j][k]; //distance vector
		if(period) { //periodic system, find dist within one cell
			for(k=0;k<3;k++) {
			redu= (xr[k]/cell[k]); //reduced coordinates
			redu= redu - round (redu); //between -0.5 and 0.5
			xr[k] = redu*cell[k]; //real coordinates
			}
		}
		r2= xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2]; //distance square
		if(r2<rc2) { //within cutoff distance
			r2i= 1/r2;
			r6i= r2i*r2i*r2i;
			energy += 4*eps*(r6i*r6i-r6i);
		}
	}
	return 0;
}


int MODEL::sample()
{
    //calculation of system temperature and kinetic energy
    int i,j,k;
	cal_eng();
    Ep += Eptail;    
    Etot=Ek+Ep;


    //calcualate stress and pressure
    double unitc=1; 
    for(i=0;i<6;i++) strs[i] *= unitc/V;  
    for(i=0;i<3;i++) strs[i] += (Ek*2*unitc/(3.0*V) + Ptail) ;
    P= (strs[0]+strs[1]+strs[2])/3.0;

    //Display current information
    char null[1024];
	sprintf(null,"Current Step %d, Time: %.3f , ",icycl,t);
    cout<<null;
    logf<<null;	
    sprintf(null,"Tset %.2lf , P %.6lf , V %.0lf , ",Tset,P,V);
    cout<<null;
    logf<<null;
    sprintf(null,"Energy: Ek %.2f Eenergy per atom %.2f Etot %.2f ",Ek,Ep/npart,Etot);
    cout<<null<<endl;
    logf<<null<<endl;  
	sprintf(null, "acceptance ratio %.2f", accpt/(double)icycl);
	cout<<null<<endl;
    logf<<null;    
    
    dump_trj();
    return 0;
}


int MODEL::dump_trj()
{
    char null[1024];
    int i;
    
    outf << "ITEM: TIMESTEP" << endl;
    sprintf(null, "%.3f",t);
    outf << null << endl;
    outf << "ITEM: NUMBER OF ATOMS" << endl;
    sprintf(null, "%d",npart);
    outf << null << endl;
    outf << "ITEM: BOX BOUNDS pp pp pp" << endl;
    sprintf(null, "%f %f", 0.0,cell[0]);
    outf << null << endl;
    sprintf(null, "%f %f", 0.0,cell[1]);
    outf << null << endl;
    sprintf(null, "%f %f", 0.0,cell[2]);
    outf << null << endl;
    outf << "ITEM: ATOMS id type mass x y z\n";
    for (int i=0; i<npart; ++i){
        sprintf(null, "%d\t %d\t %d\t %lf\t %lf\t %lf\t", i+1, 1, 1, x[i][0], x[i][1], x[i][2]);
        outf << null << endl;
    }
    return 0;
}


double MODEL::myrand01()
{
    return rand()/double(RAND_MAX);
}
        
MODEL::MODEL()
{
    
};

MODEL::~MODEL()
{
    delete [] x;
    delete [] v;    
    delete [] f;    
    delete [] m;
    outf.close();
};

