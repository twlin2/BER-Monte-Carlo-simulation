//This is a MD program for periodic system. By Shiang-Tai Lin (stlin@ntu.edu.tw)

#include <iostream>    //header for basic io
#include <cmath>       //header for math functions
#include <fstream>     //header for file io
#include <cstdlib>     //header for srand
#include <vector>
#include <string>

using namespace std;

//For polyatomic molecules, it is better to define atoms, bonds, angles, etc separately
class ATOM
{
   public:
   ATOM(){x[0]=x[1]=x[2]=v[0]=v[1]=v[2]=m=chg=Ro=eps=type=aid=0;};
   public:
   double x[3]; //position vector of atom in angstroms
   double v[3]; //velocity vector of atom in angstroms/ps
   double f[3]; //force vector of atom in in kcal/mol/A
   double chg; //partial charges in electrons
   double m;   //atomic mass g/mol
   double Ro;  //LJ Ro parameter in Angstroms
   double eps;  //LJ Do parameter in kcal/mol
   int aid;
   int type;
};

class BOND
{
   public:
   BOND(){ len=len0=Kb=eng=xr[0]=xr[1]=xr[2]=bls[0]=bls[1]=0; };
   ATOM *atm[2];
   int bls[2];
   double xr[3];
   double len; //in angstroms
   double len0; //equilib bond length in angstroms
   double Kb;   //force constant in kcal/mol A2
   double eng;   //bond energy in kcal/mol
   double cal_len(); //calculate bond length from two vectors
   double cal_eng(); //calculate bond energy
};

//This is the main driver of a MD code
class MODEL
{
   public:
   MODEL();
   ~MODEL();
   string dir_loc;
   //model description
   int Nmol;
   double frac_free;
   int Nfree; 
   int natom;             //number of particles
   int nbond;             //number of bonds
   ATOM *atom;
   BOND *bond;   
   vector <vector<double> > neigh;	// neighbor list between A and B 
   vector <vector<int> > mspar; 	// sparse matrix to keep track of bond composition
   double cell[3];      //unit cell length (in A)
   bool   period;       //periodicity 0:nonperiodic 1:periodic

   //MD parameters
   double dt;           //time step (in ps)
   double t;            //current time (in ps)   
   double Tset;         //temperature used in velocity initialization
   int nstep;           //total steps to run  
   int istep;           //current step 
   double Rcut;         //cutoff distance for vdw interaction
   double Rc2;          //Rcut*Rcut;
   double unitc;        //conver from kcal/mol/A3 to GPa

   //MC parameters
   int seed_num;
   int icycl; //currnet MC step
   int ncycl; //total MC step
   double delx; //step size in a MC move
   double clus_delx; //step size in a MC move
   double beta; //1/kT
   int accpt; //numder of accepted move
   int bercount; // number of successful ber
   int bertotal; // number of ber MC performed

   //model properties
   double T;                //temperature (in K)
   int df;                  //degree of freedom
   double Ek;               //kinetic energy in kcal/mol
   double Ep;               //potential energy in kcal/mol
   double Evdw;				//vdw energy
   double Ebnd;				//bond energy
   double Etot;             //Ek+Ep
   double P;                //pressure    (in GPa)
   double V;                //volume      (in A3)
   double Ptail;            //tail correction for pressure (GPa)
   double Eptail;           //tail correction for P.E. (kcal/mol)
   double strs[6]; //stress tensor in kcal/mol-A3
   double tmass;   //total mass
   double cmv[3];  //center of mass velocity

   //class functions
   int init();               //initialization of variables
   int force();              //force calculation
   int f_vdw();              //vdw contribution to force
   int f_bond();             //bond contribution to force
   int integrate();          //verlet integration
   int sample();             //calculation of properties
   double myrand01(); //generate a random number between 0 and 1 [0,1]
   int dump_trj();    //output trajectory data to a file
   int mcmov();         //MC move
   int ener(int, double &); //calculate change in potential energy due to particle displacement
   int clustermove();
   int clusterener(int, double &); //calculate change in potential energy due to cluster displacement
   int rotatemove();
   int write_data();
   int buildlist();
   int updatemspar();
   int bermc();
   int enerber(int, int, double &); // calcualte change in potential energy due to bond exchange
   int printbondmatrix();


   ofstream outf;     //file stream of the trajectory
   ofstream logf;     //file stream for the log
   ofstream dataf;    //file stream of data file
   ofstream bondf;    //file stream of bond list
   ofstream neighborf;	//file stream of neighbor list
};

int main(int argc, const char * argv[]){	
	
	int berrate;
	string dir_loc;
    MODEL model;       //declare a system
    model.Nmol = atoi(argv[1]);
    berrate = atoi(argv[2]);
    model.seed_num = atoi(argv[3]);
    model.dir_loc = argv[4];
    model.frac_free = atof(argv[5]);

    model.init();      //initialization
    model.sample();    //save the initial frame   
    model.buildlist();
    model.updatemspar();
    model.printbondmatrix();
	double move_prob; 
	int c_canonical, c_cluster, c_rot;
	c_canonical = c_cluster = c_rot = 0;
    while(model.icycl<model.ncycl) {  //MC loop
    	move_prob = model.myrand01()*3.0;
    	if (move_prob<1.0){
    		model.mcmov();
    		c_canonical++;
		}
    	else if (move_prob<2.0){
    		model.clustermove();
    		c_cluster++;
		}
    	else{
    		model.rotatemove();
    		c_rot++;
		}
		if(model.icycl%berrate==0) model.bermc(); 	// frequency of performing BER MC
        if(model.icycl%1000==0){
        	model.sample();  	// frequency of outputting system properties
		}
		if(model.icycl%100==0){
			model.printbondmatrix();
		}
	}

	cout << "Simulation Done\n";
	cout << "Canonial moves = " << c_canonical << ", Cluster moves = " << c_cluster << ", Rotation moves = " << c_rot << endl; 
	cout << "Successful BER = " << model.bercount << endl;
	model.write_data();
    return 0;
}

int MODEL::init()
{ 
	srand(seed_num);
    int i,j,k;
    string fname;
	double rho = 0.3;		// density 
	Nfree = Nmol*frac_free;
    istep=0;                                  //current step
    natom=Nmol*2+Nfree;                              //number of particles
    Tset=0.9;                                //target temperature
    dt=0.01;                                 //time step
    nbond=Nmol;                   
	neigh.assign(Nmol, vector<double>(Nmol+Nfree, 0));
	mspar.assign(Nmol, vector<int>(Nmol+Nfree, 0));
    atom=new ATOM [natom];
    bond=new BOND [nbond];
	
	//MC parameter init
	icycl = 0;
	ncycl = 100000;
	delx = 0.2; //
	clus_delx = 0.2;
	beta = 1.0/(Tset);
	accpt = 0;
	bercount = 0;
	bertotal = 0;

    double Ro=1.0;    // in Angstrom
    double eps=1.0;   // in kcal/mol
    
    for(i=0;i<natom;i++) {
    	atom[i].eps=eps;    // in kcal/mol
    	atom[i].Ro=Ro;    // in Angstrom
    	atom[i].chg=0;
    	atom[i].m=1.0;      // g/mol
	}
 
    for(i=0;i<nbond;i++) {
    	bond[i].len0=1.11; // angstrom
    	bond[i].Kb=1400;   // kcal/mol A2
	}
 	double box_len = pow(double(natom)/rho, 1.0/3.0);
    cell[0]=cell[1]=cell[2]=box_len;    //length of unit cell in Angstroms
    period=1;                      //flag for periodicity
    Rcut=3.0;                        //cutoff distance for vdw interaction
    Rc2=Rcut*Rcut;                 //Rcut2
	
    int nt=0;
	int testflag = 0;
	double theta, phi, lbond, vecx, vecy, vecz, dx, dy, dz;
	lbond = 1.11;

	if (period){
		for (int nb=0; nb<nbond; ++nb){
			for (int nm=0; nm<2; ++nm){
				testflag = 0;
				while(testflag==0){
					testflag = 1;
					if (nm==0){
						atom[nb].x[0] = box_len*myrand01()-0.5*box_len;
						atom[nb].x[1] = box_len*myrand01()-0.5*box_len;
						atom[nb].x[2] = box_len*myrand01()-0.5*box_len;
						atom[nb].type = 1;
						atom[nb].aid = nb;
						bond[nb].bls[nm] = nb;
						bond[nb].atm[0]=&atom[nb];
    					
						if (nb!=0){
							for (int k=0; k<(nb); ++k){
		                        dx = atom[k].x[0] - atom[nb].x[0];
		                        dy = atom[k].x[1] - atom[nb].x[1];
		                        dz = atom[k].x[2] - atom[nb].x[2];
		                        dx -= box_len*round(dx/box_len);
		                        dy -= box_len*round(dy/box_len);
		                        dz -= box_len*round(dz/box_len);
		                        if(dx*dx+dy*dy+dz*dz<1.0){
		                            testflag = 0;
		                        }
		                        
		                        dx = atom[k+Nmol].x[0] - atom[nb].x[0];
		                        dy = atom[k+Nmol].x[1] - atom[nb].x[1];
		                        dz = atom[k+Nmol].x[2] - atom[nb].x[2];
		                        dx -= box_len*round(dx/box_len);
		                        dy -= box_len*round(dy/box_len);
		                        dz -= box_len*round(dz/box_len);
		                        if(dx*dx+dy*dy+dz*dz<1.0){
		                            testflag = 0;
		                        }
		                    }
						}
					}
					else{
						testflag = 1;
	                    theta = myrand01()*2.0*3.14158; //theta ranges from 0 to 2pi
	                    phi = acos(2.0*myrand01()-1.0); //phi ranges from 0 to pi
	                    vecx = atom[nb].x[0] + lbond*cos(theta)*sin(phi);
	                    vecy = atom[nb].x[1] + lbond*sin(theta)*sin(phi);
	                    vecz = atom[nb].x[2] + lbond*cos(phi);
	                    vecx -= box_len*round(vecx/box_len);
	                    vecy -= box_len*round(vecy/box_len);
	                    vecz -= box_len*round(vecz/box_len);
						
						/// Prevent bead in the same chain from overlapping (Do we really need this step if we are using soft potential first?)
	                    for (int k=0; k<((nbond-1)+nb+nm); ++k){
	                        dx = atom[k].x[0] - vecx;
	                        dy = atom[k].x[1] - vecy;
	                        dz = atom[k].x[2] - vecz;
	                        dx -= box_len*round(dx/box_len);
	                        dy -= box_len*round(dy/box_len);
	                        dz -= box_len*round(dz/box_len);
	                        if(dx*dx+dy*dy+dz*dz<1.0){
	                            testflag = 0;
	                        }
	                        else{
	                            atom[(nbond-1)+nb+nm].x[0] = vecx;
	                            atom[(nbond-1)+nb+nm].x[1] = vecy;
	                            atom[(nbond-1)+nb+nm].x[2] = vecz;
	                            atom[(nbond-1)+nb+nm].type = 2;
	                            atom[(nbond-1)+nb+nm].aid = (nbond-1)+nb+nm;
								bond[nb].bls[nm] = (nbond-1)+nb+nm;
								bond[nb].atm[1]=&atom[(nbond-1)+nb+nm];
	                        }
	                    }
					}
				}
			}
		}
		
		for (int nf=0; nf<Nfree; ++nf){
			testflag = 0;
			while(testflag==0){
				testflag = 1;
				atom[Nmol*2+nf].x[0] = box_len*myrand01()-0.5*box_len;
				atom[Nmol*2+nf].x[1] = box_len*myrand01()-0.5*box_len;
				atom[Nmol*2+nf].x[2] = box_len*myrand01()-0.5*box_len;
				atom[Nmol*2+nf].type = 3;
				atom[Nmol*2+nf].aid = Nmol*2+nf;
				for (int k=0; k<(Nmol*2+nf); ++k){
                    dx = atom[k].x[0] - atom[Nmol*2+nf].x[0];
                    dy = atom[k].x[1] - atom[Nmol*2+nf].x[1];
                    dz = atom[k].x[2] - atom[Nmol*2+nf].x[2];
                    dx -= box_len*round(dx/box_len);
                    dy -= box_len*round(dy/box_len);
                    dz -= box_len*round(dz/box_len);
                    if(dx*dx+dy*dy+dz*dz<1.0){
                        testflag = 0;
                    }
                }
			}
			
		}
	}

	fname = "./" + dir_loc + "/seed_" + to_string(seed_num) + "_f_" + to_string(frac_free) + "_trj_ber.txt";
    outf.open(fname,ios::out);
    fname = "./" + dir_loc + "/seed_" + to_string(seed_num) + "_f_" + to_string(frac_free) + "_thermo_ber.log";
    logf.open(fname,ios::out);
    fname = "./" + dir_loc + "/seed_" + to_string(seed_num) + "_f_" + to_string(frac_free) + "_data_ber.txt";
    dataf.open(fname,ios::out);
    fname = "./" + dir_loc + "/seed_" + to_string(seed_num) + "_f_" + to_string(frac_free) + "_bondlist.txt";
    bondf.open(fname,ios::out);

   //calculate volume
    V=cell[0]*cell[1]*cell[2];   //in A^3
    
    //Tail correction for Ep, the following is true for systems containing same atoms
    double Ro_Rc3=pow(Ro/Rcut,3.0);
    double Ro_Rc6=Ro_Rc3*Ro_Rc3;
    double Ro_Rc9=Ro_Rc6*Ro_Rc3;
    double Ro3=Ro*Ro*Ro;
    Eptail = (8.0/3.0)*3.1415926*natom*(natom/V)*eps*Ro3*(Ro_Rc9/3.0-Ro_Rc3);

    unitc=1.0; //conver from kcal/mol/A3 to GPa
	Ptail = (32.0)*3.1415926*natom*natom/V/V*eps*Ro3*(Ro_Rc9/9.0-Ro_Rc3/6.0)*unitc; 
    force();
    return 0;
}

int MODEL::buildlist()
{
	double dx, dy, dz, dr;
	for (int i=0; i<Nmol; ++i){
		for (int j=0; j<(natom-Nmol); ++j){
			dx = atom[i].x[0] - atom[j+Nmol].x[0];
			dy = atom[i].x[1] - atom[j+Nmol].x[1];
			dz = atom[i].x[2] - atom[j+Nmol].x[2];
			dx -= cell[0]*round(dx/cell[0]);
            dy -= cell[1]*round(dy/cell[1]);
            dz -= cell[2]*round(dz/cell[2]);
			dr = pow(dx*dx+dy*dy+dz*dz, 0.5);
			neigh[i][j] = dr;
		}
	}
	
	return 0;
}

int MODEL::updatemspar()
{
	int a1, a2; 
	char null[1024];
	// initialisze sparse matrix
	for (int i=0; i<Nmol; ++i){
		for (int j=0; j<(natom-Nmol); ++j){
			mspar[i][j] = 0;
		}
	}	
	
	for (int nb=0; nb<nbond; ++nb){
		a1 = bond[nb].bls[0];
		a2 = bond[nb].bls[1];
		mspar[a1][a2-Nmol] = 1;
	}
	
	return 0;
}

int MODEL::printbondmatrix()
{
	int a1, a2; 
	char null[1024];
	
	for (int i=0; i<Nmol; ++i){
		for (int j=0; j<(natom-Nmol); ++j){
			sprintf(null, "%d\t", mspar[i][j]);
			bondf << null;
		}
		sprintf(null, "\n");
		bondf << null;
	}	
	sprintf(null, "\n");
	bondf << null;
	return 0;
}

int MODEL::force()
{
    //This function determines the net force on each particle

    int i,j,k;
    for(i=0;i<natom;i++) {
        for(k=0;k<3;k++) atom[i].f[k]=0;   //set forces to zero
    }
    for(i=0;i<6;i++) strs[i]=0;       //set stress to zero 
	Ep=Evdw=Ebnd=0;
		
	f_vdw();
	f_bond();
	return 0;
}

int MODEL::f_vdw()
{  
    double r2,r2i,r6i,ff,xd[3],redu, Ro,eps;
    int i, j, k;
    double p0, p1;  
    
    for(i=0;i<natom;i++) {
        for(j=i+1;j<natom;j++) {
			if(!((i%2==0)&&(j==(i+1)))){		      	
	            for(k=0;k<3;k++) xd[k]= atom[i].x[k]-atom[j].x[k];  //distance vector
	            if(period==1) { //periodic system, find distance within one cell
					for(k=0;k<3;k++) {
						redu= (xd[k]/cell[k]); //reduced coordinates
						redu= redu - int (redu); //between -1 and 1
						if (redu>0.5) redu--; //between -1 and 0.5
						if (redu<-0.5) redu++; //between -0.5 and 0.5
						xd[k] = redu*cell[k]; //real coordinates
					}
				}
	            r2= xd[0]*xd[0]+xd[1]*xd[1]+xd[2]*xd[2];   //distance squared
	            
		       	Ro=0.5*(atom[i].Ro+atom[j].Ro);
		       	eps=sqrt(atom[i].eps*atom[j].eps);
		       	p0= 12*eps*pow(Ro,12); //in unit kcal/mol*A12
    			p1= 12*eps*pow(Ro,6); //in unit kcal/mol*A6
	            
	            if(r2<Rc2 || period==0){  //within cutoff distance
		            r2i= 1/r2;
		            r6i= r2i*r2i*r2i;
		            ff = r6i*(4.0*p0*r6i-2.0*p1)*r2i;
		            for(k=0;k<3;k++) {
		                 atom[i].f[k]+= ff*xd[k];
		                 atom[j].f[k]-= ff*xd[k];  //Newton's 3rd law
		            }

					
					//the stress tensor
					strs[0]+= ff*xd[0]*xd[0]; //xx in unit of kcal/mol
					strs[1]+= ff*xd[1]*xd[1]; //yy
					strs[2]+= ff*xd[2]*xd[2]; //zz
					strs[3]+= ff*xd[0]*xd[1]; //xy
					strs[4]+= ff*xd[0]*xd[2]; //xz
					strs[5]+= ff*xd[1]*xd[2]; //yz  

	            Evdw += 4*eps*(r6i*r6i-r6i) ;
	         	}
        	}
		}	
			
    }     
    
   return 0;
}

int MODEL::f_bond()
{  
	int i, k;
	double ff;
	double r2, redu, disp[3];
	Ebnd=0;
    for(i=0;i<nbond;i++) {
    	
		for(k=0;k<3;k++){
			disp[k]=bond[i].atm[0]->x[k]-bond[i].atm[1]->x[k];	
			redu= (disp[k]/cell[k]); //reduced coordinates
			redu= redu - round (redu); //between -0.5 and 0.5
			disp[k] = redu*cell[k]; //real coordinates
		}
	    r2= disp[0]*disp[0]+disp[1]*disp[1]+disp[2]*disp[2];
	    bond[i].len = sqrt(r2);

        for(k=0;k<3;k++) {
        	ff = - bond[i].Kb*(1-bond[i].len0/bond[i].len);
	        atom[bond[i].bls[0]].f[k]+= ff*disp[k];
	        atom[bond[i].bls[1]].f[k]-= ff*disp[k];  //Newton's 3rd law
	    }
					
		//the stress tensor
		strs[0]+= ff*disp[0]*disp[0]; //xx in unit of kcal/mol
		strs[1]+= ff*disp[1]*disp[1]; //yy
		strs[2]+= ff*disp[2]*disp[2]; //zz
		strs[3]+= ff*disp[0]*disp[1]; //xy
		strs[4]+= ff*disp[0]*disp[2]; //xz
		strs[5]+= ff*disp[1]*disp[2]; //yz  

		bond[i].cal_eng();
	    Ebnd += bond[i].eng;
         	
        		
    } 
	   
    return 0;
}

int MODEL::sample()
{
	force();
    //calculation of system temperature and kinetic energy
    int i,j,k;
    T=Ek=0;
    for(i=0;i<natom;i++) {
       Ek += (atom[i].m*(atom[i].v[0]*atom[i].v[0]+atom[i].v[1]*atom[i].v[1]+atom[i].v[2]*atom[i].v[2]));
    }

    T = Ek/(df*8.314*0.1);  //0.1 is due to unit conversion
    Ek /= (2*4.184*100);    //in kcal/mol

    Ep = Evdw + Ebnd + Eptail;
    Etot=Ek+Ep;

    //calcualate pressure
    P=((Ek*2+strs[0]+strs[1]+strs[2])*unitc)/(3.0*V); //in GPa
    P+= Ptail;

    //Display current information
    char null[1024];
	sprintf(null,"Current Step %d, Time: %.3f ps, ",icycl,t);
    cout<<null;
    logf<<null;	
    sprintf(null,"T %.2f K, P %.2f MPa, V %.0f A3, ",Tset,P,V);
    cout<<null;
    logf<<null;
    sprintf(null,"E(kcal/mol) Ek %.2f Evdw %.2f Ebnd %.2f Etail %.2f Ep %.2f Et %.2f ",Ek, Evdw, Ebnd, Eptail, Ep, Etot);
    cout<<null<<endl;
    logf<<null<<endl;  
    sprintf(null,"Tail contribution %.0f%% in Ep %.0f%% in P ",Eptail/Ep*100,Ptail/P*100);
    cout<<null<<endl;
    logf<<null<<endl;
    sprintf(null, "acceptance ratio %.2f", accpt/(double)icycl);
    cout<<null<<endl;
    logf<<null<<endl; 
    sprintf(null, "BER acceptance ratio %.2f", bercount/(double)bertotal);
	cout<<null<<endl;
    logf<<null<<endl;   
    
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
    sprintf(null, "%d",natom);
    outf << null << endl;
    outf << "ITEM: BOX BOUNDS pp pp pp" << endl;
    sprintf(null, "%f %f", -0.5*cell[0],0.5*cell[0]);
    outf << null << endl;
    sprintf(null, "%f %f", -0.5*cell[1],0.5*cell[1]);
    outf << null << endl;
    sprintf(null, "%f %f", -0.5*cell[2],0.5*cell[2]);
    outf << null << endl;
    outf << "ITEM: ATOMS id type mass x y z\n";
    for (int i=0; i<natom; ++i){
        sprintf(null, "%d\t %d\t %d\t %lf\t %lf\t %lf\t", atom[i].aid+1, atom[i].type, 1, atom[i].x[0], atom[i].x[1], atom[i].x[2]);
        outf << null << endl;
    }
    return 0;
}

int MODEL::write_data()
{
    char null[1024];
    int i;
    
    dataf << "Bond and Atom Information\n\n";

    sprintf(null, "%d atoms",natom);
    dataf << null << endl;
    sprintf(null, "%d bonds",nbond);
    dataf << null << endl;
    dataf << "0 angles" << endl;
    dataf << "0 dihedrals" << endl;
    dataf << "0 impropers\n" << endl;
    
    dataf << "3 atom types" << endl;
    dataf << "2 bond types\n" << endl;
    
	sprintf(null, "%f %f xlo xhi", -0.5*cell[0],0.5*cell[0]);
    dataf << null << endl;
    sprintf(null, "%f %f ylo yhi", -0.5*cell[1],0.5*cell[1]);
    dataf << null << endl;
    sprintf(null, "%f %f zlo zhi\n", -0.5*cell[2],0.5*cell[2]);
    dataf << null << endl;
    dataf << "Masses\n\n";
    dataf << "1 1.000000\n";
    dataf << "2 1.000000\n";
    dataf << "3 1.000000\n\n";
    dataf << "Atoms	# atom-ID mol-ID atom-type x y z \n\n";
    for (int i=0; i<natom; ++i){
        sprintf(null, "%d\t %d\t %d\t %lf\t %lf\t %lf\t", atom[i].aid+1, 1, atom[i].type, atom[i].x[0], atom[i].x[1], atom[i].x[2]);
        dataf << null << endl;
    }
    
    dataf << "\n\nBonds\n\n"; 
	for (int i=0; i<nbond; ++i){
        sprintf(null, "%d\t %d\t %d %d", i+1, 1, bond[i].bls[0]+1, bond[i].bls[1]+1);
        dataf << null << endl;
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
    outf.close();
    logf.close();
};

double BOND::cal_len()
{
	int k;
	double r2;
	for(k=0;k<3;k++) xr[k]=atm[0]->x[k]-atm[1]->x[k];	
    r2= xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
    len= sqrt(r2);
    return len;
}

double BOND::cal_eng()
{
    eng = 0.5*Kb*(len-len0)*(len-len0);
    return eng;
}

int MODEL::bermc()
{
	vector<int> rxnA;
	vector<int> rxnB;
	int checkber=0;
	int rxnpair[2];
	int rxnflag, o, m_old, n_old, m_new, n_new; // o: selected molecule, m: selected type3 particle in reaction zone, n: original particle bonded to type1 particle
	double mindist, whichpair;
	double Rmax = 1.5;	// reaction zone
	double Rmin = 0.9;
	double eno, enn;
	// 1. For each type1 particle, check if there is any type3 particle in its reaction zone
	for (int i=0; i<Nmol; ++i){
		mindist = 100.0; 
		rxnflag = 0;
		for (int j=0; j<(natom-Nmol); ++j){
			if (atom[Nmol+j].type==3){	// If this atom does not belong to any dumbbell
				// If this atom is in the reaction zone and its distance to type1 particle is the shortest
				if ((neigh[i][j]>Rmin) && (neigh[i][j]<Rmax) && (neigh[i][j]<mindist)){	
					rxnflag = 1;
					rxnpair[0] = i;
					rxnpair[1] = j+Nmol;
					mindist = neigh[i][j];
				}
			}
		}
		if (rxnflag==1){
			rxnA.push_back(rxnpair[0]);
			rxnB.push_back(rxnpair[1]);
		} 
	}
	
	// 2. For each reaction pair, perform bond exchange MC
	
	for (int i=0; i<rxnA.size(); ++i){
		checkber++;
		if (atom[rxnB[i]].type != 2){
			o = rxnA[i];
			m_old = rxnB[i];
			n_old = bond[o].bls[1];

			// 3. Perform bond exchange reaction MC 
			// 3-1 Calculate energy before the bond exchange
			enerber(o, m_old, eno);
			// 3-2 Exchange bond
			m_new = n_old;
			n_new = m_old;
			bond[o].bls[1] = m_old;
			atom[m_old].type = 2;
			atom[n_old].type = 3;

			// 3-3 Calculate energy after bond exchange
			enerber(o, m_new, enn);

			
			// 3-4 Accept bond exchange MC according to energy difference
			if( myrand01() > exp(-beta*(enn-eno)) ) { //reject the move

		      bond[o].bls[1] = n_old;
				atom[m_old].type = 3;
				atom[n_old].type = 2;

		    }
		    else bercount++;
		    bertotal++;
		}
	}
	buildlist();
    updatemspar();
	return 0;
}

int MODEL::enerber(int i, int j, double &energy)
{
	double Epair, Ebond;
	double r2,r2i,r6i,xr[3], redu;
	double p00, p10, Ro0, eps0;
	double len; 
	// Epair between type 1 and type 3
	// Ebond between type 1 and type 2
	for(int k=0;k<3;k++) xr[k]= atom[bond[i].bls[0]].x[k] - atom[j].x[k]; //distance vector
	if(period) { //periodic system, find dist within one cell
		for(int k=0;k<3;k++){
			redu= (xr[k]/cell[k]); //reduced coordinates
			redu= redu - round (redu); //between -0.5 and 0.5
			xr[k] = redu*cell[k]; //real coordinates
		}	
	}
	r2= xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2]; //distance square
	if(r2<Rc2) { //within cutoff distance
		r2i= 1/r2;
		r6i= r2i*r2i*r2i;
		Ro0=0.5*(atom[i].Ro+atom[j].Ro);
	   	eps0=sqrt(atom[i].eps*atom[j].eps);
	   	p00= 12*eps0*pow(Ro0,12); //in unit kcal/mol*A12
		p10= 12*eps0*pow(Ro0,6); //in unit kcal/mol*A6
		Epair = 4*eps0*(r6i*r6i-r6i);
	}
	
	for(int k=0;k<3;k++) xr[k]= atom[bond[i].bls[0]].x[k] - atom[bond[i].bls[1]].x[k];
		for(int k=0;k<3;k++){
			redu= (xr[k]/cell[k]); //reduced coordinates
			redu= redu - round (redu); //between -0.5 and 0.5
			xr[k] = redu*cell[k]; //real coordinates
		}	
		r2= xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2]; //distance square
		len= sqrt(r2);
		Ebond = 0.5*bond[i].Kb*(len-bond[i].len0)*(len-bond[i].len0);
	energy = Epair + Ebond;
	return 0;
}

int MODEL::mcmov()
{ //attempts to displace a particle
    int k,o,n;
    double eno, enn, xo[3];
    o= int(rand()*double(natom)/(RAND_MAX+1.0)); //a random # btwn 0 and npart-1
    ener(o,eno); //calculates interaction energy of o with other particles
    for(k=0;k<3;k++) {
        xo[k]=atom[o].x[k]; //store the original coordinates
        atom[o].x[k] = atom[o].x[k] + (myrand01()-0.5)*delx; //give particle random movement
        atom[o].x[k] -= cell[k]*round(atom[o].x[k]/cell[k]);
    }
    ener(o,enn); //interaction energy of particle o with other particles
    if( myrand01() > exp(-beta*(enn-eno)) ) { //reject the move
        for(k=0;k<3;k++) atom[o].x[k] = xo[k];
    }
    else accpt++; //accpt records the number of successful attempts
    
    icycl++; //icycl records the total number of attempts
    t = icycl*dt;
    return 0;
}

int MODEL::ener(int i,double &energy)
{ //this is almost the same as the energy calculation in force()
	int j,k,nbox;
	double r2,r2i,r6i,xr[3], redu, disp[3];
	double p00, p10, Ro0, eps0;
	double len;  
	
	energy=0;
	
	/// Calculate bond energy between particle i and it's bonded neighbor
	// 1. Check if this atom is bonded to another atom
	int inbondflag=0;
	int bpairid=0;
	int bid=0;
	for (int nb=0; nb<Nmol; ++nb){
		if (bond[nb].bls[0]==i) {
			inbondflag=1; //The chosen atom is type 1
			bpairid = bond[nb].bls[1];
			bid=nb;
		}
		if (bond[nb].bls[1]==i) {
			inbondflag=2; //The chosen atom is type 2
			bpairid = bond[nb].bls[0];
			bid=nb;
		}
	}

	
	// 2. If it belongs to a dumbbell, then calculate the difference in bond energy
	if (inbondflag!=0){

		for(k=0;k<3;k++){
			xr[k]=atom[bpairid].x[k] - atom[i].x[k];
			redu= (xr[k]/cell[k]); //reduced coordinates
			redu= redu - round (redu); //between -0.5 and 0.5
			xr[k] = redu*cell[k]; //real coordinates
		}
	    r2= xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
	    len= sqrt(r2);
	    energy += 0.5*bond[bid].Kb*(len-bond[bid].len0)*(len-bond[bid].len0);
	}
		/// Calculate pair potential between particle i (chosen particle) and other particles  
	for(j=0;j<natom;j++) {
		if(j==i) continue;
		else if ((inbondflag==1) && (j==bpairid)) continue;
		
		for(k=0;k<3;k++) xr[k]= atom[i].x[k] - atom[j].x[k]; //distance vector
		if(period) { //periodic system, find dist within one cell
			for(k=0;k<3;k++) {
				redu= (xr[k]/cell[k]); //reduced coordinates
				redu= redu - round (redu); //between -0.5 and 0.5
				xr[k] = redu*cell[k]; //real coordinates	
			}
			
		}
		r2= xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2]; //distance square
		if(r2<Rc2) { //within cutoff distance
			r2i= 1/r2;
			r6i= r2i*r2i*r2i;
			Ro0=0.5*(atom[i].Ro+atom[j].Ro);
		   	eps0=sqrt(atom[i].eps*atom[j].eps);
		   	p00= 12*eps0*pow(Ro0,12); //in unit kcal/mol*A12
			p10= 12*eps0*pow(Ro0,6); //in unit kcal/mol*A6
			energy += 4*eps0*(r6i*r6i-r6i);
		}
	}
	
	
	return 0;
}

int MODEL::clustermove()
{ //attempts to displace a particle
    int k,o,n;
    double eno, enn, xo_firstatm[3], xo_secondatm[3];
    double mol_del;
    o= int(rand()*double(nbond)/(RAND_MAX+1.0)); //a random # btwn 0 and nbond-1
    clusterener(o,eno); //calculates interaction energy of o with other particles
    for(k=0;k<3;k++) {
        //store the original coordinates
        xo_firstatm[k] = atom[bond[o].bls[0]].x[k];
        xo_secondatm[k] = atom[bond[o].bls[1]].x[k];
        mol_del = (myrand01()-0.5)*clus_delx;
		atom[bond[o].bls[0]].x[k] = atom[bond[o].bls[0]].x[k] + mol_del; //give particle random movement
		atom[bond[o].bls[0]].x[k] -= cell[k]*round(atom[bond[o].bls[0]].x[k]/cell[k]);
		atom[bond[o].bls[1]].x[k] = atom[bond[o].bls[1]].x[k] + mol_del; //give particle random movement	
		atom[bond[o].bls[1]].x[k] -= cell[k]*round(atom[bond[o].bls[1]].x[k]/cell[k]);	
    }
    clusterener(o,enn); //interaction energy of particle o with other particles
    if( myrand01() > exp(-beta*(enn-eno)) ) { //reject the move
        for(k=0;k<3;k++){
        	atom[bond[o].bls[0]].x[k] = xo_firstatm[k];
			atom[bond[o].bls[1]].x[k] = xo_secondatm[k];
		}
    }
    else accpt++; //accpt records the number of successful attempts
    
    icycl++; //icycl records the total number of attempts
    t = icycl*dt;
    return 0;
}

int MODEL::clusterener(int i,double &energy) // i: bumber of bond
{ //this is almost the same as the energy calculation in force()
	int j,k,nbox, batm;
	double r2,r2i,r6i,xr[3], redu;
	double p00, p10, Ro0, eps0;
	double len;  
	
	energy=0;
	
	/// Calculate pair potential between particle i (chosen particle) and other particles  
	for(batm=0;batm<2;batm++){
		for(j=0;j<natom;j++) {
			if( (j==(bond[i].bls[0])) || (j==(bond[i].bls[1])) )  continue; // do not calculate LJ energy within dumbbell
			for(k=0;k<3;k++) xr[k]= atom[bond[i].bls[batm]].x[k] - atom[j].x[k]; //distance vector
			if(period) { //periodic system, find dist within one cell
				for(k=0;k<3;k++) {
					redu= (xr[k]/cell[k]); //reduced coordinates
					redu= redu - round (redu); //between -0.5 and 0.5
					xr[k] = redu*cell[k]; //real coordinates
				}
				
			}
			r2= xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2]; //distance square
			if(r2<Rc2) { //within cutoff distance
				r2i= 1/r2;
				r6i= r2i*r2i*r2i;
				Ro0=0.5*(atom[bond[i].bls[batm]].Ro+atom[j].Ro);
			   	eps0=sqrt(atom[bond[i].bls[batm]].eps*atom[j].eps);
			   	p00= 12*eps0*pow(Ro0,12); //in unit kcal/mol*A12
				p10= 12*eps0*pow(Ro0,6); //in unit kcal/mol*A6
				energy += 4*eps0*(r6i*r6i-r6i);
			}
		}
	}
	
	return 0;
}

int MODEL::rotatemove()
{ //attempts to displace a particle
    int k,o,n;
    double eno, enn, xo_firstatm[3], xo_secondatm[3], com[3], vec_first[3], vec_second[3];
    double mol_del;
    double alpha, beta, gamma;
    alpha = myrand01()*2.0*3.14158;		// from 0 to 2pi
    beta = acos(2.0*myrand01()-1.0);	// from 0 o pi
    gamma = myrand01()*2.0*3.14158;		// from 0 to 2pi
    double rotm[3][3];
	rotm[0][0] = cos(alpha)*cos(gamma) - cos(beta)*sin(alpha)*sin(gamma);
	rotm[0][1] = sin(alpha)*cos(gamma) + cos(beta)*cos(alpha)*sin(gamma);
	rotm[0][2] = sin(beta)*sin(gamma); 
	rotm[1][0] = -cos(alpha)*sin(gamma) - cos(beta)*sin(alpha)*cos(gamma);
	rotm[1][1] = -sin(alpha)*sin(gamma) + cos(beta)*cos(alpha)*cos(gamma);
	rotm[1][2] = sin(beta)*cos(gamma); 
	rotm[2][0] = sin(beta)*sin(alpha);
	rotm[2][1] = -sin(beta)*cos(alpha);
	rotm[2][2] = cos(beta);
    
    o= int(rand()*double(nbond)/(RAND_MAX+1.0)); //a random # btwn 0 and nbond-1
    clusterener(o,eno); //calculates interaction energy of o with other particles
	
	double dis_bond[3], atom2_uwrap[3];
	for(k=0;k<3;k++){
		//store the original coordinates
		xo_firstatm[k] = atom[bond[o].bls[0]].x[k];
        xo_secondatm[k] = atom[bond[o].bls[1]].x[k];	 
		
		// Calculate distance vector between atom1 and atom2
		dis_bond[k] = atom[bond[o].bls[0]].x[k] - atom[bond[o].bls[1]].x[k];
		// Correct this vector with pbc
		dis_bond[k] -= cell[k]*round(dis_bond[k]/cell[k]);
		// Directly add distance vector to atom1 to get atom2 position (now atom2 can be outside of box)
		atom2_uwrap[k] = atom[bond[o].bls[0]].x[k] + dis_bond[k];
		// Calculate center of mass (can also be outside of box)
		com[k] = 0.5*(atom[bond[o].bls[0]].x[k] + atom2_uwrap[k]);
		// Calculate vectors of two atoms from center of mass
		vec_first[k] = atom[bond[o].bls[0]].x[k] - com[k];
		vec_second[k] = atom2_uwrap[k] - com[k];
	}
	
	for(k=0;k<3;k++){
		// Perform rotation based on unwrapped coordinate
		atom[bond[o].bls[0]].x[k] = com[k] + rotm[k][0]*vec_first[0] + rotm[k][1]*vec_first[1] + rotm[k][2]*vec_first[2];  //give particle random movement
		atom[bond[o].bls[1]].x[k] = com[k] + rotm[k][0]*vec_second[0] + rotm[k][1]*vec_second[1] + rotm[k][2]*vec_second[2];  //give particle random movement	
		// Wrap the atom position into pbc 
		atom[bond[o].bls[0]].x[k] -= cell[k]*round(atom[bond[o].bls[0]].x[k]/cell[k]);
		atom[bond[o].bls[1]].x[k] -= cell[k]*round(atom[o*2+1].x[k]/cell[k]);
	}

	clusterener(o,enn); //interaction energy of particle o with other particles
    if( myrand01() > exp(-beta*(enn-eno)) ) { //reject the move
        for(k=0;k<3;k++){
        	atom[bond[o].bls[0]].x[k] = xo_firstatm[k];
			atom[bond[o].bls[1]].x[k] = xo_secondatm[k];
		}
    }
    else {
		accpt++; //accpt records the number of successful attempts
	}
    icycl++; //icycl records the total number of attempts
    t = icycl*dt;
    return 0;
}






