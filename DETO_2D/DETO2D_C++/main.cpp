#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <fstream>

// DEFINING OUTPUT FILES
FILE *screen;   // screen output
FILE *flog;   // output log file -- where the code tells you what it has been doing
//FILE *fwarn;   // output warnings file
FILE *fth;   // output thermo file  -- your "it" file
FILE *fdump;   // output dump file  -- your XY file, but in LAMMPS format here. Also, another file recording deformation movies at each DETO step

// DEFINING GLOBAL VARIABLES - all functions can read and edit this
// suggestion is to use local variable whenever possible and move them to here to make them global when you realise that some functions need to modify them
int nelx, nely;     // number of elements in x and y
double D, mass;           // particle diameters and masses (all the same)
int Npart;           // number of particles in the system
std::vector<double> chi;    // vector of design variables (formerly "m" in MATLAB)
std::vector<double> xi;     // vectors of initial particle positions
std::vector<double> yi;
std::vector<double> x;     // vectors of current particle positions
std::vector<double> y;
int loop;                   // number of step during DETO
int astep;                  // cumulative minimsation step during DETO to record all deformation movies in dump
double Utot;                // total interaction energy, viz. our objective function
double tol_chi;             // tolerance on chi to state convergence of DETO
std::vector<int> nn;     // per-particle number of neighbour for interactions and filtering
std::vector<int> nnf;
std::vector< std::vector<int> > Nlist;  // neighbour lists for interactions and filtering
std::vector< std::vector<int> > Nflist;
std::vector< std::vector<double> > Li;  // initial interparticle distances for interactions and filtering
std::vector< std::vector<double> > Lif;
std::vector< std::vector<double> > L;  // current interparticle distances for interactions
std::vector< std::vector<double> > dx;
std::vector< std::vector<double> > dy;
std::vector< std::vector<double> > drij;
std::vector< std::vector<double> > kpen;  // penalised k parameters of potentials
std::vector< std::vector<double> > F;  // current interaction forces
double rmin,beta;                // filtering length and beta exponent
double ppow;                // penalisation power
double kpot;                // k parameter in the interaction potential
double apot;                // a parameter in the interaction potential
std::string pottype;        // potential type
double tol_min, dmax, dtinp; //, Ftol;   // quickmin parameters
int stepAve;          // number of steps over which Utot tolerance is averaged
std::vector<int> FexID;     // ID of particles with applied external forces
std::vector<int> FeyID;
std::vector<double> Fex;    // values of applied external forces
std::vector<double> Fey;
std::vector<double> Xfix;    // vectors of supprts fixing corresponding forces
std::vector<double> Yfix;
std::vector<int> DixID;     // imposed external displacements with given rate
std::vector<int> DiyID;
std::vector<double> Dix;
std::vector<double> Diy;
std::vector<double> DixR;
std::vector<double> DiyR;
double chi_min, target_f;   // minimum chi and target volume fraction of solid
std::vector<int> thFxID;
std::vector<int> thFyID;
std::vector<int> thXID;
std::vector<int> thYID;
std::vector<double> Fx; // forces on particles
std::vector<double> Fy;
std::vector<double> sx; // stresses
std::vector<double> sy;
std::vector<double> sxy;
std::vector<double> syx;
std::vector<double> shyd;
std::vector<double> sdev;
int out_every;          // output dump and thermo every these quickmin steps
double move;    // max chi move in DETO
bool flag_Adump;    // flag to activate dumping of quickmin trajectories

// DEFINING FUNCTIONS
void read_input(std::string);   // reading input file and create initial system
void write_dump(std::string,int);            // write entry in dump file
void write_thermo(std::string,int);            // write entry in thermo file
void neighbours();          // compute neighbour lists for interactions and filtering, as well as initial lengths for those
void postprocess();         // postprocessing
void DETOrun();         // freeing arrays if any is created
void quickmin();        // energy minimization
void part_stress();     // computing per particle stresses
void read_dump(std::string,int);    // reading initial chi ditribution from input file


// MAIN PROGRAM
int main(int argc, char **argv)
{
    screen = stdout;    // setting screen output
    fprintf(screen,"\n --------------------- \n PROGRAM STARTED \n -----------------\n\n");
    
    //first argument will be the input file name. Error if not found
    std::string fname;
    if (argc < 2) {
        fprintf(screen,"ERROR: input file name (for now a dummy one) must be specified as first argument \n");
        return 0;
    }
    else {fname = argv[1];}
    
    // initiating log file
    flog = fopen ("logfile.log","w");
    fprintf(flog,"LOG FILE CREATED SUCCESFULLY \n\n");
    
    fprintf(screen,"READING INPUT FILE AND CREATING INITIAL STRUCTURE.....\n");
    read_input(fname);
    fprintf(screen,"....DONE\n\n");
    
    // initiating dump files
    fdump = fopen ("DETO.dump","w");
    fclose(fdump);
    if (flag_Adump) {
        fdump = fopen ("ALL.dump","w");
        fclose(fdump);
    }
    
    
    
    // initiating thermo file
    fth = fopen ("thermo.txt","w");
    fprintf(fth,"step chi_tot c");
    for (int i=0; i<thXID.size(); i++) fprintf(fth," x[%d]",thXID[i]);
    for (int i=0; i<thYID.size(); i++) fprintf(fth," y[%d]",thYID[i]);
    for (int i=0; i<thFxID.size(); i++) fprintf(fth," Fx[%d]",thFxID[i]);
    for (int i=0; i<thFyID.size(); i++) fprintf(fth," Fy[%d]",thFyID[i]);
    fprintf(fth,"\n");
    fclose(fth);
    fth = fopen ("thermoALL.txt","w");
    fprintf(fth,"step chi_tot c");
    for (int i=0; i<thXID.size(); i++) fprintf(fth," x[%d]",thXID[i]);
    for (int i=0; i<thYID.size(); i++) fprintf(fth," y[%d]",thYID[i]);
    for (int i=0; i<thFxID.size(); i++) fprintf(fth," Fx[%d]",thFxID[i]);
    for (int i=0; i<thFyID.size(); i++) fprintf(fth," Fy[%d]",thFyID[i]);
    fprintf(fth,"\n");
    fclose(fth);
    
    fprintf(screen,"COMPUTING NEIGHBOUR LISTS.....\n");
    neighbours();
    fprintf(screen,"....DONE\n\n");
    
    
    // Print outputs before starting run
    loop = 0;
    astep = 0;
    Utot = 0.;
    for (int i=0;i<Npart;i++) {
        x.push_back(xi[i]);
        y.push_back(yi[i]);
        Fx.push_back(0.);
        Fy.push_back(0.);
        sx.push_back(0.);
        sy.push_back(0.);
        sxy.push_back(0.);
        syx.push_back(0.);
        shyd.push_back(0.);
        sdev.push_back(0.);
    }

    write_dump("DETO.dump",loop);
    if (flag_Adump) write_dump("ALL.dump",astep);
    write_thermo("thermo.txt",loop);
    write_thermo("thermoALL.txt",astep);
    
    fprintf(screen,"RUNNING OPTIMIZATION.....\n");
    DETOrun();
    fprintf(screen,".....OPTIMIZATION DONE\n\n");
    
    fprintf(screen,"POSTPROCESSING.....\n");
    postprocess();
    fprintf(screen,".....DONE\n\n");
    
    fprintf(screen,"CLOSING OUTPUT FILES.....\n");
    fclose(flog);
    fprintf(screen,".....DONE\n");
    
    fprintf(screen,"\n---------------------------------------\n PROGRAM ENDED SUCCESFULLY \n--------------------------------------- \n");
    
    return 0;
}



void read_input(std::string fname)
{
    
    fprintf(flog,"READING INPUT FILE  %s \n--------------------------------------- \n",fname.c_str());
    
    std::string read_string,word;
    std::ifstream inFile(fname.c_str());
    
    // initialing variables to be read to default values
    nelx = -1;
    nely = -1;
    D = -1.;
    rmin = -1.;
    beta = 0.;
    ppow = 2.;
    tol_min = 1e30;
    dmax = 1e10;
    dtinp = 0.001;
    mass = -1.;
    //Ftol = 1e30;
    chi_min = 0.;
    out_every = 1000; //default
    move = 0.2;
    stepAve = 10;
    flag_Adump = false;
    
    // reading input file.
    while (!inFile.eof()){
        
        // get next line....
        std::getline (inFile, read_string);
        
        if (!read_string.empty()){
            
            // read next word in line....
            std::istringstream lss(read_string);
            bool getout = false;
            while (lss >> word && getout == false) {
                getout = false;
                // end line reading if comment found, or go through line
                if (strncmp(word.c_str(), "#", 1) == 0) break;
                else if (strcmp(word.c_str(),"geometry")==0) {
                    lss >> nelx >> nely >> D >> mass;
                    fprintf(flog,"Number of elements nelx = %d , nely = %d\n",nelx,nely);
                    fprintf(flog,"Particle diameter D = %.5f, mass = %.5f \n",D,mass);
                    // defininf geomtry accordingly
                    Npart = nelx*nely - (int)(nely/2);  // no restrictions on odd/even
                    for (int j=0; j<nely; j++){
                        for (int i=0; i<(nelx-(j%2)); i++){
                            xi.push_back( (j%2)*D/2. + D*i );    // removed symmetry around zero to simplify
                            yi.push_back( j*D*cos(M_PI/6.) );
                        }
                    }
                    fprintf(flog,"Assigned initial positions to particles\n");
                }
                else if (strcmp(word.c_str(),"init_chi")==0) {
                    lss >> word;
                    if (strcmp(word.c_str(),"uniform")==0) {
                        double fvalue;
                        lss >> fvalue;
                        for (int i=0; i<Npart; i++)  chi.push_back(fvalue);
                        fprintf(flog,"Design variable initialised with uniform value: %f\n",fvalue);
                    }
                    else if (strcmp(word.c_str(),"file")==0) {
                        int tstep;
                        lss >> word >> tstep;
                        read_dump(word.c_str(),tstep);
                    }
                    else{
                        fprintf(screen,"ERROR: unknown argument for init_chi: %s\n",word.c_str());
                        exit(0);
                    }
                }
                else if (strcmp(word.c_str(),"tol_chi")==0)  lss >> tol_chi;
                else if (strcmp(word.c_str(),"move")==0)  lss >> move;
                else if (strcmp(word.c_str(),"filtering")==0)  {
                    lss >> word;
                    if (strcmp(word.c_str(),"on")==0) lss >> rmin >> beta;
                    else break;
                }
                else if (strcmp(word.c_str(),"ppow")==0)  lss >> ppow;
                else if (strcmp(word.c_str(),"interaction")==0){
                    lss >> pottype;
                    if (strcmp(pottype.c_str(),"lin")==0) {
                        lss >> kpot;
                        getout = true;
                    }
                    else if (strcmp(pottype.c_str(),"e-x")==0 || strcmp(pottype.c_str(),"e+x")==0 || strcmp(pottype.c_str(),"sinh")==0 || strcmp(pottype.c_str(),"tanh")==0) lss >> kpot >> apot;
                    else {
                        fprintf(screen,"ERROR: unknown type of interaction potential: %s\n",pottype.c_str());
                        exit(0);
                    }
                }
                else if (strcmp(word.c_str(),"quickmin")==0)  {
                    //lss >> tol_min >> Ftol >> dmax >> dtinp >> out_every;
                    std::string dump_yn;
                    lss >> tol_min >> stepAve >> dmax >> dtinp >> out_every >> dump_yn;
                    if  (strcmp(dump_yn.c_str(),"yes")==0) flag_Adump = true;
                }
                else if (strcmp(word.c_str(),"setforce")==0)  {
                    std::string x_y;
                    lss >> x_y;
                    if (strcmp(x_y.c_str(),"x")==0)  {
                        int pid;
                        double value;
                        lss >> pid >> value;
                        FexID.push_back(pid);
                        Fex.push_back(value);
                    }
                    else if (strcmp(x_y.c_str(),"y")==0)  {
                        int pid;
                        double value;
                        lss >> pid >> value;
                        FeyID.push_back(pid);
                        Fey.push_back(value);
                    }
                    else{
                        fprintf(screen,"ERROR: imposed forces can be in direction x or y; instead I found %s\n",x_y.c_str());
                        exit(0);
                    }
                }
                else if (strcmp(word.c_str(),"fix")==0)  {
                    std::string x_y;
                    lss >> x_y;
                    if (strcmp(x_y.c_str(),"x")==0)  {
                        int pid;
                        lss >> pid;
                        Xfix.push_back(pid);
                    }
                    else if (strcmp(x_y.c_str(),"y")==0)  {
                        int pid;
                        lss >> pid;
                        Yfix.push_back(pid);
                    }
                    else{
                        fprintf(screen,"ERROR: imposed constraint can be in direction x or y; instead I found %s\n",x_y.c_str());
                        exit(0);
                    }
                }
                else if (strcmp(word.c_str(),"setdisp")==0)  {
                    std::string x_y;
                    lss >> x_y;
                    if (strcmp(x_y.c_str(),"x")==0)  {
                        int pid;
                        double dvalue, rvalue;
                        lss >> pid >> dvalue >> rvalue;
                        DixID.push_back(pid);
                        Dix.push_back(dvalue);
                        DixR.push_back(rvalue);
                    }
                    else if (strcmp(x_y.c_str(),"y")==0)  {
                        int pid;
                        double dvalue, rvalue;
                        lss >> pid >> dvalue >> rvalue;
                        DiyID.push_back(pid);
                        Diy.push_back(dvalue);
                        DiyR.push_back(rvalue);
                    }
                    else{
                        fprintf(screen,"ERROR: imposed forces can be in direction x or y; instead I found %s\n",x_y.c_str());
                        exit(0);
                    }
                }
                else if (strcmp(word.c_str(),"chi_min")==0) lss >> chi_min;
                else if (strcmp(word.c_str(),"target_f")==0) lss >> target_f;
                else if (strcmp(word.c_str(),"add_thermo")==0)  {
                    lss >> word;
                    int pid;
                    if (strcmp(word.c_str(),"x")==0) {
                        lss >> pid;
                        thXID.push_back(pid);
                    }
                    if (strcmp(word.c_str(),"y")==0) {
                        lss >> pid;
                        thYID.push_back(pid);
                    }
                    if (strcmp(word.c_str(),"Fx")==0) {
                        lss >> pid;
                        thFxID.push_back(pid);
                    }
                    if (strcmp(word.c_str(),"Fy")==0) {
                        lss >> pid;
                        thFyID.push_back(pid);
                    }
                }
                else{
                    fprintf(screen,"ERROR: unknown keyword in input file: %s\n",word.c_str());
                    exit(0);
                }
            }
        
        }
    }
    
    fprintf(flog,"Applied F in x (particle and value)\n");
    for (int i=0; i<Fex.size(); i++) fprintf(flog,"%d %e\n",FexID[i],Fex[i]);
    fprintf(flog,"\nApplied F in y (particle and value)\n");
    for (int i=0; i<Fey.size(); i++) fprintf(flog,"%d %e\n",FeyID[i],Fey[i]);
    
    fprintf(flog,"------------------------------------- \n\n");
}


void neighbours()
{
    for (int i=0; i<Npart; i++) {
        nn.push_back(0);
        nnf.push_back(0);
    }
    
    double rij;   //interparticle distance
    std::vector <int> ivecN;
    std::vector <int> fvecN;
    std::vector <double> ivecL;
    std::vector <double> fvecL;
    std::vector <double> kvec;
    std::vector <double> vecF;
    std::vector <double> vecX;
    std::vector <double> vecY;
    std::vector <double> vecR;
    
    for (int i=0; i<Npart; i++) {
        
        ivecN.clear();
        fvecN.clear();
        ivecL.clear();
        fvecL.clear();
        kvec.clear();
        vecF.clear();
        vecR.clear();
        
        // building half neighbour lists only, with j>i
        for (int j=i+1; j<Npart; j++) {
            rij = sqrt( (yi[i]-yi[j])*(yi[i]-yi[j]) + (xi[i]-xi[j])*(xi[i]-xi[j]) );
            
            // interaction list
            if (rij < 1.01*D){
                ivecN.push_back(j);
                ivecL.push_back(rij);
                kvec.push_back( pow(chi[i],ppow) * pow(chi[j],ppow) * kpot );
                vecF.push_back(0.);
                vecX.push_back(xi[i]-xi[j]);
                vecY.push_back(yi[i]-yi[j]);
                vecR.push_back(0.);   //initial bond extension assumed = 0
                nn[i]++;
            }
        }
        
        
        for (int j=0; j<Npart; j++) {
            rij = sqrt( (yi[i]-yi[j])*(yi[i]-yi[j]) + (xi[i]-xi[j])*(xi[i]-xi[j]) );
            // filtering list (full neighbour list in this case, not half)
            if (i!=j && rij < rmin){
                fvecN.push_back(j);
                fvecL.push_back(rij);
                nnf[i]++;
            }
        }
        
        Nlist.push_back(ivecN);
        Nflist.push_back(fvecN);
        Li.push_back(ivecL);
        L.push_back(ivecL);     // initialising current lengths to initial ones
        dx.push_back(vecX);
        dy.push_back(vecY);
        drij.push_back(vecR);
        Lif.push_back(fvecL);
        kpen.push_back(kvec);
        F.push_back(vecF);
    }
    fflush(flog);
}

void postprocess()
{
    
}

void write_dump(std::string dname, int dstep)
{
    fdump = fopen (dname.c_str(),"a");
    
    fprintf(fdump,"ITEM: TIMESTEP\n");
    fprintf(fdump,"%d\n",dstep);
    fprintf(fdump,"ITEM: NUMBER OF ATOMS\n");
    fprintf(fdump,"%d\n",Npart);
    fprintf(fdump,"ITEM: BOX BOUNDS ff pp ff\n");
    fprintf(fdump,"%f %f\n",0.,nelx*D);
    fprintf(fdump,"%f %f\n",0.,nely*D);
    fprintf(fdump,"%f %f\n",-D,D);
    fprintf(fdump,"ITEM: ATOMS id type x y z radius mass chi pxx pyy pxy pyx phyd pdev\n");
    
    for (int i=0; i<Npart;i++){
        fprintf(fdump,"%d %d %.6f %.6f %.6f %e %.5f %e %e %e %e %e %e %e\n",i,1,x[i],y[i],0.,D/2.,mass,chi[i],sx[i],sy[i],sxy[i],syx[i],shyd[i],sdev[i]);
    }

    fclose(fdump);
}


void write_thermo(std::string tname, int tstep)
{
    fth = fopen (tname.c_str(),"a");
    
    double totchi = 0.;
    for (int i=0; i<Npart; i++) totchi += chi[i];
        
    double Fxx, Fyy;
    std::vector<double> tempFx;
    std::vector<double> tempFy;
    
    for (int i=0; i<Npart; i++){
        tempFx.push_back(0.);
        tempFy.push_back(0.);
    }
    
    for (int i=0; i<Npart; i++){
        for (int s=0; s<nn[i]; s++){
            int j = Nlist[i][s];
            Fxx = F[i][s] * dx[i][s]/L[i][s];
            Fyy = F[i][s] * dy[i][s]/L[i][s];
            tempFx[i] += Fxx;
            tempFx[j] -= Fxx;
            tempFy[i] += Fyy;
            tempFy[j] -= Fyy;
        }
    }
    
    
    
    fprintf(fth,"%d %.6f %e",tstep,totchi,Utot);
    fflush(fth);
    for (int i=0; i<thXID.size(); i++) fprintf(fth," %f",x[thXID[i]]);
    for (int i=0; i<thYID.size(); i++) fprintf(fth," %f",y[thYID[i]]);
    for (int i=0; i<thFxID.size(); i++) fprintf(fth," %e",tempFx[thFxID[i]]);
    for (int i=0; i<thFyID.size(); i++) fprintf(fth," %e",tempFy[thFyID[i]]);
    fprintf(fth,"\n");
    
    fclose(fth);
    
}


void DETOrun()
{
    
    std::vector<double> dc; // sensitivities
    
    double dchi = 2.*tol_chi;       // same as "change" in MATLAB
    std::vector<double> ochi;   // vector of old chi values
    for (int i=0; i<Npart; i++){
        ochi.push_back(chi[i]);
        dc.push_back(0.);
    }
    
    while (dchi > tol_chi){
        loop ++;
        fprintf(screen,"DETO step %d\n",loop);
        
        // always restarting step from underformed state
        for (int i=0; i<Npart; i++){
            x[i] = xi[i];
            y[i] = yi[i];
            for (int s=0; s<nn[i]; s++) F[i][s]=0.;
            Fx[i]=0.;
            Fy[i]=0.;
            ochi[i] = chi[i];
        }
    
        // penalise kpot
        for (int i=0; i<Npart; i++){
            for (int s=0; s<nn[i]; s++){
                int j = Nlist[i][s];
                kpen[i][s] = pow(chi[i],ppow) * pow(chi[j],ppow)  * kpot;
            }
        }
        
        
        
        // write output files
        part_stress();
        if (flag_Adump)  write_dump("ALL.dump",astep);
        write_thermo("thermoALL.txt",astep);
        
        quickmin();
        
        // write output files
        part_stress();
        write_dump("DETO.dump",loop);
        write_thermo("thermo.txt",loop);
        
        // computing sensitivities (tp-to-date drij already from quickmin)
        for (int i=0; i<Npart; i++) dc[i]=0.;
        double potfac;
        
        if (strcmp(pottype.c_str(),"lin")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    potfac = 1./4.*kpot * drij[i][s] * drij[i][s];
                    int j = Nlist[i][s];
                    dc[i] -= ppow * pow(chi[i],ppow-1.)*pow(chi[j],ppow) * potfac;
                    dc[j] -= ppow * pow(chi[j],ppow-1.)*pow(chi[i],ppow) * potfac;
                }
            }
        }
        
        if (strcmp(pottype.c_str(),"e-x")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    potfac = 1./2.* ( kpot/apot * (exp(-apot*drij[i][s])/apot+drij[i][s] ) - kpot/apot/apot) ;
                    int j = Nlist[i][s];
                    dc[i] -= ppow * pow(chi[i],ppow-1.)*pow(chi[j],ppow) * potfac;
                    dc[j] -= ppow * pow(chi[j],ppow-1.)*pow(chi[i],ppow) * potfac;
                }
            }
        }
        
        if (strcmp(pottype.c_str(),"e+x")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    potfac = 1./2. * ( kpot/apot * (exp(apot*drij[i][s])/apot-drij[i][s] ) - kpot/apot/apot ) ;
                    int j = Nlist[i][s];
                    dc[i] -= ppow * pow(chi[i],ppow-1.)*pow(chi[j],ppow) * potfac;
                    dc[j] -= ppow * pow(chi[j],ppow-1.)*pow(chi[i],ppow) * potfac;
                }
            }
        }
        
        if (strcmp(pottype.c_str(),"tanh")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    potfac = 1./2. * kpot/apot/apot * log(cosh(apot*drij[i][s]));
                    int j = Nlist[i][s];
                    dc[i] -= ppow * pow(chi[i],ppow-1.)*pow(chi[j],ppow) * potfac;
                    dc[j] -= ppow * pow(chi[j],ppow-1.)*pow(chi[i],ppow) * potfac;
                }
            }
        }
        
        if (strcmp(pottype.c_str(),"sinh")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    potfac = 1./2. * ( kpot/apot/apot * cosh(apot*drij[i][s]) - kpot/apot/apot ) ;
                    int j = Nlist[i][s];
                    dc[i] -= ppow * pow(chi[i],ppow-1.)*pow(chi[j],ppow) * potfac;
                    dc[j] -= ppow * pow(chi[j],ppow-1.)*pow(chi[i],ppow) * potfac;
                }
            }
        }
 

        // Filtering of sensitivities (aka coarse graining)
        if (rmin > 0.){
            std::vector<double> dcn;
            double tot,fac;
            for (int i=0; i<Npart; i++) {
                dcn.push_back(dc[i] * pow(chi[i],beta) * rmin);
                tot = rmin;
                for (int s=0; s<nnf[i]; s++){
                    int j = Nflist[i][s];
                    fac = rmin-Lif[i][s];
                    tot += fac;
                    dcn[i] += fac*pow(chi[j],beta)*dc[j];
                }
                dcn[i] /= ( pow(chi[i],beta)*tot );
            }
            for (int i=0; i<Npart; i++) dc[i] = dcn[i];
        }
 
            
        // Updating chi while respecting constraints on volume fraction and range
        double l1 = 0.;
        double l2 = 100000.;
        double lmid;
        //double move = 0.2;
        double chi_sum;
        while (l2-l1 > 1e-8){
            chi_sum = 0.;
            lmid = 0.5*(l2+l1);
            for (int i=0; i<Npart; i++) {
                chi[i] = ochi[i] * sqrt(-dc[i]/lmid);
                if (chi[i] > ochi[i] + move) chi[i] = ochi[i] + move;
                if (chi[i] > 1.) chi[i] = 1.;
                if (chi[i] < ochi[i] - move) chi[i] = ochi[i] - move;
                if (chi[i] < chi_min) chi[i] = chi_min;
                chi_sum += chi[i];
            }
            if ( ( chi_sum - target_f*(double)Npart) > 0 ) l1 = lmid;
            else l2 = lmid;
        }
        
        dchi = 0.;
        for (int i=0; i<Npart; i++){
            if ( dchi < fabs(chi[i]-ochi[i]) ) dchi = fabs(chi[i]-ochi[i]);
        }

    }
}




void quickmin()
{
    std::vector<double> vx;
    std::vector<double> vy;
    int freeze_step = 0;    // counter of steps after freezing
    double DUtot = 2.*tol_min;  // relative change of energy to check against tolerance for convergence
    int nstep = 0;
    double vdotf;
    double Uold;
    bool disp_reached = true;
    double dt;
    double cumDU = 0.;   // cumulative of DUtot over stepAve steps
    
    for (int i=0; i<Npart; i++) {
        vx.push_back(0.);
        vy.push_back(0.);
    }
    
    astep++;
    
    if (Dix.size()>0) disp_reached=false;
    
    //double Fmax = 1.;
    int count = 0;  // counts steps to meet stepAve
    
    //while ( (DUtot > tol_min || Fmax > Ftol) || disp_reached==false){
    while ( DUtot > tol_min || disp_reached==false ){
        nstep++;
        astep++;
        dt = dtinp;
        //Fmax = 0.;
        count++;
        
        // reset particle forces
        for (int i=0; i<Npart; i++) {
            Fx[i]=0.;
            Fy[i]=0.;
        }
        
        // impose external forces
        for (int i=0; i<Fex.size(); i++) Fx[FexID[i]] = Fex[i];
        for (int i=0; i<Fey.size(); i++) Fy[FeyID[i]] = Fey[i];
        
        // update interparticle distances in neighbour list
        for (int i=0; i<Npart; i++){
            for (int s=0; s<nn[i]; s++){
                int j = Nlist[i][s];
                dx[i][s] = x[j]-x[i];
                dy[i][s] = y[j]-y[i];
                L[i][s]= sqrt( dx[i][s]*dx[i][s] + dy[i][s]*dy[i][s]);
                drij[i][s] = L[i][s] - Li[i][s];
            }
        }
        
        // compute forces
        if (strcmp(pottype.c_str(),"lin")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    F[i][s] = kpen[i][s] * drij[i][s];
                }
            }
        }
        
        if (strcmp(pottype.c_str(),"e-x")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    F[i][s] = kpen[i][s]/apot * (-exp( -apot* drij[i][s])+1.);
                }
            }
        }
        
        if (strcmp(pottype.c_str(),"e+x")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    F[i][s] = kpen[i][s]/apot * (exp( apot* drij[i][s])-1.);
                }
            }
        }
        
        if (strcmp(pottype.c_str(),"tanh")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    F[i][s] = kpen[i][s]/apot * (tanh( apot* drij[i][s]) );
                }
            }
        }
        
        if (strcmp(pottype.c_str(),"sinh")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    F[i][s] = kpen[i][s]/apot * (sinh( apot* drij[i][s]) );
                }
            }
        }
        
        double Fxx, Fyy;
        for (int i=0; i<Npart; i++){
            for (int s=0; s<nn[i]; s++){
                int j = Nlist[i][s];
                Fxx = F[i][s] * dx[i][s]/L[i][s];
                Fyy = F[i][s] * dy[i][s]/L[i][s];
                Fx[i] += Fxx;
                Fx[j] -= Fxx;
                Fy[i] += Fyy;
                Fy[j] -= Fyy;
            }
        }
        
        
        // constraints on forces, computing forces elsewhere for printing
        for (int i=0; i<Xfix.size(); i++) Fx[Xfix[i]] = 0.;
        for (int i=0; i<Yfix.size(); i++) Fy[Yfix[i]] = 0.;
        for (int i=0; i<Dix.size(); i++) Fx[DixID[i]] = 0.;
        for (int i=0; i<Diy.size(); i++) Fy[DiyID[i]] = 0.;
        
        
        //Compute scale factor
        vdotf = 0.;
        for (int i=0; i<Npart; i++) vdotf += vx[i]*Fx[i] + vy[i]*Fy[i];
        
        // Overshoot check
        if (vdotf < 0.){
            freeze_step = nstep;
            for (int i=0; i<Npart; i++){
                vx[i]=0.;
                vy[i]=0.;
            }
        }
        else{
            // Scale velocities
            double fdotf = 0.;
            double scale = 0.;
            for (int i=0; i<Npart; i++) fdotf += Fx[i]*Fx[i] + Fy[i]*Fy[i];
            if (fdotf == 0) scale = 0.;
            else scale = vdotf/fdotf;
            
            for (int i=0; i<Npart; i++){
                vx[i] = scale * Fx[i];
                vy[i] = scale * Fy[i];
            }

        }
        
        
        // Limit dmax
        double vmax = 0;
        for (int i=0; i<Npart; i++){
            if (fabs(vx[i]) > vmax)   vmax = fabs(vx[i]);
            if (fabs(vy[i]) > vmax)   vmax = fabs(vy[i]);
        }
        if (dt*vmax > dmax) dt = dmax/vmax;
    
        
        // Euler integration step
        for (int i=0; i<Npart; i++){
            x[i] += dt*vx[i];
            y[i] += dt*vy[i];
            vx[i] += dt/mass*Fx[i];
            vy[i] += dt/mass*Fy[i];
        }
    
        // Tolerance
        Uold = Utot;
        Utot = 0.;
        double potfac;
        
        if (strcmp(pottype.c_str(),"lin")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    potfac = kpen[i][s] * drij[i][s] * drij[i][s];
                    Utot += potfac;
                }
            }
            Utot *= 1./2.;
        }
        
        if (strcmp(pottype.c_str(),"e-x")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    potfac = kpen[i][s]/apot * (exp(-apot*drij[i][s])/apot+drij[i][s] ) - kpen[i][s]/apot/apot ;
                    Utot += potfac;
                }
            }
        }
        
        if (strcmp(pottype.c_str(),"e+x")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    potfac = kpen[i][s]/apot * (exp(apot*drij[i][s])/apot-drij[i][s] ) - kpen[i][s]/apot/apot ;
                    Utot += potfac;
                }
            }
        }
        
        if (strcmp(pottype.c_str(),"tanh")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    potfac = kpen[i][s]/apot/apot * log(cosh(apot*drij[i][s]));
                    Utot += potfac;
                }
            }
        }
        
        if (strcmp(pottype.c_str(),"sinh")==0) {
            for (int i=0; i<Npart; i++){
                for (int s=0; s<nn[i]; s++){
                    potfac = kpen[i][s]/apot/apot * cosh(apot*drij[i][s]) - kpen[i][s]/apot/apot;
                    Utot += potfac;
                }
            }
        }
        
        if (count == stepAve){
            if  ( (nstep - freeze_step) < 10 || nstep < 2*stepAve) DUtot = 2.*tol_min;
            else { //sum last value and average
                cumDU += (fabs(Utot-Uold))/Uold;
                DUtot = cumDU / (double)stepAve;
            }
            count = 0; // reset counter and cumulative
            cumDU = 0.;
        }
        else{
            cumDU += (fabs(Utot-Uold))/Uold;
            DUtot = 2.*tol_min;
        }
        

        // Imposing particle displacements
        disp_reached = true;
        for (int i=0; i<Dix.size(); i++){
            if ( fabs(x[DixID[i]]-xi[DixID[i]]) < fabs(Dix[i]) ){
                x[DixID[i]] += DixR[i] * dt * Dix[i]/fabs(Dix[i]);
                if ( fabs(x[DixID[i]]-xi[DixID[i]]) > fabs(Dix[i]) ){
                    x[DixID[i]] = xi[DixID[i]] + Dix[i];
                }
                else disp_reached = false;
            }
        }
        for (int i=0; i<Diy.size(); i++){
            if ( fabs(y[DiyID[i]]-yi[DiyID[i]]) < fabs(Diy[i]) ){
                y[DiyID[i]] += DiyR[i] * dt * Diy[i]/fabs(Diy[i]);
                if ( fabs(y[DiyID[i]]-yi[DiyID[i]]) > fabs(Diy[i]) ){
                    y[DiyID[i]] = yi[DiyID[i]] + Diy[i];
                }
                else disp_reached = false;
            }
        }
        
        if (nstep%out_every==0) fprintf(screen,"%d %e %e %e %d %e %e\n",nstep,Utot,Uold,DUtot,freeze_step,dt,cumDU);
        part_stress();
        if (nstep%out_every==0 && flag_Adump) write_dump("ALL.dump",astep);
        if (nstep%out_every==0) write_thermo("thermoALL.txt",astep);
    }
}

void part_stress()
{
    
    //stress per particle
    for (int i=0; i<Npart; i++){
        sx[i]=0.;
        sy[i]=0.;
        sxy[i]=0.;
        syx[i]=0.;
        shyd[i]=0.;
        sdev[i]=0.;
    }
    
    double Fx1, Fx2,Fy1,Fy2;
    double  pvol =((M_PI*D*D)/4.)/0.9069;
    
    for (int i=0; i<Npart; i++){
        for (int s=0; s<nn[i]; s++){
            int j = Nlist[i][s];
            Fx1 = F[i][s] * dx[i][s]/L[i][s];
            Fy1 = F[i][s] * dy[i][s]/L[i][s];
            Fx2 = -Fx1;
            Fy2 = -Fy1;
            sx[i] -= 0.5 * (Fx1*x[i]+Fx2*x[j]);
            sx[j] -= 0.5 * (Fx1*x[i]+Fx2*x[j]);
            sy[i] -= 0.5 * (Fy1*y[i]+Fy2*y[j]);
            sy[j] -= 0.5 * (Fy1*y[i]+Fy2*y[j]);
            sxy[i] -= 0.5* (Fy1*x[i]+Fy2*x[j]);
            sxy[j] -= 0.5* (Fy1*x[i]+Fy2*x[j]);
            syx[i] -= 0.5* (Fx1*y[i]+Fx2*y[j]);
            syx[j] -= 0.5* (Fx1*y[i]+Fx2*y[j]);
        }
    }
    
    for (int i=0; i<Npart; i++){
        shyd[i] = (sx[i]+sy[i])/3.;
        sdev[i] = sqrt(sx[i]*sx[i] + sy[i]*sy[i] + 3.*sxy[i]*sxy[i] - sx[i]*sy[i] );
        
        sx[i] /= pvol;
        sy[i] /= pvol;
        sxy[i] /= pvol;
        syx[i] /= pvol;
        shyd[i] /= pvol;
        sdev[i] /= pvol;
        
    }
}

void read_dump(std::string dname,int entry)
{
    std::string read_string, word;
    std::ifstream inFile(dname.c_str());
    bool found_entry = false;
    
    // reading input file.
    while (!inFile.eof() && !found_entry){
        
        // get next line....
        std::getline (inFile, read_string);
        
        //fprintf(screen,"%s\n",read_string.c_str());
        
        if (!read_string.empty()){
            
            // read next word in line....
            std::istringstream lss(read_string);
            bool getout = false;
            while (lss >> word && getout == false) {
                getout = false;
                if (strcmp(word.c_str(),"ITEM:")==0) {
                    lss >> word;
                    //fprintf(screen,"%s\n",word.c_str());
                    if (strcmp(word.c_str(),"TIMESTEP")==0) {
                        std::string newline;
                        std::getline(inFile, newline);
                        std::istringstream lst(newline);
                        int rtime;
                        lst >> rtime;
                        if (rtime == entry){
                            found_entry = true;
                            std::getline (inFile, newline);
                            std::getline (inFile, newline);
                            std::istringstream lst(newline);
                            int natoms;
                            lst >> natoms;
                            if (natoms != Npart){
                                fprintf(screen,"ERROR: inconsistent number of particles at requested time entry in dump file");
                                exit(0);
                            }
                            std::getline (inFile, newline);
                            std::getline (inFile, newline);
                            std::getline (inFile, newline);
                            std::getline (inFile, newline);
                            std::getline (inFile, newline);
                            std::istringstream lst1(newline);
                            lst1 >> word >> word;  //getting rid of ITEM: and ATOMS in title row in dump file
                            bool found_chi = false;
                            int chi_col = -1;
                            while (lst1 >> word && found_chi == false) {
                                chi_col++;
                                if (strcmp(word.c_str(),"chi")==0) {
                                    found_chi = true;
                                    fprintf(screen,"Found chi on dump column %d \n",chi_col);
                                }
                            }
                            if (!found_chi){
                                fprintf(screen,"ERROR: column named \"chi\" not found in dump file");
                                exit(0);
                            }
                            
                            for (int i=0; i<Npart; i++){
                                double chi_value, totrash;
                                std::getline (inFile, newline);
                                std::istringstream lst3(newline);
                                for (int j=0; j<chi_col; j++)  lst3 >> totrash;
                                lst3 >> chi_value;
                                chi.push_back(chi_value);
                            }
            
                        }
                        else getout = true;
                    }
                    else getout = true;
                }
                else getout = true;
            }
        }
    }
    
    if (!found_entry){
        fprintf(screen,"ERROR: specified time entry not found in dump file");
        exit(0);
    }
}