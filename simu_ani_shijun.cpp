#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TArrow.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TGraphErrors.h>
#include <TClonesArray.h>

#include <chealpix.h>

#include <DmpHKDSatStatus.h>

#include <time.h>
#include <math.h>
#include <cstdio>
#include <string>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <cstring>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

#define D2R 1.7453292519943295769e-2
#define R2D 5.7295779513082320877e+1

using namespace std;

/////////////////////对角度取360度的余数//////////////////////////////////////////////
double mod_360(double x) {
    long k;
    double y;

    k = (long)(x / 360.0);
    y = x - k * 360.0;
    while (y < 0.0)
        y += 360.0;
    while (y >= 360.0)
        y -= 360.0;

    return(y);
}

/////////////////////赤道坐标转黄道坐标（输入输出均以度为单位）///////////////////////////
int Equ2Ecl(double DampeTC,double l, double b, double& ra,double& dec)
{
    double x, y, z, r, e;

    double JD20000101=2451543.5; //JD of 2000/01/01 00:00:00
    double JD20130101=2456293.5; //JD of 2013/01/01 00:00:00 (When DAMPE Time Code Starts)
    double d=DampeTC/86400.+JD20130101-JD20000101; //DAMPE Time Code to Julian Day Since 2000/01/01 00:00:00
    e = 23.4393 - 3.563 * 0.0000001 * d;//黄赤交角

    l *= D2R;
    b *= D2R;
    e *= D2R;

    x = cos(b) * cos(l);
    y = cos(e) * cos(b) * sin(l) + sin(e) * sin(b);
    z = 0.0 - sin(e) * cos(b) * sin(l) + cos(e) * sin(b);

    if ((r = sqrt(x * x + y * y)))
        ra = mod_360(atan2(y, x) * R2D);
    else
        ra = 0;
    if (z)
        dec = mod_360(atan2(z, r) * R2D);
    else
        dec = 0;

    if (dec > 180) dec -= 360.;

    return 1;
}

//////////////////////////////////////////////////////////////////////////////////////
//Subroutine to Calculate the Sun Position at a Given Time
// Input: Dampe Time Code (DampeTC) in Unit of Second
//Output: Sun Position in Ecliptical Coordinate (lon,lat) in Degree
//////////////////////////////////////////////////////////////////////////////////////
int SunPosition_Hu(double DampeTC,double& lon, double& lat)
{
    int i;
    double d,w,a,e,M,oblecl,L,E,xe,ye,r,v,x,y,z,xequt,yequt,zequt,dist;

    double JD20000101=2451543.5; //JD of 2000/01/01 00:00:00
    double JD20130101=2456293.5; //JD of 2013/01/01 00:00:00 (When DAMPE Time Starts)
    d=DampeTC/86400.+JD20130101-JD20000101; //DAMPE Time Code to Julian Day Since 2000/01/01 00:00:00

    ////////////////////////////////////////////轨道根数//////////////////////////////////////////////
    w = 282.9404 + 4.70935 * 0.00001 * d;    //升交点经度
    a = 1;
    e = 0.016709 - 1.151 * 0.000000001 * d;  //偏心率
    M = 356.0470 + 0.9856002585 * d;         //平近点角
    oblecl = 23.4393 - 3.563 * 0.0000001 * d;//黄赤交角

    ////////////////////////////////////////////////////////////////////////////////////////////////
    L=w+M;//太阳的平均经度
    L=fmod(L,360);
    E=M + R2D * e * sin(M*D2R) * (1 + e * cos(M*D2R));//（开普勒方程的近似解）
    E=fmod(E,360);
    xe=cos(E*D2R) - e;//椭圆轨道上的直角坐标x
    ye=sin(E*D2R) * sqrt(1 - e*e);//椭圆轨道上的直角坐标y
    r=sqrt(xe*xe+ye*ye);//距离
    v=atan2(ye,xe)*R2D; //真近点角
    lon=fmod((v+w),360);//太阳黄经
    lat=0;

    return 1;
}

//Generate a Random Direction Uniformly Distributed On the Sphere Based on the Marsaglia Method
TVector3 RandDir()
{
    TVector3 RanDir;

    double u,v,r2;
    do{
        u=2*drand48()-1.0;
        v=2*drand48()-1.0;
        r2=u*u+v*v;
    } while (r2>=1.0 || r2==0.0);

    RanDir.SetXYZ(2*u*sqrt(1-r2),2*v*sqrt(1-r2),1.0-2*r2);

    return  RanDir;
}

//Module Transforming Equatorial Coord into Orbital Coord based on Coord Rotation
TVector3 Equ2Orb_RM(const TVector3& r_sat, const TVector3& v_sat, const TVector3& r_Particle)
{
    //The Track Direction is the Inverse Direction of the Incident Particle
    TVector3 r_Track = -r_Particle;

    //Calculate the Matrix for Coordinate Transformation Using Satellite Position and Velocity
    TMatrixD T(3,3);

    TVector3 Z_sat = -r_sat.Unit(); //Z_sat is in the Inverse Direction of Satellite Position
    TVector3 Y_sat =  v_sat.Cross(r_sat).Unit(); //Y_sat is in the Inverse Direction of Orbit Normal
    TVector3 X_sat =  Y_sat.Cross(Z_sat).Unit(); //X = Y Cross Z

    //The Row Vectors of the Rotation Matrix are X_sat, Y_sat, and Z_sat
    T(0,0) = X_sat.X();
    T(0,1) = X_sat.Y();
    T(0,2) = X_sat.Z();

    T(1,0) = Y_sat.X();
    T(1,1) = Y_sat.Y();
    T(1,2) = Y_sat.Z();

    T(2,0) = Z_sat.X();
    T(2,1) = Z_sat.Y();
    T(2,2) = Z_sat.Z();

    //Transform from Equatorial Coord to Orbit Coord by Applying the Rotation Matrix
    return T * r_Track;
}

//Module Transforming Orbital Coord into Equatorial Coord based on Coord Rotation
TVector3 Orb2Equ_RM(const TVector3& r_sat, const TVector3& v_sat, const TVector3& r_Track)
{
    //The Incident Particle Direction is the Inverse Track Direction
    TVector3 r_Particle = -r_Track;

    //Calculate the Matrix for Coordinate Transformation Using Satellite Position and Velocity
    TMatrixD T(3,3);

    TVector3 Z_sat = -r_sat.Unit(); //Z_sat is in the Inverse Direction of Satellite Position
    TVector3 Y_sat =  v_sat.Cross(r_sat).Unit(); //Y_sat is in the Inverse Direction of Orbit Normal
    TVector3 X_sat =  Y_sat.Cross(Z_sat).Unit(); //X = Y Cross Z

    //The Line Vectors of the Rotation Matrix are X_sat, Y_sat, and Z_sat
    T(0,0) = X_sat.X();
    T(1,0) = X_sat.Y();
    T(2,0) = X_sat.Z();

    T(0,1) = Y_sat.X();
    T(1,1) = Y_sat.Y();
    T(2,1) = Y_sat.Z();

    T(0,2) = Z_sat.X();
    T(1,2) = Z_sat.Y();
    T(2,2) = Z_sat.Z();

    //Transform from Orbit Coord to Equatorial Coord by Applying the Rotation Matrix
    return T * r_Particle;
}

//Normalized Phi (in Unit of Degree) Depedence of DAMPE Effective Area
double PhiDep (double Phi)
{
    return (0.94233-0.05444*cos(4*Phi*D2R)-0.00323*sin(Phi*D2R));
}

//Normalized Theta (in Unit of Degree) Depedence of DAMPE Effective Area
double ThetaDep (double Theta)
{
    if (Theta<50) {
        return (1.0-1.58766e-2*Theta+1.22454e-3*pow(Theta,2)-5.66470e-5*pow(Theta,3)+6.11793e-7*pow(Theta,4));
    }

    return (0);
}

void GetStrength(long nside,float md[],float mn[])
{
    long npix = nside2npix(nside);
    double nb[24],sb[24],eb[24],PixelTheta,PixelPhi;;
    TVector3 Dir1;

    for (int i=0; i<24; i++) {
        nb[i]=0;
        sb[i]=0;
        eb[i]=0;
    }

    for (int i=0; i<npix; i++) {
        pix2ang_ring(nside,i,&PixelTheta,&PixelPhi);
        PixelPhi*=R2D;
        if (PixelPhi<0) PixelPhi+=360;
        int j=PixelPhi/15;

        nb[j]++;
        sb[j]+=md[i];
        eb[j]+=md[i]*md[i];
    }

    for (int i=0; i<24; i++) {
        sb[i]/=nb[i];
        eb[i]=sqrt(eb[i]/nb[i]-sb[i]*sb[i]);
        cout << i*15+7.5 << "  " << nb[i] << "  " << sb[i] << "  " << eb[i] << endl;
    }
}

//Module Randomly Shuffling a Sequence of n Numbers
void RandomShuffle(long s[],long n)
{
    for(long i=0; i<n-1; i++) {
        long j=lrand48()%(n-i)+i;
        long a=s[i];
        s[i]=s[j];
        s[j]=a;
    }
}

//Main Function//////////////////////////////////////////////////////////////////////
int main()
{
    long nside,npix,ip,npixdr;
    float a[6],mr[49152],mm[49152],md[49152],mt[49152],dr[49152];
    int   mc15[49152],mc30[49152],mc45[49152],mc60[49152],mc90[49152];
    float ms15[49152],ms30[49152],ms45[49152],ms60[49152],ms90[49152];
    double EvtTime,Rx,Ry,Rz,Vx,Vy,Vz,l,b,lsun,bsun,PixelTheta,PixelPhi;
    TVector3 Dir1,Dir2,Dir3,Dir4,ParDir;

    long* rs = new long[25000000];
    double* Time = new double[25000000];
    TVector3* SatPos = new TVector3[25000000];
    TVector3* SatVel = new TVector3[25000000];
    TVector3* StkDir = new TVector3[25000000];

    //Read the DAMPE Response Map
    nside=64;
    npixdr = nside2npix(nside);
    ifstream InFile;
    InFile.open("CRPAngDisN64.dat",ios::in|ios::binary);
    InFile.read((char *)dr,(sizeof(dr[0]))*npixdr);
    InFile.close();

    //Randomize the Random Number Generator Seed Using Machine Time
    srand48((int)time(0));

    //Initialize the Maps
    nside=32;
    npix = nside2npix(nside);
    for (long i=0; i<npix; i++) {
        mr[i]=0;
        mm[i]=0;
        md[i]=0;
    }
/*
    //Read Root File
    TFile *f = new TFile("DAMPE_Proton_201601_201806_20GeV.root");
    if(!f || !f->IsOpen()) {
        cerr << "Can't open the ROOT file!" << endl;
        return 1;
    }
    TTree* ct=(TTree*)f->Get("tree");
*/
    //Read Multi Root Files from a File List Using TChain
    TChain *ct = new TChain("tree");

    ifstream ifs;
    ifs.open("DAMPE_Proton_201601_201812_20GeV", ifstream::in);
    if (ifs.is_open()) {
        std::string fname;
        ifs >> fname;
        while (!ifs.eof()) {
            if (fname[0] != '#') {
                ct->Add(fname.c_str());
            }
            ifs >> fname;
        }
    }
    else {
        std::cerr << "Error: Can NOT Open Input File." << std::endl;
        throw;
    }
    ifs.close();

    //Event Time
    ct->SetBranchAddress("EvtTime",&EvtTime);

    //House Keeping Data Rx
    ct->SetBranchAddress("px",&Rx);
    //House Keeping Data Ry
    ct->SetBranchAddress("py",&Ry);
    //House Keeping Data Rz
    ct->SetBranchAddress("pz",&Rz);

    //House Keeping Data Vx
    ct->SetBranchAddress("vx",&Vx);
    //House Keeping Data Vy
    ct->SetBranchAddress("vy",&Vy);
    //House Keeping Data Vz
    ct->SetBranchAddress("vz",&Vz);

    long TotEvt=ct->GetEntries();
    cout << "Total Proton Data: " << TotEvt << endl;

    long Events=0;
    for (long i=0; i<TotEvt; i++) {
        ct->GetEntry(i);

        //if (EvtTime<94608000 || EvtTime>157766400) continue; //Use Data During the Period of 20160101-20180101

        //Make an Event for the Current Satellite Position and Velocity
        Time[i]=EvtTime;
        SatPos[i].SetXYZ(Rx,Ry,Rz);
        SatVel[i].SetXYZ(Vx,Vy,Vz);

        bool NoEvent=true;
        while (NoEvent) {
            //Dir1.SetXYZ(0,1,0);  //Point Source at (0,1,0)
            Dir1=RandDir();  //Uniform Distribution on Sphere
            Dir2.SetXYZ(0,1,0);  //Direction of the Dipple
            if (drand48()*1.05<1.0+0.05*cos(Dir1.Angle(Dir2))) {//Add a Dipole of 0.05 Strength
                Dir2=Equ2Orb_RM(SatPos[i],SatVel[i],Dir1); //Transfore to Orbit Coord
                ang2pix_ring(64,Dir2.Theta(),Dir2.Phi(),&ip);
                if (drand48()<dr[ip]) {//DAMPE Response

                    mt[ip]++;
                    NoEvent=false;
                    StkDir[i]=Dir2;

                    //In Equatorial Coordinate
                    ang2pix_ring(nside,Dir1.Theta(),Dir1.Phi(),&ip);

                    //In Sun Fixed Ecliptic Coordinate
                    //SunPosition_Hu(Time[i],lsun,bsun);
                    //Equ2Ecl(Time[i],Dir1.Phi()*R2D,90.0-(Dir1.Theta())*R2D,l,b);
                    //ang2pix_ring(nside,(90.0-b)*D2R,(l-lsun)*D2R,&ip);

                    mr[ip]++;
                    Events++;
                }
            }
        }
    }
    cout << Events << endl;
/*
    int mtmax=0;
    for (int i=0; i<npixdr; i++) {
        if (mt[i]>mtmax) {
            mtmax=mt[i];
        }
    }
    for (int i=0; i<npixdr; i++) mt[i]/=mtmax;

    for (long i=0; i<Events; i++) rs[i]=i;

    int nround=20;
    for (int k=0; k<nround; k++) {

        //Random Shuffle the Entry
        RandomShuffle(rs,Events);

        for (long i=0; i<Events; i++) {

            //Generate Track Direction following DAMPE Response
            do {
                ParDir=RandDir();
                ang2pix_ring(64,ParDir.Theta(),ParDir.Phi(),&ip);
            } while (drand48()>dr[ip]); //Use mt[ip] Instead of dr[ip]

            //In Equatorial Coordinate
            ParDir=Orb2Equ_RM(SatPos[i],SatVel[i],ParDir); //Event Rate Method
            //ParDir=Orb2Equ_RM(SatPos[i],SatVel[i],StkDir[rs[i]]); //Random Shuffle Method
            ang2pix_ring(nside,ParDir.Theta(),ParDir.Phi(),&ip);

            //In Sun Fixed Ecliptic Coordinate
            //ParDir=Orb2Equ_RM(SatPos[i],SatVel[i],ParDir); //Event Rate Method
            //ParDir=Orb2Equ_RM(SatPos[i],SatVel[i],StkDir[rs[i]]); //Random Shuffle Method
            //Equ2Ecl(Time[i],ParDir.Phi()*R2D,90.0-(ParDir.Theta())*R2D,l,b);
            //SunPosition_Hu(Time[i],lsun,bsun);
            //ang2pix_ring(nside,(90.0-b)*D2R,(l-lsun)*D2R,&ip);

            mm[ip]++;
        }
    }
*/
    delete [] rs;
    delete [] Time;
    delete [] SatPos;
    delete [] SatVel;
    delete [] StkDir;
/*
    //Diffrentiate the Two Maps
    for (int i=0; i<npix; i++) {
        mm[i]/=nround;
        if (mr[i]>1.0 && mm[i]>1.0) {
            md[i]=mr[i]/mm[i]-1.0;
        }
    }
*/
    GetStrength(nside,mr,mr);
    //GetStrength(nside,mm,mr);
    //GetStrength(nside,md,mr);

    //Write the Real Map
    write_healpix_map(mr,nside,"MCERFMapOriDp005P3YN32.fits",0,"G");
/*
    //Write the Iso Map
    write_healpix_map(mm,nside,"MCERFMapIsoDp005P3YN32.fits",0,"G");

    //Write the Dif Map
    write_healpix_map(md,nside,"MCERFMapDifDp005P3YN32.fits",0,"G");


    //Initialize the Maps
    for (int i=0; i<npix; i++) {
        mc90[i]=0;
        mc60[i]=0;
        mc45[i]=0;
        mc30[i]=0;
        mc15[i]=0;
        ms90[i]=0;
        ms60[i]=0;
        ms45[i]=0;
        ms30[i]=0;
        ms15[i]=0;
    }

    Dir3.SetMagThetaPhi(1.0,97.4*D2R,  0.0*D2R);
    Dir4.SetMagThetaPhi(1.0,82.6*D2R,180.0*D2R);

    //Loop All Pixels
    for (int i=0; i<npix; i++) {
        pix2ang_ring(nside,i,&PixelTheta,&PixelPhi);
        Dir1.SetMagThetaPhi(1.0,PixelTheta,PixelPhi);

        //Skip the Direction Around Dir3 and Dir4
        //if (Dir1.Angle(Dir3)*R2D<45 || Dir1.Angle(Dir4)*R2D<45) continue;

        for (int j=0; j<npix; j++) {
            pix2ang_ring(nside,j,&PixelTheta,&PixelPhi);
            Dir2.SetMagThetaPhi(1.0,PixelTheta,PixelPhi);

            //Skip the Direction Around Dir3 and Dir4
            //if (Dir2.Angle(Dir3)*R2D<45 || Dir2.Angle(Dir4)*R2D<45) continue;

            if (Dir1.Angle(Dir2)*R2D<90) {
                ms90[i]+=md[j];
                mc90[i]++;
                if (Dir1.Angle(Dir2)*R2D<60) {
                    ms60[i]+=md[j];
                    mc60[i]++;
                    if (Dir1.Angle(Dir2)*R2D<45) {
                        ms45[i]+=md[j];
                        mc45[i]++;
                        if (Dir1.Angle(Dir2)*R2D<30) {
                            ms30[i]+=md[j];
                            mc30[i]++;
                            if (Dir1.Angle(Dir2)*R2D<15) {
                                ms15[i]+=md[j];
                                mc15[i]++;
                            }
                        }
                    }
                }
            }
        }
        ms90[i]/=mc90[i];
        ms60[i]/=mc60[i];
        ms45[i]/=mc45[i];
        ms30[i]/=mc30[i];
        ms15[i]/=mc15[i];
    }

    //Write the Smoothed Map
    write_healpix_map(ms90,nside,"MCERFMapDifDp005P3YN32S90.fits",0,"G");
    write_healpix_map(ms60,nside,"MCERFMapDifDp005P3YN32S60.fits",0,"G");
    write_healpix_map(ms45,nside,"MCERFMapDifDp005P3YN32S45.fits",0,"G");
    write_healpix_map(ms30,nside,"MCERFMapDifDp005P3YN32S30.fits",0,"G");
    write_healpix_map(ms15,nside,"MCERFMapDifDp005P3YN32S15.fits",0,"G");
*/
    return(0);
}
