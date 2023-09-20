#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iomanip>
#include <time.h>
#include <random>

//Equation for gate variables:
inline double x_inf(double v,double V12,double k) { return 1./(1.+exp(-(v-V12)/k)); }
inline double tau_inf(double tau,double v,double tauV12,double tauk) { return tau/(cosh((v-tauV12)/tauk)); }
inline double Gr(double v) { return 4.34e-5*exp(-0.0211539274*v); }//paper has error 4.34e5 should be e-5
inline double Gd1(double v) { return 0.075+0.043*tanh((v+20)-20); }
inline double Gv(double v) { return (10.6408- 14.6408*exp(-v/42.7671))/v; }


//NaF
const double mV12 = -43.8, mk = 6;
const double mtauV12 = -43.8, mtauk = 14;
const double hV12 = -67.5, hk = -11.8;
const double htau = 8.46, htauk = 12.8;
double gNaF = 150, htauV12 = -67.5, mtau = 0.25;
double dV12=0;

//NaP
const double mpV12 = -47.1, mpk = 3.1; 
const double mptau = 1, mptauV12 = -47.1, mptauk = 6.2;
double hpV12 = -60, hpk = -9;
double hptau = 5000, hptauV12 = -60;
double GNAPBLK = 1, hptauk = 9,FIXEDHNAP=0, hnapfixedvalue=0;
double T_HNAP_FIX=0,MNAP_FIX_onoff=0;
////for gspk[0:12]
double bursterID[100]={0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,0,0,1,0,1,0,1,1,1,1,0,0,0,0,1,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,1,1,0,1,1,0,1,1,0,0,1,0,1,1,0,0,1,1,0,1,1,0,0,0,1,0,0,1,1,0,1};
//for gspk=6
//double bursterID[100]=  {0,0,0,1,0,0,0,0,1,1,0,1,0,1,0,1,1,1,0,0,0,0,1,0,1,0,1,1,0,0,1,0,1,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,1,0,0,1,0,1,0,0,0,1,1,1,0,0,1,1,1,0,0,0,1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1};
double blockburster=0;

//CaL
const double mcV12 = -27.5, mck = 5.7;
const double mctau = 0.5;
const double hcV12 = -52.4, hck = -5.2;
const double hctau = 18;

//////////////////////////
//Na for spike height
//double gSPK =0;
const double mspkV12 = -27.5, mspkk = 1;
const double mspktau = 0.5;
const double hspkV12 = -27.5, hspkk = -1;//hspkV12 = -52.4, hspkk = -5.2;
const double hspktau = 5;//hspktau = 18;
const double h2spkV12 = -27.5, h2spkk = -1;
const double h2spktau = 1000, h2spktauV12 = -60, h2spktauk = 8.5119;
double hspkht_const=0;//spk height fast inactivation is constant at 1 if hspkht_const=1
double h2spkht_const=1;//spk height slow inactivation is constant at 1 if h2spkht_const=1
double GSPKRAMP=0,scale_gspk=1;
//////////////////////////

//////////////////////////
//IAHP for AHP Magnitude
//double gAHP =0;
const double mahpV12 = -27.5, mahpk = 1;//mahpV12 = -27.5, mahpk = 1;
const double mahptau = 5;
const double hahpV12 = -27.5, hahpk = -1;
const double hahptau = 1000, hahptauV12 = -60, hahptauk = 8.43;
double hahp_const=1;//ahp inactivation is constant at 1 if hahp_const=1
//////////////////////////


double Caout=4, alphaCa=2.5e-5,tauCa=500, Cain0=0.0000000001, Cav=0; 
double gCaL=0;
double Kpump=1e-3,Vpump=Kpump/tauCa;

//////////////////////////////////////////////////
//IP3 variables
//////////////////////////////////////////////////
//turn ER ca dynamics on (1) or off (0)
double ER_on_off = 1;

//m_ER (activation) parameters for J_ER_in
double IP3=0.0015,ki=0.001,kaa=0.0001;

// h_ER (inactivation) parameters for of J_ER_in
double lp=0.1,kd=0.0002;

//Current to concentration conversion factors (fca and Aca are the same and already declared above with alphaCa)
double fca=alphaCa,Aca=alphaCa,sigma_ca=0.185;

//Pump strength (Already declared above at 50ms)
//double tauca=500;

//SERCA pump variables: strength (gSERCA mM/ms) and activation (kSERCA mM)
double gSERCA=0.45,kSERCA=0.00005;

//J_ER_in leak conductance (gL_ER) and Maximal conductance (gER)
double gL_ER=0.1,gER=31000*2.5;

///////////////////////////////////////////////////
///////////////////////////////////////////////////

//////////////////////////////////////////////////
//Kout diffustion and glia variables
//////////////////////////////////////////////////
//#Diffusion
double taudiff=750;

//#Glia #Gmax is in units of mM/mS
double zg=10,kf=6.25,Gmax=0.5;

//#intra to extracellular ratio?
double alpha_kout=0.000105;
///////////////////////////////////////////////////
///////////////////////////////////////////////////


//Kdr
const double   nAk = 5;
const double nB = 0.17, nBk = 40;
double nA = 0.011, nAV = 44, nBV = 49;
double gKdr = 220; 
//CAN
double K_CAN = 0.74*1e-3, nc = 0.97, sigma_CAN=-0.05e-3;
double ECAN=0;
const double hcanC12 = .00015, hcank = -.0002;
const double hcantau = 100;

//ChR2
double Echr2 = 0;
double gma_chr2=0.1, eps1=0.8535, eps2=0.14;
double sig_ret=10e-20, Gd2=0.05, tau_chr2=1.3;



//Arch
double Earch = -145;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Synaptic Depression parameters: Each neuron has it's own D_syn value which represents the current level of depression for that neurons synapes
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double f_D = 0.2, D_0 = 1, tau_D = 1000; //f_D is the relative amount that the synapse is depressed with each spike, D_0 is the starting value for depression, and tau_D is the recovery rate from depression
					      //parameters are estimated from Kottick and Del Negro 2015
double pD = 0.0, pDr = 1-pD;//(not used right now)pD is the proportion not able to depress
double D_on_off = 1; //Turns depression on (1) or off (0)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double c=36;
double I_scale=1, syn_scale=1, INaP_scale=1,mINaP=1;
double NPM = 100;//number of pacemakers or Sst-
double TEMP=27,DTDt=0;
//double ICm;
//concentrations
double Nain=15, Naout = 120;
double imposed_Nain=15;
double Kin = 125, kbath=4;                    
double breakdown=1;
double alphaNa=5e-5;
double Nain0=15, tauNain=1000;
double dynamicNain=0; //set to 1 for dynamic Nain

//Pump
double Rpump=10;

//double ENa=26.54*log(Naout/Nain);

const double Esyn=-10;

//flag for distribiting gleak and gnap in a grid
double IBFIND=0;

double rnd() {return double(rand())/RAND_MAX; }



double tausyn=5;

double IAPP = 0;
double stim_ON=0, stim_W=0, cells[2]={75,85};

//Ca2+ squarewave 
double square_ON_OFF = 0;//Turning the ca2+ square wave transient on/off: 0=off 1=on
double Iapp_ca=0;

//tonic inhib
double ginh = 0;
double Einh=-70;

//OPIOIDS
double gmu=0;//conductance of mu-opioid activated potassium channel
double musyndep=1;//depression of synaptic strength by opioids

//Values for leak dist. exp((kbath-A_gl)/B_gl);
double A_gl = (8.5-((8.5-3.425)/4.05)*5);//2.2345679;//changed (7/13/2023) 3.425;
double B_gl = 5;//changed (7/13/2023) 4.05;


//////////
//hypoxia
//////////
double hypox=1;
double hypox_spk=0;
double hypox_ahp=0;

//////////
//for 2D plots in Fig 1 gTonic vs gSPK or gAHP
//////////
double uniformdist=0;


//////////////////////////////////////
//synaptic weights
//////////////////////////////////////
//nonPM to nonPM(2->2)
double Wnpm2npm = 0.0063;//Weight of nonPM to nonPM
double Pnpm2npm = 0.02;//Prob of nonPM to nonPM
//PM to nonPM (1->2)
double Wpm2npm = 0.000175;//Weight of PM to nonPM
double Ppm2npm = 0.3;//0.25;//Prob of PM to nonPM
//nonPM to PM(2->1)
double Wnpm2pm = 0.05*5;//Weight of PM to nonPM
//double Pnpm2pm = prob;//This is defined below because the parameter prob is passed through the code.
//////////////////////////////////////
//synaptic weights
//////////////////////////////////////

//////////////////////////////////////
//gL, gNaP dist. parameters
//////////////////////////////////////
double rho=0.8;//covariance
double sL_pct=0.05;//gleak dist pct of mean
//////////////////////////////////////
//gL, gNaP dist. parameters
//////////////////////////////////////


//Global var for imposed PM spike freq.
int imposedPMfr=0;
double PMfr = 0, NPMfr = 0, Dsyn_NPM=1;
double SPIKETHRESH=-35;

class Neuron
{
	public:

// parameters	

	double gCa,gNaP,gCAN,gChR2,gArch,gSPK,gAHP;
	double g_L;
	int ID;
	double stimIDs;
	double tlast_AP;
	double drivestart;

// dynamical variables

	double v;
	double m,h; // INaF
	double mp,hp; // INaP
	double mc,hc; // ICa
	double n; // IKdr
	double mspk,hspk,h2spk; // Ispk for spk heigh manipulations
	double mahp,hahp; //IAHP for AHP manipulations
	double Cain;
	double Cav;
	double gsyn; 
	double Nain_test;
	double hcan;
	double OP1,OP2,CL1, CL2, Pchr2; //IChR2
	double Carch,Oarch,Iarch;//arch channel states
	double m_ER,h_ER,J_ER_in, J_ER_out, caER, catot;
	double EK,Kout;
	double D_syn;
	double E_leak;
	double ENa;
	double ICm;//Current induced by changes in capacitance
	double Ipump;//Na/K current
	double FIXMNAP;
	

	Neuron();
	void init(double,double,double,double,double,double);
	void step(double dt,double drive,double fcan, double fsynCa, double irr);
	int spike(double vth,double dt,double drive,double fcan, double fsynCa, double irr);
};

struct connection
{
	int source,target;
	double weight;
};

class Population
{
	public:

	Neuron* net;
	int size;

	connection *w;
	int sex;

	double vth;
	
	Population(int n,double gca0,double gca1,double gl0,double gl1,double gnap0,double gnap1,
		double gcan0,double gcan1,double gchr20,double gchr21,double garch0,double garch1,double prob,double wght, double stimnum,double gspk0,double gspk1,double gahp0,double gahp1,double seed)
	{
		vth=SPIKETHRESH;//-35;//-35;
		size=n;
		net=new Neuron [size];
		std::cout << "The value of seed is: " << seed << std::endl;
		std::default_random_engine generator;//random number generator for gnap/gleak dist. needs to be before loop
		if(seed!=1){ generator.seed(seed);}//added the line above so that I don't change from the default distribution unless the seed is something other than 1

		for(int i=0;i<size;i++)
		{	
			
                     	net[i].ID=i;
                     	net[i].stimIDs=0;
                     	net[i].FIXMNAP=0;
                     	net[i].catot=0.001;
                     	net[i].Cain=0.0000001;
			net[i].gCa=gca0+(gca1-gca0)*rnd();
			net[i].gCAN=gcan0+(gcan1-gcan0)*rnd();
			net[i].gNaP=gnap0+(gnap1-gnap0)*rnd();
			net[i].g_L=gl0+(gl1-gl0)*rnd();
			net[i].gChR2=0.0;//rnd();
			net[i].gArch=0.0;//rnd();//Distrabution of I_Arch conductance
			net[i].gChR2=gchr20+(gchr21-gchr20)*rnd();//rnd();
			net[i].gSPK=gspk0+(gspk1-gspk0)*rnd();
			net[i].gAHP=gahp0+(gahp1-gahp0)*rnd();
			net[i].drivestart=0;
			net[i].Nain_test=15;
			

			
			//assign gnap,gleak and gcan values to rhythm generating population
			if (i<NPM && size>2){		
			double mu_nap=3.33;
			double s_nap=0.75;
			std::normal_distribution<double> dist_gnap_PM (mu_nap,s_nap);
			net[i].gNaP=dist_gnap_PM(generator);		
			double M_gl=exp((kbath-A_gl)/B_gl);
			double s_gl=sL_pct*M_gl;
			double mu_gl=M_gl+rho*(s_gl/s_nap)*(net[i].gNaP-mu_nap);
			double sigma_gl=sqrt((1-rho*rho)*s_gl*s_gl);
			
			std::normal_distribution<double> dist_gl_PM (mu_gl,sigma_gl);
			net[i].g_L=dist_gl_PM(generator);	
			net[i].gCAN=0.0;
			}
			
			//assign gnap,gleak and gcan values to pattern generating population
			if (i>=NPM && size>2){
			net[i].g_L=gl0+(gl1-gl0)*rnd();
			std::normal_distribution<double> dist_gnap_nPM (net[i].g_L*0.28,net[i].g_L*0.03);
			double mu_nap=1.5;
			double s_nap=0.25;
			std::normal_distribution<double> dist_gnap_PM (mu_nap,s_nap);
			net[i].gNaP=dist_gnap_PM(generator);
			
			double M_gl=exp((kbath-A_gl)/B_gl);
			double s_gl=0.025*M_gl;
			double mu_gl=M_gl+rho*(s_gl/s_nap)*(net[i].gNaP-mu_nap);
			double sigma_gl=sqrt((1-rho*rho)*s_gl*s_gl);
			std::normal_distribution<double> dist_gl_PM (mu_gl,sigma_gl);
			net[i].g_L=dist_gl_PM(generator);
			
			std::normal_distribution<double> dist_gcan(gcan0,gcan1);
//			net[i].gCAN=gcan0+(gcan1-gcan0)*rnd();//dist_gcan(generator);//changed 08/10/22 to try and get distributed recruitment of follower neurons
			net[i].gCAN=dist_gcan(generator);
			}
			
//			//2 NEURON NETWORK
			if (i==0 && size==2)
			{
			net[i].gCAN=0.0;
			net[i].gNaP=3.33;
			net[i].g_L==gl0;
			}
			if (i==1 && size==2)
			{
			net[i].gCAN=0;
			net[i].gNaP=3.33;
			net[i].g_L==gl1;
			}
//			
			//distribute gL and gNaP in grid for identification of burst capable neurons
			if(IBFIND==1)
			{
				int D1=20;
				//loop for assigning NaP
				for(int n=0;n<D1;n++)
				{
					//loop fro assigning Leak
					for(int m=0;m<D1;m++)
					{
						net[D1*n+m].gNaP=gnap0+n*(gnap1-gnap0)/(D1-1);
						net[D1*n+m].g_L=gl0+m*(gl1-gl0)/(D1-1);
//						net[i].drivestart= 0.16*net[D1*n+m].g_L +(-0.225) - net[D1*n+m].gNaP*(0.3/5.0);//Equation for assing a starting drive value as a fn of gnap and gleak
//						if(net[i].drivestart<0){net[i].drivestart=0;}
					}
				}
			
			}
			//Uniformly distribute gSPK or gAHP and ramp drive
			if (uniformdist==1 && size>1){
						net[i].gNaP=gnap0;
						net[i].g_L=gl0;
//						net[i].drivestart=0.4+i*((0.15-0.4)/(size-1));
//						if(gnap0==gnap1){net[i].drivestart=-0.05*gnap0+0.4;}
//						if(gspk0==gspk1 ||gahp0==gahp1){net[i].drivestart=0;}
						net[i].gSPK=gspk0+i*((gspk1-gspk0)/(size-1));
						net[i].gAHP=gahp0+i*((gahp1-gahp0)/(size-1));
						}
			//Uniformly distribute gSPK or gAHP and ramp drive
			if (size==1){
						net[i].gNaP=gnap0;
						net[i].g_L=gl0;
						net[i].gSPK=gspk0;
						net[i].gAHP=gahp0;
						}
		}
		
		w=new connection [size*size];
		sex=0;
		for(int i=0;i<size;i++) for(int j=0;j<size;j++)
		{
			/////////////////////
			/////////////////////
			//Neurons 0-24=inhib;25-99=PM; 100-end=NPM
			////////////////////
			////////////////////
			
			double conprob = prob;
			double tmpweight = wght;
			
			
			
			//nonPM to PM
			double Pnpm2pm = prob;//Prob of PM to nonPM
			
			if(size>2){
			//nonPM to nonPM
			if(i>=NPM && j>=NPM) {tmpweight = Wnpm2npm; conprob=Pnpm2npm;}
			//PM to nonPM
			if(i<NPM && j>=NPM) {tmpweight = Wpm2npm; conprob=Ppm2npm;}
			//nonPM to PM
			if(i>=NPM && j<NPM) {tmpweight = Wnpm2pm; conprob=Pnpm2pm; }
			//PM to PM
			if(i<NPM && j<NPM) {tmpweight = wght; conprob=prob;}
			}

			//std::default_random_engine generator;
			if(i==j || rnd()>conprob) continue;
			w[sex].source=i;
			w[sex].target=j;
			w[sex].weight=tmpweight*rnd();
			if(size==2){w[sex].weight=tmpweight;}
			sex++;
		}
	}

	~Population() { delete w; delete net; }

	int step(double dt,double drive,int* spk,double fcan,double fw, double fsynCa, double irr)
	{
		int sp=0;
		for(int i=0;i<size;i++) { spk[i]=net[i].spike(vth,dt,drive,fcan,fsynCa,irr); sp+=spk[i]; }
		for(int i=0;i<sex;i++){	 
					 if(spk[w[i].source]&& w[i].source<NPM) {net[w[i].target].gsyn+=w[i].weight*fw*musyndep*net[w[i].source].D_syn;}//if source is a PM neuron (includes opioid reduction of synapses)
					 else if(spk[w[i].source])              {net[w[i].target].gsyn+=w[i].weight*fw*net[w[i].source].D_syn;}//if source is a follower neuron
					 
					}
		//If synaptic depression is on
		if(D_on_off==1){
				for(int i=0;i<size;i++){ if(spk[i]) {net[i].D_syn += -f_D*net[i].D_syn; net[i].tlast_AP=0;}} //if neuron i spikes then update D_syn
				}
		return sp;
	}
};

//initial conditons
Neuron::Neuron()
{
		v=-65+1*rnd(); m=.1; h=.1; mp=.1; hp=.4*rnd(); n=.1; Cain=.000001; Cav=5*1e-5; mc=.1; hc=.1;mspk=.1;hspk=1;h2spk=1;mahp=.1;hahp=1;Nain_test=15;
		gCa=.003;
		gNaP=5;
		gCAN=1;
		gChR2=0.4;
		gArch=0.4;
		E_leak=-68; g_L=2.5;
		hcan=1;
		OP1 =0.00; OP2=0.00; CL1=0.99; CL2=0.01; Pchr2=0.1;
		Carch = 0.9; Oarch = 0.0;
		m_ER=0.1; h_ER=0.99;
		Kout = 4;
		D_syn=D_0;
		caER=0.002; //catot=0.0005;
}

int Neuron::spike(double vth,double dt,double drive,double fcan, double fsynCa, double irr)
{
		double vpre=v;
		step(dt,drive,fcan,fsynCa, irr);
		return (vpre<vth && v>=vth);
}

void Neuron::init(double gca,double glk,double gnap,double gcan,double gchr2,double garch)
{
	gCa=gca;
	g_L=glk;
	gNaP=gnap;
	gCAN=gcan;
	gChR2=gchr2;
	gArch=garch;
	m=0.1*rnd();
	h=0.1*rnd();
	mp=0.1*rnd();
	hp=0.1*rnd();
	mc=0.1*rnd();
	hc=rnd();
	n=rnd();
	Cain=0.000001;
	caER=0.002; 
	Cav=5e-5*rnd();
	OP1=0.00; OP2=0.00; CL1=0.99; CL2=0.01; Pchr2=0.1;
	Carch = 0.9; Oarch = 0.0;
	
}

int flag_hcan=0;

void Neuron::step(double dt,double drive,double fcan,double fsynCa, double irr)
{
	double Q10=pow(1.5,((27-TEMP)/10));
	double deltaT=(TEMP-27);
	double RTF=8.314*(308+deltaT)/96.485;//initial value 26.54
	//double EQ10=0.86*(21-TEMP)
//	if(ID<NPM && imposedPMfr==1)	{v=-65;
//					tlast_AP+=dt;
////					if((PMfr+1.0*Dsyn_NPM*NPMfr)*dt/1000>=rnd()){v=0.0;tlast_AP=0;}
//					if((PMfr)*dt/1000>=rnd()){v=0.0;tlast_AP=0;D_syn=1;}
////					if((PMfr)*dt/1000>=rnd()){v=0.0;tlast_AP=0;}
//					}
	if(ID>=NPM || imposedPMfr==0)
	{
	
	double ECa=26.54*log(Caout/Cain)/2;
	EK=breakdown*26.54*log(Kout/Kin);
	ENa=breakdown*26.54*log(Naout/Nain_test);
	double pna = 1, pk=42;
	E_leak = (-26.54)*log((pna*Nain_test + pk*Kin)/(pna*Naout + pk*Kout));
	double INaF = gNaF*m*m*m*h;
	double INaP = gNaP*mp*hp;
	if(blockburster==1 && bursterID[ID]==1)
		{
		INaP = gNaP*mp*hp*GNAPBLK;
		}
	if(blockburster==2 && bursterID[ID]==0)
		{
		INaP = gNaP*mp*hp*GNAPBLK;
		}
	if(blockburster==0)
		{
		INaP = gNaP*mp*hp*GNAPBLK;
		}
	
	double ICa = gCa*mc*hc+gCaL;
	double IChR2 = gChR2*Gv(v)*(OP1+gma_chr2*OP2)*(v-Echr2);
	double IArch = gArch*Oarch*(v-Earch);
	double k0 = v+nAV;
	double k1 = nA*k0/(1-exp(-k0/nAk));
	double k2 = nB*exp(-(v+nBV)/nBk);
	double taun_inf = 1/(k1+k2);
	double n_inf = k1/(k1+k2);
	double IKdr = gKdr*n*n*n*n;
	double ICAN = gCAN*fcan/(1.+pow(K_CAN/Cain,nc))*hcan; // Tporikova	
	double Istim = 0;
//	ICm=(v-130)*c*1.003*DTDt;//Pinto et al biophysical reviews 2022
	ICm=(v-130)*c*0.003*DTDt;//Equation and the 130mV is from Pinto et al biophysical reviews 2022, alpha (0.003%) is from plaksin Physical Review X 2018
	Ipump=Rpump*( 1/(1+pow(Nain0/Nain_test,3.0)) - 0.5);
	////////////////////////////////////////////////////
	//*New currents for spike height & AHP manipulations
	////////////////////////////////////////////////////
	double ISPK = scale_gspk*hypox*(gSPK+hypox_spk)*mspk*hspk*h2spk*(v-ENa);
	double IAHP = hypox*(gAHP+hypox_ahp)*mahp*hahp*(v-EK);
	
	
	if(stimIDs==1){Istim =stim_W;}

	double Iinh=ginh*(v-Einh);
	if(ID<=NPM){Iinh=0;}
	
	//opioid DEP POTASSIUM CURRENT Imu
	double Imu=0;
	if(ID<NPM){Imu=gmu*(v-EK);}
	

	
	v+=(-I_scale*INaF*(v-ENa)-INaP_scale*INaP*(v-ENa)-I_scale*IKdr*(v-EK)-I_scale*ICAN*(v-ECAN)-I_scale*ICa*(v-ECa)-I_scale*IChR2-I_scale*IArch-I_scale*g_L*(v-E_leak)-ICm - I_scale*ISPK- I_scale*IAHP-(drive+syn_scale*gsyn)*(v-Esyn)-Iinh - Imu + IAPP + Istim)/c*dt;

	//INa
	m+=(x_inf(v,(mV12+dV12),mk)-m)*(1.-exp(-dt/(Q10*tau_inf(mtau,v,mtauV12,mtauk))));
	h+=(x_inf(v,(hV12+dV12),hk)-h)*(1.-exp(-dt/(Q10*tau_inf(htau,v,htauV12,htauk))));
	//INaP
	if(MNAP_FIX_onoff==1 && tlast_AP>1000 ){FIXMNAP=1;}// && FIXMNAP==0 && tlast_AP>20
	if(FIXMNAP==0){mp+=(x_inf(v,mpV12,mpk)-mp)*(1.-exp(-dt/(Q10*tau_inf(mptau,v,mptauV12,mptauk))));}
	if(FIXEDHNAP==0)hp+=(x_inf(v,hpV12,hpk)-hp)*(1.-exp(-dt/(Q10*tau_inf(hptau,v,hptauV12,hptauk))));
	if(FIXEDHNAP==1 || FIXEDHNAP==2 || FIXEDHNAP==3){hp=hnapfixedvalue;}
//	if(FIXEDHNAP==3){hp=x_inf(v,hpV12,hpk);}
	//Ica 
	mc+=(x_inf(v,mcV12,mck)-mc)*(1.-exp(-dt/(Q10*mctau)));
	hc+=(x_inf(v,hcV12,hck)-hc)*(1.-exp(-dt/(Q10*hctau)));
	//IK
	n+=(n_inf-n)*(1.-exp(-dt/(Q10*taun_inf)));
	//IP3
	m_ER=(Cain*IP3)/((Cain+kaa)*(IP3+ki));
	h_ER=h_ER+(lp*(kd-(Cain+kd)*h_ER))*dt;
	J_ER_in =(gL_ER + gER*m_ER*m_ER*m_ER*h_ER*h_ER*h_ER)*(caER-Cain);
	J_ER_out = gSERCA*Cain*Cain/(kSERCA*kSERCA+Cain*Cain);
	
	if(ER_on_off==0){m_ER=0;J_ER_in=0;J_ER_out=0;}
	
	////////////////////////////////////////////////////
	//*(In)activation for new currents for spike height & AHP manipulations
	////////////////////////////////////////////////////
	//Ispk
	mspk+=(x_inf(v,mspkV12+dV12,mspkk)-mspk)*(1.-exp(-dt/mspktau));
	if(hspkht_const==0){
			hspk+=(x_inf(v,hspkV12+dV12,hspkk)-hspk)*(1.-exp(-dt/hspktau));
		}else{
			hspk==1;
		}
	if(h2spkht_const==0){
			h2spk+=(x_inf(v,h2spkV12+dV12,h2spkk)-h2spk)*(1.-exp(-dt/tau_inf(h2spktau,v,h2spktauV12,h2spktauk)));
		}else{
			h2spk==1;
		}
	
	//IAHP
	mahp+=(x_inf(v,mahpV12,mahpk)-mahp)*(1.-exp(-dt/mahptau));
	if(hahp_const==0){
//			h2spk+=(x_inf(v,h2spkV12,h2spkk)-h2spk)*(1.-exp(-dt/h2spktau));
//			hahp+=(x_inf(v,hahpV12,hahpk)-hahp)*(1.-exp(-dt/hahptau));
			hahp+=(x_inf(v,hahpV12,hahpk)-hahp)*(1.-exp(-dt/tau_inf(hahptau,v,hahptauV12,hahptauk)));
//			h2spk+=(x_inf(v,h2spkV12,h2spkk)-h2spk)*(1.-exp(-dt/tau_inf(h2spktau,v,h2spktauV12,h2spktauk)));
		}else{
			hahp==1;
		}
	

	//Ca2+ dynamics
	Cain+=(fca*(J_ER_in-J_ER_out) - alphaCa*ICa*(v-ECa) - alphaCa*(gsyn)*fsynCa*(v-ECa) + 0*alphaCa*fsynCa*Istim +alphaCa*Iapp_ca + (Cain0-Cain)/tauCa)*dt;
	catot+=(-alphaCa*ICa*(v-ECa) - alphaCa*(gsyn)*fsynCa*(v-ECa) + 0*alphaCa*fsynCa*Istim  +alphaCa*Iapp_ca + (Cain0-Cain)/tauCa)*dt;
	caER=(catot-Cain)/sigma_ca;
	
	//Nain dynamics
	double xna=(EK-E_leak)/(EK-ENa);
	if(dynamicNain==1){Nain_test+=(-alphaNa*INaF*(v-ENa) -alphaNa*xna*g_L*(v-E_leak)-alphaNa*INaP*(v-ENa) + (Nain0-Nain_test)/tauNain)*dt;}
	if(dynamicNain==2){Nain_test+=(-alphaNa*INaF*(v-ENa)-alphaNa*INaP*(v-ENa) - 3*alphaNa*Ipump)*dt;}
	if(dynamicNain==0){Nain_test=Nain;}
	if(dynamicNain==3 || dynamicNain==4){Nain_test=imposed_Nain;}
	
	//Make sure that concentrations cant go below 0;
	if(Cain<0){Cain=0;}
	if(catot<0){catot=0;}
	
	Kout=kbath;
	gsyn*=exp(-dt/(Q10*tausyn));
	
	}//end of if statment for imposed rhythm
	tlast_AP+=dt;
	if(D_on_off==1){D_syn+=((D_0-D_syn)/(Q10*tau_D))*dt;}//Exponential decay back to D_syn=1.
}//End of "Neurons::step"

using namespace std;

ostream& operator <<(ostream& os,Neuron& N)
{
	return (os<<N.v<<'\t'<<N.m<<'\t'<<N.h<<'\t'<<N.mp<<'\t'<<N.hp<<'\t'<<N.n<<'\t'<<N.Cain);
}

istream& operator >>(istream& is,Neuron& N)
{
    return (is>>N.v>>N.m>>N.h>>N.mp>>N.hp>>N.n>>N.Cain);
}

ostream& operator <<(ostream& os,Population& p)
{
	for(int i=0;i<p.size;i++) os<<p.net[i]<<'\t';
	return os;
}

istream& operator >>(istream& is,Population& p)
{
	for(int i=0;i<p.size;i++) is>>p.net[i];
	return is;
}

int main(int argc,char** argv)
{
	double dt=.025,T=10000;
	double DT=20;
	int size=100;
	double gca[2]={0.0000065,0.0000065},gleak[2]={2.75,2.75},gnap[2]={3,5},gcan[2]={2,1},prob=.1,wght=.1;
	double dr[2]={0,1},fc[2]={1,1},fcw[2]={1,1}, fsca[2]={0.01,0.01}, tauDrug = -5000, pctblk = 0, irrpower[2]={0.0,0.0}, gchr2[2] = {0.0,0.0}, garch[2] = {0.0,0.0}, gspk[2] = {0.0,0.0}, gahp[2] = {0.0,0.0};
	double pulseON=0, pulsefreq=0.0, pulsedur=100,pfmax =.25, pfmin=.25, dl=0.1, numdl=0, stimnum= 100, seed=1, max_gnap_blk=0.5, APPramp[2] ={0.0,0.0};
	double stim_time=2000,strength = 25, Nstim=0,stimcount=0, stimseed =0, GINHIB=0,stim_delay=0, stim_advance=0;
	double caTon=250, caP=5000, caToff=caP-caTon, caAmp=0,CRC=0,nain_ramp[2]={15,15},hypox_ramp[2]={1,1},hypox_ramp2[2]={1,1},hypox_spk_ramp[2]={0,0},fignum=0,TEMP_ramp[2]={27,27};
	
	//vars for imposed PM rhythm; AMP=spks/s, all other params are in ms.
	double AMP = 30, period =5000, t_up=200, t_down=400, dFR_up=AMP/t_up, dFR_down=-AMP/t_down, tlast=0,IBIamp[2]={0,0},refract=0.0,Cm_ramp[2]={36,36},TEMPsteps=-1,TBATH=27,DBATHTDt=0, DCmDt=0, dV12_ramp[2]={0.0,0.0}, dV12_ramp2[2]={0.0,0.0};
	double tauNapump_ramp[2]={1000,1000},inap_blk_ramp[2]={1,1};
	double T_DSYN_FIX=0,T_SYN_OFF=0, T_SYN_MOD=0, syn_mod_amp=1,T_MNAP_FIX=0.0,RAMPSTART=20000,DEVMODE=0,SYNCUT=0;
	int ini=0, block_type = -9, Nholostim = 0, RASTONLY=0;
	char fname[256]="dat";
	int output=0;
	for(int i=1;i<argc;i++)
	{
		//Flags for changing model parameters
		if(strcmp(argv[i],"-DT")==0) DT=atof(argv[++i]);
		else if(strcmp(argv[i],"-i")==0) ini=1;
		else if(strcmp(argv[i],"-hcan")==0) flag_hcan=1;
		else if(strcmp(argv[i],"-o")==0) output=1;
		else if(strcmp(argv[i],"-dt")==0) dt=atof(argv[++i]);
		else if(strcmp(argv[i],"-DT")==0) DT=atof(argv[++i]);
		else if(strcmp(argv[i],"-s")==0) size=atoi(argv[++i]);
		else if(strcmp(argv[i],"-T")==0) T=atof(argv[++i])*1000;
		else if(strcmp(argv[i],"-d")==0) { dr[0]=atof(argv[++i]); dr[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-w")==0) wght=atof(argv[++i]);
		else if(strcmp(argv[i],"-tgs")==0) tausyn=atof(argv[++i]);
		else if(strcmp(argv[i],"-nap")==0) { gnap[0]=atof(argv[++i]); gnap[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-na")==0) gNaF=atof(argv[++i]);
		else if(strcmp(argv[i],"-kdr")==0) gKdr=atof(argv[++i]);
		else if(strcmp(argv[i],"-ca")==0) { gca[0]=atof(argv[++i]); gca[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-gleak")==0) { gleak[0]=atof(argv[++i]); gleak[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-cal")==0) gCaL=atof(argv[++i]);
		else if(strcmp(argv[i],"-tca")==0) tauCa=atof(argv[++i]);
		else if(strcmp(argv[i],"-kp")==0) Kpump=atof(argv[++i]);
		else if(strcmp(argv[i],"-vp")==0) Vpump=atof(argv[++i]);
		else if(strcmp(argv[i],"-kbath")==0) kbath=atof(argv[++i]);
		else if(strcmp(argv[i],"-prob")==0) prob=atof(argv[++i]);
		else if(strcmp(argv[i],"-can")==0) { gcan[0]=atof(argv[++i]); gcan[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-fc")==0) { fc[0]=atof(argv[++i]); fc[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-fw")==0) { fcw[0]=atof(argv[++i]); fcw[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-fsca")==0) { fsca[0]=atof(argv[++i]); fsca[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-ca0")==0) Cain0=atof(argv[++i]);
		else if(strcmp(argv[i],"-exp")==0) block_type=atof(argv[++i]); //1=can 2=TRPC3 
		else if(strcmp(argv[i],"-tauDrg")==0) tauDrug=atof(argv[++i]);//tau drug
		else if(strcmp(argv[i],"-sd")==0) seed=atof(argv[++i]);
		else if(strcmp(argv[i],"-scalegNaK")==0) {gNaF=gNaF*atof(argv[++i]); gKdr=gKdr*atof(argv[++i]);}
		else if(strcmp(argv[i],"-na_h12")==0) htauV12=atof(argv[++i]);
		else if(strcmp(argv[i],"-mtau")==0) mtau=atof(argv[++i]);
		else if(strcmp(argv[i],"-ktau")==0) nA=atof(argv[++i]);
		else if(strcmp(argv[i],"-Cm")==0) c=atof(argv[++i]);
		else if(strcmp(argv[i],"-K_v12")==0) {nAV=atof(argv[++i]); nBV=atof(argv[++i]);}
		else if(strcmp(argv[i],"-hp12")==0) hpV12=atof(argv[++i]);
		else if(strcmp(argv[i],"-iapp")==0) {APPramp[0]=atof(argv[++i]); APPramp[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-hptauk")==0) hptauk=atof(argv[++i]);
		else if(strcmp(argv[i],"-ttxblk")==0) GNAPBLK=atof(argv[++i]);
		else if(strcmp(argv[i],"-inap_blk_ramp")==0) {inap_blk_ramp[0]=atof(argv[++i]); inap_blk_ramp[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-dsyn")==0) D_on_off=atof(argv[++i]);
		else if(strcmp(argv[i],"-ERoff")==0) ER_on_off=atof(argv[++i]);
		else if(strcmp(argv[i],"-CRC")==0) CRC=atof(argv[++i]);
		else if(strcmp(argv[i],"-spotstimIDs")==0) {cells[0]=atof(argv[++i]); cells[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-spotstimW")==0) strength=atof(argv[++i]);
		else if(strcmp(argv[i],"-spotstimT")==0) stim_time=atof(argv[++i]);
		else if(strcmp(argv[i],"-NUMspotstim")==0) Nstim=atof(argv[++i]);
		else if(strcmp(argv[i],"-ca_sqONOFF")==0) square_ON_OFF=atof(argv[++i]);
		else if(strcmp(argv[i],"-ca_sqtONtOFF")==0) {caTon=atof(argv[++i]); caToff=atof(argv[++i]);}
		else if(strcmp(argv[i],"-ca_sqAMP")==0) caAmp=atof(argv[++i]);
		else if(strcmp(argv[i],"-ca_sqP")==0) caP=atof(argv[++i]);
		else if(strcmp(argv[i],"-gmu")==0) gmu=atof(argv[++i]);
		else if(strcmp(argv[i],"-Agl")==0) A_gl=atof(argv[++i]);
		else if(strcmp(argv[i],"-Bgl")==0) B_gl=atof(argv[++i]);
		else if(strcmp(argv[i],"-w12")==0) Wpm2npm=atof(argv[++i]);
		else if(strcmp(argv[i],"-p22")==0) Pnpm2npm=atof(argv[++i]);
		else if(strcmp(argv[i],"-w22")==0) Wnpm2npm=atof(argv[++i]);
		else if(strcmp(argv[i],"-w21")==0) Wnpm2pm=atof(argv[++i]);
		else if(strcmp(argv[i],"-syn_step")==0) f_D=atof(argv[++i]);
		else if(strcmp(argv[i],"-tau_syn")==0) tau_D=atof(argv[++i]);
		else if(strcmp(argv[i],"-GSERCA")==0) gSERCA=atof(argv[++i]);
		else if(strcmp(argv[i],"-GER")==0) gER=atof(argv[++i]);
		else if(strcmp(argv[i],"-gldist")==0) sL_pct=atof(argv[++i]);
		else if(strcmp(argv[i],"-GINHIB")==0) GINHIB=atof(argv[++i]);
		else if(strcmp(argv[i],"-kSERCA")==0) kSERCA=atof(argv[++i]);
		else if(strcmp(argv[i],"-gL_ER")==0) gL_ER=atof(argv[++i]);
		else if(strcmp(argv[i],"-Nholo")==0) Nholostim=atof(argv[++i]);
		else if(strcmp(argv[i],"-stimseed")==0) stimseed=atof(argv[++i]);
		else if(strcmp(argv[i],"-stimdelay")==0) stim_delay=atof(argv[++i]);
		else if(strcmp(argv[i],"-stimadvance")==0) stim_advance=atof(argv[++i]);
		else if(strcmp(argv[i],"-musyndep")==0) musyndep=atof(argv[++i]);//#should be between 1 and 0.5
		else if(strcmp(argv[i],"-gSPK")==0) {gspk[0]=atof(argv[++i]); gspk[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-gAHP")==0) {gahp[0]=atof(argv[++i]); gahp[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-fixedhnap")==0) FIXEDHNAP=atof(argv[++i]);
		else if(strcmp(argv[i],"-rastonly")==0) RASTONLY=atof(argv[++i]);
		else if(strcmp(argv[i],"-nain")==0) {nain_ramp[0]=atof(argv[++i]); nain_ramp[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-IBfind")==0) IBFIND=atof(argv[++i]);
		else if(strcmp(argv[i],"-HYPOX")==0) {hypox_ramp[0]=atof(argv[++i]); hypox_ramp[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-fig")==0) fignum=atof(argv[++i]);
		else if(strcmp(argv[i],"-imposedPMfr")==0) imposedPMfr=atof(argv[++i]);
		else if(strcmp(argv[i],"-Blfreq")==0) period=1000/(atof(argv[++i]));
		else if(strcmp(argv[i],"-Blamp")==0) AMP=(atof(argv[++i]));
		else if(strcmp(argv[i],"-IBIamp")==0) {IBIamp[0]=atof(argv[++i]); IBIamp[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-refract")==0) refract=(atof(argv[++i]));
		else if(strcmp(argv[i],"-uniformdist")==0) uniformdist=atof(argv[++i]);
		else if(strcmp(argv[i],"-dahp")==0) hahp_const=atof(argv[++i]);
		else if(strcmp(argv[i],"-dspk")==0) h2spkht_const=atof(argv[++i]);
		else if(strcmp(argv[i],"-CAP")==0) c=atof(argv[++i]);
		else if(strcmp(argv[i],"-HYPOX2")==0) {hypox_ramp2[0]=atof(argv[++i]); hypox_ramp2[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-HYPOXspk")==0) {hypox_spk_ramp[0]=atof(argv[++i]); hypox_spk_ramp[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-CAPramp")==0) {Cm_ramp[0]=atof(argv[++i]); Cm_ramp[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-TEMP")==0) TEMP=atof(argv[++i]);
		else if(strcmp(argv[i],"-TEMPsteps")==0) TEMPsteps=atof(argv[++i]);
		else if(strcmp(argv[i],"-TEMPramp")==0) {TEMP_ramp[0]=atof(argv[++i]); TEMP_ramp[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-INa_dV12")==0) dV12=atof(argv[++i]);
		else if(strcmp(argv[i],"-INa_dV12_ramp")==0) {dV12_ramp[0]=atof(argv[++i]); dV12_ramp[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-INa_dV12_ramp2")==0) {dV12_ramp2[0]=atof(argv[++i]); dV12_ramp2[1]=atof(argv[++i]);}		
		else if(strcmp(argv[i],"-tauNapump")==0) tauNain=atof(argv[++i]);
		else if(strcmp(argv[i],"-tauNapump_ramp")==0) {tauNapump_ramp[0]=atof(argv[++i]); tauNapump_ramp[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-dynamicNain")==0) dynamicNain=atof(argv[++i]);
		else if(strcmp(argv[i],"-t_hnap_fix")==0) T_HNAP_FIX=atof(argv[++i]);
		else if(strcmp(argv[i],"-t_dsyn_fix")==0) T_DSYN_FIX=atof(argv[++i]);
		else if(strcmp(argv[i],"-t_syn_off")==0) T_SYN_OFF=atof(argv[++i]);
		else if(strcmp(argv[i],"-t_syn_mod")==0) {T_SYN_MOD=atof(argv[++i]); syn_mod_amp=atof(argv[++i]);}
		else if(strcmp(argv[i],"-t_mnap_fix")==0) T_MNAP_FIX=atof(argv[++i]);
		else if(strcmp(argv[i],"-rampstart")==0) RAMPSTART=atof(argv[++i]);
		else if(strcmp(argv[i],"-devmode")==0) DEVMODE=atof(argv[++i]);	
		else if(strcmp(argv[i],"-scale")==0) {I_scale=atof(argv[++i]); syn_scale=atof(argv[++i]); mINaP=atof(argv[++i]);}
		else if(strcmp(argv[i],"-syncut")==0) SYNCUT=atof(argv[++i]);
		else if(strcmp(argv[i],"-gspkramp")==0) GSPKRAMP=atof(argv[++i]);
		else if(strcmp(argv[i],"-BLKbursters")==0) blockburster=atof(argv[++i]);
		else if(strcmp(argv[i],"-SPIKETHRESH")==0) SPIKETHRESH=atof(argv[++i]);
		else if(strcmp(argv[i],"-kegLdep")==0) {B_gl=atof(argv[++i]);  A_gl=(8.5-((8.5-3.425)/4.05)*B_gl);}
		else strcpy(fname,argv[i]);
	}
	


//Calculating the scaling of for INaP
INaP_scale=mINaP*(I_scale-1)+1;



	if(seed !=0) srand (seed);
	int freq=int(DT/dt);


	pulseON = (1000/pulsefreq);
	Population pop(size,gca[0],gca[1],gleak[0],gleak[1],gnap[0],gnap[1],gcan[0],gcan[1],gchr2[0],gchr2[1],garch[0],garch[1],prob,wght,stimnum,gspk[0],gspk[1],gahp[0],gahp[1],seed);
    	if(ini) { ifstream is("ini"); is>>pop; }
	ofstream out(fname);
	out << setprecision(7);
	
	//////////////
	//////////////
	int sp=0;	
	int h[size];
	
	double id_ton_sil=0;	
	double N_silent=0;
	double N_tonic=0;
	int tonID[size];
	for(int i=0;i<size;i++) {tonID[i]=0;}
	int tonic_sp=0;
	
	int silentID[size];
	for(int i=0;i<size;i++) {silentID[i]=0;}
	int silent_sp=0;
	
	double last_sp_t[size];
	for(int i=0;i<size;i++) last_sp_t[i]=0;
	double irr=0;
        double PM_sp=0,nPM_sp=0;
        double t_last_stim=-5000;
        double t_Ca_last=0;
        double t_ca_on=0;
        double t_after_B;//time after last burst
        double pop_count[size];// array for counting the neurons that spike in a given time step
        for(int i=0;i<size;i++) {pop_count[i]=0;}
        double N_active=0;//Number of active neurons
        double tlast_pop[size];// array for recording the last spike time fore each neuron
        for(int i=0;i<size;i++) {tlast_pop[i]=0;}
        double mnap_fix_count=0;
        //////////////
	//////////////
        
	caToff=caP-caTon;
	//set seed for stimulation
	if(stimseed !=0) srand (stimseed);
	int count=0;
	while(count<Nholostim){
				pop.net[(rand() % size)].stimIDs=1;
				count=0;
				//count the number of stimulated neurons
				for(int i=0;i<size;i++){if(pop.net[i].stimIDs==1){count++;}}
				}
//	cerr<<Nholostim<<'\t'<<stimseed<<'\t'<<endl;
//	for(int i=0;i<size;i++){if(pop.net[i].stimIDs==1){cerr<<i<<'\t'<<pop.net[i].stimIDs<<'\t'<<endl;}}
	//set seed back to what it was
	if(seed !=0) srand (seed);

//	//write gLeak,gNaP,CellID, .ca file if t=0
	if(T==0 && fignum==0)	{
			for(int i=0;i<size;i++)
				{
				cerr<<pop.net[i].g_L<<'\t'<<pop.net[i].gNaP<<'\t'<<i<<'\t'<<gleak[0]<<'\t'<<gleak[1]<<endl;
				}
			}

//	//write Source ID and Target ID and synaptic weight to .ca file if t=0
	if(T==0 && fignum==6)	{
				cerr<<"Source,Target,Weight,gLeak,gNaP"<<endl;
				cout<<"id,label,gLeak,gNaP,popID"<<endl;
				//loop through each source
				for(int i=0;i<size*size;i++)
				{
				if(pop.w[i].weight>0){cerr<<pop.w[i].source<<','<<pop.w[i].target<<','<<pop.w[i].weight<<','<<pop.net[pop.w[i].source].g_L<<','<<pop.net[pop.w[i].source].gNaP<<endl;}
				}
				double popid=1;
				for(int i=0;i<size;i++)
				{
				if(i>99){popid=2;}
				cout<<i<<','<<popid<<','<<pop.net[i].g_L<<','<<pop.net[i].gNaP<<','<<popid<<endl;
				}
			}
	//write gLeak,gNaP,CellID,gspk .ca file if t=0
	if(T==0 && fignum==4)	{
			for(int i=0;i<size;i++)
				{
				cerr<<pop.net[i].g_L<<'\t'<<pop.net[i].gNaP<<'\t'<<i<<'\t'<<pop.net[i].gSPK<<'\t'<<seed<<endl;
				}
			}

	///////////////////////////////////////////////////////
	//Main Simulation Loop
	//////////////////////////////////////////////////////
	dFR_up=dFR_up=(AMP-IBIamp[1])/t_up;
	
	if(DEVMODE>1){GNAPBLK=DEVMODE*((gspk[0]+gspk[1])/2-6)+1;}
	for(double t=0;t<=T;t+=dt)
	{
		
	///////////////////
	//Imposed PM pop firing rate
	///////////////////
	if(imposedPMfr==1)	{
				tlast+=dt;
				//if out of IBI period and before middle of burst
				if(tlast>=(period-t_up-t_down) && tlast<(period-t_down)){PMfr+=dt*dFR_up;}
				//if out of IBI period and past middle of burst
				if(tlast>(period-t_down) && tlast<period){PMfr+=dt*dFR_down;}
				//If end of period
				if(tlast>=(period)){tlast=0; PMfr=IBIamp[0];}
				//if in percolation period of IBI
				if(tlast<((period-t_up-t_down)*refract)){PMfr=0;}
				//if in percolation period of IBI
				if(tlast>=((period-t_up-t_down)*refract) && tlast<(period-t_up-t_down)){ PMfr=IBIamp[0]+((tlast-(period-t_up-t_down)*refract)/((period-t_up-t_down)*(1-refract)))*(IBIamp[1]-IBIamp[0]);}
			  	}	
	
	////////////////////
	//END: Imposed PM pop firing rate
	///////////////////
		
		
		
		
		//Parameters that can vary during simulations
//		double drive=dr[0]+t/T*(dr[1]-dr[0]);//synaptic drive

		if(T_HNAP_FIX>0 && t>=T_HNAP_FIX){FIXEDHNAP=4;}
		if(mnap_fix_count>0){MNAP_FIX_onoff=0;}//Fixing mnap only needs to be set at one time step, so this line turns this back off after it has been turned on once
		if(T_MNAP_FIX>0 && t>=T_MNAP_FIX && mnap_fix_count<1){MNAP_FIX_onoff=1;mnap_fix_count=1;}//This line says that if t is >=T_HNAP_Fix and MNAP_FIX_onoff hasn't been turned on yet then turn it on
		if(T_DSYN_FIX>0 && t>=T_DSYN_FIX){D_on_off=0;}
		
		


		double drive=I_scale*dr[0];
		if (t>=20000){ drive=I_scale*(dr[0]+(t-20000)/(T-20000)*(dr[1]-dr[0]));}//synaptic drive with 20s delay
		double fcan= fc[0]+t/T*(fc[1]-fc[0]);//Can conducatance
		
		if (GSPKRAMP==1){scale_gspk=(0);}
		if (GSPKRAMP==1 && t>=50000){ scale_gspk=(0+(t-50000)/(T-50000));}//ramp in gspk w/ 20sec delay
		
		double fw=0;
		if(RAMPSTART==20000){fw=fcw[0];}
		if (t>=RAMPSTART){ fw=fcw[0]+(t-RAMPSTART)/(T-RAMPSTART)*(fcw[1]-fcw[0]);}//synaptic weights ramp with delay defined by RAMPSTART
		if(T_SYN_OFF>0 && t>=T_SYN_OFF){fw=0;}
		if(T_SYN_MOD>0 && t>=T_SYN_MOD){fw=syn_mod_amp;}
		
		//Synaptic cut experiments
		if(SYNCUT==1 && t<30000){fw=0;}
		if(SYNCUT==1 && t>=30000){fw=1;}

		
		double fsynCa=fsca[0]+t/T*(fsca[1]-fsca[0]);//Psynca
		IAPP=APPramp[0]+(t)/(T)*(APPramp[1]-APPramp[0]);//Ramp in applied current
		
		if (dV12_ramp[0]==dV12_ramp[1] && dV12_ramp[0]>0){dV12 = dV12_ramp[0];}
		
		if (t<20000  && dV12_ramp[0]!=dV12_ramp[1]){dV12 = dV12_ramp[0];}
		if (t>=20000 && t<=60000 && dV12_ramp[0]!=dV12_ramp[1] ){dV12 = dV12_ramp[0]+(t-20000)/(40000)*(dV12_ramp[1]-dV12_ramp[0]);}//ramp in INa V12 parameters
		if (dV12_ramp2[0]!=dV12_ramp2[1]){dV12 =0+dV12_ramp2[1]/(1+exp((45000-t)/4000));}
//		if (t>=20000 && t<=60000 && dV12_ramp[0]!=dV12_ramp[1] ){dV12 = dV12_ramp[0]+(t-20000)/(60000)*(dV12_ramp[1]-dV12_ramp[0]);}//ramp in INa V12 parameters
		
		if (inap_blk_ramp[0]!=inap_blk_ramp[1] && t>=100000){ GNAPBLK=inap_blk_ramp[0]+(t-100000)/(T-100000)*(inap_blk_ramp[1]-inap_blk_ramp[0]);}//synaptic drive with 20s delay
		
		Nain= nain_ramp[0];
		if (t>=20000){Nain = nain_ramp[0]+(t-20000)/(T-20000)*(nain_ramp[1]-nain_ramp[0]);}//ramp in Nain
//		if (t>=20000 && tauNapump_ramp[0]!=tauNapump_ramp[1]){tauNain = tauNapump_ramp[0]+(t-20000)/(T-20000)*(tauNapump_ramp[1]-tauNapump_ramp[0]);}//ramp in Nain
		if(dynamicNain==3)	{imposed_Nain=15+ (47.5-15)/(1+exp((72500-t)/5000));
//					dV12=dV12_ramp[0]+dV12_ramp[1]/(1+exp((40000-t)/4000));
					}
		
		if(dynamicNain==4)	{imposed_Nain=15 + ((47.5-15)/(1+exp((72500-t)/5000)) - (47.5-15)/(1+exp((100000-t)/5000)));
					}
		
		if (tauNapump_ramp[0]!=tauNapump_ramp[1]){tauNain = tauNapump_ramp[0]+ (tauNapump_ramp[1]-tauNapump_ramp[0])/(1+exp((60000-t)/1000));}
		
		
		if (TEMPsteps==0 && t>=20000 && t<=120000){TEMP = TEMP_ramp[0]+(t-20000)/(T-40000)*(TEMP_ramp[1]-TEMP_ramp[0]);}//ramp in TEMP
		if (TEMPsteps==1)
				{
				TEMP=21;
				DTDt=0;
				if(t>=50000 && t<70000){TEMP = 21+(t-50000)/(20000)*6; DTDt=6./20000;}
				if(t>=70000 && t<110000){TEMP = 27;DTDt=0;}
				if(t>=110000 && t<145000){TEMP = 27+(t-110000)/(35000)*10; DTDt=10./20000;}
				if(t>=145000 && t<185000){TEMP = 37;DTDt=0;}
				if(t>=185000 && t<195000){TEMP = 37+(t-185000)/(10000)*3; DTDt=3./20000;}
				if(t>=195000){TEMP = 40;DTDt=0;}
				}
		if (TEMPsteps==2)
				{
				DTDt=0.0001*(TBATH-TEMP);
				TEMP=TEMP+dt*DTDt;
				if(t>=50000 && t<70000){TBATH = 21+(t-50000)/(20000)*6;}
				if(t>=70000 && t<110000){TBATH = 27;}
				if(t>=110000 && t<145000){TBATH = 27+(t-110000)/(35000)*10; }
				if(t>=145000 && t<185000){TBATH = 37;}
				if(t>=185000 && t<195000){TBATH = 37+(t-185000)/(10000)*3;}
				if(t>=195000){TBATH = 40;}
				}
		if (TEMPsteps==3)
				{
				DTDt=0.0001*(TBATH-TEMP);
				TEMP=TEMP+dt*DTDt;
				if(t>=00000 && t<7000){TBATH = 27;}
				if(t>=70000 && t<100000){TBATH = 27+(t-70000)/(50000)*10; }
				if(t>=100000 && t<180000){TBATH = 37;}
				if(t>=180000 && t<210000){TBATH = 37-(t-180000)/(50000)*10; }
				if(t>=210000){TBATH = 27;}
				}
		if (TEMPsteps==4)
				{
//				DTDt=0.0000025*(TBATH-TEMP)/dt;//the 0.0000025 is the rate constant determining how fast it decays to the bath temp
				DTDt=0.0008*(TBATH-TEMP)/dt;//the 0.0000025 is the rate constant determining how fast it decays to the bath temp
				double TBATH_last=TBATH;
				double Cm_last=c;
				TEMP=TEMP+dt*DTDt;
				if(t<200000){TBATH = 27+ 10/(1+exp((120000-t)/5500));}//DTDt=0*(1+exp((85000-t)/3000)-1+exp((85000-(t-dt))/3000))/(1000*dt);}
				if(t>=200000){TBATH = 37- 10/(1+exp((270000-t)/5500));}//DTDt= -0*(10/(1+exp((195000-t)/3000))-10/(1+exp((195000-(t-dt))/3000)))/(1000*dt);}
				DBATHTDt=(TBATH-TBATH_last)/(dt);
				DCmDt=(c-Cm_last)/(dt);
//				TEMP=TBATH;
				}
				
		if (TEMPsteps==5)
				{
//				DTDt=0.0000025*(TBATH-TEMP)/dt;//the 0.0000025 is the rate constant determining how fast it decays to the bath temp
				DTDt=0.0008*(TBATH-TEMP)/dt;//the 0.0000025 is the rate constant determining how fast it decays to the bath temp
				double TBATH_last=TBATH;
				double Cm_last=c;
				TEMP=TEMP+dt*DTDt;
				if(t<200000){TBATH = 27+ 10/(1+exp((120000-t)/10000));}//DTDt=0*(1+exp((85000-t)/3000)-1+exp((85000-(t-dt))/3000))/(1000*dt);}if(t>=200000){TBATH = 37- 10/(1+exp((270000-t)/5500));}//DTDt= -0*(10/(1+exp((195000-t)/3000))-10/(1+exp((195000-(t-dt))/3000)))/(1000*dt);}
				DBATHTDt=(TBATH-TBATH_last)/(dt);
				DCmDt=(c-Cm_last)/(dt);
//				TEMP=TBATH;
				}
				
				
		//Set temp dep of capacitance: represents a 0.3% increase per degree C adapted by Plaksin et al 2018 Physical review X
		if(TEMP!=27){c=.108*TEMP+33.08;} //1.0015*TEMP+8.96;
		
		
//		if (t>=20000){breakdown = hypox_ramp2[0]+(t-20000)/(T-20000)*(hypox_ramp2[1]-hypox_ramp2[0]);}//ramp in Nain
//		if(t>=20000 && t<=70000 && (hypox_ramp[0]!=1 || hypox_ramp[1]!=1)) hypox=hypox_ramp[0]+(t-20000)/(T-(70000-20000))*(hypox_ramp[1]-hypox_ramp[0]);//decrease ISPK and IAHP during simulated hypoxia
		if(t>=44000 && t<=68000 &&  (hypox_ramp[0]!=1 || hypox_ramp[1]!=1)) hypox=hypox_ramp[0]+(t-44000)/((24000))*(hypox_ramp[1]-hypox_ramp[0]);
		if(t>=92000 && t<=116000 && (hypox_ramp[0]!=1 || hypox_ramp[1]!=1)) hypox=hypox_ramp[1]-(t-92000)/((24000))*(hypox_ramp[1]-hypox_ramp[0]);		
		if (t>=5000){hypox_spk = hypox_spk_ramp[0]+(t-5000)/(T-5000)*(hypox_spk_ramp[1]-hypox_spk_ramp[0]);}//{hypox_spk=hypox_spk_ramp[1];}
		if (t>=20000 && Cm_ramp[0]!=Cm_ramp[1]){c = Cm_ramp[0] +(t-20000)/(T-20000)*(Cm_ramp[1]-Cm_ramp[0]);}//ramp in Cap
		
		if(FIXEDHNAP==1) {hnapfixedvalue=0.425/(1+exp((drive-.31)/.042))+.095;}//{hnapfixedvalue=-2.*drive+0.9;}
//		if(FIXEDHNAP==2) {hnapfixedvalue=-1.875*drive+0.9;}
		if(FIXEDHNAP==2) {hnapfixedvalue=-1.4521*drive+0.772035;}
		if(FIXEDHNAP==3) {hnapfixedvalue=-1.72374*drive+0.826245;}
		//Ca square wave 
		if(square_ON_OFF==1)
			{	
			Iapp_ca=0;
			if(t_Ca_last>caToff){Iapp_ca=caAmp; t_ca_on+=dt;}
			if(t_ca_on>caTon){t_Ca_last=0; t_ca_on=0; Iapp_ca=0;}
			t_Ca_last+=dt;
			}
		
		//spot stimulation
//		if(t_last_stim>stim_time && t_after_B>=stim_delay && t_after_B<(stim_delay+dt) && stimcount<Nstim){
//							     stim_ON=1;stim_W=strength;stimcount+=1;t_last_stim=0;
//							     stim_delay = stim_delay+stim_advance;
////							     for(int i=0;i<(cells[1]-cells[0]);i++) cerr<<(t/1000)<<'\t'<<cells[0]+i<<'\t'<<-1<<'\t'<<endl;
//							     for(int i=0;i<size;i++){if(pop.net[i].stimIDs==1){cerr<<(t/1000)<<'\t'<<i<<'\t'<<-1<<'\t'<<endl;}}
//							    }					    
		t_after_B+=dt;
		stim_W*=exp(-dt/100);
		t_last_stim+=dt;
		//END spot stimulation


		//Array for cells that spike during timestep
		int spk[size];
	

		//Find silent and tonic neurons
		if(SYNCUT==1 && t>=20000 && id_ton_sil==0)	{
								id_ton_sil=1;
								for(int i=0;i<size;i++)	{
											if((t-tlast_pop[i])>10000){silentID[i]=1; N_silent+=1;}//if the neuron hasn't spiked in the last 10sec it's silent
											if((t-tlast_pop[i])<=10000){tonID[i]=1; N_tonic+=1;}//if the neuron has spiked in the last 10sec it's tonic
											}
								}

		//CUT SYN FROM TONIC TO SILENT NEURONS	
		
		if(SYNCUT==1 && t>=80000)
		{
		//loop through each silent neuron
		for(int i=0;i<size*size;i++)for(int n=0;n<size;n++)
							{
//							//if source=silent && target =tonic - i.e. cut synapses from silent back to tonic neurons
							if( silentID[pop.w[i].source]==1 &&  tonID[pop.w[i].target]==1){pop.w[i].weight=0;}
							}
		SYNCUT=2;//increasing this to 2 will make it so the simulation does go through this look every time. 

		}
		if(SYNCUT==2 && t>=130000)
		{
		//loop through each silent neuron
		for(int i=0;i<size*size;i++)for(int n=0;n<size;n++)
							{
//							//if source ==tonic && target=silent -- i.e. cut synapes from tonic neurons to silent neurons
							if( tonID[pop.w[i].source]==1 &&  silentID[pop.w[i].target]==1){pop.w[i].weight=0;}
							}
		SYNCUT=3;//increasing this to 2 will make it so the simulation does go through this look every time. 

		}
		
		

		sp+=pop.step(dt,drive,spk,fcan,fw,fsynCa,irr);

		if(T>0)
		{
		for(int i=0;i<size;i++)
			{
				h[i]+=spk[i];
				if(h[i]>=1)
				{
					pop_count[i]=1;
				
					cout<<(t/1000)<<'\t'<<i<<'\t'<<pop.net[i].drivestart+drive<<'\t'<<pop.net[i].hp<<'\t'<<pop.net[i].gNaP<<'\t'<<pop.net[i].g_L<<'\t'<<pop.net[i].gSPK<<'\t'<<pop.net[i].gAHP<<'\t'<<Nain<<'\t'<<(t-tlast_pop[i])<<'\t'<<pop.net[i].D_syn<<'\t'<<kbath<<'\t'<<pop.net[0].EK<<'\t'<<pop.net[0].E_leak<<'\t'<<c<<'\t'<<TEMP<<'\t'<<DCmDt<<endl;//uncomment this line to write spike times to .hst file for plotting raster plot
					tlast_pop[i]=t;
					if(tonID[i]==1){tonic_sp+=1;}
					if(silentID[i]==1){silent_sp+=1;}
				}

			}
		}
		
		

		
		for(int i=0;i<size;i++) h[i]=0;
		


		if(int(t/dt)%freq==0) 
		{	
			//count the number of neurons that spike in this bin.
			for(int i=0;i<size;i++){N_active+=pop_count[i];}

			if(size>2)	{
					double meanvm=0;
					double meanhnap=0;
					double meanmnap=0;
					double meansynd=0;
					double meanICm=0;
					double meanNain=0;
					double meanPump=0;
					double meantlast_ap=0;
					double meanFIXMNAP=0;
					
					double meanINAP=0;
					double meanISYN=0;
					
					
					
					
					for(int i=0;i<size;i++){meanINAP+=(pop.net[i].gNaP*pop.net[i].mp*pop.net[i].hp*(pop.net[i].v-pop.net[i].ENa))/size; meanISYN+=(pop.net[i].gsyn*(pop.net[i].v-Esyn))/size; meanvm+=pop.net[i].v/size;meanhnap+=pop.net[i].hp/size;meanmnap+=pop.net[i].mp/size;meanICm+=pop.net[i].ICm/size;meanNain+=pop.net[i].Nain_test/size;meanPump+=pop.net[i].Ipump/size;meansynd+=pop.net[i].D_syn/size;meantlast_ap+=pop.net[i].tlast_AP/size; meanFIXMNAP+=pop.net[i].FIXMNAP/size;}
					if(RASTONLY==0 && fignum!=6 && size>2){out<<(t/1000)<<'\t'<<sp/(.001*DT*size)<<'\t'<<meanvm<<'\t'<<meanhnap<<'\t'<<drive<<'\t'<<fw<<'\t'<<gahp[0]<<'\t'<<Nain<<'\t'<<pop.net[0].ENa<<'\t'<<dV12<<'\t'<<N_active<<'\t'<<TEMP<<'\t'<<gNaF<<'\t'<<meanICm<<'\t'<<DTDt<<'\t'<<TBATH<<'\t'<<c<<'\t'<<DBATHTDt<<'\t'<<DCmDt<<'\t'<<meanNain<<'\t'<<meanPump<<'\t'<<meansynd<<'\t'<<meanmnap<<'\t'<<GNAPBLK<<'\t'<<meantlast_ap<<'\t'<<mnap_fix_count<<'\t'<<meanFIXMNAP<<'\t'<<MNAP_FIX_onoff<<'\t'<<silent_sp/(.001*DT*N_silent)<<'\t'<<tonic_sp/(.001*DT*N_tonic)<<'\t'<<meanINAP<<'\t'<<meanISYN<<endl;}
					if(RASTONLY==0 && fignum==6){out<<(t/1000)<<'\t'<<sp/(.001*DT*size)<<'\t'<<meanvm<<'\t'<<meanhnap<<'\t'<<drive<<'\t'<<gspk[0]<<'\t'<<gahp[0]<<'\t'<<Nain<<'\t'<<pop.net[0].ENa<<'\t'<<hypox<<'\t'<<N_active<<endl;}
					
					}
			if(size==2){out<<(t/1000)<<'\t'<<pop.net[0].v<<'\t'<<pop.net[1].v<<'\t'<<IAPP<<'\t'<<endl;}
			if(size==1){ out<<(t/1000)<<'\t'<<pop.net[0].v<<'\t'<<drive<<'\t'<<pop.net[0].Cain<<'\t'<<pop.net[0].caER<<'\t'<<kbath<<'\t'<<pop.net[0].EK<<'\t'<<pop.net[0].E_leak<<'\t'<<Nain<<'\t'<<pop.net[0].ENa<<'\t'<<c<<'\t'<<TEMP<<'\t'<<gNaF<<endl;}
			sp=0;
			tonic_sp=0;
			silent_sp=0;
			N_active=0;
			for(int i=0;i<size;i++){pop_count[i]=0;}
	
		}
		
		//Writing individual neurons to file
		if(int(t/dt)%int(0.1/dt)==0 && t>=20000 && t<=200000 && RASTONLY==0 && fignum==6) 
			{
			cerr<<(t/1000)<<'\t';
			for(int n=30;n<40;n++) {cerr<<pop.net[n].v<<'\t';}
			for(int n=130;n<140;n++) {cerr<<pop.net[n].v<<'\t';}
			for(int n=30;n<40;n++) {cerr<<pop.net[n].gsyn<<'\t';}
			for(int n=130;n<140;n++) {cerr<<pop.net[n].gsyn<<'\t';}
			cerr<<endl;
			}
		if(int(t/dt)%int(0.1/dt)==0 && t>=20000  && fignum==5.8)
			{
			cerr<<(t/1000)<<'\t';
			for(int n=0;n<20;n++) {cerr<<pop.net[n].v<<'\t';}
			cerr<<endl;
			}
			
		
		if(int(t/dt)%int(0.5/dt)==0 && t>=00000  && fignum==4 && size>2)
			{
			cerr<<(t/1000)<<'\t';
			for(int n=0;n<20;n++) {cerr<<pop.net[n].v<<'\t';}
			for(int n=0;n<20;n++) {cerr<<pop.net[n].hp<<'\t';}
			cerr<<endl;
			}
		if(int(t/dt)%int(0.5/dt)==0 && t>=00000  && fignum==4 && size==2)
			{
			cerr<<(t/1000)<<'\t';
			for(int n=0;n<2;n++) {cerr<<pop.net[n].v<<'\t';}
			for(int n=0;n<2;n++) {cerr<<pop.net[n].hp<<'\t';}
			cerr<<endl;
			}
		if(int(t/dt)%int(0.5/dt)==0 && t>=00000  && fignum==6.1 && size>2)
			{
			cerr<<(t/1000)<<'\t';
			for(int n=0;n<20;n++) {cerr<<pop.net[n].v<<'\t';}
//			for(int n=0;n<20;n++) {cerr<<pop.net[n].hp<<'\t';}
			cerr<<endl;
			}

	}
   	ofstream os("ini"); os<<pop;
	return 0;
}

