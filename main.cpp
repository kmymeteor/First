#include<iostream>
#include<fstream>
#include<stdio.h>
#include<time.h>
#include<cmath>
#include<stdlib.h>
#include <deque>
#include<queue>
#include<float.h>
#include<random>
#include<algorithm>

using namespace std;

#define INFTY DBL_MAX

long long int Ntotal=5000000,NNtotal=10*Ntotal,randn,mod = 1;
int arrivestate,servicestate;
double bur;

//double arrTon=20,arrToff=10,serTon=10,serToff=0;
double arrTon=0,arrToff=20,serTon=20,serToff=10;
double arriverenewT,servicenewT;
double t;
double buflength;
double sumlength;
double meanlength;
double meandelay;
//double lambda0=10.125,lambda1=0.125;
//double mu0=8.125,mu1=3.125;
//double sizeB=0.01;
double lambda0=11,lambda1=5;
double mu0=10,mu1=0;
double sizeB=0;
double sizeBd=20;
double aver_mu=(serTon*mu0+serToff*mu1)/(serTon+serToff);

void init()
{
    int i;
    randn = time(NULL);
    for(i=0;i<31;i++)
    mod *= 2;
    mod = mod - 1;
    //srand(time(NULL));
}

long long int ran(long long int mod)
{
    return randn = 7*7*7*7*7* randn% mod;
}

double exprand(double para)             //exponential distribution random variable
{
    return -1/para*log((double)ran(mod)/(mod-1));
}
double prtrand(double xm,double a)
{
    return xm/(pow((double)ran(mod)/(mod-1),1/a));
}

int main()
{
    init();
    ofstream outfile("meandelay.txt",ios::out);
    //for(bur=0.002;bur<0.011;bur+=0.002)
    //{

     //bur=3;
        //arrTon=10*bur,arrToff=10*bur;
    arriverenewT=exprand(1.0/arrToff);servicenewT=exprand(1.0/serToff);
    arrivestate=0;servicestate=0;
    sumlength=0;buflength=0;meanlength=0;meandelay=0;t=0;
    //std::default_random_engine generator;
    //std::exponential_distribution<double> arroff_exp(1/arrToff);
    //std::exponential_distribution<double> arron_exp(1/arrTon);
    //std::exponential_distribution<double> seroff_exp(1/serToff);
    //std::exponential_distribution<double> seron_exp(1/serTon);
    //arriverenewT=arroff_exp(generator);servicenewT=seroff_exp(generator);
    for(int long long i=0;i<5*NNtotal;i++){

        if(arrivestate==0){

            //arrivestate=1;
            if(servicestate==0){//00
                //servicestate=1;
                if((arriverenewT-t)<=(servicenewT-t)){
                    if(mu1<lambda1){
                        sumlength+=(arriverenewT-t)*(2*buflength+(arriverenewT-t)*(lambda1-mu1))/2;
                        buflength+=(arriverenewT-t)*(lambda1-mu1);
                    }
                    else{
                        if(buflength<(arriverenewT-t)*(mu1-lambda1)){
                            sumlength+=buflength*buflength/(2*(mu1-lambda1));
                            buflength=0;
                        }
                        else{
                            sumlength+=(arriverenewT-t)*(2*buflength-(arriverenewT-t)*(mu1-lambda1))/2;
                            buflength=buflength-(arriverenewT-t)*(mu1-lambda1);
                        }
                    }
                    t=arriverenewT;
                    arrivestate=1;
                    arriverenewT=t+exprand(1.0/arrTon);
                    //arriverenewT=t+10;
                }
                else{
                    if(mu1<lambda1){
                        sumlength+=(servicenewT-t)*(2*buflength+(servicenewT-t)*(lambda1-mu1))/2;
                        buflength+=(servicenewT-t)*(lambda1-mu1);
                    }
                    else{
                        if(buflength<(servicenewT-t)*(mu1-lambda1)){
                            sumlength+=buflength*buflength/(2*(mu1-lambda1));
                            buflength=0;
                        }
                        else{
                            sumlength+=(servicenewT-t)*(2*buflength-(servicenewT-t)*(mu1-lambda1))/2;
                            buflength=buflength-(servicenewT-t)*(mu1-lambda1);
                        }
                    }
                    t=servicenewT;
                    servicestate=1;
                    servicenewT=t+exprand(1.0/serTon);
                    //servicenewT=t+30;
                }
                //servicestate=1;
            }
            else{//01
                //servicestate=0;
                if((arriverenewT-t)<=(servicenewT-t)){
                    if(buflength<(arriverenewT-t)*(mu0-lambda1)){
                        sumlength+=buflength*buflength/(2*(mu0-lambda1));
                        buflength=0;
                    }
                    else{
                        sumlength+=(arriverenewT-t)*(2*buflength-(arriverenewT-t)*(mu0-lambda1))/2;
                        buflength=buflength-(arriverenewT-t)*(mu0-lambda1);
                    }
                    t=arriverenewT;
                    arrivestate=1;
                    arriverenewT=t+exprand(1.0/arrTon);
                    //arriverenewT=t+10;
                }
                else{
                    if(buflength<(servicenewT-t)*(mu0-lambda1)){
                        sumlength+=buflength*buflength/(2*(mu0-lambda1));
                        buflength=0;
                    }
                    else{
                        sumlength+=(servicenewT-t)*(2*buflength-(servicenewT-t)*(mu0-lambda1))/2;
                        buflength=buflength-(servicenewT-t)*(mu0-lambda1);
                    }
                    t=servicenewT;
                    servicestate=0;
                    servicenewT=t+exprand(1.0/serToff);
                    //servicenewT=t+10;
                }
                //servicestate=0;
            }
            //arrivestate=1;
        }
        else{
            //arrivestate=0;
            if(servicestate==0){
                    //servicestate=1;
               //10
                if((arriverenewT-t)<=(servicenewT-t)){
                    sumlength+=(arriverenewT-t)*(2*buflength+(arriverenewT-t)*(lambda0-mu1))/2;
                    buflength+=(arriverenewT-t)*(lambda0-mu1);
                    t=arriverenewT;
                    arrivestate=0;
                    arriverenewT=t+exprand(1.0/arrToff);
                    //arriverenewT=t+10;
                }
                else{
                    sumlength+=(servicenewT-t)*(2*buflength+(servicenewT-t)*(lambda0-mu1))/2;
                    buflength+=(servicenewT-t)*(lambda0-mu1);
                    t=servicenewT;
                    servicestate=1;
                    servicenewT=t+exprand(1.0/serTon);
                    //servicenewT=t+30;
                }
               // servicestate=1;
            }
            else{
                //servicestate=0;
                //11
                if((arriverenewT-t)<=(servicenewT-t)){
                    sumlength+=(arriverenewT-t)*(2*buflength+(arriverenewT-t)*(lambda0-mu0))/2;
                    buflength+=(arriverenewT-t)*(lambda0-mu0);
                    t=arriverenewT;
                    arrivestate=0;
                    arriverenewT=t+exprand(1.0/arrToff);
                    //arriverenewT=t+10;
                }
                else{
                    sumlength+=(servicenewT-t)*(2*buflength+(servicenewT-t)*(lambda0-mu0))/2;
                    buflength+=(servicenewT-t)*(lambda0-mu0);
                    t=servicenewT;
                    servicestate=0;
                    servicenewT=t+exprand(1.0/serToff);
                    //servicenewT=t+10;
                }
               // servicestate=0;
            }
            //arrivestate=0;
        }
    }
    meanlength=sumlength/t;
    //meandelay=1000*(meanlength-sizeB)/aver_mu;
    meandelay=1000*meanlength/aver_mu;
    cout<<"sumlength="<<sumlength<<"  t="<<t<<"  meanlength="<<meanlength<<"  meandelay= "<<meandelay<<endl;
    outfile<<meandelay<<endl;
    //}
}
