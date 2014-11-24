//#include <QtGui/QApplication>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include "device.h"
#include <usb.h>

long m_DeviceHandle;
unsigned int l_Time=0;
unsigned short m_NrPixels;
double m_pLambda[MAX_NR_PIXELS];
double l_pSpectrum[MAX_NR_PIXELS];

int t_exp=100;
int n_avg=10;
int dark=0;
int smoth=0;
int interact(usb_dev_handle *dev)
{
    char comm[100],name[100];
    int inp,cnt;
    inp=0;
    cnt=1;
    while (cnt>=1) {
        cnt+=1;
        if (fgets(comm,100,stdin)==NULL) return -1;
        printf("# got %s",comm);
        if (sscanf(comm,"EXP %i\n",&t_exp)==1) printf("set exposure\n");
        else if (sscanf(comm,"AVG %i\n",&n_avg)==1) printf("set averaging\n");
        else if (sscanf(comm,"SMO %i\n",&smoth)==1) printf("set smoothing\n");
        else if (sscanf(comm,"DYD %i\n",&dark)==1) printf("set dyn. dark\n");
        else if (sscanf(comm,"MEAS %s\n",name)==1) {
            measure(dev,t_exp,1,n_avg,smoth,dark);
            save(name,inp);
            //strcpy(comm,"");
        }
        else if (strcmp(comm,"END\n")==0) break;
        else cnt-=1;
        //fflush(stdout);
    }
    return cnt;
}

int main(int argc, char *argv[])
{
	//QAppliction app(argc, argv);
	//Qtdemo *w = new Qtdemo;
    usb_dev_handle *dev=NULL;
    float fit[10];
	init(dev,fit);//show();
    interact(dev);
    //if (argc>1) save(argv[1],0);
    //else save("data1.txt",0);
    usb_close(dev);
    //w->on_CloseCommBtn_clicked();
	//app.connect(&app, SIGNAL(lastWindowClosed()), &app, SLOT(quit()));
	return 1;//app.exec();
}
