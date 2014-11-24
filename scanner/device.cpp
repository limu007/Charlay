#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include "device.h"
#include <usb.h>


usb_dev_handle *open_dev(void);

usb_dev_handle *open_dev(void)
{
	struct usb_bus *bus;
	struct usb_device *dev;

	for(bus = usb_get_busses(); bus; bus = bus->next) 
	{
		for(dev = bus->devices; dev; dev = dev->next) 
		{
			if(dev->descriptor.idVendor == USB_VENDOR
			&& dev->descriptor.idProduct == USB_PRODUCT)
			{
				return usb_open(dev);
			}
		}
	}
	return NULL;
}

float swapsingle(float floatin)
{
	union s
	{
		char sa[4];
		float res;
	} temp;
	temp.res = floatin;
	s floatout;
	for (int teller=0; teller<4; ++teller){
		floatout.sa[teller]=temp.sa[3-teller];
	}
    return floatout.res;
}

uint32 swap32(uint32 uint32in)
{
	union s
	{
		char sa[4];
		uint32 res;
	} temp;
	temp.res = uint32in;
	s uint32out;
	for (int teller=0; teller<4; ++teller){
		uint32out.sa[teller]=temp.sa[3-teller];
	}
    return uint32out.res;
}

uint16 swap16(uint16 uint16in)
{
	union s
	{
		char sa[2];
		uint16 res;
	} temp;
	temp.res = uint16in;
	s uint16out;
	for (int teller=0; teller<2; ++teller){
		uint16out.sa[teller]=temp.sa[1-teller];
	}
    return uint16out.res;
}

int init(usb_dev_handle *dev,float *fit) {
    //float fit[NR_WAVELEN_POL_COEF];
	unsigned char bufout[6];
	unsigned char bufout2[6];
	usb_init();
	usb_find_busses();
	usb_find_devices();
	// struct usb_device *devstruc;

	if (!(dev = open_dev()))
		return -1;//if (doqt>0) QMessageBox::critical(this, tr("QT4 Demo"),tr("Device not found"),QMessageBox::Ok, QMessageBox::Ok);
	if (usb_set_configuration(dev, 1) < 0){
		//if (doqt>0) QMessageBox::critical(this, tr("QT4 Demo"),tr("Error in usb_set_configuration"),QMessageBox::Ok, QMessageBox::Ok);
		usb_close(dev);
	}
	if (usb_claim_interface(dev, 0) < 0){
		//if (doqt>0) QMessageBox::critical(this, tr("QT4 Demo"),	tr("Error in usb_claim_interface"), QMessageBox::Ok, QMessageBox::Ok);
		usb_close(dev);
	}
	bufout[00]=0x20;
	bufout[01]=0x00;
	bufout[02]=0x02;  // length of the command
	bufout[03]=0x00;
	bufout[04]=0x13;  // get_ident
	bufout[05]=0x00;

	if (usb_bulk_write(dev,0x02,(char *) bufout,sizeof(bufout),5000) < 0)
		return -2;//if (doqt>0) QMessageBox::critical(this, tr("QT4 Demo"),tr("Bulk write failed"),QMessageBox::Ok, QMessageBox::Ok);
	AvsIdentityType Avs_Id;
	if (usb_bulk_read(dev,0x86,(char*) &Avs_Id,sizeof(Avs_Id),5000) < 0) { 
		return -3;//if (doqt>0) QMessageBox::critical(this, tr("QT4 Demo"),tr("Bulk read failed"),QMessageBox::Ok, QMessageBox::Ok);
	} else {
        //if (doqt>0) SerialNum->setText(QString(Avs_Id.SerialNumber).left(9));
        printf("device Id: %s\n",Avs_Id.SerialNumber);
    }
	bufout2[00]=0x20;
	bufout2[01]=0x00;
	bufout2[02]=0x02;  // length of the command
	bufout2[03]=0x00;
	bufout2[04]=0x01;  // get_device_configuration
	bufout2[05]=0x00;
	if (usb_bulk_write(dev,0x02,(char *)bufout2,sizeof(bufout2),5000) < 0)
		return -4;//if (doqt>0) QMessageBox::critical(this, tr("QT4 Demo"),tr("Bulk write failed"),QMessageBox::Ok, QMessageBox::Ok);
	DeviceConfigType devcon;
	if (usb_bulk_read(dev,0x86,(char*) &devcon,sizeof(devcon),5000) < 0) {
		return -5;//if (doqt>0) QMessageBox::critical(this, tr("QT4 Demo"),tr("Bulk read failed"),QMessageBox::Ok, QMessageBox::Ok);
	} else {
		for (int teller=0;teller<NR_WAVELEN_POL_COEF;++teller) {
		fit[teller]=swapsingle(devcon.m_Detector.m_aFit[teller]);
	}
	extern unsigned short m_NrPixels;
	extern double m_pLambda[MAX_NR_PIXELS];
	m_NrPixels=swap16(devcon.m_Detector.m_NrPixels);
    printf("Init OK (%i pixels)\n",m_NrPixels);
	for (int teller=0;teller<m_NrPixels;++teller) {
		m_pLambda[teller]=	fit[0] +
						fit[1]*teller*1.0 +
						fit[2]*teller*teller*1.0 +
                        fit[3]*teller*teller*teller*1.0 +
                        fit[4]*teller*teller*teller*teller*1.0;
	  }
	}
	usb_release_interface(dev,0);
	//usb_close(dev);
	usb_reset(dev);
	//is_open=1;
    return 0;
};


int measure(usb_dev_handle *dev,double inttime=100,short l_NrOfScans=1,uint32 numavg=10,int smooth=0,int dark=0) {
 
	unsigned char bufin[63483];
	unsigned char bufout2[8];
	unsigned char bufout3[6];
	//uint32 size=0;
	//bool ok;
	SendMeasConfigType l_PrepareMeasData;
	sony_single_measdatatype sony_single_meas;
	sony_multi_measdatatype sony_multi_meas;

	l_PrepareMeasData.prefix[00]=0x20;
	l_PrepareMeasData.prefix[01]=0x00;
	l_PrepareMeasData.prefix[02]=0x2B;   // length of the command
	l_PrepareMeasData.prefix[03]=0x00;
	l_PrepareMeasData.prefix[04]=0x05;   // prepare_measurement
	l_PrepareMeasData.prefix[05]=0x00;
	l_PrepareMeasData.m_Meas.m_StartPixel = 0;
	extern unsigned short m_NrPixels;
	l_PrepareMeasData.m_Meas.m_StopPixel = swap16(m_NrPixels-1);
	//QLocale::setDefault(QLocale::Dutch);
	l_PrepareMeasData.m_Meas.m_IntegrationTime = swapsingle(inttime);
	l_PrepareMeasData.m_Meas.m_IntegrationDelay = 0;
	l_PrepareMeasData.m_Meas.m_NrAverages = swap32(numavg);
    if (dark>0) {
    	l_PrepareMeasData.m_Meas.m_CorDynDark.m_Enable = 1;
	    l_PrepareMeasData.m_Meas.m_CorDynDark.m_ForgetPercentage = dark;
    } else l_PrepareMeasData.m_Meas.m_CorDynDark.m_Enable = 0;
	l_PrepareMeasData.m_Meas.m_Smoothing.m_SmoothPix = smooth;
	l_PrepareMeasData.m_Meas.m_Smoothing.m_SmoothModel = 0;
	l_PrepareMeasData.m_Meas.m_SaturationDetection = 0;
	l_PrepareMeasData.m_Meas.m_Trigger.m_Mode = 0;
	l_PrepareMeasData.m_Meas.m_Trigger.m_Source = 0;
	l_PrepareMeasData.m_Meas.m_Trigger.m_SourceType = 0;
	l_PrepareMeasData.m_Meas.m_Control.m_StrobeControl = 0;
	l_PrepareMeasData.m_Meas.m_Control.m_LaserDelay = 0;
	l_PrepareMeasData.m_Meas.m_Control.m_LaserWidth = 0;
	l_PrepareMeasData.m_Meas.m_Control.m_LaserWaveLength = 0;
	l_PrepareMeasData.m_Meas.m_Control.m_StoreToRam = 0;

	usb_init();
	usb_find_busses();
	usb_find_devices();

    // usb_dev_handle *dev = NULL; /* the device handle */
	if (!(dev = open_dev()))
		return -1; //QMessageBox::critical(this, tr("QT4 Demo"),tr("Device not found"),QMessageBox::Ok, QMessageBox::Ok);
	if (usb_set_configuration(dev, 1) < 0){
    	usb_close(dev);
		return -2; //QMessageBox::critical(this, tr("QT4 Demo"),tr("Error in usb_set_configuration"),QMessageBox::Ok, QMessageBox::Ok);
	}
    printf("data prepared\n");
    fflush(stdout);
	if (usb_claim_interface(dev, 0) < 0){
		//QMessageBox::critical(this, tr("QT4 Demo"),tr("Error in usb_claim_interface"),QMessageBox::Ok, QMessageBox::Ok);
    	usb_close(dev);
	}
	int pmsize = sizeof(l_PrepareMeasData);
	if (usb_bulk_write(dev,0x02,(char *)&l_PrepareMeasData,pmsize,5000) != pmsize)
		return -3; //QMessageBox::critical(this, tr("QT4 Demo"),tr("Bulk write failed"),QMessageBox::Ok, QMessageBox::Ok);
	int retval;
	retval = usb_bulk_read(dev,0x86,(char*) bufin,6,5000);
	if (retval < 0) {
			return -4; //QMessageBox::critical(this, tr("QT4 Demo"),tr("Bulk read failed (prepare_measurement)"),QMessageBox::Ok, QMessageBox::Ok);
	} else {
		if ((retval!= 6 ) || (bufin[4]!=0x85))
			return -5;//QMessageBox::critical(this, tr("QT4 Demo"),tr("Error in prepare_measurement"),QMessageBox::Ok, QMessageBox::Ok);
	}
	bufout2[00]=0x20;
	bufout2[01]=0x00;
	bufout2[02]=0x04;   // length of the command
	bufout2[03]=0x00;
	bufout2[04]=0x06;   // start_measurement
	bufout2[05]=0x00;
	bufout2[06]=HIBYTE(l_NrOfScans);
	bufout2[07]=LOBYTE(l_NrOfScans);

	if (usb_bulk_write(dev,0x02,(char *)bufout2,sizeof(bufout2),5000) < 0)
		return -6;//QMessageBox::critical(this, tr("QT4 Demo"),tr("Bulk write failed"),QMessageBox::Ok, QMessageBox::Ok);
	retval = usb_bulk_read(dev,0x86,(char *)bufin,6,5000);
	if (retval < 0)
		return -7;//QMessageBox::critical(this, tr("QT4 Demo"),tr("Bulk read failed (start_measurement)"),QMessageBox::Ok, QMessageBox::Ok);
	else {
		if ((retval!= 6 ) || (bufin[4]!=0x86))
			return -8;//QMessageBox::critical(this, tr("QT4 Demo"),tr("Error in Start_measurement"),QMessageBox::Ok, QMessageBox::Ok);
	}
    printf("starting measurement\n");
    fflush(stdout);
	extern unsigned int l_Time;
	extern double l_pSpectrum[MAX_NR_PIXELS];
	int request_size;
	if (numavg <=1)
		request_size = sizeof(sony_single_meas);
	else
		request_size = sizeof(sony_multi_meas);
	int retval2;
	int measnr=0;
	while (measnr < l_NrOfScans) {
		//printf("starting %i scan",measnr);
		retval2 = usb_bulk_read(dev,0x86,(char *) &bufin,request_size,5000);
		if (retval2 < 0) {
			return -9;//QMessageBox::critical(this, tr("QT4 Demo"),tr("Bulk read failed (measurement)"),QMessageBox::Ok, QMessageBox::Ok);
		}
		else switch (bufin[4])
			{
				case 0xB0:
					if (retval2 == sizeof(sony_single_meas)) {
						memcpy(&sony_single_meas,bufin,sizeof(sony_single_meas));
						l_Time = swap32(sony_single_meas.timestamp);
						for (int teller=0;teller<m_NrPixels;++teller) 
							l_pSpectrum[teller]=swap16(sony_single_meas.pixels[teller]); 
					}
					break;
				case 0xB1:
					memcpy(&sony_multi_meas,bufin,sizeof(sony_multi_meas));
					l_Time = swap32(sony_multi_meas.timestamp);
					if (numavg != swap16(sony_multi_meas.averages)) printf("wrong nb. of averages\n");
						//QMessageBox::critical(this, tr("QT4 Demo"),tr("Error in Number of Averages"),QMessageBox::Ok, QMessageBox::Ok);
					for (int teller=0;teller<m_NrPixels;++teller) 
						l_pSpectrum[teller]=swap32(sony_multi_meas.pixels[teller])/numavg; 
					break;
				printf("passed %i scan",measnr);
			}
		++measnr;
		//printf("Measurement %i\n",measnr);
		if (measnr<l_NrOfScans){
			bufout3[00]=0x21;
			bufout3[01]=0x00;
			bufout3[02]=0x02;  // length of the command
			bufout3[03]=0x00;
			bufout3[04]=0xC0;  // acknowledge
			bufout3[05]=0x00;
    		if (usb_bulk_write(dev,0x02,(char *) bufout3,sizeof(bufout3),5000)< 0) printf("bad measurement\n");
			//QMessageBox::critical(this, tr("QT4 Demo"),tr("Writing acknowledgement to COM1 failed"),QMessageBox::Ok, QMessageBox::Ok);
		}
	} // while
    double spec_sum=0;
    for (int teller=0;teller<m_NrPixels;++teller) spec_sum+=l_pSpectrum[teller];
    spec_sum/=m_NrPixels;
    printf("..closing measurement (%.4f)\n",spec_sum);
    fflush(stdout);
	usb_release_interface(dev,0);
	//usb_close(dev);
	usb_reset(dev);
    return 0;
}


int save(char *name,int tlen=0) {
    FILE *file;
 	extern unsigned short m_NrPixels;
	extern double m_pLambda[MAX_NR_PIXELS];
	extern double l_pSpectrum[MAX_NR_PIXELS];
    if (name==NULL) file=stdout;
    else {
        file=fopen(name,"w");
        if (file==NULL) return -1;
    }
    fprintf(file,"# wavelength\tabsolute\n");
    if (tlen==0) tlen=m_NrPixels;
    for (int t=0;t<tlen;++t) fprintf(file,"%.5f\t%.5f\n",m_pLambda[t],l_pSpectrum[t]);
    fflush(file);
    if (name==NULL) fclose(file); 
    return 0; 
}
