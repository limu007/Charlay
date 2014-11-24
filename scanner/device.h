#include <usb.h>

#ifndef DEV_H
#define DEV_H

#define LOBYTE(w)		(w & 0xff)
#define HIBYTE(w)		(w >> 8)

#define HANDLE_PRAGMA_PACK_PUSH_POP 1
#pragma pack(push,1)

typedef signed char     int8;
typedef unsigned char   uint8;
typedef signed short    int16;
typedef unsigned short  uint16;
typedef unsigned int    uint32;

uint8 const     USER_ID_LEN             = 64;
uint8 const     NR_WAVELEN_POL_COEF     = 5;
uint8 const     NR_NONLIN_POL_COEF      = 8;
uint8 const     NR_DEFECTIVE_PIXELS     = 30;
uint16 const    MAX_NR_PIXELS           = 4096;
uint8 const     NR_TEMP_POL_COEF        = 5;
uint8 const     MAX_TEMP_SENSORS        = 3;
uint8 const     ROOT_NAME_LEN           = 6;
uint8 const     AVS_SERIAL_LEN          = 10;
uint16 const    MAX_PIXEL_VALUE         = 0xFFFC;
uint8 const     MAX_VIDEO_CHANNELS      = 2;
uint16 const    MAX_LASER_WIDTH         = 0xFFFF;
uint8 const     HW_TRIGGER_MODE		    = 1;
uint8 const 	SW_TRIGGER_MODE	    	= 0;
uint8 const 	EDGE_TRIGGER_SOURCE  	= 0;
uint8 const 	LEVEL_TRIGGER_SOURCE	= 1;
uint8 const     MAX_TRIGGER_MODE        = 1;
uint8 const     MAX_TRIGGER_SOURCE      = 1;
uint8 const     MAX_TRIGGER_SOURCE_TYPE = 1;
uint32 const    MAX_INTEGRATION_TIME    = 600000;    // 600 seconds
uint8 const     SAT_DISABLE_DET         = 0;
uint8 const     SAT_ENABLE_DET          = 1;
uint8 const     SAT_PEAK_INVERSION      = 2;
uint8 const     NR_DAC_POL_COEF         = 2;
uint8 const     MAX_NR_SMOOTHING        = 51;

unsigned short const USB_VENDOR = 0x1992;   // Avantes
unsigned short const USB_PRODUCT = 0x0667;  // AS5216

float const wave_min=390;
float const wave_max=930;

int const smoothing=19;

typedef struct
{
    uint16                  m_StrobeControl;
    uint32                  m_LaserDelay;
    uint32                  m_LaserWidth;
    float                   m_LaserWaveLength;
    uint16                  m_StoreToRam;
} ControlSettingsType; 

typedef struct
{
    uint8                   m_Enable;
    uint8                   m_ForgetPercentage;
} DarkCorrectionType;

typedef uint8 SensorType;

typedef struct
{
    SensorType              m_SensorType;
    uint16                  m_NrPixels;
    float                   m_aFit[NR_WAVELEN_POL_COEF];
    bool                    m_NLEnable;
    double                  m_aNLCorrect[NR_NONLIN_POL_COEF];
    double                  m_aLowNLCounts;
    double                  m_aHighNLCounts;
    float                   m_Gain[MAX_VIDEO_CHANNELS];
    float                   m_Reserved;
    float                   m_Offset[MAX_VIDEO_CHANNELS];
    float                   m_ExtOffset;
    uint16                  m_DefectivePixels[NR_DEFECTIVE_PIXELS];
} DetectorType;

typedef struct
{
    uint16                  m_SmoothPix;
    uint8                   m_SmoothModel;
} SmoothingType;

typedef struct
{
    SmoothingType           m_Smoothing;
    float                   m_CalInttime;
    float                   m_aCalibConvers[MAX_NR_PIXELS];
} SpectrumCalibrationType;

typedef struct
{
    SpectrumCalibrationType m_IntensityCalib;
    uint8                   m_CalibrationType;
    uint32                  m_FiberDiameter;
} IrradianceType;

typedef struct
{
    uint8                   m_Mode;
    uint8                   m_Source;
    uint8                   m_SourceType;
} TriggerType;

typedef struct
{
    uint16                  m_StartPixel;
    uint16                  m_StopPixel;
    float                   m_IntegrationTime;
    uint32                  m_IntegrationDelay;
    uint32                  m_NrAverages;
    DarkCorrectionType      m_CorDynDark;
    SmoothingType           m_Smoothing;
    uint8                   m_SaturationDetection;
    TriggerType             m_Trigger;
    ControlSettingsType     m_Control;
} MeasConfigType;

// bugfix 09-09-08: MeasConfigType is nested in StandAloneType
// do not add prefix, define a new type

typedef struct
{
    uint8                   prefix[6];
    MeasConfigType          m_Meas;
} SendMeasConfigType;

typedef struct
{
    uint16                  m_Date;
    uint16                  m_Time;
} TimeStampType;

typedef struct
{
    bool                    m_Enable;
    uint8                   m_SpectrumType;
    char                    m_aFileRootName[ROOT_NAME_LEN];
    TimeStampType           m_TimeStamp;
} SDCardType;

typedef struct
{
    float                   m_aSpectrumCorrect[MAX_NR_PIXELS];
} SpectrumCorrectionType; 

typedef struct
{
    bool                    m_Enable;
    MeasConfigType          m_Meas;
    int16                   m_Nmsr;
    SDCardType              m_SDCard;
} StandAloneType;

typedef struct
{
    float                   m_aFit[NR_TEMP_POL_COEF];
} TempSensorType;

typedef struct
{
    bool                    m_Enable;
    float                   m_Setpoint;     // [degree Celsius]
    float                   m_aFit[NR_DAC_POL_COEF];
} TecControlType;

typedef struct
{
    float                   AnalogLow[2];
    float                   AnalogHigh[2];
    float                   DigitalLow[10];
    float                   DigitalHigh[10];
} ProcessControlType;

uint16 const    SETTINGS_RESERVED_LEN   = ((62*1024) -  sizeof(uint32) -
                                                        (sizeof(uint16) +   // m_Len
                                                         sizeof(uint16) +  // m_ConfigVersion
                                                         USER_ID_LEN +
                                                         sizeof(DetectorType) +
                                                         sizeof(IrradianceType) +
                                                         sizeof(SpectrumCalibrationType) +
                                                         sizeof(SpectrumCorrectionType) +
                                                         sizeof(StandAloneType) +
                                                        (sizeof(TempSensorType)*MAX_TEMP_SENSORS) +
                                                         sizeof(TecControlType) +
                                                         sizeof(ProcessControlType)
                                                        )
                                           );


typedef struct
{
    char                    prefix[6];
	uint16                  m_Len;
    uint16                  m_ConfigVersion;
    char                    m_aUserFriendlyId[USER_ID_LEN];
    DetectorType            m_Detector;
    IrradianceType          m_Irradiance;
    SpectrumCalibrationType m_Reflectance;
    SpectrumCorrectionType  m_SpectrumCorrect;
    StandAloneType          m_StandAlone;
    TempSensorType          m_aTemperature[MAX_TEMP_SENSORS];
    TecControlType          m_TecControl;
    ProcessControlType      m_ProcessControl;
    uint8                   m_aReserved[SETTINGS_RESERVED_LEN];
} DeviceConfigType;

typedef struct
{
    char            prefix[6];
	uint32          version1;
	uint32          version2;
	uint32          version3;
	char            SerialNumber[AVS_SERIAL_LEN];
    char            UserFriendlyName[USER_ID_LEN];
    unsigned char   Status;
} AvsIdentityType;

typedef struct
{
	char	prefix[6];
	uint32	timestamp;
	uint16	deadpix[18];
	uint16	pixels[2048];
} sony_single_measdatatype;

typedef struct
{
	char	prefix[6];
	uint32	timestamp;
	uint16	averages;
	uint32	deadpix[18];
	uint32	pixels[2048];
} sony_multi_measdatatype;


int init(usb_dev_handle *dev,float *fit);

int measure(usb_dev_handle *dev,double inttime,short l_NrOfScans,uint32 numavg,int smooth,int dark);

int save(char *name,int tlen);

#pragma pack(pop)

#endif
