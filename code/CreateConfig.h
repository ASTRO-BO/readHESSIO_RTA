/***************************************************************************
 CreateConfig.h  -  description
 -------------------
 copyright            : (C) 2015 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software for non commercial purpose              *
 *   and for public research institutes; you can redistribute it and/or    *
 *   modify it under the terms of the GNU General Public License.          *
 *   For commercial purpose see appropriate license terms                  *
 *                                                                         *
 ***************************************************************************/

#ifndef _CREATECONFIG_H
#define _CREATECONFIG_H

#include <rtautils/OutputFileFITS.h>
#include <iostream>
#include <stdlib.h>

using namespace std;

class CreateConfig {
        
    private:
        qlbase::OutputFileFITS conf_file;
        
        
    public:
        /// It writes the CTA configuration file in FITS format
        /// \param confInputFileName The input file name of the CTA configuration file
        CreateConfig(const string& confInputFileName);

        void writeConfig_L0( int conf_Nrows_L0,
			vector<int64_t> vecL0ID,
        		vector<int16_t> vecTelID,
        		vector<int64_t> vecTelType,
        		vector<float> vecTelX,
        		vector<float> vecTelY,
        		vector<float> vecTelZ,
        		vector<float> vecFL,
        		vector<float> vecFOV,
        		vector<float> vecCameraScaleFactor,
        		vector<float> vecCameraCentreOffset,
        		vector<float> vecCameraRotation,
        		vector<int16_t> vecNPixel,
        		vector<int16_t> vecNPixel_active,
        		vector<int16_t> vecNSamples,
        		vector<float> vecSample_time_slice,
        		vector<int16_t> vecNGains,
        		vector<float> vecHiLoScale,
        		vector<int16_t> vecHiLoThreshold,
        		vector<float> vecHiLoOffset,
        		vector<int16_t> vecNTubesOFF,
        		vector<int16_t> vecNMirrors,
        		vector<float> vecMirrorArea);
        
	void writeConfig_L1( int conf_Nrows_L1,
			vector<int64_t> vecL1ID,
        		vector<int64_t> vecL0ID_L1,
        		vector<int16_t> vecPixelID,
        		vector<float> vecXTubeMM,
        		vector<float> vecYTubeMM,
        		vector<float> vecRTubeMM,
        		vector<float> vecXTubeDeg,
        		vector<float> vecYTubeDeg,
        		vector<float> vecRTubeDeg,
        		vector<float> vecAreaTube_m2,
        		vector<int16_t> vecTubeOFF);
        
        ~CreateConfig();

        
        
        /// Number of telescope types
        int NTelType;

};

#endif
