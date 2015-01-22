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

#include <ctautils/OutputFileFITS.h>
#include <iostream>
#include <stdlib.h>


class CreateConfig {
        
    private:
        CTAUtils::OutputFileFITS conf_file;
        
        
    public:
        /// It writes the CTA configuration file in FITS format
        /// \param confInputFileName The input file name of the CTA configuration file
        CreateConfig(const std::string& confInputFileName);

        void writeConfig_L0( int conf_Nrows_L0,
			std::vector<int64_t> vecL0ID,
        		std::vector<int16_t> vecTelID,
        		std::vector<int64_t> vecTelType,
        		std::vector<float> vecTelX,
        		std::vector<float> vecTelY,
        		std::vector<float> vecTelZ,
        		std::vector<float> vecFL,
        		std::vector<float> vecFOV,
        		std::vector<float> vecCameraScaleFactor,
        		std::vector<float> vecCameraCentreOffset,
        		std::vector<float> vecCameraRotation,
        		std::vector<int16_t> vecNPixel,
        		std::vector<int16_t> vecNPixel_active,
        		std::vector<int16_t> vecNSamples,
        		std::vector<float> vecSample_time_slice,
        		std::vector<int16_t> vecNGains,
        		std::vector<float> vecHiLoScale,
        		std::vector<int16_t> vecHiLoThreshold,
        		std::vector<float> vecHiLoOffset,
        		std::vector<int16_t> vecNTubesOFF,
        		std::vector<int16_t> vecNMirrors,
        		std::vector<float> vecMirrorArea);
        
	void writeConfig_L1( int conf_Nrows_L1,
			std::vector<int64_t> vecL1ID,
        		std::vector<int64_t> vecL0ID_L1,
        		std::vector<int16_t> vecPixelID,
        		std::vector<float> vecXTubeMM,
        		std::vector<float> vecYTubeMM,
        		std::vector<float> vecRTubeMM,
        		std::vector<float> vecXTubeDeg,
        		std::vector<float> vecYTubeDeg,
        		std::vector<float> vecRTubeDeg,
        		std::vector<float> vecAreaTube_m2,
        		std::vector<int16_t> vecTubeOFF);
        
        ~CreateConfig();

        
        
        /// Number of telescope types
        int NTelType;

};

#endif
