/***************************************************************************
 CreateConfig.cpp  -  description
 -------------------
 copyright            : (C) 2015 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software for non commercial purpose              *
 *   modify it under the terms of the GNU General Public License.          *
 *   For commercial purpose see appropriate license terms                  *
 *                                                                         *
 ***************************************************************************/

#include <algorithm>

#include "CreateConfig.h"
#include <rtautils/OutputFileFITS.h>

#define CONF_L0HEADER 0
#define CONF_L1HEADER 1


CreateConfig::CreateConfig(const string& confInputFileName) {
	
	try {
		
		if (confInputFileName != "") {

			/// if exists overwrite
			//checkExistsOverWrite(confInputFileName);
			/// create FITS file
			conf_file.create("!"+confInputFileName);            
		}
        
	}
	catch (qlbase::IOException& e) {
		cout << "ERROR: File "<< confInputFileName <<" does not exist. Error code: " << e.getErrorCode() << endl;
		exit(1);
	}
}

void CreateConfig::writeConfig_L0(int conf_Nrows_L0,
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
        		vector<float> vecMirrorArea) {

			/// Moving to the Telescope level Header L0
			conf_file.moveToHeader(CONF_L0HEADER);

			// create the field name vector
			std::vector<string> vec_conf_fields_L0;			
			vec_conf_fields_L0.push_back("L0ID");
			vec_conf_fields_L0.push_back("TelID");
			vec_conf_fields_L0.push_back("TelType");
			vec_conf_fields_L0.push_back("TelX");
			vec_conf_fields_L0.push_back("TelY");
			vec_conf_fields_L0.push_back("TelZ");
			vec_conf_fields_L0.push_back("FL");
			vec_conf_fields_L0.push_back("FOV");
			vec_conf_fields_L0.push_back("CameraScaleFactor");
			vec_conf_fields_L0.push_back("CameraCentreOffset");
			vec_conf_fields_L0.push_back("CameraRotation");
			vec_conf_fields_L0.push_back("NPixel");
			vec_conf_fields_L0.push_back("NPixel_active");
			vec_conf_fields_L0.push_back("NSamples");
			vec_conf_fields_L0.push_back("Sample_time_slice");
			vec_conf_fields_L0.push_back("NGains");
			vec_conf_fields_L0.push_back("HiLoScale");
			vec_conf_fields_L0.push_back("HiLoThreshold");
			vec_conf_fields_L0.push_back("HiLoOffset");
			vec_conf_fields_L0.push_back("NTubesOFF");
			vec_conf_fields_L0.push_back("NMirrors");
			vec_conf_fields_L0.push_back("MirrorArea");

			// writing a new binary table on first header should raise an error
			std::vector<qlbase::field> fields_L0(vec_conf_fields_L0.size());

			// create the fields structure
			for(unsigned int i=0; i<fields_L0.size(); i++)
			{
				fields_L0[i].name = vec_conf_fields_L0[i];
				fields_L0[i].vsize = 1;
			}
            
            		/// assign the type for each field

			fields_L0[0].type = qlbase::INT64;
			fields_L0[0].unit = "/";
			fields_L0[1].type = qlbase::INT16;
			fields_L0[1].unit = "/";
			fields_L0[2].type = qlbase::INT64;
			fields_L0[2].unit = "/";
			fields_L0[3].type = qlbase::FLOAT;
			fields_L0[3].unit = "m";
			fields_L0[4].type = qlbase::FLOAT;
			fields_L0[4].unit = "m";
			fields_L0[5].type = qlbase::FLOAT;
			fields_L0[5].unit = "m";
			fields_L0[6].type = qlbase::FLOAT;
			fields_L0[6].unit = "m";
			fields_L0[7].type = qlbase::FLOAT;
			fields_L0[7].unit = "";
			fields_L0[8].type = qlbase::FLOAT;
			fields_L0[8].unit = "/";
			fields_L0[9].type = qlbase::FLOAT;
			fields_L0[9].unit = "/";
			fields_L0[10].type = qlbase::FLOAT;
			fields_L0[10].unit = "deg";
			fields_L0[11].type = qlbase::INT16;
			fields_L0[11].unit = "/";
			fields_L0[12].type = qlbase::INT16;
			fields_L0[12].unit = "/";
			fields_L0[13].type = qlbase::INT16;
			fields_L0[13].unit = "/";
			fields_L0[14].type = qlbase::FLOAT;
			fields_L0[14].unit = "ns";
			fields_L0[15].type = qlbase::INT16;
			fields_L0[15].unit = "/";
			fields_L0[16].type = qlbase::FLOAT;
			fields_L0[16].unit = "/";
			fields_L0[17].type = qlbase::INT16;
			fields_L0[17].unit = "/";
			fields_L0[18].type = qlbase::FLOAT;
			fields_L0[18].unit = "/";
			fields_L0[19].type = qlbase::INT16;
			fields_L0[19].unit = "/";
			fields_L0[20].type = qlbase::INT16;
			fields_L0[20].unit = "/";
			fields_L0[21].type = qlbase::FLOAT;
			fields_L0[21].unit = "m^2";

			// Create binary table
			conf_file.createTable("TELESCOPE_LEVEL0", fields_L0);

			// Write the columns
			conf_file.write64i(0, vecL0ID, 0, conf_Nrows_L0-1);
			conf_file.write16i(1, vecTelID, 0, conf_Nrows_L0-1);
			conf_file.write64i(2, vecTelType, 0, conf_Nrows_L0-1);
			conf_file.write32f(3, vecTelX, 0, conf_Nrows_L0-1);
			conf_file.write32f(4, vecTelY, 0, conf_Nrows_L0-1);
			conf_file.write32f(5, vecTelZ, 0, conf_Nrows_L0-1);
			conf_file.write32f(6, vecFL, 0, conf_Nrows_L0-1);
			conf_file.write32f(7, vecFOV, 0, conf_Nrows_L0-1);
			conf_file.write32f(8, vecCameraScaleFactor, 0, conf_Nrows_L0-1);
			conf_file.write32f(9, vecCameraCentreOffset, 0, conf_Nrows_L0-1);
			conf_file.write32f(10, vecCameraRotation, 0, conf_Nrows_L0-1);
			conf_file.write16i(11, vecNPixel, 0, conf_Nrows_L0-1);
			conf_file.write16i(12, vecNPixel_active, 0, conf_Nrows_L0-1);
			conf_file.write16i(13, vecNSamples, 0, conf_Nrows_L0-1);
			conf_file.write32f(14, vecSample_time_slice, 0, conf_Nrows_L0-1);
			conf_file.write16i(15, vecNGains, 0, conf_Nrows_L0-1);
			conf_file.write32f(16, vecHiLoScale, 0, conf_Nrows_L0-1);
			conf_file.write16i(17, vecHiLoThreshold, 0, conf_Nrows_L0-1);
			conf_file.write32f(18, vecHiLoOffset, 0, conf_Nrows_L0-1);
			conf_file.write16i(19, vecNTubesOFF, 0, conf_Nrows_L0-1);
			conf_file.write16i(20, vecNMirrors, 0, conf_Nrows_L0-1);
			conf_file.write32f(21, vecMirrorArea, 0, conf_Nrows_L0-1);
            		

}

void CreateConfig::writeConfig_L1( int conf_Nrows_L1,
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
        		vector<int16_t> vecTubeOFF){

			/// Moving to the Telescope level Header L1
			conf_file.moveToHeader(CONF_L1HEADER);

			// create the field name vector
			std::vector<string> vec_conf_fields_L1;	
			vec_conf_fields_L1.push_back("L1ID");
			vec_conf_fields_L1.push_back("L0ID");
			vec_conf_fields_L1.push_back("PixelID");
			vec_conf_fields_L1.push_back("XTubeMM");
			vec_conf_fields_L1.push_back("YTubeMM");
			vec_conf_fields_L1.push_back("RTubeMM");
			vec_conf_fields_L1.push_back("XTubeDeg");
			vec_conf_fields_L1.push_back("YTubeDeg");
			vec_conf_fields_L1.push_back("RTubeDeg");
			vec_conf_fields_L1.push_back("AreaTube_m2");
			vec_conf_fields_L1.push_back("TubeOFF");

			// writing a new binary table on first header should raise an error
			std::vector<qlbase::field> fields_L1(vec_conf_fields_L1.size());

			// create the fields structure
			for(unsigned int i=0; i<fields_L1.size(); i++)
			{
				fields_L1[i].name = vec_conf_fields_L1[i];
				fields_L1[i].vsize = 1;
			}            

            		/// assign the type for each field

			fields_L1[0].type = qlbase::INT64;
			fields_L1[0].unit = "/";
			fields_L1[1].type = qlbase::INT64;
			fields_L1[1].unit = "/";
			fields_L1[2].type = qlbase::INT16;
			fields_L1[2].unit = "/";
			fields_L1[3].type = qlbase::FLOAT;
			fields_L1[3].unit = "mm";
			fields_L1[4].type = qlbase::FLOAT;
			fields_L1[4].unit = "mm";
			fields_L1[5].type = qlbase::FLOAT;
			fields_L1[5].unit = "mm";
			fields_L1[6].type = qlbase::FLOAT;
			fields_L1[6].unit = "deg";
			fields_L1[7].type = qlbase::FLOAT;
			fields_L1[7].unit = "deg";
			fields_L1[8].type = qlbase::FLOAT;
			fields_L1[8].unit = "deg";
			fields_L1[9].type = qlbase::FLOAT;
			fields_L1[9].unit = "m^2";
			fields_L1[10].type = qlbase::INT16;
			fields_L1[10].unit = "";

			// Create binary table
			conf_file.createTable("PIXEL_LEVEL1", fields_L1);

			// Write the columns
			conf_file.write64i(0, vecL1ID, 0, conf_Nrows_L1-1);
			conf_file.write64i(1, vecL0ID_L1, 0, conf_Nrows_L1-1);
			conf_file.write16i(2, vecPixelID, 0, conf_Nrows_L1-1);
			conf_file.write32f(3, vecXTubeMM, 0, conf_Nrows_L1-1);
			conf_file.write32f(4, vecYTubeMM, 0, conf_Nrows_L1-1);
			conf_file.write32f(5, vecRTubeMM, 0, conf_Nrows_L1-1);
			conf_file.write32f(6, vecXTubeDeg, 0, conf_Nrows_L1-1);
			conf_file.write32f(7, vecYTubeDeg, 0, conf_Nrows_L1-1);
			conf_file.write32f(8, vecRTubeDeg, 0, conf_Nrows_L1-1);
			conf_file.write32f(9, vecAreaTube_m2, 0, conf_Nrows_L1-1);
			conf_file.write16i(10, vecTubeOFF, 0, conf_Nrows_L1-1);

}

CreateConfig::~CreateConfig() {
	conf_file.writeKeyword("AUTHOR", "Valentina Fioretti", "INAF/IASF Bologna");
	conf_file.close();
    
}



