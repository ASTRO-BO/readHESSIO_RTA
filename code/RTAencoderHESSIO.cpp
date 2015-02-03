/***************************************************************************
                          RTAencoderHESSIO.cpp  -  description
                             -------------------
    copyright            : (C) 2013 Andrea Bulgarelli
                               2013 Andrea Zoli
                               2014 Valentina Fioretti
    email                : bulgarelli@iasfbo.inaf.it
                           zoli@iasfbo.inaf.it
                           fioretti@iasfbo.inaf.it
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define MATH_PI 3.14159265359

#define CTA_PROD1 1
// This flag sets the presence of the array configuration file
#define ARRAYCONFIG 0

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <math.h>

#include <time.h>
#include <vector>

#include <packet/PacketLibDefinition.h>
#include <packet/PacketExceptionIO.h>
#include <packet/OutputPacketStream.h>
#include <packet/ByteStream.h>
#include <packet/Output.h>
#include <packet/OutputFile.h>

#include <io_hess.h>
#include <EventIO.hh>

#include "CreateConfig.h"

using namespace std;


int FillTelConfig( AllHessData* hsdata , const string& fits_filename){
	
	if( !hsdata )
	{
		return 0;
	}
	
	int NTel = hsdata->run_header.ntel;
	cout << "Number of telescopes: " << NTel << endl;
	

	// From G. Maier

	bool   fApplyCameraScaling = true;   // apply camera plate scaling according for DC telescopes

	// HARDCODED: all large telescopes are parabolic
	//            all telescopes with a mirror area larger than this value are
	//            considered to be parabolic (in m^2)
	double fParabolic_mirrorArea = 380.;
	// HARDCODED: all telescopes with 2 mirrors only are SC telescopes
	int    fSC_number_of_mirrors = 2;
	cout << "Info:" << endl;
	cout << " assume that all telescopes with mirror area larger than " << fParabolic_mirrorArea;
	cout << " m^2 are of parabolic type" << endl;
	cout << " assume that all telescopes with " << fSC_number_of_mirrors;
	cout << " mirrors are Schwarzschild-Couder telescopes" << endl;

	int16_t check_NTel = 0;
        vector<int64_t> vecL0ID;
        vector<int16_t> vecTelID;
	int64_t TelType;
        vector<int64_t> vecTelType;
        vector<float> vecTelX;
        vector<float> vecTelY;
        vector<float> vecTelZ;
	float FL;
        vector<float> vecFL;
        vector<float> vecFOV;
	float CameraScaleFactor;
        vector<float> vecCameraScaleFactor;
        vector<float> vecCameraCentreOffset;
        vector<float> vecCameraRotation;
	int16_t NPixel;
        vector<int16_t> vecNPixel;
        vector<int16_t> vecNPixel_active;
        vector<int16_t> vecNSamples;
        vector<float> vecSample_time_slice;
        vector<int16_t> vecNGains;
        vector<float> vecHiLoScale;
        vector<int16_t> vecHiLoThreshold;
        vector<float> vecHiLoOffset;
	int16_t nDisabled;
        vector<int16_t> vecNTubesOFF;
	int16_t NMirrors;
        vector<int16_t> vecNMirrors;
	float MirrorArea;
        vector<float> vecMirrorArea;
        
	int64_t counterL1IDrow = 0;
	int64_t L1ID = 0;
        vector<int64_t> vecL1ID;
        vector<int64_t> vecL0ID_L1;
        vector<int16_t> vecPixelID;
        vector<float> vecXTubeMM;
        vector<float> vecYTubeMM;
        vector<float> vecRTubeMM;
	float XTubeDeg;
        vector<float> vecXTubeDeg;
	float YTubeDeg;
        vector<float> vecYTubeDeg;
	float RTubeDeg;
        vector<float> vecRTubeDeg;
	vector<float> vecAreaTube_m2; 
	int16_t TubeOFF;
	vector<int16_t> vecTubeOFF;


	// fill the parameters
	for( int itel = 0; itel <  NTel; itel++ )
	{
		int16_t NPixel_active = 0;


		check_NTel++;
		if (check_NTel > NTel){
			return 0;
		}

		// Level 0 ID
		vecL0ID.push_back(itel);

		// Telescope ID
		vecTelID.push_back(hsdata->run_header.tel_id[itel]);

		// Telescope position
		// From Konrad:
 		// x,y,z positions of the telescopes [m]. x is counted from array reference 
		// position towards North, y towards West, z upwards.
		// ---------------------------------------------------------------------
		// From G. Maier:
		// Observe: transform to VERITAS coordinate system
		// x: East, y: North
		// VF: Unit still [m]			
		vecTelX.push_back(-1.*hsdata->run_header.tel_pos[itel][1]);
		vecTelY.push_back(hsdata->run_header.tel_pos[itel][0]);
		vecTelZ.push_back(hsdata->run_header.tel_pos[itel][2]);

		// Focal length [m]
		FL = hsdata->camera_set[itel].flen;
		vecFL.push_back(hsdata->camera_set[itel].flen);		

		// Camera conf (added by G. Maier)
		CameraScaleFactor = 1.;
		vecCameraCentreOffset.push_back(0.);

		// Rotation angle of camera (counter-clock-wise from back side for prime focus camera).
		// in [deg] as changed by G. Maier
		#ifdef CTA_PROD2		
			vecCameraRotation.push_back( -1.*hsdata->camera_set[itel].cam_rot * (180.0 / MATH_PI));
		#else
			vecCameraRotation.push_back( 0.);
		#endif

		// Number of pixels
		NPixel = hsdata->camera_set[itel].num_pixels;
		vecNPixel.push_back(hsdata->camera_set[itel].num_pixels);
		
		// Pixels disabled in HV
		nDisabled = hsdata->pixel_disabled[itel].num_HV_disabled;

		// Number of samples
		vecNSamples.push_back(hsdata->event.teldata[itel].raw->num_samples);

		//  Width of readout time slice (i.e. one sample) [ns]
		vecSample_time_slice.push_back(hsdata->pixel_set[itel].time_slice);

		// Number of gains
		vecNGains.push_back(hsdata->event.teldata[itel].raw->num_gains);

		// The scale factor (denominator) in shrinking h-g data
		vecHiLoScale.push_back(hsdata->event.teldata[itel].raw->scale_hg8);

		// Threshold (in high gain) for recording low-gain data
		vecHiLoThreshold.push_back(hsdata->event.teldata[itel].raw->threshold);

		// The offset to be used in shrinking high-gain data
		vecHiLoOffset.push_back(hsdata->event.teldata[itel].raw->offset_hg8);

		// Number of mirror tiles
		NMirrors = hsdata->camera_set[itel].num_mirrors;
		vecNMirrors.push_back(hsdata->camera_set[itel].num_mirrors);

		// Total area of individual mirrors corrected for inclination [m^2].
		MirrorArea = hsdata->camera_set[itel].mirror_area;
		vecMirrorArea.push_back(MirrorArea);

		float maxPix_dist = 0.;
		float pix_size = 0.;

		for( int p = 0; p < NPixel; p++ )
		{
			// Count the rows
			counterL1IDrow++;

			// Level 1 ID
			vecL1ID.push_back(L1ID++);

			// Level 0 in LEVEL 1 ID
			vecL0ID_L1.push_back(itel);

			// Pixel ID added by VF
			vecPixelID.push_back(p);

			// Pixel x,y position in camera [mm].			
			vecXTubeMM.push_back(hsdata->camera_set[itel].xpix[p] * 1.e3);
			vecYTubeMM.push_back(hsdata->camera_set[itel].ypix[p] * 1.e3);
			// use as size the radius of the active area of the tube
			vecRTubeMM.push_back(sqrt( hsdata->camera_set[itel].area[p] / MATH_PI ) * 1.e3);
			// Pixel active area ([m^2])
			vecAreaTube_m2.push_back(hsdata->camera_set[itel].area[p]);
			// camera_set.size = Pixel diameter (flat-to-flat, [m]).
			if( p == 0 )
			{
				pix_size = atan2( ( double )hsdata->camera_set[itel].size[p], ( double )FL ) * (180.0 / MATH_PI);
			}
					
			// Pixel x position in camera.
			// mm -> deg
			XTubeDeg = atan2( ( double )hsdata->camera_set[itel].xpix[p] / 1000., ( double )FL ) * (180.0 / MATH_PI);
			vecXTubeDeg.push_back(atan2( ( double )hsdata->camera_set[itel].xpix[p] / 1000., ( double )FL ) * (180.0 / MATH_PI));
			YTubeDeg = atan2( ( double )hsdata->camera_set[itel].ypix[p] / 1000., ( double )FL ) * (180.0 / MATH_PI);
			vecYTubeDeg.push_back(atan2( ( double )hsdata->camera_set[itel].ypix[p] / 1000., ( double )FL ) * (180.0 / MATH_PI));
			RTubeDeg = atan2( ( double )sqrt( hsdata->camera_set[itel].area[p] / MATH_PI ), ( double )FL ) * (180.0 / MATH_PI);
			vecRTubeDeg.push_back(atan2( ( double )sqrt( hsdata->camera_set[itel].area[p] / MATH_PI ), ( double )FL ) * (180.0 / MATH_PI));
			
					
			// FoV computed by G. Maier
			// Given the pixel position in degrees, the FoV is computed for that pixel distance. 
			// The maximum distance, i.e. the farthest pixel, gives the FoV
			float x2 = XTubeDeg * XTubeDeg;
			float y2 = YTubeDeg * YTubeDeg;
			float pix_dist = sqrt( x2 + y2 ) * 2.;
			float half_pix_dist = sqrt( x2 + y2 );

			if( pix_dist > maxPix_dist )
			{
				maxPix_dist = sqrt( x2 + y2 ) * 2.;
			}
			// disable pixels which are too far out
			// This check requires reading a configuration file for the array
			
			TubeOFF = hsdata->pixel_disabled[itel].HV_disabled[p];
			#ifdef ARRAYCONFIG
				// TODO reading the camera radius from the array configuration file
				float Camera_Radius = 0;		
				if( (half_pix_dist - RTubeDeg) > Camera_Radius )
				{
					TubeOFF = 2;
					nDisabled++;
				}			
			#endif

			// Pixels disabled in HV
			vecTubeOFF.push_back(TubeOFF);
			
			if( TubeOFF == 0 )
			{
				NPixel_active++;
			}
				
		}
		
		// FoV
		vecFOV.push_back(maxPix_dist);

		// telescope types
		TelType  = round( pix_size * 100. );
		TelType += round( maxPix_dist * 10. ) * 100;
		TelType += round( MirrorArea ) * 100 * 10 * 100;
		
		// Telescope type and Camera Scale Factor as defined by G. maier 
		// all large telescopes are parabolic, all others are Davies-Cotton (hardwired)
		if( MirrorArea > fParabolic_mirrorArea )
		{
			TelType += 100000000;
		}
		// Schwarzschild-Couder: check number of mirrors
		// (assumption is that SC telescope has 2 mirrors only)
		else if( NMirrors == fSC_number_of_mirrors )
		{
			TelType += 200000000;
		}
		else
		{
			// from raytracing by G.Hughes and optics note by S.Fegan and V.Vassiliev (2009)
			if( fApplyCameraScaling )
			{
				// assume that CTA MST is of intermediate-1.2 design
				if( abs( MirrorArea - 100. ) < 5. )
				{
					// G.Hughes (2011/10/21): 2.79+-0.06 %
					CameraScaleFactor = 1. - 0.0279;
				}
				// eq 1 from optics note by S.Fegan and V.Vassiliev (2009)
				else
				{
					// need f/D; unknown?
					//                 fCameraScaleFactor = 1.-1./(TMath::Power(2.,5.)*f_D*f_D) + 1./(TMath::Power(2.,8.) * TMath::Power(f_D,4.) );
					//
					CameraScaleFactor = 1.;
				}
			}
		}
		
		// Camera conf
		vecCameraScaleFactor.push_back(CameraScaleFactor);

		// Number of active pixels
		vecNPixel_active.push_back(NPixel_active);
		
		// Pixels disabled in HV
		vecNTubesOFF.push_back(nDisabled);

		// Telescope type
		vecTelType.push_back(TelType);		

	}

	cout << "Number of L1 rows" << counterL1IDrow << endl;


        CreateConfig telconfig( fits_filename );  
	telconfig.writeConfig_L0(       NTel,
					vecL0ID,
        				vecTelID,
        				vecTelType,
        				vecTelX,
        				vecTelY,
        				vecTelZ,
        				vecFL,
        				vecFOV,
        				vecCameraScaleFactor,
        				vecCameraCentreOffset,
        				vecCameraRotation,
        				vecNPixel,
        				vecNPixel_active,
        				vecNSamples,
        				vecSample_time_slice,
        				vecNGains,
        				vecHiLoScale,
        				vecHiLoThreshold,
        				vecHiLoOffset,
        				vecNTubesOFF,
        				vecNMirrors,
        				vecMirrorArea);
	
	telconfig.writeConfig_L1(       counterL1IDrow,
					vecL1ID,
        				vecL0ID_L1,
        				vecPixelID,
        				vecXTubeMM,
        				vecYTubeMM,
        				vecRTubeMM,
        				vecXTubeDeg,
        				vecYTubeDeg,
        				vecRTubeDeg,
        				vecAreaTube_m2,
        				vecTubeOFF);

	return 1;
}

/*
int FillPacketFADC(AllHessData* hsdata , const string& packet_filename){

	//Objects that manages a stream of byte
	OutputPacketStream* ops = 0;
	try
	{

		//create output packet stream - configFileName is the name of the xml file
    		ops = new OutputPacketStream("./conf/rta_fadc1.xml");

    		// Create and open an output device: file
    		Output* out = (Output*) new OutputFile(ops->isBigEndian());
    		vector<string> param;
    		param.push_back(argv[2]+"_FADC.raw");
    		out->openDevice(param);

   		// connect the output packet stream with the output
    		ops->setOutput(out)

		//get a packet to encode the data of a camera that manage 40 samples for each pixel
		Packet* p = ops->getPacketType("triggered_telescope1");

		// number of events		
		int numberOfEvent = 100;


		//get the number of samples managed by the packet
		word numPixelSamples = p->getPacketSourceDataField()->getFieldValue("Number of samples");

		//set the number of pixel managed by the packet
		word numberOfCameraPixels = 2048;
		p->getPacketSourceDataField()->setNumberOfBlocks(numberOfCameraPixels);

		//time of trigger
		unsigned long timetrigger = 100000;

		for(int i=0; i<nevents; i++) {

			//store some informations on the headers
		    	p->getPacketDataFieldHeader()->setFieldValue_32ui("Ttime:secs", timetrigger++);
		    	p->getPacketDataFieldHeader()->setFieldValue_32ui("Ttime:nsecs", rand());
		    	p->getPacketSourceDataField()->setFieldValue_32ui("eventNumber", i);

		    	//store the information of the pixels

 		    	//set each single field
		    	for(word pixel=0; pixel<numberOfCameraPixels; pixel++){
		        	for(word sample=0; sample<numPixelSamples; sample++) {
 		           		//word sampleValue = rand() % 100 + 50;
 		           		word sampleValue = sample;
 		           		p->getPacketSourceDataField()->getBlock(pixel)->setFieldValue(sample, sampleValue);
		        	}
			}

	  		//encode the packet
 		   	p->encode();

 		   	//compress the data
 		   	if(compress) p->compressData(LZ4, compress);

 		   	//write the encoded packet to output
 		   	ops->writePacket(p);

 		   	//get the size of the packet (only for measurement of performances)
 		   	dword packetSize = p->size();
 		   	totbytes += packetSize;

		}

        ///Number of events
        //int numberOfEvent = dst_tree->GetEntriesFast();
        int numberOfEvent = 100;
        cout << numberOfEvent << endl;
        /// Looping in the triggered events
        srand(0);
	long counts = 0;
	PacketLib::word SSC_index;
	vector<int16_t> SSC_array;
	SSC_array.resize(300);
		
	unsigned short ssc = 0;
        for(int evtindex = 0; evtindex<numberOfEvent; evtindex++) {
            //cout << "--" << evtindex << endl;

            //for each triggere telescope, generate a telemetry packet
            for(int telindex = 0; telindex<300; telindex++) {
				
				

            }		

        }
        
        t = clock() - t;
        //printf ("It took me %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        cout << "END " << counts << endl;
        
        return 0;

	}
	catch (PacketException* e)
	{
    		cout << "Error on stream constructor: ";
    		cout << e->geterror() << endl;
	}

	return 1;



}
*/

/// Writing the Packet
int main(int argc, char *argv[])
{
    try
    {
        clock_t t;
        
        string ctarta;
        const char* home = getenv("CTARTA");
        if(argc > 3) {
        	/// The Packet containing the FADC value of each triggered telescope
         	if (!home)
        	{
        	   std::cerr << "CTARTA environment variable is not defined." << std::endl;
        	   return 0;
        	}

        	ctarta = home;
        } else {            
            if(argc == 1){
        	     cerr << "Please, provide the names of the input and output files" << endl;
        	     return 0;
            }
            if(argc == 2){
        	     cerr << "Please, provide the name of the .raw file" << endl;
        	     return 0;
            }
            if(argc == 3){
        	     cerr << "Please, provide the name of the .fits file" << endl;
        	     return 0;
            }
        }

	AllHessData* hsdata;
	hsdata = NULL;
	PacketLib::word nev = 0, ntrg = 0;
	long int item_type = 0;   /// ID of the last item
	int itel = 0;
	int rc = 0;
	PacketLib::word verbose = 0, ignore = 1, quiet = 1;
	PacketLib::word showdata = 0;
	PacketLib::word events = 0, max_events = 0;
	PacketLib::dword tel_id;
	PacketLib::word ntel_trg = 0, min_tel_trg = 0;
        PacketLib::dword find_block, read_block; ///< Wether the file has more blocks to read

	// VF parameters
	PacketLib::word NTel;
	
 	// Reading the input hessio file
	eventio::EventIO iobuf;  //eventio::EventIO _input;
	iobuf.OpenInput(argv[1]);
	hsdata = new AllHessData;

	// ------------> preparing the packetlib output
	//Objects that manages a stream of byte
	PacketLib::OutputPacketStream* ops = 0;

	//create output packet stream - configFileName is the name of the xml file
    	ops = new PacketLib::OutputPacketStream("./conf/rta_fadc1.xml");

    	// Create and open an output device: file
    	PacketLib::Output* out = (PacketLib::Output*) new PacketLib::OutputFile(ops->isBigEndian());
    	vector<string> param;
	char name_fadc_out[100];   // array to hold the packet name.

	strcpy(name_fadc_out,argv[2]);
	strcat(name_fadc_out,"_FADC.raw"); 
    	param.push_back(name_fadc_out);
    	out->openDevice(param);

   	// connect the output packet stream with the output
    	ops->setOutput(out);

	//get a packet to encode the data of a camera that manage 40 samples for each pixel
	PacketLib::Packet* p = ops->getPacketType("triggered_telescope1");

        /// Looping in the triggered events
        srand(0);
	long counts = 0;
	unsigned short SSC_index;
	vector<int16_t> SSC_array;		
	PacketLib::word ssc = 0;

	//get a packet to encode the data of a camera that manage 40 samples for each pixel
	PacketLib::Packet* p = ops->getPacketType("triggered_telescope1");

	// number of events		
	//int numberOfEvent = 100;


	// from Gernot/Konrad/Etienne
	/////////////////////////////////////////////
	// Loop over all data in the input data file
		for( ;; )
		{
			if (!hsdata){
				cout << "Open the file before reading it" << endl;
				break;
			}

			find_block = iobuf.Find();

			if( find_block != 0 ){
				break;
			}

			item_type = iobuf.ItemType();
			cout << "item_type " << item_type << endl;
			read_block = iobuf.Read();

			if( read_block != 0 )
			{
				break;
			}
			
			int ret_code = 0;
			
			printf( "Found I/O block of type %ld\n", item_type );
	
			//////////////////////////////////////////
			// check header types
			switch( item_type )
			{
				/* =================================================== */
				case IO_TYPE_HESS_RUNHEADER:

					cout << "Activating IO_TYPE_HESS_RUNHEADER" << endl;			
					rc = read_hess_runheader( iobuf.Buffer(), &hsdata->run_header);
					cout << "rc: " << rc << endl;
					if( rc < 0 )
					{
						cout << "Reading run header failed." << endl;
						exit( 1 );
					}
            
					NTel = hsdata->run_header.ntel;
					for( itel = 0; itel < NTel; itel++ )
					{
						tel_id = hsdata->run_header.tel_id[itel]; // Telescope ID
						// Assign tel_id to array of telescopes configuration structures
						//cout << "Telescope ID " << tel_id << " is telescope # " << itel << endl;
						hsdata->camera_set[itel].tel_id = tel_id; // camera set
						hsdata->camera_org[itel].tel_id = tel_id; // camera_org
						hsdata->pixel_set[itel].tel_id = tel_id; // pixel_set
						hsdata->pixel_disabled[itel].tel_id = tel_id; // pixel_disabled
						hsdata->cam_soft_set[itel].tel_id = tel_id; //camera soft(ware?) set
						hsdata->tracking_set[itel].tel_id = tel_id; // tracking set
						hsdata->point_cor[itel].tel_id = tel_id; // pointing correction
						
						// event data set-up
						hsdata->event.num_tel = NTel;
						hsdata->event.teldata[itel].tel_id = tel_id;
						hsdata->event.trackdata[itel].tel_id = tel_id;

                    				hsdata->event.teldata[itel].raw           = new AdcData;
                    				hsdata->event.teldata[itel].raw->known    = 0;
                    				hsdata->event.teldata[itel].raw->tel_id   = tel_id;
                    				hsdata->event.teldata[itel].pixtm         = new PixelTiming;
                    				hsdata->event.teldata[itel].pixtm->tel_id = tel_id;
                    				hsdata->event.teldata[itel].img           = new ImgData[2];
                    				hsdata->event.teldata[itel].max_image_sets = 2;
                    				hsdata->event.teldata[itel].img[0].tel_id = tel_id;
                    				hsdata->event.teldata[itel].img[1].tel_id = tel_id;
                    				hsdata->tel_moni[itel].tel_id             = tel_id;
                    				hsdata->tel_lascal[itel].tel_id           = tel_id;
					}

					// From Etienne
					/*
                			_pointing_dir[0] = _hess_data->run_header.direction[0] + _hess_data->run_header.offset_fov[0];// + M_PI/2;
                			_pointing_dir[1] = _hess_data->run_header.direction[1] + _hess_data->run_header.offset_fov[1];
                			for (int32 i=0;i<_hess_data->run_header.ntel;i++)
                			{
                    				_tels_pos[_hess_data->run_header.tel_id[i]][0] = _hess_data->run_header.tel_pos[i][0];
                    				_tels_pos[_hess_data->run_header.tel_id[i]][1] = _hess_data->run_header.tel_pos[i][1];
                    				_tels_pos[_hess_data->run_header.tel_id[i]][2] = _hess_data->run_header.tel_pos[i][2];
                			}*/
					break;
				case IO_TYPE_HESS_CAMSETTINGS:
					tel_id = iobuf.ItemIdent(); // Telescope ID is in the header
					itel = find_tel_idx( tel_id );
					rc = read_hess_camsettings( iobuf.Buffer(), &hsdata->camera_set[itel] );

					break;
				case IO_TYPE_HESS_MCRUNHEADER:
					// MC run header info
					cout << "Activating IO_TYPE_HESS_MCRUNHEADER" << endl;	
					rc = read_hess_mcrunheader( iobuf.Buffer(), &hsdata->mc_run_header );
		
					break;
				case IO_TYPE_MC_INPUTCFG:
					// Corsika inputs info
					cout << "Activating IO_TYPE_MC_INPUTCFG" << endl;			
				case IO_TYPE_HESS_CAMORGAN:
					// Camera organization info
					cout << "Activating IO_TYPE_HESS_CAMORGAN" << endl;			
					tel_id = iobuf.ItemIdent(); // Telescope ID is in the header
					rc = read_hess_camorgan( iobuf.Buffer(), &hsdata->camera_org[itel] );

					break;
				case IO_TYPE_HESS_PIXELSET:
					// Pixel set info. The time slice is here.
					cout << "Activating IO_TYPE_HESS_PIXELSET" << endl;			
					tel_id = iobuf.ItemIdent(); // Telescope ID is in the header
					rc = read_hess_pixelset( iobuf.Buffer(), &hsdata->pixel_set[itel] );

					break;
				case IO_TYPE_HESS_PIXELDISABLE:
					// Pixel disable set info.
					cout << "Activating IO_TYPE_HESS_PIXELDISABLE" << endl;			
					tel_id = iobuf.ItemIdent(); // Telescope ID is in the header
					rc = read_hess_pixeldis( iobuf.Buffer(), &hsdata->pixel_disabled[itel] );

					break;
				case IO_TYPE_HESS_CAMSOFTSET:
					//Camera software settings.
					cout << "Activating IO_TYPE_HESS_CAMSOFTSET" << endl;			
					tel_id = iobuf.ItemIdent(); // Telescope ID is in the header
					rc = read_hess_camsoftset( iobuf.Buffer(), &hsdata->cam_soft_set[itel] );

					break;
				case IO_TYPE_HESS_POINTINGCOR:
					//Pointing correction info
					cout << "Activating IO_TYPE_HESS_POINTINGCOR" << endl;			
					tel_id = iobuf.ItemIdent(); // Telescope ID is in the header
					rc = read_hess_pointingcor( iobuf.Buffer(), &hsdata->point_cor[itel] );

					break;
				case IO_TYPE_HESS_TRACKSET:
					//Tracking settings
					cout << "Activating IO_TYPE_HESS_TRACKSET" << endl;			
					tel_id = iobuf.ItemIdent(); // Telescope ID is in the header
					rc = read_hess_trackset( iobuf.Buffer(), &hsdata->tracking_set[itel] );

					break;
				case IO_TYPE_HESS_EVENT:
					cout << "Activating IO_TYPE_HESS_EVENT" << endl;			
					rc = read_hess_event( iobuf.Buffer(), &hsdata->event, -1 );

					cout << "EventID: " << events << endl;
					cout << "Number of triggered telescopes: " << hsdata->event.num_teldata << endl;
					events++;

					//Count number of telescopes (still) present in data and triggered 
					ntel_trg = 0;
					for( itel = 0; itel < hsdata->run_header.ntel; itel++ )
						if( hsdata->event.teldata[itel].known )
						{
							// If non-triggered telescopes record data (like HEGRA),
							//   we may have to check the central trigger bit as well,
							//   but ignore this for now. 
							ntel_trg++;
						}
					if( hsdata->event.shower.known )
					{
						hsdata->event.shower.num_trg = ntel_trg;
					}
					if( ntel_trg < min_tel_trg )
					{
						continue;
					}
					ntrg++;

					// Writing the FADC packet
					SSC_array.resize(NTel);
		
				        //for each triggere telescope, generate a telemetry packet
            				for(int telindex = 0; telindex<ntel_trg; telindex++) {
				
				                //set the header of the tm packet
                				trtel->header.setAPID(vecTelID[telindex]); 	//the data generator (for now, the telescope)

				                for (int j = 0; j < vecTelID.size(); j++){
               					    if (TelescopeId == vectorTelID[j]){
                 				     	SSC_index = j;
                    					break;
                  				    }
          					}
                
                trtel->header.setSSC(ssc=	SSC_array[SSC_index]);	//a unique counter of packets
                cout << "ssc " << ssc << endl;



						//set the number of samples managed by the packet
						word numPixelSamples = p->getPacketSourceDataField()->getFieldValue("Number of samples");

						//set the number of pixel managed by the packet
						word numberOfCameraPixels = 2048;
						p->getPacketSourceDataField()->setNumberOfBlocks(numberOfCameraPixels);
	
						//time of trigger
						unsigned long timetrigger = 100000;
				

            				}		


			
					break;
					
				case IO_TYPE_HESS_CALIBEVENT:
					{
						cout << "Activating IO_TYPE_HESS_CALIBEVENT" << endl;			

						int type = -1;
						rc = read_hess_calib_event( iobuf.Buffer(), &hsdata->event, -1, &type );
					}
					break;
				
				case IO_TYPE_HESS_MC_SHOWER:
					cout << "Activating IO_TYPE_HESS_MC_SHOWER" << endl;			
					rc = read_hess_mc_shower( iobuf.Buffer(), &hsdata->mc_shower );
					break;
					
				case IO_TYPE_HESS_MC_EVENT:
					cout << "Activating IO_TYPE_HESS_MC_EVENT" << endl;			
					rc = read_hess_mc_event( iobuf.Buffer(), &hsdata->mc_event );
					
					break;
					
				case IO_TYPE_MC_TELARRAY:
					cout << "Activating IO_TYPE_MC_TELARRAY" << endl;			
					if( hsdata && hsdata->run_header.ntel > 0 )
					{
						rc = read_hess_mc_phot( iobuf.Buffer(), &hsdata->mc_event );
					}
					break;
					
				case IO_TYPE_HESS_MC_PE_SUM:
					cout << "Activating IO_TYPE_HESS_MC_PE_SUM" << endl;			
					rc = read_hess_mc_pe_sum( iobuf.Buffer(), &hsdata->mc_event.mc_pesum );
					break;
					
				case IO_TYPE_HESS_TEL_MONI:
					cout << "Activating IO_TYPE_HESS_TEL_MONI" << endl;			
					// Telescope ID among others in the header
					tel_id = ( iobuf.ItemIdent() & 0xff ) |
							 ( ( iobuf.ItemIdent() & 0x3f000000 ) >> 16 );

					rc = read_hess_tel_monitor( iobuf.Buffer(), &hsdata->tel_moni[itel] );

					break;

				case IO_TYPE_HESS_LASCAL:
					cout << "Activating IO_TYPE_HESS_LASCAL" << endl;			
					tel_id = iobuf.ItemIdent(); // Telescope ID is in the header
					rc = read_hess_laser_calib( iobuf.Buffer(), &hsdata->tel_lascal[itel] );

					break;					
				case IO_TYPE_HESS_RUNSTAT:
					cout << "Activating IO_TYPE_HESS_RUNSTAT" << endl;			
					rc = read_hess_run_stat( iobuf.Buffer(), &hsdata->run_stat );
					break;
					
				case IO_TYPE_HESS_MC_RUNSTAT:
					cout << "Activating IO_TYPE_HESS_MC_RUNSTAT" << endl;			
					rc = read_hess_mc_run_stat( iobuf.Buffer(), &hsdata->mc_run_stat );
					break;
					
				default:
					if( !ignore )
					{
						cout << "Ignoring data block type "<< item_type << endl;
					}

		}
		
		
		if( nev > 0 )
		{
			printf( "%d of %d events triggered\n", ntrg, nev );
		}
		
		if( hsdata != NULL )
		{
			hsdata->run_header.run = 0;
		}

	}
	// end while loop over all input files
	////////////////////////////////////////////////////

	iobuf.CloseInput();


	// Fill the telconfig parameters
	FillTelConfig(hsdata, argv[3]);

	// Fill the telconfig parameters
	//FillPacketFADC(hsdata, argv[2]);

    }
    catch(PacketLib::PacketExceptionIO* e)
    {
        cout << e->geterror() << endl;;
    }
    catch(PacketLib::PacketException* e)
    {
        cout << e->geterror() << endl;
    }
     	    
	return 1;
}


