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

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include<vector>
#include "packet/PacketLibDefinition.h"
#include "packet/PacketExceptionIO.h"
#include "io_hess.h"
#include "EventIO.hh"
//#include <map>

using namespace std;

/// Writing the Packet
int main(int argc, char *argv[])
{
    try
    {
        clock_t t;
        
        string ctarta;
        const char* home = getenv("CTARTA");
        if(argc > 2) {
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
        }

	//static AllHessData* hsdata;
	hess_all_data_struct*     hsdata;
	hsdata = NULL;
	uint16_t nev = 0, ntrg = 0;
	uint64_t item_type = 0;   /// ID of the last item
	uint32_t itel = 0;
	int32_t rc = 0;
	uint16_t verbose = 0, ignore = 0, quiet = 0;
	uint16_t showdata = 0;
	size_t events = 0, max_events = 0;
	uint32_t tel_id;
	uint16_t ntel_trg = 0, min_tel_trg = 0;
        int32_t _find_result; ///< Wether the file has more blocks to read

	// VF parameters
	uint16_t NTel;

	// Reading the input hessio file
	eventio::EventIO iobuf;  //eventio::EventIO _input;
	iobuf.OpenInput(argv[1]);
	hsdata = new AllHessData;


	// from Maier/Konrad
	/////////////////////////////////////////////
	// Loop over all data in the input data file
		for( ;; )
		{
			_find_result = iobuf.Find();
			item_type = iobuf.ItemType();
			cout << "item_type " << item_type << endl;
			iobuf.Read();
			
			if( _find_result != 0 )
			{
				break;
			}
			
			printf( "Found I/O block of type %ld\n", item_type );
	
			//////////////////////////////////////////
			// check header types
			switch( item_type )
			{
				/* =================================================== */
				case IO_TYPE_HESS_RUNHEADER:
					/* Summary of a preceding run in the same file ? */
					if( nev > 0 )
					{
						printf( "%d of %d events triggered.\n", ntrg, nev );
					}
					nev = ntrg = 0;
					//hsdata = ( AllHessData* ) calloc( 1, sizeof( AllHessData ) );
					if( ( rc = read_hess_runheader( iobuf.Buffer(), &hsdata->run_header ) ) < 0 )
					{
						cout << "Reading run header failed." << endl;
						exit( 1 );
					}
					if( showdata )
					{
						print_hess_runheader( iobuf.Buffer() );
					}
					
					NTel = hsdata->run_header.ntel;
					cout << "Number of telescopes: " << NTel << endl;
					for( itel = 0; itel < NTel; itel++ )
					{
						tel_id = hsdata->run_header.tel_id[itel]; // Telescope ID
						// Assign tel_id to array of telescopes configuration structures
						cout << "Telescope ID " << tel_id << " is telescope # " << itel << endl;
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

						cout << "DEBUG2" << endl;
					}
					break;
				case IO_TYPE_HESS_MCRUNHEADER:
				case IO_TYPE_MC_INPUTCFG:
				case IO_TYPE_HESS_CAMORGAN:
				case IO_TYPE_HESS_PIXELSET:
				case IO_TYPE_HESS_PIXELDISABLE:
				case IO_TYPE_HESS_CAMSOFTSET:
				case IO_TYPE_HESS_POINTINGCOR:
				case IO_TYPE_HESS_TRACKSET:
				case IO_TYPE_HESS_EVENT:
					rc = read_hess_event( iobuf.Buffer(), &hsdata->event, -1 );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_event(), rc = %d\n", rc );
					}
					events++;
					if( showdata )
					{
						print_hess_event( iobuf.Buffer() );
					}
					/* Count number of telescopes (still) present in data and triggered */
					ntel_trg = 0;
					for( itel = 0; itel < hsdata->run_header.ntel; itel++ )
						if( hsdata->event.teldata[itel].known )
						{
							/* If non-triggered telescopes record data (like HEGRA),
							   we may have to check the central trigger bit as well,
							   but ignore this for now. */
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

					
					break;
					
				/* =================================================== */
				case IO_TYPE_HESS_CALIBEVENT:
				{
					int type = -1;
					rc = read_hess_calib_event( iobuf.Buffer(), &hsdata->event, -1, &type );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_calib_event(), rc = %d, type=%d\n", rc, type );
					}
				}
				break;
				
				/* =================================================== */
				case IO_TYPE_HESS_MC_SHOWER:
					rc = read_hess_mc_shower( iobuf.Buffer(), &hsdata->mc_shower );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_mc_shower(), rc = %d\n", rc );
					}
					if( showdata )
					{
						print_hess_mc_shower( iobuf.Buffer() );
					}
					break;
					
				/* =================================================== */
				case IO_TYPE_HESS_MC_EVENT:
					rc = read_hess_mc_event( iobuf.Buffer(), &hsdata->mc_event );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_mc_event(), rc = %d\n", rc );
					}
					if( showdata )
					{
						print_hess_mc_event( iobuf.Buffer() );
					}

					
					break;
					
				/* =================================================== */
				case IO_TYPE_MC_TELARRAY:
					if( hsdata && hsdata->run_header.ntel > 0 )
					{
						rc = read_hess_mc_phot( iobuf.Buffer(), &hsdata->mc_event );
						if( verbose || rc != 0 )
						{
							printf( "read_hess_mc_phot(), rc = %d\n", rc );
						}
					}
					break;
					
				/* =================================================== */
				case IO_TYPE_HESS_MC_PE_SUM:
					rc = read_hess_mc_pe_sum( iobuf.Buffer(), &hsdata->mc_event.mc_pesum );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_mc_pe_sum(), rc = %d\n", rc );
					}
					if( showdata )
					{
						print_hess_mc_pe_sum( iobuf.Buffer() );
					}
					break;
					
				case IO_TYPE_HESS_TEL_MONI:
				case IO_TYPE_HESS_LASCAL:

					
				/* =================================================== */
				case IO_TYPE_HESS_RUNSTAT:
					rc = read_hess_run_stat( iobuf.Buffer(), &hsdata->run_stat );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_run_stat(), rc = %d\n", rc );
					}
					if( showdata )
					{
						print_hess_run_stat( iobuf.Buffer() );
					}
					break;
					
				/* =================================================== */
				case IO_TYPE_HESS_MC_RUNSTAT:
					rc = read_hess_mc_run_stat( iobuf.Buffer(), &hsdata->mc_run_stat );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_mc_run_stat(), rc = %d\n", rc );
					}
					if( showdata )
					{
						print_hess_mc_run_stat( iobuf.Buffer() );
					}
					break;
					
				default:
					if( !ignore )
					{
						fprintf( stderr, "Ignoring data block type %ld\n", item_type );
					}
			//}
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
