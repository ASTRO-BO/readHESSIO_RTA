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
#include <map>

using namespace std;

map< unsigned int, unsigned int > fTelescopeSimTelList;


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

	static AllHessData* hsdata;
	int nev = 0, ntrg = 0;
	double wsum_all = 0.;
	uint64_t item_type;   ///< Id type of the latest read item
	int32_t itel, rc = 0;
	int verbose = 0, ignore = 0, quiet = 0;
	int showdata = 0, showhistory = 0;
	size_t events = 0, max_events = 0;
	int tel_id;
	int ntel_trg = 0, min_tel_trg = 0;

	IO_ITEM_HEADER item_header;
	eventio::EventIO iobuf;  //eventio::EventIO _input;
	iobuf.OpenInput(argv[1]);
	

	// from Maier/Konrad
	/////////////////////////////////////////////
	// Loop over all data in the input data file
		for( ;; )
		{
			cout << "Find()" << iobuf.Find() << endl;
			/* Find and read the next block of data. */
			/* In case of problems with the data, just give up. */
			if( iobuf.Find() != 0 )
			{
				break;
			}
			if( iobuf.Read() != 0 )
			{
				break;
			}
			
			item_type = iobuf.ItemType();

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
					wsum_all = 0.;
					/* Structures might be allocated from previous run */
					if( hsdata != NULL )
					{
						/* Free memory allocated inside ... */
						for( itel = 0; itel < hsdata->run_header.ntel; itel++ )
						{
							if( hsdata->event.teldata[itel].raw != NULL )
							{
								free( hsdata->event.teldata[itel].raw );
								hsdata->event.teldata[itel].raw = NULL;
							}
							if( hsdata->event.teldata[itel].pixtm != NULL )
							{
								free( hsdata->event.teldata[itel].pixtm );
								hsdata->event.teldata[itel].pixtm = NULL;
							}
							if( hsdata->event.teldata[itel].img != NULL )
							{
								free( hsdata->event.teldata[itel].img );
								hsdata->event.teldata[itel].img = NULL;
							}
						}
						/* Free main structure */
						free( hsdata );
						hsdata = NULL;
						
						/* Perhaps some cleaning needed in ROOT as well ... */
						
					}
					hsdata = ( AllHessData* ) calloc( 1, sizeof( AllHessData ) );
					if( ( rc = read_hess_runheader( iobuf.Buffer(), &hsdata->run_header ) ) < 0 )
					{
						cout << "Reading run header failed." << endl;
						exit( 1 );
					}
					if( !quiet )
					{
						printf( "Reading simulated data for %d telescope(s)\n", hsdata->run_header.ntel );
					}
					if( verbose || rc != 0 )
					{
						printf( "read_hess_runheader(), rc = %d\n", rc );
					}
					if( showdata )
					{
						print_hess_runheader( iobuf.Buffer() );
					}
					
					for( itel = 0; itel < hsdata->run_header.ntel; itel++ )
					{
						tel_id = hsdata->run_header.tel_id[itel];
						hsdata->camera_set[itel].tel_id = tel_id;
						fTelescopeSimTelList[tel_id] = itel;
						cout << "Telescope ID " << tel_id << " is telescope # " << itel << endl;
						hsdata->camera_org[itel].tel_id = tel_id;
						hsdata->pixel_set[itel].tel_id = tel_id;
						hsdata->pixel_disabled[itel].tel_id = tel_id;
						hsdata->cam_soft_set[itel].tel_id = tel_id;
						hsdata->tracking_set[itel].tel_id = tel_id;
						hsdata->point_cor[itel].tel_id = tel_id;
						hsdata->event.num_tel = hsdata->run_header.ntel;
						hsdata->event.teldata[itel].tel_id = tel_id;
						hsdata->event.trackdata[itel].tel_id = tel_id;
						if( ( hsdata->event.teldata[itel].raw =
									( AdcData* ) calloc( 1, sizeof( AdcData ) ) ) == NULL )
						{
							cout << "Not enough memory" << endl;
							exit( 1 );
						}
						hsdata->event.teldata[itel].raw->tel_id = tel_id;
						if( ( hsdata->event.teldata[itel].pixtm =
									( PixelTiming* ) calloc( 1, sizeof( PixelTiming ) ) ) == NULL )
						{
							cout << "Not enough memory" << endl;
							exit( 1 );
						}
						hsdata->event.teldata[itel].pixtm->tel_id = tel_id;
						if( ( hsdata->event.teldata[itel].img =
									( ImgData* ) calloc( 2, sizeof( ImgData ) ) ) == NULL )
						{
							cout << "Not enough memory" << endl;
							exit( 1 );
						}
						hsdata->event.teldata[itel].max_image_sets = 2;
						hsdata->event.teldata[itel].img[0].tel_id = tel_id;
						hsdata->event.teldata[itel].img[1].tel_id = tel_id;
						hsdata->tel_moni[itel].tel_id = tel_id;
						hsdata->tel_lascal[itel].tel_id = tel_id;
					}
					break;
				/* =================================================== */
				case IO_TYPE_HESS_MCRUNHEADER:
					rc = read_hess_mcrunheader( iobuf.Buffer(), &hsdata->mc_run_header );
					
					if( verbose || rc != 0 )
					{
						printf( "read_hess_mcrunheader(), rc = %d\n", rc );
					}
					//            if ( showdata )
					print_hess_mcrunheader( iobuf.Buffer() );
					break;
					
				/* =================================================== */
				case IO_TYPE_MC_INPUTCFG:
				{
					struct linked_string corsika_inputs;
					corsika_inputs.text = NULL;
					corsika_inputs.next = NULL;
					read_input_lines( iobuf.Buffer(), &corsika_inputs );
					if( corsika_inputs.text != NULL )
					{
						struct linked_string* xl = NULL, *xln = NULL;
						if( ! quiet )
						{
							printf( "\nCORSIKA was run with the following input lines:\n" );
						}
						for( xl = &corsika_inputs; xl != NULL; xl = xln )
						{
							if( ! quiet )
							{
								printf( "   %s\n", xl->text );
							}
							free( xl->text );
							xl->text = NULL;
							xln = xl->next;
							xl->next = NULL;
							if( xl != &corsika_inputs )
							{
								free( xl );
							}
						}
					}
				}
				break;
					
				/* =================================================== */
				case IO_TYPE_HESS_CAMSETTINGS:
					tel_id = item_header.ident; // Telescope ID is in the header
					if( ( itel = find_tel_idx( tel_id ) ) < 0 )
					{
						char msg[256];
						snprintf( msg, sizeof( msg ) - 1,
								  "Camera settings for unknown telescope %d.", tel_id );
						cout << msg << endl;
						exit( 1 );
					}
					rc = read_hess_camsettings( iobuf.Buffer(), &hsdata->camera_set[itel] );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_camsettings(), rc = %d\n", rc );
					}
					break;
					
				/* =================================================== */
				case IO_TYPE_HESS_CAMORGAN:
					tel_id = item_header.ident; // Telescope ID is in the header
					if( ( itel = find_tel_idx( tel_id ) ) < 0 )
					{
						char msg[256];
						snprintf( msg, sizeof( msg ) - 1,
								  "Camera organisation for unknown telescope %d.", tel_id );
						cout << msg << endl;
						exit( 1 );
					}
					rc = read_hess_camorgan( iobuf.Buffer(), &hsdata->camera_org[itel] );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_camorgan(), rc = %d\n", rc );
					}
					if( showdata )
					{
						print_hess_camorgan( iobuf.Buffer() );
					}
					break;
					
				/* =================================================== */
				case IO_TYPE_HESS_PIXELSET:
					tel_id = item_header.ident; // Telescope ID is in the header
					if( ( itel = find_tel_idx( tel_id ) ) < 0 )
					{
						char msg[256];
						snprintf( msg, sizeof( msg ) - 1,
								  "Pixel settings for unknown telescope %d.", tel_id );
						cout << msg << endl;
						exit( 1 );
					}
					rc = read_hess_pixelset( iobuf.Buffer(), &hsdata->pixel_set[itel] );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_pixelset(), rc = %d\n", rc );
					}
					if( showdata )
					{
						print_hess_pixelset( iobuf.Buffer() );
					}
					break;
					
				/* =================================================== */
				case IO_TYPE_HESS_PIXELDISABLE:
					tel_id = item_header.ident; // Telescope ID is in the header
					if( ( itel = find_tel_idx( tel_id ) ) < 0 )
					{
						char msg[256];
						snprintf( msg, sizeof( msg ) - 1,
								  "Pixel disable block for unknown telescope %d.", tel_id );
						cout << msg << endl;
						exit( 1 );
					}
					rc = read_hess_pixeldis( iobuf.Buffer(), &hsdata->pixel_disabled[itel] );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_pixeldis(), rc = %d\n", rc );
					}
					break;
					
				/* =================================================== */
				case IO_TYPE_HESS_CAMSOFTSET:
					tel_id = item_header.ident; // Telescope ID is in the header
					if( ( itel = find_tel_idx( tel_id ) ) < 0 )
					{
						char msg[256];
						snprintf( msg, sizeof( msg ) - 1,
								  "Camera software settings for unknown telescope %d.", tel_id );
						cout << msg << endl;
						exit( 1 );
					}
					rc = read_hess_camsoftset( iobuf.Buffer(), &hsdata->cam_soft_set[itel] );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_camsoftset(), rc = %d\n", rc );
					}
					break;
					
				/* =================================================== */
				case IO_TYPE_HESS_POINTINGCOR:
					tel_id = item_header.ident; // Telescope ID is in the header
					if( ( itel = find_tel_idx( tel_id ) ) < 0 )
					{
						char msg[256];
						snprintf( msg, sizeof( msg ) - 1,
								  "Pointing correction for unknown telescope %d.", tel_id );
						cout << msg << endl;
						exit( 1 );
					}
					rc = read_hess_pointingcor( iobuf.Buffer(), &hsdata->point_cor[itel] );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_pointingco(), rc = %d\n", rc );
					}
					break;
					
				/* =================================================== */
				case IO_TYPE_HESS_TRACKSET:
					tel_id = item_header.ident; // Telescope ID is in the header
					if( ( itel = find_tel_idx( tel_id ) ) < 0 )
					{
						char msg[256];
						snprintf( msg, sizeof( msg ) - 1,
								  "Tracking settings for unknown telescope %d.", tel_id );
						cout << msg << endl;
						exit( 1 );
					}
					rc = read_hess_trackset( iobuf.Buffer(), &hsdata->tracking_set[itel] );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_trackset(), rc = %d\n", rc );
					}
					break;
					
				/* =================================================== */
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
					
				/* =================================================== */
				case IO_TYPE_HESS_TEL_MONI:
					// Telescope ID among others in the header
					tel_id = ( item_header.ident & 0xff ) |
							 ( ( item_header.ident & 0x3f000000 ) >> 16 );
					if( ( itel = find_tel_idx( tel_id ) ) < 0 )
					{
						char msg[256];
						snprintf( msg, sizeof( msg ) - 1,
								  "Telescope monitor block for unknown telescope %d.", tel_id );
						cout << msg << endl;
						exit( 1 );
					}
					rc = read_hess_tel_monitor( iobuf.Buffer(), &hsdata->tel_moni[itel] );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_tel_monitor(), rc = %d\n", rc );
					}
					if( showdata )
					{
						print_hess_tel_monitor( iobuf.Buffer() );
					}
					break;
					
				/* =================================================== */
				case IO_TYPE_HESS_LASCAL:
					tel_id = item_header.ident; // Telescope ID is in the header
					if( ( itel = find_tel_idx( tel_id ) ) < 0 )
					{
						char msg[256];
						snprintf( msg, sizeof( msg ) - 1,
								  "Laser/LED calibration for unknown telescope %d.", tel_id );
						cout << msg << endl;
						exit( 1 );
					}
					rc = read_hess_laser_calib( iobuf.Buffer(), &hsdata->tel_lascal[itel] );
					if( verbose || rc != 0 )
					{
						printf( "read_hess_laser_calib(), rc = %d\n", rc );
					}
					if( showdata )
					{
						print_hess_laser_calib( iobuf.Buffer() );
					}
					break;
					
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
		

		reset_io_block( iobuf.Buffer() );
		
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
