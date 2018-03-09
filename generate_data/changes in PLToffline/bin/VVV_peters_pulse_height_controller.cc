#include <iostream>
#include <string>
#include <map>

#include "PLTEvent_peter.h"
#include "PLTU.h"

#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TFile.h"


// FUNCTION DEFINITIONS HERE
int PulseHeights (std::string const, std::string const);



float Average (std::vector<float>& V)
{
  double Sum = 0;
  for (std::vector<float>::iterator it = V.begin(); it != V.end(); ++it) {
    Sum += *it;
  }

  return Sum / (float) V.size();
}

std::string split_last_part(std::string s)
{
  std::vector<char> ret_string_v;

  for (int i = s.size()-5; i>-1; i--)
  {
    if (s[i]=='/')
    {
      break;
    }
    else
    {
      ret_string_v.push_back(s[i]);
    }
  }

  std::string ret_string;
  for (int i=ret_string_v.size(); i>-1; i--)
  {
    ret_string += ret_string_v[i];
  }
  

  return ret_string;

}

// CODE BELOW




int PulseHeights (std::string const DataFileName, std::string const GainCalFileName, std::string const  fill_numbers)
{
  PLTU::SetStyle();
  gStyle->SetOptStat(111111);

  // Grab the plt event reader
  PLTEvent Event(DataFileName, GainCalFileName);
  Event.SetPlaneClustering(PLTPlane::kClustering_Seed_5x5, PLTPlane::kFiducialRegion_All);
  Event.SetTrackingAlgorithm(PLTTracking::kTrackingAlgorithm_NoTracking);
  //  Event.SetPlaneFiducialRegion(PLTPlane::kFiducialRegion_m2_m2);

  // Map for all ROC hists and canvas
  std::map<int, std::vector<TGraphErrors*> > gClEnTimeMap;
  std::map<int, TH1F*>    hClusterSizeMap;
  std::map<int, TCanvas*> cClusterSizeMap;
  std::map<int, std::vector<TH1F*> > hMap;
  std::map<int, std::vector<TH1F*> > hMap_ADC; //p
  std::map<int, TCanvas*>            cMap;
  std::map<int, TH2F* >              hMap2D;
  std::map<int, TCanvas*>            cMap2D;

  //double Avg2D[250][PLTU::NCOL][PLTU::NROW];
  std::map<int, std::vector< std::vector<double> > > Avg2D;
  std::map<int, std::vector< std::vector<int> > > N2D;

  std::map<int, std::vector< std::vector<double> > > Avg2D_ADC; //p
  std::map<int, std::vector< std::vector<int> > >    N2D_ADC;   //p
  //int      N2D[250][PLTU::NCOL][PLTU::NROW];

  // Bins and max for pulse height plots
  int   const NBins =     60;
  float const XMin  =  -1000;
  float const XMax  =  50000;

  // Time width in events for energy time dep plots
  // This is the time width in ms
  // const unsigned int TimeWidth = 1000 * (60 * 1);
  const unsigned int TimeWidth = 1000;
  std::map<int, std::vector< std::vector<float> > > ChargeHits;
  std::map<int, std::vector< std::vector<float> > > ADCHits; //p
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over all events in file
  int NGraphPoints = 0;
  int ientry = 0;
  for ( ; Event.GetNextEvent() >= 0; ++ientry) {
    if (ientry % 10000 == 0) {
      std::cout << "Processing event: " << ientry << " at " << Event.ReadableTime() << std::endl;
    }



    // First event time
    //static uint32_t const StartTime = Event.Time();
    //uint32_t const ThisTime = Event.Time();
    static uint32_t const StartTime = 0;
    uint32_t static ThisTime = 0;
    ++ThisTime;

    if (ientry == 300000) {
      std::cout << "Reached target of 300000 events; stopping..." << std::endl;
      break;
    }

    while (ThisTime - (StartTime + NGraphPoints * TimeWidth) > TimeWidth) {
      // make point(s)
      for (std::map<int, std::vector<TGraphErrors*> >::iterator mit = gClEnTimeMap.begin(); mit != gClEnTimeMap.end(); ++mit) {
        int const id = mit->first;
        for (size_t ig = 0; ig != mit->second.size(); ++ig) {
          TGraphErrors* g = (mit->second)[ig];

          if (g->GetN() != NGraphPoints) {
            // Play some catchup
            g->Set(NGraphPoints);
            for (int i = 0; i > NGraphPoints; ++i) {
              g->SetPoint(i, i * TimeWidth, 0);
            }
          }

          g->Set( NGraphPoints + 1 );
          if (ChargeHits[id][ig].size() != 0) {
            float const Avg = PLTU::Average(ChargeHits[id][ig]);
            g->SetPoint(NGraphPoints, NGraphPoints * TimeWidth, Avg);
            g->SetPointError( NGraphPoints, 0, Avg/sqrt((float) ChargeHits[id][ig].size()));
            ChargeHits[id][ig].clear();
            ChargeHits[id][ig].reserve(10000);
          } else {
            g->SetPoint(NGraphPoints , NGraphPoints * TimeWidth, 0);
            g->SetPointError( NGraphPoints , 0, 0 );
          }
        }
      }
      ++NGraphPoints;

//      std::cout << NGraphPoints << std::endl;

    }



    for (size_t iTelescope = 0; iTelescope != Event.NTelescopes(); ++iTelescope) {
      PLTTelescope* Telescope = Event.Telescope(iTelescope);

      for (size_t iPlane = 0; iPlane != Telescope->NPlanes(); ++iPlane) {
        PLTPlane* Plane = Telescope->Plane(iPlane);

        int Channel = Plane->Channel();
        int ROC = Plane->ROC();


        if (ROC > 2) {
          std::cerr << "WARNING: ROC > 2 found: " << ROC << std::endl;
          continue;
        }
        if (Channel > 99) {
          std::cerr << "WARNING: Channel > 99 found: " << Channel << std::endl;
          continue;
        }

        // ID the plane and roc by 3 digit number
        int const id = 10 * Channel + ROC;

        if (!Avg2D.count(id)) {
          Avg2D[id].resize(PLTU::NCOL);
          N2D[id].resize(PLTU::NCOL);
          for (int icol = 0; icol != PLTU::NCOL; ++icol) {
            Avg2D[id][icol].resize(PLTU::NROW);
            N2D[id][icol].resize(PLTU::NROW);
          }
        }

        if (!hMap.count(id)) {
          hMap[id].push_back( new TH1F( TString::Format("Pulse Height for Ch %02i ROC %1i Pixels All %s", Channel, ROC, fill_numbers.c_str()),
                TString::Format("PulseHeight_Ch%02i_ROC%1i_All", Channel, ROC), NBins, XMin, XMax) );
            hMap2D[id] = new TH2F( TString::Format("Avg Charge Ch %02i ROC %1i Pixels All %s", Channel, ROC, fill_numbers.c_str()),
              TString::Format("PixelCharge_Ch%02i_ROC%1i_All", Channel, ROC), PLTU::NCOL, PLTU::FIRSTCOL, PLTU::LASTCOL, PLTU::NROW, PLTU::FIRSTROW, PLTU::LASTROW);
          for (size_t ih = 1; ih != 4; ++ih) {
            hMap[id].push_back( new TH1F( TString::Format("Pulse Height for Ch %02i ROC %1i Pixels %i %s", Channel, ROC, (int) ih, fill_numbers.c_str()),
                   TString::Format("PulseHeight_Ch%02i_ROC%1i_Pixels%i", Channel, ROC, (int) ih), NBins, XMin, XMax) );
          }

         

          // If we're making a new hist I'd say there's a 1 in 3 chance we'll need a canvas for it
          if (!cMap.count(Channel)) {
            // Create canvas with given name
            cMap[Channel] = new TCanvas( TString::Format("PulseHeight_Ch%02i", Channel), TString::Format("PulseHeight_Ch%02i", Channel), 900, 900);
            cMap[Channel]->Divide(3, 3);
          }
        }

        if (!hMap_ADC.count(id)) {
            hMap_ADC[id].push_back( new TH1F( TString::Format("Pulse Height ADC for Ch %02i ROC %1i Pixels All %s", Channel, ROC, fill_numbers.c_str()),
                TString::Format("PulseHeightADC_Ch%02i_ROC%1i_All", Channel, ROC), 100, 0-0.5, 1200-0.5) );
          
            hMap_ADC[id].push_back( new TH1F( TString::Format("Pulse Height ADC for Ch %02i ROC %1i Pixels %i %s", Channel, ROC, 1, fill_numbers.c_str()),
                    TString::Format("PulseHeight_Ch%02i_ROC%1i_Pixels%i", Channel, ROC,  1), 256, 0-0.5, 256-0.5) );

            hMap_ADC[id].push_back( new TH1F( TString::Format("Pulse Height ADC for Ch %02i ROC %1i Pixels %i %s", Channel, ROC, 2, fill_numbers.c_str()),
                    TString::Format("PulseHeight_Ch%02i_ROC%1i_Pixels%i", Channel, ROC, 2), 256, 0-0.5, 256-0.5) );

            hMap_ADC[id].push_back( new TH1F( TString::Format("Pulse Height ADC for Ch %02i ROC %1i Pixels %i %s", Channel, ROC, 3, fill_numbers.c_str()),
                    TString::Format("PulseHeight_Ch%02i_ROC%1i_Pixels%i", Channel, ROC, 3), 256, 0-0.5, 256-0.5) );
            hMap_ADC[id].push_back( new TH1F( TString::Format("Pulse Height ADC for Ch %02i ROC %1i AllPixels %s", Channel, ROC, fill_numbers.c_str()),
                    TString::Format("PulseHeight_Ch%02i_ROC%1i_AllPixels", Channel, ROC), 256, 0-0.5, 256-0.5) );
          
        }

        if (!gClEnTimeMap.count(id)) {
          gClEnTimeMap[id].resize(4);
          for (int ig = 0; ig != 4; ++ig) {
            TString const Name = TString::Format("TimeAvgGraph_id%d_Cl%d", id, ig);
            gClEnTimeMap[id][ig] = new TGraphErrors();
            gClEnTimeMap[id][ig]->SetName(Name);
          }
        }


        if (!ChargeHits.count(id)) {
          ChargeHits[id].resize(4);
          ChargeHits[id][0].reserve(10000);
          ChargeHits[id][1].reserve(10000);
          ChargeHits[id][2].reserve(10000);
          ChargeHits[id][3].reserve(10000);
        }

        // If this id doesn't exist in the cluster size map, make the hist and possibly canvas for this channel
        if (!hClusterSizeMap.count(id)) {
          hClusterSizeMap[id] = new TH1F( TString::Format("ClusterSize_Ch%02i_ROC%i %s", Channel, ROC, fill_numbers.c_str()), TString::Format("ClusterSize_Ch%02i_ROC%i", Channel, ROC), 10, 0, 10);
          hClusterSizeMap[id]->SetXTitle("Number of pixels in Cluster");

          // One in three chance you'll need a new canvas for thnat =)
          if (!cClusterSizeMap.count(Channel)) {
            cClusterSizeMap[Channel] = new TCanvas( TString::Format("ClusterSize_Ch%02i %s", Channel, fill_numbers.c_str()), TString::Format("ClusterSize_Ch%02i", Channel), 900, 300);
            cClusterSizeMap[Channel]->Divide(3, 1);
          }
        }


        // Loop over all clusters on this plane
        for (size_t iCluster = 0; iCluster != Plane->NClusters(); ++iCluster) {
          PLTCluster* Cluster = Plane->Cluster(iCluster);

          //if (Cluster->NHits() != 1) continue;
          //if (Cluster->Hit(0)->Column() != 31 || Cluster->Hit(0)->Row() != 55) continue;

          // Get number of hits in this cluster
          size_t NHits = Cluster->NHits();

          int const col = PLTGainCal::ColIndex(Cluster->SeedHit()->Column());
          int const row = PLTGainCal::RowIndex(Cluster->SeedHit()->Row());

          // Call it once.. it's faster.
          float const ThisClusterCharge = Cluster->Charge();
          std::vector<float> ThisClusterADC    = Cluster->ADC();

          if (ThisClusterCharge < 100000 && ThisClusterCharge >= 0) {
            Avg2D[id][col][row] = Avg2D[id][col][row] * ((double) N2D[id][col][row] / ((double) N2D[id][col][row] + 1.)) + ThisClusterCharge / ((double) N2D[id][col][row] + 1.);
            ++N2D[id][col][row];
          }


          // Fill cluster size
          hClusterSizeMap[id]->Fill(NHits);

          hMap[id][0]->Fill( ThisClusterCharge );
          if (NHits == 1) {
            hMap[id][1]->Fill( ThisClusterCharge );
          } else if (NHits == 2) {
            hMap[id][2]->Fill( ThisClusterCharge );
          } else if (NHits >= 3) {
            hMap[id][3]->Fill( ThisClusterCharge );
          }

          if (ThisClusterCharge < 200000) {
            ChargeHits[id][0].push_back( ThisClusterCharge );
            if (NHits == 1) {
              ChargeHits[id][1].push_back( ThisClusterCharge );
            } else if (NHits == 2) {
              ChargeHits[id][2].push_back( ThisClusterCharge );
            } else if (NHits >= 3) {
              ChargeHits[id][3].push_back( ThisClusterCharge );
            }
          }


          float sumADC = 0;

          for (uint i2=0; i2<ThisClusterADC.size(); i2++)
          {
              hMap_ADC[id][4]->Fill( ThisClusterADC[i2] );
              sumADC += ThisClusterADC[i2];
          }

          hMap_ADC[id][0]->Fill( sumADC );
          if (NHits == 1) {
              hMap_ADC[id][1]->Fill( ThisClusterADC[0] );
          } else if (NHits == 2) {
              hMap_ADC[id][2]->Fill( ThisClusterADC[0] );
              hMap_ADC[id][2]->Fill( ThisClusterADC[1] );
          } else if (NHits >= 3) {
              for (uint i2=0; i2<ThisClusterADC.size(); i2++)
              {
                  hMap_ADC[id][3]->Fill( ThisClusterADC[i2] );
              }
          }
          

          


        }
      }
    }




  }
  std::cout << "Events read: " << ientry+1 << std::endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TFile* output_root = TFile::Open(TString("plots/pulse_heights_output_" + fill_numbers + ".root"),"RECREATE");
  std::cout << "File created:    " << TString("plots/pulse_heights_output_" + fill_numbers + ".root") << std::endl;

	std::vector<TDirectory*>	plottingFolders;

  std::vector<TH1F*> h_summaPulse_calib;
  std::vector<TH1F*> h_summaPulse_ADC;

  std::vector<TH1F*> h_singlePulse_calib;
  std::vector<TH1F*> h_doublePulse_calib;
  std::vector<TH1F*> h_triplePulse_calib;
  std::vector<TH1F*> h_singlePulse_ADC;
  std::vector<TH1F*> h_doublePulse_ADC;
  std::vector<TH1F*> h_triplePulse_ADC;
  std::vector<TH1F*> h_allPixelPulse_ADC;
  std::vector<TH1F*> h_clusterSize;

  std::vector<TH2F*> h_spatialChargeDistribution;


  //std::map<int, std::vector<TH1F*> > hMap_ADC = hMap;

  // Loop over all histograms and draw them in the root file
  for (std::map<int, std::vector<TH1F*> >::iterator it = hMap.begin(); it != hMap.end(); ++it) {

    // Decode the ID
    int const Channel = it->first / 10;
    int const ROC     = it->first % 10;
    int const id      = it->first;

    printf("Drawing hists for Channel %2i ROC %i\n", Channel, ROC);

    //form output folder, insert relevant histograms
    output_root->cd();
    std::string temp_name = "Ch" + std::to_string(Channel) + "_ROC" + std::to_string(ROC);
    plottingFolders.push_back(output_root->mkdir(temp_name.c_str()));

    
    plottingFolders[plottingFolders.size()-1]->cd();

   
    //fill calibrated pulse height histograms
    TH1F* h_tempTH1F;

    h_tempTH1F = (TH1F*)hMap[id][0]->Clone();
    h_tempTH1F->SetName("h_summaPulse_calib");
    h_summaPulse_calib.push_back(h_tempTH1F);

    h_tempTH1F = (TH1F*)hMap[id][1]->Clone();
    h_tempTH1F->SetName("h_singlePulse_calib");
    h_singlePulse_calib.push_back(h_tempTH1F);

    h_tempTH1F = (TH1F*)hMap[id][2]->Clone();
    h_tempTH1F->SetName("h_doublePulse_calib");
    h_doublePulse_calib.push_back(h_tempTH1F);

    h_tempTH1F = (TH1F*)hMap[id][3]->Clone();
    h_tempTH1F->SetName("h_triplePulse_calib");
    h_triplePulse_calib.push_back(h_tempTH1F);

    //fill ADC height spectra
    h_tempTH1F = (TH1F*)hMap_ADC[id][0]->Clone();
    h_tempTH1F->SetName("h_clusterSumPulse_ADC");
    h_summaPulse_ADC.push_back(h_tempTH1F);

    h_tempTH1F = (TH1F*)hMap_ADC[id][1]->Clone();
    h_tempTH1F->SetName("h_singlePulse_ADC");
    h_singlePulse_ADC.push_back(h_tempTH1F);

    h_tempTH1F = (TH1F*)hMap_ADC[id][2]->Clone();
    h_tempTH1F->SetName("h_doublePulse_ADC");
    h_doublePulse_ADC.push_back(h_tempTH1F);

    h_tempTH1F = (TH1F*)hMap_ADC[id][3]->Clone();
    h_tempTH1F->SetName("h_triplePulse_ADC");
    h_triplePulse_ADC.push_back(h_tempTH1F);

    h_tempTH1F = (TH1F*)hMap_ADC[id][4]->Clone();
    h_tempTH1F->SetName("h_allPixelPulse_ADC");
    h_allPixelPulse_ADC.push_back(h_tempTH1F);

    //fill cluster size histo
    h_tempTH1F = (TH1F*)hClusterSizeMap[id]->Clone();
    h_tempTH1F->SetName("h_clusterSize");
    h_clusterSize.push_back(h_tempTH1F);

    // //fill 2D histo
    // TH2F* h_tempTH2F = (TH2F*)hMap2D[id]->Clone();
    // h_tempTH2F->SetName("h_spatialChargeDistribution");
    // h_spatialChargeDistribution.push_back(h_tempTH2F);

  }


  gROOT->SetBatch(kTRUE);
	output_root->Write();
	gROOT->GetListOfFiles()->Remove(output_root);

  output_root->Close();


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  
  return 0;
}
  


int main (int argc, char* argv[])
{


  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " [DataFileName] [GainCalFileName] [fill_ID]" << std::endl;
    return 1;
  }

  std::string const DataFileName = argv[1];
  std::string const GainCalFileName = argv[2];
  std::string const fill_ID = argv[3];
  std::cout << "DataFileName:    " << DataFileName << std::endl;
  std::cout << "GainCalFileName: " << GainCalFileName << std::endl;
  std::cout << "fill_ID:         " << fill_ID << std::endl;

  PulseHeights(DataFileName, GainCalFileName, fill_ID);


    

   return 0;
}
