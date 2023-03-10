/**
 *  @copyright Copyright 2020 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file EventCategorizer.cpp
 */

#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>
#include "EventCategorizerTools.h"
#include "EventCategorizer.h"
#include <iostream>
#include <tuple>
#include "TEfficiency.h"


using namespace jpet_options_tools;
using namespace std;

EventCategorizer::EventCategorizer(const char* name): JPetUserTask(name) {}

EventCategorizer::~EventCategorizer() {}

bool EventCategorizer::init()
{
  INFO("Event categorization started.");
  // Parameter for back to back categorization
  if (isOptionSet(fParams.getOptions(), kBack2BackSlotThetaDiffParamKey)){
    fB2BSlotThetaDiff = getOptionAsFloat(fParams.getOptions(), kBack2BackSlotThetaDiffParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kBack2BackSlotThetaDiffParamKey.c_str(), fB2BSlotThetaDiff
    ));
  }
  // Parameter for scattering determination
  if (isOptionSet(fParams.getOptions(), kScatterTOFTimeDiffParamKey)) {
    fScatterTOFTimeDiff = getOptionAsFloat(fParams.getOptions(), kScatterTOFTimeDiffParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kScatterTOFTimeDiffParamKey.c_str(), fScatterTOFTimeDiff
    ));
  }
  // Parameters for deexcitation TOT cut
  if (isOptionSet(fParams.getOptions(), kDeexTOTCutMinParamKey)) {
    fDeexTOTCutMin = getOptionAsFloat(fParams.getOptions(), kDeexTOTCutMinParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kDeexTOTCutMinParamKey.c_str(), fDeexTOTCutMin
    ));
  }
  if (isOptionSet(fParams.getOptions(), kDeexTOTCutMaxParamKey)) {
    fDeexTOTCutMax = getOptionAsFloat(fParams.getOptions(), kDeexTOTCutMaxParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kDeexTOTCutMaxParamKey.c_str(), fDeexTOTCutMax
    ));
  }
  if (isOptionSet(fParams.getOptions(), kMaxTimeDiffParamKey)) {
    fMaxTimeDiff = getOptionAsFloat(fParams.getOptions(), kMaxTimeDiffParamKey);
  } else {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kMaxTimeDiffParamKey.c_str(), fMaxTimeDiff));
  }
  // Getting bool for saving histograms
  if (isOptionSet(fParams.getOptions(), kSaveControlHistosParamKey)) {
    fSaveControlHistos = getOptionAsBool(fParams.getOptions(), kSaveControlHistosParamKey);
  }
  if (isOptionSet(fParams.getOptions(), kTOTCalculationType)) {
    fTOTCalculationType = getOptionAsString(fParams.getOptions(), kTOTCalculationType);
  } else {
    WARNING("No TOT calculation option given by the user. Using standard sum.");
  }

  Event_Eff = new TEfficiency("event_eff","Event_efficiency;methods;#epsilon",20,0,10);
  //  Hit_Eff = new TEfficiency("hit_eff","cut_efficiency;cuts;#epsilon",20,0,10);


  // Input events type
  fOutputEvents = new JPetTimeWindow("JPetEvent");
  // Initialise hisotgrams
  if(fSaveControlHistos){
    initialiseHistograms();
  }
  return true;
}

bool EventCategorizer::exec()
{
  if (auto& timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {
    vector<JPetEvent> events;
       
    for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++) {
      const auto& event = dynamic_cast<const JPetEvent&>(timeWindow->operator[](i));
      TotalInitialcut.totalNumber++;

      /*      std::pair<TEfficiency*, bool> isInit = EventCategorizerTools::initialCut(
       event, getStatistics(), fSaveControlHistos);
      Hit_Eff = std::get<0>(isInit);
      bool isInitialCut = std::get<1>(isInit);*/
      bool isInitialCut = EventCategorizerTools::initialCut(              
       event, getStatistics(), fSaveControlHistos, fHitCounter);
      Event_Eff->Fill(isInitialCut,1);
      bool isAnnihilation, is2Gamma, is3Gamma, isPrompt,isNeighbourHits;
      bool isScattered = EventCategorizerTools::checkForScatter(
        event, getStatistics(), fSaveControlHistos, fScatterTOFTimeDiff, fTOTCalculationType);
      Event_Eff->Fill(isScattered,7);
      
      double sum_tot=0.0;
      double sum_tot_2g= 0.0;
      double sum_tot_3g=0.0;
      double sum_tot_scatter=0.0;
      double sum_tot_prompt=0.0;
      double sum_tot_ann=0.0;
      double sum_tot_ann_prompt= 0.0;
      double sum_tot_3gann_prompt=0.0;

      

      if (isInitialCut){TotalInitialcut.totalAccepted++;
	fAnnihilation.totalNumber++;
	
	
	isAnnihilation = EventCategorizerTools::checkForAnnihilation(
       event, getStatistics(), fSaveControlHistos);
	
	isPrompt = EventCategorizerTools::checkForPrompt(
        event, getStatistics(), fSaveControlHistos, fDeexTOTCutMin, fDeexTOTCutMax, fTOTCalculationType);

	Event_Eff->Fill(isPrompt,4);

       if(isAnnihilation){
	 fAnnihilation.totalAccepted++;
	 is2Gamma = EventCategorizerTools::checkFor2Gamma(
                    event, getStatistics(), fSaveControlHistos, fB2BSlotThetaDiff, fMaxTimeDiff
                     );
         is3Gamma = EventCategorizerTools::checkFor3Gamma(
                    event, getStatistics(), fSaveControlHistos
                    );
	 isNeighbourHits = EventCategorizerTools::removeNeighbourhits(
        event, getStatistics(),fSaveControlHistos, fTOTCalculationType
      );

	   Event_Eff->Fill(isAnnihilation,2);
	   Event_Eff->Fill(is2Gamma,5);
	   Event_Eff->Fill(is3Gamma,6);
	   Event_Eff->Fill(isNeighbourHits,3);
       }
      }
      
       if(isNeighbourHits);
       
      JPetEvent newEvent = event;
      if(is2Gamma){  newEvent.addEventType(JPetEventType::k2Gamma); }
      if(is3Gamma){  newEvent.addEventType(JPetEventType::k3Gamma); }
      if(isPrompt) newEvent.addEventType(JPetEventType::kPrompt);
      if(isScattered) newEvent.addEventType(JPetEventType::kScattered);
      
      
      if(fSaveControlHistos){
       for(auto hit : event.getHits()){
	  double tot = hit.getEnergy();
	  double hit_pos = hit.getPos().Mag();
	  double hit_time = hit.getTime();
	  double hit_timediff = fabs(hit_time - hit_pos/ kLightVelocity_cm_ps);
	  int scID = hit.getScintillator().getID();
	  double Z = hit.getPosZ();

	  sum_tot += tot;
	  if(is2Gamma) sum_tot_2g += tot;
	  if(is3Gamma) sum_tot_3g += tot;	  
	  if(is2Gamma || is3Gamma) sum_tot_ann += tot;      
	  if(isPrompt) sum_tot_prompt += tot;
	  if(isScattered) sum_tot_scatter += tot;
	  if(!isScattered) sum_tot_ann_prompt += tot;
	  if(is3Gamma || isPrompt) sum_tot_3gann_prompt +=tot;

	     //Filling of sum TOT
            getStatistics().fillHistogram("sum_TOT_2g", sum_tot_2g);
            getStatistics().fillHistogram("sum_TOT_3g", sum_tot_3g);
            getStatistics().fillHistogram("sum_TOT_ann", sum_tot_ann);
	    getStatistics().fillHistogram("sum_TOT_prompt", sum_tot_prompt);
	    getStatistics().fillHistogram("sum_TOT_scatter", sum_tot_scatter);
	    getStatistics().fillHistogram("sum_TOT_ann_prompt", sum_tot_ann_prompt);
	    getStatistics().fillHistogram("sum_TOT_3gann_prompt", sum_tot_3gann_prompt);
            getStatistics().fillHistogram("sum_TOT", sum_tot);

	    ////////
            getStatistics().fillHistogram("All_XYpos", hit.getPosX(), hit.getPosY());
            getStatistics().fillHistogram("Hit_timedifference", hit_timediff);
	    getStatistics().fillHistogram("Z Vs scID", Z, scID);	    
        }
      }
      events.push_back(newEvent);
      
    }
   saveEvents(events);
  }
  else { return false; }
  return true;

}

bool EventCategorizer::terminate()
{
  INFO("Event categorization completed.");
  INFO("Total hits:" + std::to_string(fHitCounter.totalNumber));
  INFO("Total accepted hits:" + std::to_string(fHitCounter.totalAccepted));
  INFO("Ratio of accepted hits: " + std::to_string(fHitCounter.getRatio()));

  /* // std::cout <<float(totalPrompts)/totalEvents <<std::endl;
  // std::cout <<float(totalScattered)/totalEvents <<std::endl;
  auto file1 = TFile::Open("efficiency_hit.root", "recreate");
  // if (file) std::cout << "file was created"  <<std::endl;
  //Event_Eff->SetDirectory(gDirectory);
  Hit_Eff->SetDirectory(gDirectory);
  file1->Write();
  file1->Close();
   Hit_Eff= nullptr;
   delete Hit_Eff;
  auto file2 = TFile::Open("efficiency_event.root", "recreate");
  Event_Eff->SetDirectory(gDirectory);
  file2->Write();
  file2->Close();
  Event_Eff = nullptr;
  delete Event_Eff;*/
  return true;
  
}

void EventCategorizer::saveEvents(const vector<JPetEvent>& events)
{
  for (const auto& event : events) { fOutputEvents->add<JPetEvent>(event); }
}

void EventCategorizer::initialiseHistograms(){

  // General histograms
  getStatistics().createHistogramWithAxes(
                                          new TH2D("All_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
                                          "Hit X position [cm]", "Hit Y position [cm]"
                                           );
  
  getStatistics().createHistogramWithAxes(
                                          new TH2D("Z Vs scID", "Z Vs scID", 100, -50, 50, 200, 0, 200),
                                          "Scintillator ID", "Hit Z position [cm]"
                                           );

  // Histograms for sum TOT
  
  getStatistics().createHistogramWithAxes(
  					  new TH1D("sum_TOT_2g","Sum of all 2gamma TOT", 200, 0, 100000),
  					  "SUM_TOT_2g(ps)","counts");
  getStatistics().createHistogramWithAxes(
                                          new TH1D("sum_TOT_3g","Sum of all 3gamma TOT", 200, 0, 100000),
                                          "SUM_TOT_3g(ps)","counts");
  getStatistics().createHistogramWithAxes(
                                          new TH1D("sum_TOT_ann","Sum of all annhilations TOT(2Gamma and 3Gamma)", 400, 0, 200000),
                                          "SUM_TOT_annhilation(ps)","counts");
  getStatistics().createHistogramWithAxes(
                                          new TH1D("sum_TOT_prompt","Sum of all ptompt TOT", 200, 0, 200000),
                                          "SUM_TOT_prompt(ps)","counts");
  getStatistics().createHistogramWithAxes(
                                          new TH1D("sum_TOT","Sum of all TOT", 400, 0, 200000),
                                          "SUM_TOT(ps)","counts");
  getStatistics().createHistogramWithAxes(
                                          new TH1D("sum_TOT_scatter","Sum of all scattered TOT", 400, 0, 200000),
                                          "SUM_TOT_scatter(ps)","counts");
  getStatistics().createHistogramWithAxes(
                                          new TH1D("sum_TOT_ann_prompt","Sum of all annhilation and prompt TOT", 400, 0, 200000),
                                          "SUM_TOT_ann_prompt(ps)","counts");
  getStatistics().createHistogramWithAxes(
                                          new TH1D("sum_TOT_3gann_prompt","Sum of all 3gamma and prompt TOT", 400, 0, 200000),
                                          "SUM_TOT_3g_prompt(ps)","counts");

  //For mew category checkForAnnihilation

 getStatistics().createHistogramWithAxes(
    new TH1D("Ann_TOT_before_cut", "TOT of all hits before cut", 200, -100.0, 200000.0),
    "TOT [ps]", "Number of Hits"
  );
 getStatistics().createHistogramWithAxes(
    new TH1D("simple_tot", "TOT of all hits after cut", 200, -100.0, 200000.0),
    "TOT [ps]", "Number of Hits"
   );


  
  getStatistics().createHistogramWithAxes(
    new TH1D("Ann_TOT", "TOT of all hits after cut", 200, -100.0, 200000.0),
    "TOT [ps]", "Number of Hits"
   );

getStatistics().createHistogramWithAxes(
                                          new TH1D("Hit_multiplicity_0","Multiplicity of hits(no cut)", 10, 0.5, 10.5),
                                          "No. of hits","counts");
getStatistics().createHistogramWithAxes(
                                          new TH1D("Hit_multiplicity_cut1","Multiplicity of hits(1st cut)", 10, 0.5, 10.5),
                                          "No. of hits","counts");
getStatistics().createHistogramWithAxes(
                                          new TH1D("Hit_multiplicity_cut2","Multiplicity of hits(2nd cut)", 10, 0.5, 10.5),
                                          "No. of hits","counts");
 getStatistics().createHistogramWithAxes(
                                          new TH1D("Hit_multiplicity_ann0","Multiplicity of hits(3rd cut)", 10, 0.5, 10.5),
                                          "No. of hits","counts");
getStatistics().createHistogramWithAxes(
                                          new TH1D("Hit_multiplicity_ann1","Multiplicity of hits(3rd cut)", 10, 0.5, 10.5),
                                          "No. of hits","counts");
  
 

  // Histograms for 2Gamma category
  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_Zpos", "Z-axis position of 2 gamma hits", 201, -50.25, 50.25),
    "Z axis position [cm]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_DLOR", "Delta LOR distance", 100, -0.25, 49.25),
    "Delta LOR [cm]", "Counts"
  );
  

  //Newly added
  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_absThetaDiff","Absolute angle difference of 2 gamma hits ", 180, 0, 180),
    "Hits theta diff [deg]", "Counts"
   );

  getStatistics().createHistogramWithAxes(
 new TH1D("Hit_multiplicity_2g","Multiplicity of hits", 10, 0.5, 10.5),
  "No. of hits","counts");

  /////////
  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_ThetaDiff", "Angle difference of 2 gamma hits ", 180, -0.5, 180.5),
    "Hits theta diff [deg]", "Counts"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_TimeDiff", "Time difference of 2 gamma hits", 6, -50.0, 1600.0),
    "Time Difference [ps]", "Number of Hit Pairs"
  );
  //Newly added
  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_ThetaDiff_after", "Angle difference of 2 gamma hits after cut ", 180, -0.5, 180.5),
    "Hits theta diff [deg]", "Counts"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_TimeDiff_after", "Time difference of 2 gamma hits after cut", 6, -50.0, 1600.0),
    "Time Difference [ps]", "Number of Hit Pairs"
  );
  //////

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_Dist", "B2B hits distance", 150, -0.5, 149.5),
    "Distance [cm]", "Number of Hit Pairs"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Annih_TOF", "Annihilation pairs Time of Flight", 201, -3015.0, 3015.0),
    "Time of Flight [ps]", "Number of Annihilation Pairs"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("AnnihPoint_XY", "XY position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "X position [cm]", "Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("AnnihPoint_ZX", "ZX position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Z position [cm]", "X position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("AnnihPoint_ZY", "ZY position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Z position [cm]", "Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Annih_DLOR", "Delta LOR distance of annihilation photons", 100, -0.25, 49.25),
    "Delta LOR [cm]", "Counts"
  );

  // Histograms for 3Gamama category
  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_Angles", "Relative angles - transformed", 250, -0.5, 249.5, 20, -0.5, 199.5),
    "Relative angle 1-2", "Relative angle 2-3"
  );

  getStatistics().createHistogramWithAxes(
                                          new TH1D("Hit_multiplicity_3g","Multiplicity of hits", 10, 0.5, 10.5),
                                          "No. of hits","counts");


  // Histograms for scattering category
  getStatistics().createHistogramWithAxes(
    new TH1D("ScatterTOF_TimeDiff", "Difference of Scatter TOF and hits time difference",
    3.0*fScatterTOFTimeDiff, -0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Scat_TOF - time diff [ps]", "Number of Hit Pairs"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("ScatterAngle_PrimaryTOT", "Angle of scattering vs. TOT of primary hits",
    200, -0.5, 199.5, 200, -100.0, 39900.0),
    "Scattering Angle", "TOT of primary hit [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("ScatterAngle_ScatterTOT", "Angle of scattering vs. TOT of scattered hits",
    200, -0.5, 199.5, 200, -100.0, 39900.0),
    "Scattering Angle", "TOT of scattered hit [ps]"
  );

  // Histograms for deexcitation
  //New
  getStatistics().createHistogramWithAxes(
    new TH1D("SYNC_TOT", "TOT of all hits at event level", 200, -100.0, 200000.0),
    "TOT [ps]", "Number of Hits"
  );
  ///
  getStatistics().createHistogramWithAxes(
    new TH1D("Deex_TOT_cut", "TOT of all hits with deex cut", 200, 85000.0, 140000.0),
    "TOT [ps]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
                                          new TH1D("Hit_multiplicity_prompt","Multiplicity of hits", 10, 0.5, 10.5),
                                          "No. of hits","counts");


  // Histogram for time differences between hit time and calculated hit time
  getStatistics().createHistogramWithAxes(
    new TH1D("Hit_timedifference", "Time difference of hit", 10000, -100.0, 21000000.0),
    "Time Difference [ps]", "counts"
  );

  //histogram for openingangles and scintillator corellation
  getStatistics().createHistogramWithAxes(
    new TH1D("Opening_angle_before", "opening angles between 2 hits before cut", 180, 0, 180),
    "Angle", "Number of Hit"
  );

 getStatistics().createHistogramWithAxes(
    new TH1D("Scintillator_before", "Difference in Scintillator ID before cut", 200, 0, 200),
    "Difference", "number of hits"
  );

getStatistics().createHistogramWithAxes(
    new TH1D("Opening_angle", "opening angles between 2 hits", 180, 0, 360),
    "Angle", "Number of Hit"
  );

 getStatistics().createHistogramWithAxes(
    new TH1D("Scintillator", "Difference in Scintillator ID", 200, 0, 200),
    "Difference", "number of hits"
  );


 getStatistics().createHistogramWithAxes(
    new TH2D("opening_angle_Vs_distance_before", "Opening Angle of hits vs. Distance between hits",
    360, 0, 360, 200, 0, 200),
    "Opening Angle", "Distance [cm]"
   );
   
 getStatistics().createHistogramWithAxes(
    new TH2D("opening_angle_Vs_distance", "Opening Angle of hits vs. Distance between hits",
    360, 0, 360, 200, 0, 200),
    "Opening angle", "Distance [cm]"
  );

 getStatistics().createHistogramWithAxes(
    new TH2D("opening_angle_Vs_z_before", "Opening Angle of hits vs. z",
    200, 0, 200, 160,-80, 80),
    "Opening angle", "Difference in z"
  );

 getStatistics().createHistogramWithAxes(
    new TH2D("opening_angle_Vs_z", "Opening Angle of hits vs. z",
    200, 0, 200, 160, -80, 80),
    "Opening angle", "Difference in z"
  );

 
 getStatistics().createHistogramWithAxes(
					 new TH1D("delta_time_before", "Difference in time between 2 hits before cut",500, 0, 30000),
    "#Delta_t[ps]", "counts"
  );

 getStatistics().createHistogramWithAxes(
					 new TH1D("delta_time", "Difference in time between 2 hits after cut",500, 0, 30000),
    "#Delta_t[ps]", "counts"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("opening_angle_Vs_scinID1", "Opening Angle of hits vs. scinID1",
    180, 0, 180, 200, 0, 200),
    "Opening angle", "scinID1"
  );
 getStatistics().createHistogramWithAxes(
    new TH2D("opening_angle_Vs_scinID2", "Opening Angle of hits vs. scinID2",
    180, 0, 180, 200, 0, 200),
    "Opening angle", "scinID1"
  );

 getStatistics().createHistogramWithAxes(
    new TH1D("Hit_multiplicity_n2","Multiplicity of hits(2 hITS)", 10, 0.5, 10.5),
    "No. of hits","counts"
  );
 
 getStatistics().createHistogramWithAxes(
   new TH1D("Hit_multiplicity_n3","Multiplicity of hits(3 HITS)", 10, 0.5, 10.5),
  "No. of hits","counts"
   );

 getStatistics().createHistogramWithAxes(
    new TH1D("Hit_multiplicity_n5","Multiplicity of hits(5 HITS)", 10, 0.5, 10.5),
    "No. of hits","counts"
  );


 getStatistics().createHistogramWithAxes(
   new TH3D("scID Distribution", "Angular distribution of scID",
   200, 0, 200, 200, 0, 200, 180, 0, 180), "scID1","scID2","opening angles"
   );

  getStatistics().createHistogramWithAxes(
    new TH2D("scinID1_Vs_scinID2", "scinID1 vs. scinID2",
    200, 0, 200, 200, 0, 200),
    "scinID1", "scinID2"
  );

  getStatistics().createHistogramWithAxes(
   new TH3D("scID Distribution_after", "Angular distribution of scID",
   200, 0, 200, 200, 0, 200, 180, 0, 180), "scID1","scID2","opening angles"
   );

  getStatistics().createHistogramWithAxes(
    new TH2D("distance_vs_time_diff_before", "distance_vs_time_diff",
	     150, 0, 150, 500, 0, 5000),"Distance[cm]","#Delta_t[ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("distance_vs_time_diff_before2(dis)", "distance_vs_time_diff",
             150, 0, 150,500, 0, 10000),
    "Distance[cm]", "#Delta_t[ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("distance_vs_time_diff_before3(dt)", "distance_vs_time_diff",
             150, 0, 150,500, 0, 5000),
    "Distance[cm]", "#Delta_t[ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("distance_vs_time_diff", "distance_vs_time_diff",
             150, 0, 150,500, 0, 5000),
    "Distance[cm]", "#Delta_t[ps]"
  );

  
  getStatistics().createHistogramWithAxes(
    new TH2D("sum_diff_angle_dist", "sum and difference of 2 smallest relative angles",
             250, 0, 250, 200, 0, 200),
   "#theta_{1}+#theta_{2}[degree]" , "#theta_{1}-#theta_{2}[degree]"
  );


  getStatistics().createHistogramWithAxes(
    new TH2D("xyz1", "distance_vs_time_diff(pair12)",
             150, 0, 150,250, 0, 2000000),
    "Distance[cm]","#Delta_t[ps]"
  );
    
getStatistics().createHistogramWithAxes(
    new TH2D("xyz2", "distance_vs_time_diff(pair13)",
             150, 0, 150,250, 0, 2000000),
    "Distance[cm]","#Delta_t[ps]"
  );

 getStatistics().createHistogramWithAxes(
    new TH2D("xyz3", "distance_vs_time_diff(pair23)",
             150,0,150,250, 0, 2000000),
    "Distance[cm]","#Delta_t[ps]"
  );

  getStatistics().createHistogramWithAxes(
   new TH1D("t1","t1", 250, 0, 5000000),
  "#Delta_t(12)[ps]","counts"
   );
     getStatistics().createHistogramWithAxes(
   new TH1D("t2","t2", 250, 0, 5000000),
  "Delta_t(13)[ps]","counts"
   );
     getStatistics().createHistogramWithAxes(
   new TH1D("t3","t3", 250, 0, 5000000),
  "Delta_t(23)[ps]","counts"
   );

     /*     getStatistics().createHistogramWithAxes(
   new TH2D("xy","scatter_test", 250, 0, 1700000,250, 0, 1700000),"1","2"
   );*/

     
     getStatistics().createHistogramWithAxes(
    new TH1D("hit_order1","time(t3-t2)", 100, 0, 15000),
    "#deltat_between_hits","counts" );

     getStatistics().createHistogramWithAxes(
   new TH1D("hit_order2","time(t2-t1)", 100, 0, 15000),
  "#deltat_between_hits","counts"
   );
getStatistics().createHistogramWithAxes(
   new TH1D("hit_order3","time(t3-t1)", 100, 0, 15000),
  "#deltat_between_hits[cm]","counts"
   );

 getStatistics().createHistogramWithAxes(
   new TH1D("time_difference","dt", 250, 0, 1500000),
  "#dt(12)[ps]","counts"
   );

 getStatistics().createHistogramWithAxes(
   new TH1D("distance","distance", 150, 0,150),
  "distance[cm]","counts"
   );


 //y for(auto i = 1; i<=10; i++)
   {
     getStatistics().createHistogramWithAxes(
	new TH2D("D_vs_dt", "distance_vs_dt", 150, 0, 150,250, 0, 1500000),
        "Distance[cm]","#Delta_t[ps]");

     getStatistics().createHistogramWithAxes(
        new TH2D("dt_vs_dt", "dt_vs_dt", 250, 0, 5000000,250, 0, 5000000),
        "#Delta_t[ps]","#Delta_t[ps]");

   }
 getStatistics().createHistogramWithAxes(new TH1D("3hit_tot", "TOT for 3 hits", 200, -100, 200000), "TOT[ps]", "counts");
 getStatistics().createHistogramWithAxes(new TH1D("5hit_tot", "TOT for 5 hits", 200, -100, 200000), "TOT[ps]", "counts");
 getStatistics().createHistogramWithAxes(new TH1D("efficiency",  "TOT efficiency",  200, -100.0, 200000.0), "TOT[ps]", "Efficiency");
	    
}
