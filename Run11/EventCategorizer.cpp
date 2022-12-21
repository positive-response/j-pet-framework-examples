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

using namespace jpet_options_tools;
using namespace std;

int a, counta = 0;
int b, countb = 0;

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
      b++;
      // Check types of current event
      bool is2Gamma = EventCategorizerTools::checkFor2Gamma(
        event, getStatistics(), fSaveControlHistos, fB2BSlotThetaDiff, fMaxTimeDiff
      );
      bool is3Gamma = EventCategorizerTools::checkFor3Gamma(
        event, getStatistics(), fSaveControlHistos
      );
      bool isPrompt = EventCategorizerTools::checkForPrompt(
        event, getStatistics(), fSaveControlHistos, fDeexTOTCutMin, fDeexTOTCutMax, fTOTCalculationType
      );
      bool isScattered = EventCategorizerTools::checkForScatter(
        event, getStatistics(), fSaveControlHistos, fScatterTOFTimeDiff, fTOTCalculationType
      );

      bool isAnnihilation = EventCategorizerTools::checkForAnnihilation(
	event, getStatistics(), fSaveControlHistos, fTOTCalculationType
      );

      bool isNeighbourHits = EventCategorizerTools::removeNeighbourhits(
	event, getStatistics(),fSaveControlHistos, fTOTCalculationType
      );


      double sum_tot=0.0;
      double sum_tot_2g= 0.0;
      double sum_tot_3g=0.0;
      double sum_tot_scatter=0.0;
      double sum_tot_prompt=0.0;
      double sum_tot_ann=0.0;
      double sum_tot_ann_prompt= 0.0;
      double sum_tot_3gann_prompt=0.0;
      
      JPetEvent newEvent = event;
      if(is2Gamma) newEvent.addEventType(JPetEventType::k2Gamma);
      if(is3Gamma) newEvent.addEventType(JPetEventType::k3Gamma);
      if(isPrompt) newEvent.addEventType(JPetEventType::kPrompt);
      if(isScattered) newEvent.addEventType(JPetEventType::kScattered);
      //      if(isAnnihilation) newEvent.addEventType(JPetEventType::kAnnihilation);
      if(isAnnihilation) a++;
      if(isNeighbourHits);

      if(fSaveControlHistos){
        for(auto hit : event.getHits()){
	  double tot = hit.getEnergy();
	  double hit_pos = hit.getPos().Mag();
	  double hit_time = hit.getTime();
	  double hit_timediff = fabs(hit_time - hit_pos/ kLightVelocity_cm_ps);
	  
	  sum_tot = sum_tot+tot;
	  
	  if(is2Gamma)
	    {
	      sum_tot_2g += tot;
	    }
	  
	  if(is3Gamma)
            {
              sum_tot_3g += tot;
            }
	  
	  if(is2Gamma || is3Gamma)
            {
              sum_tot_ann += tot;
            }

	  if(isPrompt)
            {
              sum_tot_prompt += tot;
            }

	  if(isScattered)
            {
              sum_tot_scatter += tot;
            }

	   if(!isScattered)
            {
              sum_tot_ann_prompt += tot;
            }

	   if(is3Gamma || isPrompt)
            {
              sum_tot_3gann_prompt +=tot;
            }

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

	    
        }
      }
      events.push_back(newEvent);
      
    }
    countb=countb+b;
    counta=counta+a;
    b=0; a=0;
    saveEvents(events);
  }
  else { return false; }
  return true;

}

bool EventCategorizer::terminate()
{
  INFO("Event categorization completed.");
  cout<<"number of events before categorization:"<<countb<<endl;
  cout<<"number of events after categorization:"<<counta<<endl;
  //  cout<<"hit_before:"<<hcountb<<endl;
  // cout<<"hit_after:"<<hcounta<<endl;
 
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
    new TH1D("Ann_TOT", "TOT of all hits after cut", 200, -100.0, 200000.0),
    "TOT [ps]", "Number of Hits"
   );

  
 

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

  // Histogram for time differences between hit time and calculated hit time
  getStatistics().createHistogramWithAxes(
    new TH1D("Hit_timedifference", "Time difference of hit", 100, -50.0, 3000.0),
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
    new TH1D("Opening_angle", "opening angles between 2 hits", 180, 0, 180),
    "Angle", "Number of Hit"
  );

 getStatistics().createHistogramWithAxes(
    new TH1D("Scintillator", "Difference in Scintillator ID", 200, 0, 200),
    "Difference", "number of hits"
  );


}
