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
#include <TFile.h>
#include <vector>
#include <TMath.h>
#include <TSystem.h>
#include <ctime>
#include <JPetMCHit/JPetMCHit.h>


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
      
      bool isInitialCut = EventCategorizerTools::initialCut(              
       event, getStatistics(), fSaveControlHistos, fHitCounter);
      
      bool isScattered = false;
      bool isAnnihilation = false;
      bool is2Gamma= false;
      bool is3Gamma= false;
      bool isPrompt= false;
      bool isNeighbourHits= false;
      bool isAdditionalCuts = false;


      const int atLeastNAnihilationHits = 2;
      const double TOT_Cut = 65000;
      const int atleastNprompt = 1;
           

      if (isInitialCut){
	TotalInitialcut.totalAccepted++;
	fAnnihilation.totalNumber++;
	fPrompt.totalNumber++;
	

       	isNeighbourHits = EventCategorizerTools::removeNeighbourhits(
       	event, getStatistics(),fSaveControlHistos, fTOTCalculationType);
	
	isScattered = EventCategorizerTools::checkForScatter(
	event, getStatistics(), fSaveControlHistos, fScatterTOFTimeDiff, fTOTCalculationType);

	if(!isScattered)
	  {
	    fEventAnnPrompt.totalNumber++;
	    fEventNoAnnPrompt.totalNumber++;
	    fEventNoPromptAnn.totalNumber++;
            fEventNoAnnNoPrompt.totalNumber++;

	    isPrompt = EventCategorizerTools::checkForPrompt(
	    event, getStatistics(), fSaveControlHistos, fDeexTOTCutMin, fDeexTOTCutMax, fTOTCalculationType, atleastNprompt);

            isAnnihilation = EventCategorizerTools::checkForAnnihilation(
            event, getStatistics(), fSaveControlHistos, atLeastNAnihilationHits, TOT_Cut);

	if (isAnnihilation) { fAnnihilation.totalAccepted++;}
	if(isPrompt){fPrompt.totalAccepted++;}
	
	/* Categories of events left after initial cut and scatter test */
		
	if (isAnnihilation && isPrompt ){
	  fEventAnnPrompt.totalAccepted++;
	  getStatistics().fillHistogram("Hit_multiplicity_case1", event.getHits().size());
	  for(auto hit : event.getHits()){
          double tot1 = hit.getEnergy();
	  double angle_z = TMath::RadToDeg() * hit.getPos().Theta();
	  getStatistics().fillHistogram("theta_case2", angle_z);
	  getStatistics().fillHistogram("TOT_case1", tot1);
	  }
	  
	}

	if (isAnnihilation && !isPrompt ){
          fEventNoPromptAnn.totalAccepted++;
	  getStatistics().fillHistogram("Hit_multiplicity_case2", event.getHits().size());
          isAdditionalCuts = EventCategorizerTools::additionalCuts(event, getStatistics(), fSaveControlHistos, fTOTCalculationType);
	  for(auto hit : event.getHits()){
          double tot2 = hit.getEnergy();
	//  double angle_z = TMath::RadToDeg() * hit.getPos().Theta();
	  getStatistics().fillHistogram("XYZ", hit.getPosX(), hit.getPosY(), hit.getPosZ());
	//  getStatistics().fillHistogram("theta_case2", angle_z);
          getStatistics().fillHistogram("TOT_case2", tot2);
          }
        }
	
        if (!isAnnihilation && isPrompt ){
	  fEventNoAnnPrompt.totalAccepted++;
	}
	
        if (!isAnnihilation && !isPrompt ){
	  fEventNoAnnNoPrompt.totalAccepted++;
	    }
	  }

      }
      
       
      JPetEvent newEvent = event;
      if(is2Gamma){  newEvent.addEventType(JPetEventType::k2Gamma); }
      if(is3Gamma){  newEvent.addEventType(JPetEventType::k3Gamma); }
      if(isPrompt) newEvent.addEventType(JPetEventType::kPrompt);
      // if(isScattered) newEvent.addEventType(JPetEventType::kScattered);
      
      
      if(fSaveControlHistos){
       for(auto hit : event.getHits()){
	  double tot = hit.getEnergy();
	  double hit_pos = hit.getPos().Mag();
	  double hit_time = hit.getTime();
	  double hit_timediff = fabs(hit_time - hit_pos/ kLightVelocity_cm_ps);
	  int scID = hit.getScintillator().getID();
	  double Z = hit.getPosZ();
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
  // INFO("Total hits:" + std::to_string(fHitCounter.totalNumber));
  // INFO("Total accepted hits:" + std::to_string(fHitCounter.totalAccepted));
  // INFO("Ratio of accepted hits: " + std::to_string(fHitCounter.getRatio()));
  INFO("Efficiency of 1st cut: "+ std::to_string(TotalInitialcut.getRatio()));
  INFO("Fraction of prompt events: "+ std::to_string(fPrompt.getRatio()));
  INFO("Fraction of Annihilation events: "+ std::to_string(fAnnihilation.getRatio()));
  //  INFO("Efficiency of 1st Case: "+ std::to_string(case1.getRatio()));
  INFO("Ratio of accepted events(EventNoPromptAnn): " + std::to_string(fEventNoPromptAnn.getRatio()));
  INFO("Ratio of accepted events(EventNoAnnPrompt): " + std::to_string(fEventNoAnnPrompt.getRatio()));
  INFO("Ratio of accepted events(EventAnnPrompt): " + std::to_string(fEventAnnPrompt.getRatio()));
  INFO("Ratio of accepted events(EventNoAnnNoPrompt):" + std::to_string(fEventNoAnnNoPrompt.getRatio()));
  fEventPromptCounter, fEventAnihilationCounter, fEventScatteredCounter, fEventNoPromptAnn, fEventNoAnnPrompt, fEventAnnPrompt, fEventNoAnnNoPrompt = {};

  /*
  auto file1 = TFile::Open("efficiency_hit.root", "recreate");
  // if (file) std::cout << "file was created"  <<std::endl;
  //Event_Eff->SetDirectory(gDirectory);
  Hit_Eff->SetDirectory(gDirectory);
  file1->Write();
  file1->Close();*/
  return true;
  
}

void EventCategorizer::saveEvents(const vector<JPetEvent>& events)
{
  for (const auto& event : events) { fOutputEvents->add<JPetEvent>(event); }
}

void EventCategorizer::initialiseHistograms(){

  /*******************************Eventcategorizer********************/
 getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_case1","Multiplicity of hits(ann&prompt)", 10, 0.5, 10.5), "No. of hits","counts");

 getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_case2","Multiplicity of hits(ann&_no_prompt)", 10, 0.5, 10.5), "No. of hits","counts");

 getStatistics().createHistogramWithAxes( new TH1D("TOT_case1", "TOT for anni&prompt hits", 200, -100, 200000), "TOT[ps]", "counts");

 getStatistics().createHistogramWithAxes(new TH1D("TOT_case2", "TOT for anni&no_prompt hits", 200, -100, 200000), "TOT[ps]", "counts");
  
 getStatistics().createHistogramWithAxes(new TH1D("theta_case2","theta_case", 200, 0, 360), "theta(degree)","counts");

getStatistics().createHistogramWithAxes( new TH1D("Hit_timedifference", "Time difference of hit", 10000, -100.0, 21000000.0), "Time Difference [ps]", "counts");
getStatistics().createHistogramWithAxes( new TH1D("lifetime", "Difference in time between 1st hit & last hit",500, 0, 250000), "#Delta_t[ps]", "counts");
 
for(int i = 0; i< 5; i++)
{
 getStatistics().createHistogramWithAxes(  new TH1D(Form("Individual_tot%d",i),"TOT", 200, 0, 200000), "tot(ps)","counts");

 //getStatistics().createHistogramWithAxes( new TH1D("Hit_timedifference%d", "Time difference of hit", 10000, -100.0, 21000000.0), "Time Difference [ps]", "counts");

// getStatistics().createHistogramWithAxes( new TH1D("lifetime%d", "Difference in time between 1st hit & last hit",500, 0, 250000), "#Delta_t[ps]", "counts");

 getStatistics().createHistogramWithAxes(new TH2D(Form("tot_vs_positionX%d",i), "TOT vs positionX", 200, 0, 20000, 200, -60, 60), "TOT(ps)", "position");

 getStatistics().createHistogramWithAxes(new TH2D(Form("tot_vs_positionY%d",i), "TOT vs positionY", 200, 0, 20000, 200, -60, 60), "TOT(ps)", "position");

 getStatistics().createHistogramWithAxes(new TH2D(Form("tot_vs_positionZ%d",i), "TOT vs positionZ", 200, 0, 20000, 200, -60, 60),  "TOT(ps)", "position");

}

/***************************************************/ 
  //histogram for openingangles and scintillator corellation
     /*     getStatistics().createHistogramWithAxes(new TH2D("xy","scatter_test", 250, 0, 1700000,250, 0, 1700000),"1","2");*/

     /*   
     getStatistics().createHistogramWithAxes(new TH1D("hit_order1","time(t3-t2)", 100, 0, 15000),"#deltat_between_hits","counts" );

     getStatistics().createHistogramWithAxes(new TH1D("hit_order2","time(t2-t1)", 100, 0, 15000),"#deltat_between_hits","counts");

getStatistics().createHistogramWithAxes(new TH1D("hit_order3","time(t3-t1)", 100, 0, 15000),"#deltat_between_hits[cm]","counts");*/


getStatistics().createHistogramWithAxes(new TH1D("Distance_between_1st and last hit","Distance_between_1st and last hit", 200, 0, 200), "distance","counts");

getStatistics().createHistogramWithAxes(new TH2D("X_5HITS","X coordinate of 1st and last hit", 240, -60, 60, 240, -60, 60), "X 1st hit","X 5th hit");

getStatistics().createHistogramWithAxes(new TH2D("Y_5HITS","Y coordinate of 1st and last hit", 240, -60, 60, 240, -60, 60), "Y 1st hit","Y 5th hit");

getStatistics().createHistogramWithAxes(new TH2D("Z_5HITS","Z coordinate of 1st and last hit", 100, -40, 40, 100, -40, 40), "Z 1st hit","Z 5th hit");
 
getStatistics().createHistogramWithAxes(new TH1D("opening angles between 5 hits","opening angles between 5 hits", 360, 0, 720), "theta(degree)","counts");

getStatistics().createHistogramWithAxes(new TH1D("sum_TOT_5","SUM tot of 5 hits", 200, 0, 20500), "TOT(ps)","counts");

getStatistics().createHistogramWithAxes(new TH1D("Angle_5","opening angles between 5 hits", 200, 0, 200),  "theta(degree)","counts");

getStatistics().createHistogramWithAxes(new TH1D("Angle_1","opening angles between 5 hits", 200, 0, 200),"theta(degree)","counts");

getStatistics().createHistogramWithAxes(new TH1D("Angle_2","opening angles between 5 hits", 200, 0, 200), "theta(degree)","counts");

getStatistics().createHistogramWithAxes(new TH1D("Angle_3","opening angles between 5 hits", 200, 0, 200), "theta(degree)","counts");

getStatistics().createHistogramWithAxes( new TH1D("Angle_4","opening angles between 5 hits", 200, 0, 200), "theta(degree)","counts");

getStatistics().createHistogramWithAxes(new TH1D("SAngle_1","opening angles between 5 hits", 200, 0, 200),  "theta(degree)","counts");

getStatistics().createHistogramWithAxes(new TH1D("SAngle_5","opening angles between 5 hits", 200, 0, 200), "theta(degree)","counts");
  
getStatistics().createHistogramWithAxes(new TH3D("XYZ", "XYZ", 240, -60.25, 59.75, 240, -60.25, 59.75, 240, -60.25, 59.75), "Hit X position [cm]", "Hit Y position [cm]",  "Hit Y position [cm]");
  
getStatistics().createHistogramWithAxes(new TH2D("All_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),"Hit X position [cm]", "Hit Y position [cm]");
  
getStatistics().createHistogramWithAxes(new TH2D("Z Vs scID", "Z Vs scID", 100, -50, 50, 200, 0, 200),"Scintillator ID", "Hit Z position [cm]");

getStatistics().createHistogramWithAxes(new TH2D("Anglevstot","Angle vs tot", 100, 0, 720, 100, 0,  2000000),"angle(deg)","tot(ps)");

/**********************check for 2 gamma******************/
getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_2g","Multiplicity of hits", 10, 0.5, 10.5),"No. of hits","counts");
 
 getStatistics().createHistogramWithAxes(new TH1D("2Gamma_Zpos", "Z-axis position of 2 gamma hits", 201, -50.25, 50.25),"Z axis position [cm]", "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("2Gamma_DLOR", "Delta LOR distance", 100, -0.25, 49.25),"Delta LOR [cm]", "Counts");
  
  getStatistics().createHistogramWithAxes(new TH1D("2Gamma_absThetaDiff","Absolute angle difference of 2 gamma hits ", 180, 0, 180),"Hits theta diff [deg]", "Counts");
  
  getStatistics().createHistogramWithAxes(new TH1D("2Gamma_ThetaDiff", "Angle difference of 2 gamma hits ", 180, -0.5, 180.5),"Hits theta diff [deg]", "Counts");

  getStatistics().createHistogramWithAxes(new TH1D("2Gamma_TimeDiff", "Time difference of 2 gamma hits", 6, -50.0, 1600.0),"Time Difference [ps]", "Number of Hit Pairs");
 
  getStatistics().createHistogramWithAxes(new TH1D("2Gamma_ThetaDiff_after", "Angle difference of 2 gamma hits after cut ", 180, -0.5, 180.5),"Hits theta diff [deg]", "Counts");

  getStatistics().createHistogramWithAxes(new TH1D("2Gamma_TimeDiff_after", "Time difference of 2 gamma hits after cut", 6, -50.0, 1600.0),"Time Difference [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("2Gamma_Dist", "B2B hits distance", 150, -0.5, 149.5),"Distance [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("Annih_TOF", "Annihilation pairs Time of Flight", 201, -3015.0, 3015.0),"Time of Flight [ps]", "Number of Annihilation Pairs");

  getStatistics().createHistogramWithAxes(new TH2D("AnnihPoint_XY", "XY position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75),"X position [cm]", "Y position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("AnnihPoint_ZX", "ZX position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75),"Z position [cm]", "X position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("AnnihPoint_ZY", "ZY position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75), "Z position [cm]", "Y position [cm]");

  getStatistics().createHistogramWithAxes(new TH1D("Annih_DLOR", "Delta LOR distance of annihilation photons", 100, -0.25, 49.25),"Delta LOR [cm]", "Counts");

/******************************************************/
/*********************** Histograms for 3Gamama category*********************/
getStatistics().createHistogramWithAxes(new TH2D("3Gamma_Angles", "Relative angles - transformed", 250, -0.5, 249.5, 20, -0.5, 199.5),"Relative angle 1-2", "Relative angle 2-3");

getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_3g","Multiplicity of hits", 10, 0.5, 10.5),"No. of hits","counts");

getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_prompt0","Multiplicity of hits before prompt cut", 10, 0.5, 10.5), "No. of hits","counts");

getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_prompt1","Multiplicity of hits after prompt cut", 10, 0.5, 10.5),"No. of hits","counts");


/************************************************************************/
/********************Histograms for prompt****************************/
  
getStatistics().createHistogramWithAxes(new TH1D("SYNC_TOT", "TOT of all hits at event level", 200, -100.0, 200000.0),"TOT [ps]", "Number of Hits");
  
getStatistics().createHistogramWithAxes(new TH1D("Deex_TOT_cut", "TOT of all hits with deex cut", 200, 85000.0, 140000.0),"TOT [ps]", "Number of Hits");

/**************************************************************/
/**************Additional cuts*****************/

getStatistics().createHistogramWithAxes(new TH1D("Residual_time", " Residual_time", 10, 0, 0.01), "Residual time", "counts");

getStatistics().createHistogramWithAxes(new TH1D("TOF_all_hits", "TOF_all_hits", 20000, -21000, 0), "TOF(ns)", "counts");

getStatistics().createHistogramWithAxes(new TH1D("meantime", "meantime of all hits", 20000, -21000, 100), "mean time(ns)", "counts");

/* 
for(int i = 0; i< 10 ; i++)
{

getStatistics().createHistogramWithAxes(new TH1D(Form("Hit_Multi%d",i),"Multiplicity of hits(ann&_no_prompt)", 10, 0.5, 10.5),"No. of hits","counts");

}
*/
getStatistics().createHistogramWithAxes(new TH1D("Residual_time4", " Residual_time", 10, 0, 0.01), "Residual time", "counts");

getStatistics().createHistogramWithAxes(new TH1D("TOF_4_hits", "TOF_4_hits", 20000, -21000, 0), "TOF(ns)", "counts");

getStatistics().createHistogramWithAxes(new TH1D("meantime_4hits", "meantime of 4 hits", 20000, -10000, 100), "mean time(ns)", "counts");

/*************************************/

/**********************check for Annhilation*************/

 getStatistics().createHistogramWithAxes(new TH1D("Ann_TOT_before_cut", "TOT of all hits before cut", 200, -100.0, 200000.0), "TOT [ps]", "Number of Hits");

 getStatistics().createHistogramWithAxes(new TH1D("Ann_TOT", "TOT of all hits after cut", 200, -100.0, 200000.0), "TOT [ps]", "Number of Hits");

 getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_ann0","Multiplicity of hits(3rd cut)", 10, 0.5, 10.5),"No. of hits","counts");

 getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_ann1","Multiplicity of hits(3rd cut)", 10, 0.5, 10.5), "No. of hits","counts");
 
/************************************************/ 
/***********************Initial cuts***************/

getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_0","Multiplicity of hits(no cut)", 10, 0.5, 10.5),"No. of hits","counts");

getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_cut1","Multiplicity of hits(1st cut)", 10, 0.5, 10.5),"No. of hits","counts");

getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_cut2","Multiplicity of hits(2nd cut)", 10, 0.5, 10.5),"No. of hits","counts");

 getStatistics().createHistogramWithAxes(new TH1D("Hit_Z_POS_before","Hit position Z ", 80, -40, 40), "Z_Hit_Pos[cm]","counts");

 getStatistics().createHistogramWithAxes(new TH1D("Hit_Z_POS","Hit position Z ", 80, -40, 40), "Z_Hit_Pos[cm]","counts");

/*********************************************/

/*****************Remove neighbouring hits*********************/

 getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_n2","Multiplicity of hits(2 hITS)", 10, 0.5, 10.5), "No. of hits","counts");
 
 getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_n3","Multiplicity of hits(3 HITS)", 10, 0.5, 10.5), "No. of hits","counts");

 getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_n5","Multiplicity of hits(5 HITS)", 10, 0.5, 10.5), "No. of hits","counts");

 getStatistics().createHistogramWithAxes(new TH1D("Opening_angle_before", "opening angles between 2 hits before cut", 180, 0, 180), "Angle", "Number of Hit");

 getStatistics().createHistogramWithAxes(new TH1D("Scintillator_before", "Difference in Scintillator ID before cut", 200, 0, 200), "Difference", "number of hits");

 getStatistics().createHistogramWithAxes(new TH1D("Opening_angle", "opening angles between 2 hits", 180, 0, 360), "Angle", "Number of Hit");

 getStatistics().createHistogramWithAxes(new TH1D("Scintillator", "Difference in Scintillator ID", 200, 0, 200), "Difference", "number of hits");

 getStatistics().createHistogramWithAxes(new TH2D("opening_angle_Vs_distance_before", "Opening Angle of hits vs. Distance between hits", 360, 0, 360, 200, 0, 200),"Opening Angle", "Distance [cm]");
   
 getStatistics().createHistogramWithAxes(new TH2D("opening_angle_Vs_distance", "Opening Angle of hits vs. Distance between hits",360, 0, 360, 200, 0, 200),"Opening angle", "Distance [cm]");

 getStatistics().createHistogramWithAxes(new TH2D("opening_angle_Vs_z_before", "Opening Angle of hits vs. z",200, 0, 200, 160,-80, 80),"Opening angle", "Difference in z");

 getStatistics().createHistogramWithAxes(new TH2D("opening_angle_Vs_z", "Opening Angle of hits vs. z",200, 0, 200, 160, -80, 80),"Opening angle", "Difference in z");

 getStatistics().createHistogramWithAxes(new TH1D("delta_time_before", "Difference in time between 2 hits before cut",500, 0, 30000),"#Delta_t[ps]", "counts");

 getStatistics().createHistogramWithAxes(new TH1D("delta_time", "Difference in time between 2 hits after cut",500, 0, 30000), "#Delta_t[ps]", "counts");

 getStatistics().createHistogramWithAxes(new TH2D("opening_angle_Vs_scinID1", "Opening Angle of hits vs. scinID1",180, 0, 180, 200, 0, 200), "Opening angle", "scinID1");
 
 getStatistics().createHistogramWithAxes(new TH2D("opening_angle_Vs_scinID2", "Opening Angle of hits vs. scinID2",180, 0, 180, 200, 0, 200),"Opening angle", "scinID1");

 getStatistics().createHistogramWithAxes(new TH3D("scID Distribution", "Angular distribution of scID",  200, 0, 200, 200, 0, 200, 180, 0, 180), "scID1","scID2","opening angles");

 getStatistics().createHistogramWithAxes(new TH2D("scinID1_Vs_scinID2", "scinID1 vs. scinID2",200, 0, 200, 200, 0, 200),"scinID1", "scinID2");

 getStatistics().createHistogramWithAxes(new TH3D("scID Distribution_after", "Angular distribution of scID",200, 0, 200, 200, 0, 200, 180, 0, 180),"scID1","scID2","opening angles");

 getStatistics().createHistogramWithAxes(new TH2D("distance_vs_time_diff_before", "distance_vs_time_diff",150, 0, 150, 500, 0, 5000),"Distance[cm]","#Delta_t[ps]");

 getStatistics().createHistogramWithAxes(new TH2D("distance_vs_time_diff_before2(dis)", "distance_vs_time_diff", 150, 0, 150,500, 0, 10000),"Distance[cm]", "#Delta_t[ps]");

 getStatistics().createHistogramWithAxes(new TH2D("distance_vs_time_diff_before3(dt)", "distance_vs_time_diff", 150, 0, 150,500, 0, 5000), "Distance[cm]", "#Delta_t[ps]");

 getStatistics().createHistogramWithAxes(new TH2D("distance_vs_time_diff", "distance_vs_time_diff",150, 0, 150,500, 0, 5000), "Distance[cm]", "#Delta_t[ps]");

 getStatistics().createHistogramWithAxes(new TH2D("sum_diff_angle_dist", "sum and difference of 2 smallest relative angles", 250, 0, 250, 200, 0, 200),"#theta_{1}+#theta_{2}[degree]" , "#theta_{1}-#theta_{2}[degree]");

 getStatistics().createHistogramWithAxes(new TH2D("dist_vs_M_3h", "distance_vs_time_diff", 150, 0, 150,25, -0.5, 3.0*fScatterTOFTimeDiff-0.5),"Distance[cm]","#Delta_t[ps]");
  
 getStatistics().createHistogramWithAxes(new TH1D("time_difference_3hits","t1", 25, -0.5, 4.0*fScatterTOFTimeDiff-0.5), "#Delta_t(12)[ps]","counts");
 
 getStatistics().createHistogramWithAxes(new TH1D("time_difference_2hits","dt", 26, -0.5, 4.0*fScatterTOFTimeDiff-0.5),"#dt(12)[ps]","counts");

 getStatistics().createHistogramWithAxes(new TH1D("distance","distance", 150, 0,150),"distance[cm]","counts");

 getStatistics().createHistogramWithAxes(new TH2D("D_vs_dt", "distance_vs_dt(5hits)", 150, 0, 150, 26, -0.5, 3.0*fScatterTOFTimeDiff-0.5),"Distance[cm]","#Delta_t[ps]");

 getStatistics().createHistogramWithAxes(new TH2D("dt_vs_dt", "dt_vs_dt", 250, 0, 5000000,250, 0, 1500000), "#Delta_t[ps]","#Delta_t[ps]");
 
 getStatistics().createHistogramWithAxes(new TH1D("3hit_tot", "TOT for 3 hits", 200, -100, 200000),"TOT[ps]", "counts");
 
 getStatistics().createHistogramWithAxes(new TH1D("5hit_tot", "TOT for 5 hits", 200, -100, 200000), "TOT[ps]", "counts");
 
 getStatistics().createHistogramWithAxes(new TH1D("dt_5_hits", "Time _diff_for 5 hits",26, -0.5, 4.0*fScatterTOFTimeDiff-0.5), "Time_diff[ps]", "counts");

/*************************************************/

/*************************scatter_test_check**********************/
 getStatistics().createHistogramWithAxes(new TH2D("Scatt_Dist_vs_M_before", "Scatt_Dist_vs_M_before", 150, 0, 150, 26, -0.5, 3.0*fScatterTOFTimeDiff-0.5),"Distance[cm]","#Delta_t[ps]");
 
 getStatistics().createHistogramWithAxes(new TH2D("Scatt_Dist_vs_M", "Scatt_Dist_vs_M_before", 150, 0, 150, 26, -0.5, 3.0*fScatterTOFTimeDiff-0.5), "Distance[cm]","#Delta_t[ps]");
 
 getStatistics().createHistogramWithAxes(new TH1D("Scatt_time_diff_before","dt", 26, -0.5, 3.0*fScatterTOFTimeDiff-0.5),"#dt(12)[ps]","counts");
 
 getStatistics().createHistogramWithAxes(new TH1D("Scatt_time_diff","dt", 26, -0.5, 3.0*fScatterTOFTimeDiff-0.5),"#dt(12)[ps]","counts");

 getStatistics().createHistogramWithAxes(new TH2D("ScatterAngle_PrimaryTOT_before", "Angle of scattering vs. TOT of primary hits_before_cut",200, -0.5, 199.5, 200, -100.0, 39900.0),"Scattering Angle", "TOT of primary hit [ps]");

getStatistics().createHistogramWithAxes(new TH2D("ScatterAngle_ScatterTOT_before", "Angle of scattering vs. TOT of scattered hits_before_cut",200, -0.5, 199.5, 200, -100.0, 39900.0),"Scattering Angle", "TOT of scattered hit [ps]");
  
getStatistics().createHistogramWithAxes(new TH1D("ScatterTOF_TimeDiff", "Difference of Scatter TOF and hits time difference", 25, -0.5, 3.0*fScatterTOFTimeDiff-0.5),"Scat_TOF - time diff [ps]", "Number of Hit Pairs");

getStatistics().createHistogramWithAxes(new TH2D("ScatterAngle_PrimaryTOT", "Angle of scattering vs. TOT of primary hits",200, -0.5, 199.5, 200, -100.0, 39900.0), "Scattering Angle", "TOT of primary hit [ps]");

getStatistics().createHistogramWithAxes(new TH2D("ScatterAngle_ScatterTOT", "Angle of scattering vs. TOT of scattered hits",200, -0.5, 199.5, 200, -100.0, 39900.0), "Scattering Angle", "TOT of scattered hit [ps]");

getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_scatt0","Multiplicity of hits(scattered)", 10, 0.5, 10.5), "No. of hits","counts");

getStatistics().createHistogramWithAxes(new TH1D("Hit_multiplicity_scatt1","Multiplicity of hits(scattered)", 10, 0.5, 10.5),"No. of hits","counts");

/*****************************************************/
	    
}
