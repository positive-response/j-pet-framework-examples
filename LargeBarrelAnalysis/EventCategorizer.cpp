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
#include <JPetMCHit/JPetMCHit.h>
#include "HitFinderTools.h"

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
  if(fSaveControlHistos) initialiseHistograms();
  return true;
}

bool EventCategorizer::exec()
{ 	
	if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {
	auto timeWindowMC = dynamic_cast<const JPetTimeWindowMC* const>(fEvent);
	vector<JPetEvent> events;
    
	for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++) {
	  	const auto& event = dynamic_cast<const JPetEvent&>(timeWindow->operator[](i));
		bool isMLParams = false;

		bool isInitialCutPassed = EventCategorizerTools::checkForInitialCuts(
		event, getStatistics(), fSaveControlHistos, flowEnergyCut, fAnnihilationEnergyCut
		);

		if(isInitialCutPassed)
		{
			if(timeWindowMC)
			{
				isMLParams = EventCategorizerTools::calculateMLParamsBefore(event, getStatistics(), fSaveControlHistos); 
			}

			if(event.getHits().size() == 5)
			{
				isMLParams = EventCategorizerTools::calculateMLParamsAfter(event, getStatistics(), fSaveControlHistos); 
			}

		}

		for(const auto & hit : event.getHits())
		{
			if(timeWindowMC)
			{
				auto mcHit = timeWindowMC->getMCHit<JPetMCHit>(hit.getMCindex());
			        getStatistics().fillHistogram("True Energy", mcHit.getEnergy());
			        getStatistics().fillHistogram("Gengamma_multi_all", mcHit.getGenGammaMultiplicity());
			}
		}
		getStatistics().fillHistogram("Multiplicity_all", event.getHits().size());
		

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

      JPetEvent newEvent = event;
      if(is2Gamma) newEvent.addEventType(JPetEventType::k2Gamma);
      if(is3Gamma) newEvent.addEventType(JPetEventType::k3Gamma);
      if(isPrompt) newEvent.addEventType(JPetEventType::kPrompt);
      if(isScattered) newEvent.addEventType(JPetEventType::kScattered);

      if(fSaveControlHistos){
        for(auto hit : event.getHits()){
          getStatistics().fillHistogram("All_XYpos", hit.getPosX(), hit.getPosY());
		  getStatistics().fillHistogram("Energy_all", hit.getEnergy());
		  

        }
      }
      events.push_back(newEvent);
    }
    saveEvents(events);
  } else { return false; }
  return true;
}

bool EventCategorizer::terminate()
{
  INFO("Event categorization completed.");
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
////////////////for after
 getStatistics().createHistogramWithAxes(new TH1D("M_ij", "Difference in hit time and time calculated by d/c",3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5 ), "M_ij", "counts");
  
  getStatistics().createHistogramWithAxes(new TH1D("TOF", "Time of flight",3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5 ), "TOF", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("3DOpenAngle", "3D opening angles",180, 0, 180 ), "Angle", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("2DOpenAngle", "2D opening angle", 360, 0, 360), "Angle", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("Theta", "Theta",100, 50, 150 ), "Angle", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("Phi", "Phi", 400, -200, 200),  "Angle", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("SmallestEnergy", "Smallest Energy",500, 0, 600 ), "smallest Energy", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("LargestEnergy", "Largest Energy", 500, 0, 1500), "Largest Energy", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("EnergySum", "Sum of energies of hits", 500, 0, 1500), "EnergySum", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("EnergySum-largest", "EnergySum-largest", 500, 0, 1500), "EnergySum-largest", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("EnergySum-Smallest", "EnergySum-Smallest",  500, 0, 1500), "EnergySum-Smallest", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("Hit_Z_POS(before)", "Hit_Z_POS(before)", 100, -50, 50), "Hit_Z_POS(before)", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("Energy(before)", "Energy(before)", 750, 0, 1200), "Energy(before)", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("Hit_Z_POS(after)", "Hit_Z_POS(after)", 100, -50, 50), "Hit_Z_POS(after)", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("Energy(after)", "Energy(after)", 500, 0, 800), "Energy(after)", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("Multiplicity_all","Multiplicity of all hits(without cut)", 10, 0.5, 10.5), "No. of hits","counts");
  getStatistics().createHistogramWithAxes( new TH1D("Energy_all", "Energy of hits", 750, 0, 1200), "Energy(keV)", "counts");
  getStatistics().createHistogramWithAxes( new TH1D("True Energy", "Energy of hits(true)", 750, 0, 1200), "Energy(keV)", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("Gengamma_multi_all", "Generated gamma multiplicity_all", 1000, 0.5, 10000.5), "No. of hits","counts");


  //2d plots

  getStatistics().createHistogramWithAxes(new TH2D("TOF_vs_distance", "TOF_vs_distance", 3.0*fScatterTOFTimeDiff, -0.5, 3.0*fScatterTOFTimeDiff-0.5, 150, 0, 150), "TOF[ps]", "Distance[cm]");
 getStatistics().createHistogramWithAxes(new TH2D("TOF_vs_timeDifference", "TOF_vs_timeDifference", 5.0*fScatterTOFTimeDiff, -0.5*fScatterTOFTimeDiff, 5.0*fScatterTOFTimeDiff-0.5, 5.0*fScatterTOFTimeDiff, -0.5*fScatterTOFTimeDiff, 5.0*fScatterTOFTimeDiff-0.5), "TOF[ps]", "Time Difference[ps]" );
  getStatistics().createHistogramWithAxes(new TH2D("EnergySum-Smallest_vs_EnergySum-largest", "EnergySum-Smallest_vs_EnergySum-largest", 500, 0, 1500, 500, 0, 1500), "Energy_sum -smallestE(keV)", "Energy_sum -largestE(keV)");

  ///////////////////////////////////after/ends here

 getStatistics().createHistogramWithAxes(new TH1D("M_ijB", "Difference in hit time and time calculated by d/c",3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5 ), "M_ij", "counts");
  
  getStatistics().createHistogramWithAxes(new TH1D("TOFB", "Time of flight",3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5 ), "TOF", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("3DOpenAngleB", "3D opening angles",180, 0, 180 ), "Angle", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("2DOpenAngleB", "2D opening angle", 360, 0, 360), "Angle", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("ThetaB", "Theta",100, 50, 150 ), "Angle", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("PhiB", "Phi", 400, -200, 200),  "Angle", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("SmallestEnergyB", "Smallest Energy",500, 0, 600 ), "smallest Energy", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("LargestEnergyB", "Largest Energy", 500, 0, 1500), "Largest Energy", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("EnergySumB", "Sum of energies of hits", 500, 0, 1500), "EnergySum", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("EnergySum-largestB", "EnergySum-largest", 500, 0, 1500), "EnergySum-largest", "counts");
  getStatistics().createHistogramWithAxes(new TH1D("EnergySum-SmallestB", "EnergySum-Smallest",  500, 0, 1500), "EnergySum-Smallest", "counts");

  //2d plots

  getStatistics().createHistogramWithAxes(new TH2D("TOF_vs_distanceB", "TOF_vs_distance", 3.0*fScatterTOFTimeDiff, -0.5, 3.0*fScatterTOFTimeDiff-0.5, 150, 0, 150), "TOF[ps]", "Distance[cm]");
 getStatistics().createHistogramWithAxes(new TH2D("TOF_vs_timeDifferenceB", "TOF_vs_timeDifference", 5.0*fScatterTOFTimeDiff, -0.5*fScatterTOFTimeDiff, 5.0*fScatterTOFTimeDiff-0.5, 5.0*fScatterTOFTimeDiff, -0.5*fScatterTOFTimeDiff, 5.0*fScatterTOFTimeDiff-0.5), "TOF[ps]", "Time Difference[ps]" );
  getStatistics().createHistogramWithAxes(new TH2D("EnergySum-Smallest_vs_EnergySum-largestB", "EnergySum-Smallest_vs_EnergySum-largest", 500, 0, 1500, 500, 0, 1500), "Energy_sum -smallestE(keV)", "Energy_sum -largestE(keV)");
  // Histograms for 2Gamma category
  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_Zpos", "Z-axis position of 2 gamma hits", 201, -50.25, 50.25),
    "Z axis position [cm]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_DLOR", "Delta LOR distance", 100, -0.25, 49.25),
    "Delta LOR [cm]", "Counts"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_ThetaDiff", "Angle difference of 2 gamma hits ", 181, -0.5, 180.5),
    "Hits theta diff [deg]", "Counts"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_TimeDiff", "Time difference of 2 gamma hits", 200, -10100.0, 99900.0),
    "Time Difference [ps]", "Number of Hit Pairs"
  );

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
  getStatistics().createHistogramWithAxes(
    new TH1D("Deex_TOT_cut", "TOT of all hits with deex cut (30,50) ns", 200, 24950.0, 54950.0),
    "TOT [ps]", "Number of Hits"
  );
}
