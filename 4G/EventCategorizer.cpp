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

const double EventCategorizer::kUnknownEventType = 66666666;

EventCategorizer::EventCategorizer(const char* name): JPetUserTask(name) {}

EventCategorizer::~EventCategorizer() {}

bool EventCategorizer::init()
{
  INFO("Event categorization started.");

pTree22 = new TTree("T","flat tree");
pTree22->Branch("eventNumber",&fEventNumber,"eventNumber/I");
pTree22->Branch("numberOfHits",&fNumberOfHits,"numberOfHits/I");
pTree22->Branch("x",&fPosX);
pTree22->Branch("y",&fPosY);
pTree22->Branch("z",&fPosZ);
pTree22->Branch("energy",&fEnergy);
pTree22->Branch("time",&fTime);
pTree22->Branch("hitType",&fHitType);
pTree22->Branch("vtxIndex",&fVtxIndex);


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
if (isOptionSet(fParams.getOptions(), kDataType)) {
    kDType = getOptionAsBool(fParams.getOptions(), kDataType);
  } else {
    WARNING("No data type option given by the user. Using default data type.");
  }

  // Input events type
  fOutputEvents = new JPetTimeWindow("JPetEvent");
  // Initialise hisotgrams
  if(fSaveControlHistos) initialiseHistograms();
  return true;
}

bool EventCategorizer::exec()
{

fEventNumber = 0;
fNumberOfHits = 0; 

JPetTimeWindowMC *timeWindowMC = nullptr;
 if (timeWindowMC = dynamic_cast<JPetTimeWindowMC* const>(fEvent)) {
    fIsMC = true;
    INFO("The input file is MC.");
  } else {
    INFO("The input file is DATA.");
  }

if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {

//auto timeWindowMC = dynamic_cast<const JPetTimeWindowMC* const>(fEvent);
vector<JPetEvent> events;
  
for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++) {
fEnergy.clear();
fTime.clear();
fPosX.clear();
fPosY.clear();
fPosZ.clear();
fHitType.clear();
fVtxIndex.clear();

const auto& event = dynamic_cast<const JPetEvent &>(timeWindow->operator[](i));		

bool cut1 = false;
bool cut2 = false;
bool cut3 = false;
bool cut4 = false;
bool cut5 = false;
bool cut5A = false;  //combined cut
bool cut6 = false;
bool scatter = false;

EventCategorizerTools::Plots1hit(event, getStatistics(), "noCut");
EventCategorizerTools::Plots2hits(event, getStatistics(), "noCut");
EventCategorizerTools::Plots3hits(event, getStatistics(), "noCut");

JPetEvent sel1Event;  //for z
JPetEvent sel2Event;  // for  same scintillator check
JPetEvent sel3Event;  //for tot
JPetEvent sel5Event;
JPetEvent sel6Event;
JPetEvent STEvents;

cut1 = EventCategorizerTools::checkForCut1(event, getStatistics(), sel1Event); //z cut
cut3 = EventCategorizerTools::checkForCut3(sel1Event, getStatistics(), sel2Event); //same scintillator check
cut2 = EventCategorizerTools::checkForCut2(sel2Event, getStatistics(), sel3Event, 2);  //tot cut

if(sel3Event.getHits().size() > 3){
cut5 = EventCategorizerTools::checkForCut5(sel3Event, getStatistics(), sel5Event);  //ETS CHECK
cut6 = EventCategorizerTools::checkForCut6(sel5Event, getStatistics(), sel6Event);   //theta sum cut
scatter = EventCategorizerTools::checkForScat( sel6Event, getStatistics(), STEvents);
}


fEventNumber = i;
fNumberOfHits = event.getHits().size();

for(int i = 0; i < fNumberOfHits; i++)
{
	const JPetHit &hit = event.getHits().at(i);
  const JPetMCHit & mcHit = timeWindowMC->getMCHit<JPetMCHit>(hit.getMCindex());

	fPosX.push_back(hit.getPosX());
	fPosY.push_back(hit.getPosY());
	fPosZ.push_back(hit.getPosZ());
	fEnergy.push_back(hit.getEnergy()); 
	fTime.push_back(hit.getTime());
  getStatistics().fillHistogram("Energy_all", hit.getEnergy()); 
	if(fIsMC)
		{
			fHitType.push_back(mcHit.getGenGammaMultiplicity());
			fVtxIndex.push_back(mcHit.getMCVtxIndex());

      //if(mcHit.getGenGammaMultiplicity() == 4)
      {
			    getStatistics().fillHistogram("True Energy", mcHit.getEnergy());		        			
		   //	getStatistics().fillHistogram("Gengamma_multi_all", mcHit.getGenGammaMultiplicity());		
      }			
		}/* else {
			fHitType.push_back(kUnknownEventType);	
			fVtxIndex.push_back(0);
		} */ 
}

pTree22->Fill();
fEnergy.clear();
fTime.clear();
fPosX.clear();
fPosY.clear();
fPosZ.clear();
fHitType.clear();
fVtxIndex.clear();
/*	
	if ( std::equal(fVtxIndex.begin() + 1, fVtxIndex.end(), fVtxIndex.begin()  ))
	{	
		fIsAcc = kFALSE;
	} else {
		fIsAcc = kTRUE;
	}
	for(auto i : fHitType){
	auto number = i%10;
	auto n = 1;
	if(i%100!=0 && number!=0 && i>10)
		fIsScattered=kTRUE;
	if(number==1)
		fContainsPrompt=kTRUE;
	else if(number==2)
		fIsPickOff=kTRUE;
	else if(number==3)
		fIsOPs=kTRUE;
	else if(number==0){
		fIsSecondary=kTRUE;	
		while(number==0&&i>0){
		n = n*10;
		number = (i/n)%10;
		if(i%n*100 == 0 && number!=0 && i/n>10)
			fIsScattered=kTRUE;
		if(number==1)
			fContainsPrompt=kTRUE;
		else if(number==2)
			fIsPickOff=kTRUE;
		else if(number==3)
			fIsOPs=kTRUE;
		}
	}
}
	fIsSecondary = kFALSE;
	fIsScattered = kFALSE;
	fContainsPrompt = kFALSE;
	fIsPickOff = kFALSE;
	fIsOPs = kFALSE;
*/
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

    /*  if(fSaveControlHistos){
  		for(auto hit : event.getHits()){	
       // getStatistics().fillHistogram("All_XYpos", hit.getPosX(), hit.getPosY());
       // getStatistics().fillHistogram("Energy_all", hit.getEnergy());  
        }
  }*/
      events.push_back(newEvent);
    }
    saveEvents(events);
  } else { 
    return false; 
    }
  return true;
}

bool EventCategorizer::terminate()
{
  INFO("FROM PARAMS: "+ std::to_string(kDType) );
  INFO("Event categorization completed.");
/*  if(!pTree22) std::cout << "wtf??" << std::endl;
  file->cd();
  pTree22->Write();
  file->Write("",TObject::kWriteDelete);
  file->Close();*/
  return true;
}

void EventCategorizer::saveEvents(const vector<JPetEvent>& events)
{
  for (const auto& event : events) { fOutputEvents->add<JPetEvent>(event); }
}

void EventCategorizer::initialiseHistograms(){

  auto minScinID = 0;
  auto maxScinID = 200;
  bool MC = kDType;
  bool Data = !MC;

const std::string fix[9] = {"noCut", "zCut", "totCut", "SCCut", "ETSCut",   "ttCut", "STCut", "STCutRej", "CombCut"};
  for(int i = 0; i < 9; i++)
  {
    const std::string postfix = fix[i];
    getStatistics().createHistogramWithAxes(new TH2D(
    ("ap_XY " + postfix).c_str( ), " Annihilation point XY", 
    120, -60, 60, 120, -60, 60), "  X [cm]", "Y [cm]");

    //getStatistics().createHistogramWithAxes(new TH2D(
    //("ap_XZ " + postfix).c_str( ), " Annihilation point XZ", 
    //120, -60, 60, 120, -60, 60), "  X [cm]", "Z [cm]");

    //getStatistics().createHistogramWithAxes(new TH2D(
    //("ap_YZ " + postfix).c_str( ), " Annihilation point YZ", 
    //120, -60, 60, 120, -60, 60), "  Y [cm]", "Z [cm]");

    getStatistics().createHistogramWithAxes(new TH1D(
    ("HM " + postfix).c_str(), "Number of hits", 
    40, -0.5, 39.5), "Number of hits", "Counts");

    
    
    //getStatistics().createHistogramWithAxes(new TH1D(
    //("TOTSum " + postfix).c_str(), "TOT gamma", 
    //201, 0.0, 3*fToTHistoUpperLimit), "Time over Threshold [ps]", "Counts");

    getStatistics().createHistogramWithAxes(new TH1D(
    ("Z " + postfix).c_str(), "Zpos", 
    100, -50.0, 50.0), "z position of hits [cm]", "Counts");

    getStatistics().createHistogramWithAxes(new TH1D(
    ("planeCentre " + postfix).c_str(), "planeCentre", 
    65, 0, 65.0), "distance of plane from centre [cm]", "Counts");

    getStatistics().createHistogramWithAxes(new TH1D(
    ("ETS " + postfix).c_str(), "Emission time spread", 
    201,  0, 2*fMaxTimeDiff), "Emission time spread  [ps]", "Counts");

    getStatistics().createHistogramWithAxes(new TH2D(
    ("ETS_vs_DOP " + postfix).c_str(), "ETS VS PlaneCentre",
     30,  0, fMaxTimeDiff, 60, 0, 60.0), "Emission time spread [ps]", "distance of plane from centre [cm]");


   // getStatistics().createHistogramWithAxes(new TH1D(
   //("Phi " + postfix).c_str(), "Phi",
   //50, -200, 200),  "#phi Angle", "counts");

   // getStatistics().createHistogramWithAxes(new TH1D(
   //("Theta " + postfix).c_str(), "Theta",
   // 200, 0, 200 ), "#theta Angle", "counts");

    //getStatistics().createHistogramWithAxes(new TH2D(
    //("hit_XY " + postfix).c_str( ), " Hit position", 
    //120, -60, 60, 120, -60, 60), "  X [cm]", "Y [cm]");

    /*getStatistics().createHistogramWithAxes(new TH2D(
    ("scinID_VS_TOT " + postfix).c_str(), "scintillator Id vs tot", 
    maxScinID - minScinID + 1, minScinID - 0.5, maxScinID + 0.5, 201, -1.0, 200.0 ), "Scintillator ID", "TOT [ps]");

    getStatistics().createHistogramWithAxes(new TH1D(
    ("SmallestEnergy " + postfix).c_str(), "Smallest Energy",
    500, 0, 200 ), "smallest Energy", "counts");

    getStatistics().createHistogramWithAxes(new TH1D(
    ("LargestEnergy " + postfix).c_str(), "Largest Energy", 
    500, 0, 200), "Largest Energy", "counts");
*/

if(MC) {
  double uRange = 1.0;

  // Creating the histogram for timeWindowMC with energy ranges
   /* getStatistics().createHistogramWithAxes(new TH2D(
    ("smallest_vs_largest " + postfix).c_str(), "EnergySum-Smallest_vs_EnergySum-largest", 
    100, 0, uRange, 100, 0, uRange), "smallestE [MeV]", "largestE [MeV]");

    getStatistics().createHistogramWithAxes(new TH2D(
    ("EnergySum-Smallest_vs_EnergySum-largest "+ postfix).c_str(),  "EnergySum-Smallest_vs_EnergySum-largest", 
    100, 0, 1.2, 100, 0, 1.2), "Energy_sum -smallestE [MeV]", "Energy_sum -largestE [MeV]"); */

    getStatistics().createHistogramWithAxes(new TH1D(
    ("EnergySum " + postfix).c_str(),  "Sum of energies of hits",
    100, 0, 2), "EnergySum [MeV]", "counts");

    getStatistics().createHistogramWithAxes(new TH1D(
    ("TOT " + postfix).c_str(), "Energy of hits", 
    100, 0, 1.5), "Energy [MeV]", "Counts");

} else if(Data) {
  /*getStatistics().createHistogramWithAxes(new TH2D(
    ("smallest_vs_largest" + postfix).c_str(), 
    "SmallestTOT_vs_LargestTOT", 
    310, -10, 300, 310, -10, 300), "sTOT [ns]", "lTOT [ns]");

    getStatistics().createHistogramWithAxes(new TH2D(
    ("EnergySum-Smallest_vs_EnergySum-largest "+ postfix).c_str(),  "TOTSum-STOT_vs_TOTSum-LTOT", 
    600, 0, 600, 600, 0, 600), "TOTSum-STOT [ns]", "TOTSum-LTOT [ns]");*/

    getStatistics().createHistogramWithAxes(new TH1D(
    ("EnergySum " + postfix).c_str(),  "Sum of TOT of hits",
    600, 0, 600), "TOTSum [ns]", "counts");

    getStatistics().createHistogramWithAxes(new TH1D(
    ("TOT " + postfix).c_str(), "TOT gamma", 
    201, -1.0, 200.0), "Time over Threshold [ns]", "Counts");

}

/*    
    getStatistics().createHistogramWithAxes(new TH1D(
    ("EnergySum-Smallest " + postfix).c_str(),"EnergySum-Smallest",
      600, 0, 200), "EnergySum-Smallest [ns]", "counts");

    getStatistics().createHistogramWithAxes(new TH1D(
    ("EnergySum-largest " + postfix).c_str(), "EnergySum-largest", 500, 0, 200),
    "EnergySum-largest [ns]", "counts");
*/
  }
  
  const std::string fix2[9] = {"noCut", "zCut", "totCut", "SCCut",  "ETSCut", "ttCut", "STCutRej", "STCut", "CombCut"};

  for (int i = 0; i < 9; i++)
  {
    const std::string pfix = fix2[i];

    /*getStatistics().createHistogramWithAxes(new TH1D( 
    ("timeDiff " + pfix).c_str(), "hit time difference", 
    201,  -100, fMaxTimeDiff), "Time difference [ps]", "Counts");
*/
    getStatistics().createHistogramWithAxes(new TH1D(
    ("TOF " + pfix).c_str(),"Time of flight",
    200, -5*fScatterTOFTimeDiff-0.5, 5*fScatterTOFTimeDiff-0.5 ), "TOF", "counts");

    getStatistics().createHistogramWithAxes(new TH1D(
    ("Distance " + pfix).c_str(), "Distance between hits", 
    150, 0, 150), "Distance between hits [cm]", "counts");

    getStatistics().createHistogramWithAxes(new TH1D(
    ("openAngles " + pfix).c_str(), "open Angles", 
    180, 0, 180 ), "#theta Open Angle", "counts");
    
    //2d PLots
    /*getStatistics().createHistogramWithAxes(new TH2D(
      ("timeDiff_vs_angle " + pfix).c_str(), "time difference vs angle between 2 hits", 
      201,  -fMaxTimeDiff, 2*fMaxTimeDiff, 180, 0, 180 ), "Time difference [ps]", "#theta Open Angle");
    
    getStatistics().createHistogramWithAxes(new TH2D(
    ("Distance_vs_angle " + pfix).c_str(), "Distance vs angle between 2 hits", 
    150, 0, 150, 180, 0, 180 ),"Distance [cm]", "#theta Open Angle");
    
    getStatistics().createHistogramWithAxes(new TH2D(
    ("timeDiff_vs_Distance " + pfix).c_str(), "time difference vs Distance between 2 hits", 
    201,  -fMaxTimeDiff, 2*fMaxTimeDiff, 150, 0, 150 ), "Time difference [ps]", "Distance [cm]"); 
    */
    //scattering plots just for inspection

    getStatistics().createHistogramWithAxes(new TH1D(
    ("M_ij " + pfix).c_str(), "Difference in hit time and time calculated by d/c(minimum)",
    2.0*fScatterTOFTimeDiff, -10*fScatterTOFTimeDiff-0.5, 10*fScatterTOFTimeDiff-0.5 ), "M_ij [ps]", "counts");


   /* getStatistics().createHistogramWithAxes(new TH2D(
    ("MDist_vs_angle " + pfix).c_str(), "Distance vs angle between 2 hits(Mij)", 
    500, 0, 500, 180, 0, 180 ), "Distance [cm]", "#theta Open Angle");
*/
    getStatistics().createHistogramWithAxes(new TH2D(
    ("Mdt_vs_angle " + pfix).c_str(), "time difference vs angle between 2 hits(Mij)", 
    2.0*fScatterTOFTimeDiff, -10*fScatterTOFTimeDiff-0.5, 10*fScatterTOFTimeDiff-0.5, 180, 0, 180 ), "Scatter test [ps]", "scattering angle"); 
/*
    getStatistics().createHistogramWithAxes(new TH1D(
    ("M_ij_minimum " + pfix).c_str(), "MInimum of difference between hit time and time calculated by d/c(minimum)",
    201, -fScatterTOFTimeDiff-0.5, fScatterTOFTimeDiff-0.5 ), "M_ij_min", "counts");

    getStatistics().createHistogramWithAxes(new TH1D(
    ("M_ij_maximum " + pfix).c_str(), "Maximum of difference between hit time and time calculated by d/c(minimum)",
    201, -3.5*fScatterTOFTimeDiff-0.5, fScatterTOFTimeDiff-0.5 ), "M_ij_max", "counts");                                       
    */
    //getStatistics().createHistogramWithAxes(new TH2D(
    //("TOF_vs_distance " + pfix).c_str(), "TOF_vs_distance", 
    //2.0*fScatterTOFTimeDiff, -10*fScatterTOFTimeDiff-0.5, 10*fScatterTOFTimeDiff-0.5, 150, 0, 150), "TOF[ps]", "Distance[cm]");

    //getStatistics().createHistogramWithAxes(new TH2D(
    //("TOF_vs_timeDifference " + pfix).c_str(), "TOF_vs_timeDifference",
    //2.0*fScatterTOFTimeDiff, -10*fScatterTOFTimeDiff, 10*fScatterTOFTimeDiff-0.5, 2*fScatterTOFTimeDiff, -10*fScatterTOFTimeDiff, 10*fScatterTOFTimeDiff-0.5), "TOF[ps]", "Time Difference[ps]" );

    getStatistics().createHistogramWithAxes(new TH2D(
    ("3g_rel_angles " + pfix).c_str(), "Sum vs. difference of two smallest relative angles in 3 gamma event", 
    250, 0.0, 250, 200, 0.0, 200.0), "ang1+ang2 [deg]", "ang2-ang1 [deg]");

    getStatistics().createHistogramWithAxes(new TH2D(
    ("relAngles12_23 " + pfix).c_str(), "relative angles between 2 hits", 
    250, 0.0, 250, 200, 0.0, 200.0), "#theta_12 [deg]", "#theta_23 [deg]");

  }



// HIstograms for general and initial cuts

  getStatistics().createHistogramWithAxes( new TH1D(
  "True Energy", "Energy of hits(true)",  
  750, 0, 1200), "Energy(keV)", "counts");

  getStatistics().createHistogramWithAxes( new TH1D(
  "Energy_all", "Energy of hits(reco)",  
  750, 0, 1200), "Energy(keV)", "counts");

  
/*
  getStatistics().createHistogramWithAxes(new TH1D(
  "Gengamma_multi_all", "Generated gamma multiplicity_all",  
   1000, 0.5, 10000.5), "No. of hits","counts");

*/

  // Histograms for 2Gamma category
 /* getStatistics().createHistogramWithAxes(
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
*/
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


}
