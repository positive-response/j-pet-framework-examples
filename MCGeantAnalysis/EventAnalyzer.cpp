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
 *  @file EventAnalyzer.cpp
 */

#include <JPetMCHit/JPetMCHit.h>
#include <JPetOptionsTools/JPetOptionsTools.h>
#include "EventAnalyzer.h"
#include "../LargeBarrelAnalysis/EventCategorizerTools.h"

using namespace jpet_options_tools;
using namespace std;

const double EventAnalyzer::kUnknownEventType = 66666666;

EventAnalyzer::EventAnalyzer(const char *name) : JPetUserTask(name) {}

EventAnalyzer::~EventAnalyzer() {}

bool EventAnalyzer::init() {
  INFO("Event analysis started.");

  pTree22 = new TTree("T","flat tree");
  //pTree22->Branch("position", &fPos);
  pTree22->Branch("x",&fPosX);
  pTree22->Branch("y",&fPosY);
  pTree22->Branch("z",&fPosZ);
  pTree22->Branch("energy",&fEnergy);
  pTree22->Branch("time",&fTime);
  pTree22->Branch("hitType",&fHitType);
 // pTree22->Branch("vtxIndex",&fVtxIndex);

  /*getStatistics().createHistogramWithAxes(
      new TH1D("z_res", "Resolution along Z", 301, -15.05, 15.05),
      "Z_{REC}-Z_{MC} [cm]");

  getStatistics().createHistogramWithAxes(
      new TH1D("Edep_res", "Resolution of deposited energy", 201, -201., 201.),
      "E_{REC}-E_{MC} [keV]");*/

  getStatistics().createHistogramWithAxes( new TH1D(
  "Energy_all", "Energy of hits(reco)",  
  750, 0, 1200), "Energy(keV)", "counts");

  getStatistics().createHistogramWithAxes(new TH1D(
    "HM_reco ", "Number of hits", 
    40, -0.5, 39.5), "Number of hits", "Counts");

      getStatistics().createHistogramWithAxes(new TH1D(
    "HM_true ", "Number of hits", 
    20, -0.5, 19.5), "Number of hits", "Counts");
   

     getStatistics().createHistogramWithAxes(new TH1D(
    "getGenGammaMultiplicity", "gen gamma multiplicity", 
    10, -0.5, 9.5), "Number of hits", "Counts");

  getStatistics().createHistogramWithAxes(new TH1D(
   "Phi ", "Phi",
   100, -200, 200),  "#phi Angle", "counts");

    getStatistics().createHistogramWithAxes(new TH1D(
   "Theta ", "SumTheta",
    200, 0, 200 ), "#theta Angle", "counts");

    getStatistics().createHistogramWithAxes(new TH1D(
   "SumTheta ", "SumTheta",
    360, 0, 720 ), "#theta Angle", "counts");
    

    getStatistics().createHistogramWithAxes(new TH2D(
    "Theta_vs_phi", "Theta_vs_phi",
     200, 0, 200 , 100, -200, 200), "#theta Angle", "#phi Angle");

    getStatistics().createHistogramWithAxes(new TH2D(
    "Theta_vs_energy", "Theta_vs_energy",
     360, 0, 360 , 500, 0, 500), "#theta Angle", "Energy [keV]");

    getStatistics().createHistogramWithAxes(new TH2D(
    "phi_vs_energy" , "phi_vs_energy",
     100, -200, 200 , 500, 0, 500), "#theta Angle", "Energy [keV]");

  fOutputEvents = new JPetTimeWindow("JPetEvent");

  return true;
}

bool EventAnalyzer::exec() {

  // Identify whether the input events are MC or DATA.
  // In case of MC, store the pointer to the TimeWindowMC object
  // which contains "true MC" information about the generated events.
  JPetTimeWindowMC *time_window_mc = nullptr;
  if (time_window_mc = dynamic_cast<JPetTimeWindowMC *const>(fEvent)) {
    fIsMC = true;
    INFO("The input file is MC.");
  } else {
    INFO("The input file is DATA.");
  }

 bool trueEvents = false;
 JPetEvent sel1Event;
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow *const>(fEvent)) {

    for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++) {
      const auto &event =
          dynamic_cast<const JPetEvent &>(timeWindow->operator[](i));
        
      getStatistics().fillHistogram("HM_reco ", event.getHits().size());
      // if the input is MC, we fill resolution histograms
      // to check if MC smearing works fine
      if (fIsMC) {

        if(event.getHits().size()==3){
          for(const auto & hit : event.getHits())
                {
                  getStatistics().fillHistogram("Energy_all", hit.getEnergy());
                  fPosX = hit.getPosX();
                  fPosY = hit.getPosY();
                  fPosZ = hit.getPosZ();
                  fEnergy = hit.getEnergy(); 
                  fTime = hit.getTime();
                  if(time_window_mc)
                  {
                    auto mcHit = time_window_mc->getMCHit<JPetMCHit>(hit.getMCindex());
                    fHitType = mcHit.getGenGammaMultiplicity();
                    //fVtxIndex = mcHit.getMCVtxIndex();			
                  } else {
                          fHitType = kUnknownEventType;	
                          //fVtxIndex= 0;
                          }  
                          pTree22->Fill();
                }

          //bool EventSvd = saveEvt(event, time_window_mc);
          //pTree22->Fill();
          /*fEnergy.clear();
          fTime.clear();
          fPosX.clear();
          fPosY.clear();
          fPosZ.clear();
          fHitType.clear();
          fVtxIndex.clear();*/

        }

        fillResolutionHistograms(event, time_window_mc);
        trueEvents = hitParser(event, time_window_mc, sel1Event);
        
      }
    }

if(trueEvents){

  

  //Twohit(sel1Event);
  Fourhit(sel1Event);
    getStatistics().fillHistogram("HM_true ", sel1Event.getHits().size());
    double sumTheta = 0.0;
    for(int i = 0; i < sel1Event.getHits().size(); i++)
    {
      auto rhit1 = sel1Event.getHits().at(i);
      auto hit1 = time_window_mc->getMCHit<JPetMCHit>(rhit1.getMCindex());
     sumTheta += hit1.getPos().Theta();
    
     // getStatistics().fillHistogram("Energy_all", hit1.getEnergy());
      getStatistics().fillHistogram("Phi ",TMath::RadToDeg() * hit1.getPos().Phi());
      getStatistics().fillHistogram("Theta ", TMath::RadToDeg() * hit1.getPos().Theta());
      getStatistics().fillHistogram("Theta_vs_phi", TMath::RadToDeg() * hit1.getPos().Theta(), TMath::RadToDeg() * hit1.getPos().Phi());
    //  getStatistics().fillHistogram("Theta_vs_energy", TMath::RadToDeg() * hit1.getPos().Theta(), hit1.getEnergy());
      getStatistics().fillHistogram("phi_vs_energy", TMath::RadToDeg() * hit1.getPos().Phi(), hit1.getEnergy()); 
     
    }
   // getStatistics().fillHistogram("SumTheta ", TMath::RadToDeg() * sumTheta);
    sumTheta = 0.0;
}

  } else {
    return false;
  }

  return true;
}

bool EventAnalyzer::terminate() {
  INFO("Event analysis completed.");
  if(!pTree22) std::cout << "wtf??" << std::endl;
  file->cd();
  pTree22->Write();
  //file->Write("",TObject::kWriteDelete);
  file->Close();
  return true;
}

void EventAnalyzer::fillResolutionHistograms(const JPetEvent &event,
                                             const JPetTimeWindowMC *tw) {

  int hits_number = event.getHits().size();
  for (int k = 0; k < hits_number; ++k) {
    const JPetHit &reconstructed_hit = event.getHits().at(k);
    // for each reconstructed hit, we access the corresponding
    // "true MC" hit
    const JPetMCHit &mc_hit =
        tw->getMCHit<JPetMCHit>(reconstructed_hit.getMCindex());

   /* getStatistics().fillHistogram("z_res", reconstructed_hit.getPos().Z() -
                                               mc_hit.getPos().Z());
    getStatistics().fillHistogram("Edep_res", reconstructed_hit.getEnergy() -
                                                  mc_hit.getEnergy());*/
  }
}

bool EventAnalyzer::hitParser(const JPetEvent &event,const JPetTimeWindowMC *tw, JPetEvent& sel1Event){
  sel1Event = JPetEvent();
for(int i = 0; i < event.getHits().size(); i++){
  const JPetHit& reco_hit = event.getHits().at(i);
  const JPetMCHit& mc_hit = tw->getMCHit<JPetMCHit>(reco_hit.getMCindex());
  if(mc_hit.getGenGammaMultiplicity() == 3){
    getStatistics().fillHistogram("getGenGammaMultiplicity", mc_hit.getGenGammaMultiplicity());
    sel1Event.addHit(reco_hit);
   // sel1Event.addHit(mc_hit);

  }
}
if(sel1Event.getHits().size() != 0) return true;
 }


 double EventAnalyzer::calcAngles(const JPetHit& hit1, const JPetHit& hit2)
{
   return (TMath::RadToDeg() * hit1.getPos().Angle(hit2.getPos()));
}

bool EventAnalyzer::Twohit(const JPetEvent &event){
  if(event.getHits().size() < 2) return false;
  double sumTheta = 0.0;
for(int i = 0; i < event.getHits().size(); i++){
  for(int j = i; j < event.getHits().size(); j++){

      auto rhit1 = event.getHits().at(i);
      auto rhit2 = event.getHits().at(j);

      double openAng = calcAngles(rhit1, rhit2);
      sumTheta += openAng;
      getStatistics().fillHistogram("Theta_vs_energy", openAng, rhit1.getEnergy());
      getStatistics().fillHistogram("Theta_vs_energy", openAng, rhit2.getEnergy());
      
  }
}
if(sumTheta > 0.0){
getStatistics().fillHistogram("SumTheta ", sumTheta);
sumTheta = 0.0;
}
return true;
 }
 
 bool EventAnalyzer::Fourhit(const JPetEvent &event){
  if(event.getHits().size() < 4) return false;
  double sumTheta = 0.0;
for(int i = 0; i < event.getHits().size(); i++){
  for(int j = i; j < event.getHits().size(); j++){

      auto rhit1 = event.getHits().at(i);
      auto rhit2 = event.getHits().at(j);

      double openAng = calcAngles(rhit1, rhit2);
      sumTheta += openAng;

      getStatistics().fillHistogram("Theta_vs_energy", openAng, rhit1.getEnergy());
      getStatistics().fillHistogram("Theta_vs_energy", openAng, rhit2.getEnergy());
  }
}
  if(sumTheta > 0.0){
  getStatistics().fillHistogram("SumTheta ", sumTheta);
  sumTheta = 0.0;
  }
return true;
 }
/*
 bool EventAnalyzer::saveEvt(const JPetEvent &event, const JPetTimeWindowMC *tw){
  
  fEnergy.clear();
  fTime.clear();
  fPosX.clear();
  fPosY.clear();
  fPosZ.clear();
  fHitType.clear();
  fVtxIndex.clear();

  for(const auto & hit : event.getHits())
  {
    fPosX.push_back(hit.getPosX());
    fPosY.push_back(hit.getPosY());
    fPosZ.push_back(hit.getPosZ());
    fEnergy.push_back(hit.getEnergy()); 
    fTime.push_back(hit.getTime());
    if(tw)
    {
      auto mcHit = tw->getMCHit<JPetMCHit>(hit.getMCindex());
      fHitType.push_back(mcHit.getGenGammaMultiplicity());
      fVtxIndex.push_back(mcHit.getMCVtxIndex());			
    } else {
            fHitType.push_back(kUnknownEventType);	
            fVtxIndex.push_back(0);
            }  
  }
    
 }*/
