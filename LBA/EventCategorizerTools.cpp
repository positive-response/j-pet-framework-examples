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
 *  @file EventCategorizerTools.cpp
 */

#include "EventCategorizerTools.h"
#include "HitFinderTools.h"
#include <TMath.h>
#include <vector>

using namespace std;

bool EventCategorizerTools::Plots3hits(const JPetEvent& event, JPetStatistics& stats, const std::string& postfix){
  if(event.getHits().size() < 3) return false;
for (uint i = 0; i < event.getHits().size(); i++)
    {
      for (uint j = i + 1; j < event.getHits().size(); j++)
      {
        for (uint k = j + 1; k < event.getHits().size(); k++)
        {
          auto firstHit = event.getHits().at(i);
	        auto secondHit = event.getHits().at(j);
	        auto thirdHit = event.getHits().at(k);

          std::vector<double> relativeAngles = {calcAngles(firstHit, secondHit), calcAngles(secondHit, thirdHit), calcAngles(thirdHit, firstHit)};

        //  stats.fillHistogram(("relAngles12_23 " + postfix).c_str(), fabs(calcAngles(firstHit, secondHit)), fabs(calcAngles(secondHit, thirdHit)));  
          std::sort(relativeAngles.begin(), relativeAngles.end());

          double transformedX = relativeAngles.at(1) + relativeAngles.at(0);
          double transformedY = relativeAngles.at(1) - relativeAngles.at(0);

          std::vector<double> eTime = {calcET(firstHit), calcET(secondHit), calcET(thirdHit) };
          std::sort(eTime.begin(), eTime.end());
          double fTimeSpread = eTime.at(2) - eTime.at(0);

          double planeC = EventCategorizerTools::calculatePlaneCenterDistance(firstHit, secondHit, thirdHit);
          stats.fillHistogram(("planeCentre " + postfix).c_str(), planeC);
          stats.fillHistogram(("ETS " + postfix).c_str(), fTimeSpread);
          stats.fillHistogram(("3g_rel_angles " + postfix).c_str(), transformedX, transformedY);
          stats.fillHistogram(("ETS_vs_DOP " + postfix).c_str(), fTimeSpread/750, planeC/2.0);

        }
      }
    }
  return true;
}

bool EventCategorizerTools::Plots2hits(const JPetEvent& event, JPetStatistics& stats, const std::string& postfix){
  for(int i = 0; i < event.getHits().size(); i++){
      auto firstHit = event.getHits().at(i);
    for(int j = i+1; j < event.getHits().size(); j++){
      auto secondHit = event.getHits().at(j);

      TVector3 ap = calculateAnnihilationPoint(firstHit, secondHit);
      double fHitTimeDiff = (secondHit.getTime() - firstHit.getTime() );
      double fTOF = calculateTOFByConvention(firstHit, secondHit);
      double fDistance = calculateDistance(firstHit, secondHit);
      double openAng = calcAngles(firstHit, secondHit);
      double fTimeDiff = fTOF - fDistance/kLightVelocity_cm_ps;
      double scatterAngle = calculateScatteringAngle(firstHit, secondHit);
     

     // stats.fillHistogram(("timeDiff " + postfix).c_str(), fHitTimeDiff);
    //  stats.fillHistogram(("TOF " + postfix).c_str(), fTOF); 
     // stats.fillHistogram(("Distance " + postfix).c_str(), fDistance);
     // stats.fillHistogram(("openAngles " + postfix).c_str(), openAng);

	   // stats.fillHistogram(("timeDiff_vs_angle " + postfix).c_str(), fHitTimeDiff, openAng);
     // stats.fillHistogram(("Distance_vs_angle " + postfix).c_str(), fDistance, openAng);
      //stats.fillHistogram(("timeDiff_vs_Distance " + postfix).c_str(), fHitTimeDiff, fDistance);
      
      //for checking scattering
      stats.fillHistogram(("M_ij " + postfix).c_str(), fTimeDiff);
    //  stats.fillHistogram(("Mdt_vs_angle " + postfix).c_str(), fTimeDiff, scatterAngle);

     // stats.fillHistogram(("ap_XY " + postfix).c_str(), ap.X(), ap.Y() );
      //stats.fillHistogram(("ap_XZ " + postfix).c_str(), ap.X(), ap.Z() );
      //stats.fillHistogram(("ap_YZ " + postfix).c_str(), ap.Y(), ap.Z() );
     // stats.fillHistogram(("TOF_vs_distance " + postfix).c_str(), fTOF, fDistance);
     // stats.fillHistogram(("TOF_vs_timeDifference " + postfix).c_str(), fTOF,fTimeDiff);   
      
    }   
  }
return true;
}

bool EventCategorizerTools::Plots1hit(const JPetEvent& event, JPetStatistics& stats, const std::string& postfix)
{
  double fSumEnergy = 0.0;
  double smallestEnergy = std::numeric_limits<double>::max();
  double largestEnergy = std::numeric_limits<double>::lowest();

stats.fillHistogram( ("HM " + postfix).c_str(), event.getHits().size());
  for(int i = 0; i < event.getHits().size(); i++)
  {
    auto hit = event.getHits().at(i);
    auto energy =  hit.getEnergy()/1000.0;

    int scinID = hit.getScintillator().getID();

    if(energy < smallestEnergy){smallestEnergy = energy;}
	  if(energy > largestEnergy){largestEnergy = energy; }
    
    fSumEnergy += energy;

    stats.fillHistogram(("TOT " + postfix).c_str(), energy);
    stats.fillHistogram(("Z " + postfix).c_str(), hit.getPos().Z());
   // stats.fillHistogram(("Phi " + postfix).c_str(),TMath::RadToDeg() * hit->getPos().Phi());
   // stats.fillHistogram(("Theta " + postfix).c_str(), TMath::RadToDeg() * hit->getPos().Theta());
    //stats.fillHistogram(("hit_XY " + postfix).c_str(), hit->getPos().X(), hit->getPos().Y() );
    //stats.fillHistogram(("scinID_VS_TOT " + postfix).c_str(), scinID, energy);
  }
   // stats.fillHistogram(("SmallestEnergy " + postfix).c_str(), smallestEnergy/1000);
   // stats.fillHistogram(("LargestEnergy " + postfix).c_str(), largestEnergy/1000);
  if(fSumEnergy != 0){
    //stats.fillHistogram(("smallest_vs_largest " + postfix).c_str(), smallestEnergy, largestEnergy);
   // stats.fillHistogram(("EnergySum " + postfix).c_str(), fSumEnergy);
    //stats.fillHistogram(("EnergySum-Smallest_vs_EnergySum-largest "+ postfix).c_str(), (fSumEnergy - smallestEnergy), (fSumEnergy - largestEnergy));
  }
    fSumEnergy = 0.0;
    smallestEnergy = std::numeric_limits<double>::max();
    largestEnergy = std::numeric_limits<double>::lowest();
  return true;
}

/****************************************************** */
std::pair<bool, int> EventCategorizerTools::ETSTEST(const JPetEvent& event, JPetEvent& Event1, JPetEvent& Event2, JPetEvent& Event3){
  if(event.getHits().size() < 3) return std::make_pair(false,0);

  double sdop = std::numeric_limits<double>::max();
  double sets = std::numeric_limits<double>::max();
  double scomb = std::numeric_limits<double>::max();

  JPetEvent DOPEvent;
  JPetEvent etsEvent;
  JPetEvent cobEvent;

  for (uint i = 0; i < event.getHits().size(); i++)
    {
      for (uint j = i + 1; j < event.getHits().size(); j++)
	    {
	      for (uint k = j + 1; k < event.getHits().size(); k++)
        {
          auto firstHit = event.getHits().at(i);
          auto secondHit = event.getHits().at(j);
          auto thirdHit = event.getHits().at(k);

          double planeC = EventCategorizerTools::calculatePlaneCenterDistance(firstHit, secondHit, thirdHit);
          if(planeC < sdop){
            DOPEvent = JPetEvent();
            sdop = planeC;
            DOPEvent.addHit(firstHit);
            DOPEvent.addHit(secondHit);
            DOPEvent.addHit(thirdHit);
          }

          std::vector<double> eTime = {calcET(firstHit), calcET(secondHit), calcET(thirdHit) };
          std::sort(eTime.begin(), eTime.end());
          double fTimeSpread = eTime.at(2) - eTime.at(0);

          if(fTimeSpread < sets){
            etsEvent = JPetEvent();
            sets = fTimeSpread;
            etsEvent.addHit(firstHit);
            etsEvent.addHit(secondHit);
            etsEvent.addHit(thirdHit);
          }

          double combVal = TMath::Sqrt((fTimeSpread * fTimeSpread) + (planeC * planeC)); 
          if(combVal < scomb){
            cobEvent = JPetEvent();
            scomb = combVal;

            cobEvent.addHit(firstHit);
            cobEvent.addHit(secondHit);
            cobEvent.addHit(thirdHit);
          }
        }
	    }
    }
    Event1 = DOPEvent;
    Event2 = etsEvent;
    Event3 = cobEvent;

    return std::make_pair(true, event.getHits().size());
}
bool EventCategorizerTools::zTest(const JPetEvent& event, JPetEvent& zEvent, double& fZcut){
 bool isZValid = false;
  zEvent =  JPetEvent();

  for(auto & hit: event.getHits())
  {
    if(abs(hit.getPos().Z()) <= fZcut)
    {
      zEvent.addHit(hit);
    }
    else
      continue;
  }
  
  if(zEvent.getHits().size() > 0)
  {
    isZValid = true;
  }
  return isZValid;
}

bool EventCategorizerTools::totCheck(const JPetEvent& event, JPetEvent& totEvent, double totCutMin, double totCutMax){
 totEvent = JPetEvent();

 for(auto & hit: event.getHits()){
     auto tot =  hit.getEnergy()/1000.0;
     if(tot >= totCutMin && tot <= totCutMax ){ 
      totEvent.addHit(hit);
      }
   }
   if(totEvent.getHits().size() != 0) return true;
 return false;
}

bool EventCategorizerTools::ScIDCheck(const JPetEvent& event, JPetEvent& scEvent)
{
  std::set<int> addedHits;
  scEvent = JPetEvent();

  for(int i = 0; i < event.getHits().size(); i++)
  {
    for(int j = i+1; j < event.getHits().size(); j++)
    {
        auto hit1 = event.getHits().at(i);
        auto hit2 = event.getHits().at(j);

      auto scid1 = hit1.getScintillator().getID();
      auto scid2 = hit2.getScintillator().getID();

      if(scid1 == scid2){
        if(hit1.getTime() < hit2.getTime()) 
          {
            if (addedHits.find(scid1) == addedHits.end()) 
              {
                scEvent.addHit(hit1);
                addedHits.insert(scid1);            
              }
          }else if(hit1.getTime() > hit2.getTime()) 
          {
            if (addedHits.find(scid2) == addedHits.end()) 
            {
            scEvent.addHit(hit2);
            addedHits.insert(scid2); 
            }
          }
      }
      else{
        if (addedHits.find(scid1) == addedHits.end()) 
        {
          scEvent.addHit(hit1);
          addedHits.insert(scid1);
        }
         if (addedHits.find(scid2) == addedHits.end())
          {
          scEvent.addHit(hit2);
          addedHits.insert(scid2);
        }
      }
      }
    }
    return true;
  }

bool EventCategorizerTools::coplanarCheck(const JPetEvent& event, JPetEvent& CoPEvent, double &fDOP){
if(event.getHits().size() < 3) return false;
CoPEvent = JPetEvent();
double planeDist = std::numeric_limits<double>::max();
JPetEvent CEvent;

for (uint i = 0; i < event.getHits().size(); i++)
    {
      for (uint j = i + 1; j < event.getHits().size(); j++)
      {
        for (uint k = j + 1; k < event.getHits().size(); k++)
        {
          auto firstHit = event.getHits().at(i);
          auto secondHit = event.getHits().at(j);
          auto thirdHit = event.getHits().at(k);

          double planeC = EventCategorizerTools::calculatePlaneCenterDistance(firstHit, secondHit, thirdHit);

          if(abs(planeC) <= fDOP )
          {
            if(planeC < planeDist)
            {
              CEvent = JPetEvent();
              planeDist = planeC;
              CEvent.addHit(firstHit);
              CEvent.addHit(secondHit);
              CEvent.addHit(thirdHit);
            }
          }
        }
      }
    }
    if(CEvent.getHits().size() > 0) {
    CoPEvent = CEvent;
    return true;
  }
  return false;
}

bool EventCategorizerTools::emissionTimeSpreadCheck(const JPetEvent& event, JPetEvent& etsEvent, double &fETS){
if(event.getHits().size() < 3) return false;
etsEvent = JPetEvent();

double sETS = std::numeric_limits<double>::max();
JPetEvent eEvent; 

for (uint i = 0; i < event.getHits().size(); i++)
    {
      for (uint j = i + 1; j < event.getHits().size(); j++)
      {
        for (uint k = j + 1; k < event.getHits().size(); k++)
        {
          auto firstHit = event.getHits().at(i);
	        auto secondHit = event.getHits().at(j);
	        auto thirdHit = event.getHits().at(k);

          std::vector<double> eTime = {calcET(firstHit), calcET(secondHit), calcET(thirdHit) };
          std::sort(eTime.begin(), eTime.end());
          double fTimeSpread = eTime.at(2) - eTime.at(0);

          if( fTimeSpread <= fETS ){
            if(fTimeSpread < sETS)
            {
              eEvent = JPetEvent();
              sETS = fTimeSpread;
              eEvent.addHit(firstHit);
              eEvent.addHit(secondHit);
              eEvent.addHit(thirdHit);
            } 
          }
        }
      }
    }

   if(eEvent.getHits().size() > 0) {
    etsEvent = eEvent;
    return true;
  }
  return false;
}

bool EventCategorizerTools::combinedCheck(const JPetEvent& event, JPetEvent& cmbEvent, double &val){
if(event.getHits().size() < 3) return false;
cmbEvent = JPetEvent();
JPetEvent eEvent; 

for (uint i = 0; i < event.getHits().size(); i++)
    {
      for (uint j = i + 1; j < event.getHits().size(); j++)
      {
        for (uint k = j + 1; k < event.getHits().size(); k++)
        {
          auto firstHit = event.getHits().at(i);
          auto secondHit = event.getHits().at(j);
          auto thirdHit = event.getHits().at(k);

          std::vector<double> eTime = {calcET(firstHit), calcET(secondHit), calcET(thirdHit) };
          std::sort(eTime.begin(), eTime.end());
          double fTimeSpread = eTime.at(2) - eTime.at(0);

          double planeC = EventCategorizerTools::calculatePlaneCenterDistance(firstHit, secondHit, thirdHit);
	  double valToCheck = TMath::Sqrt((fTimeSpread * fTimeSpread)/(750*750) + (planeC * planeC)/(2.0 * 2.0));
            
            if( valToCheck  <= val)
            {
              eEvent = JPetEvent();
              eEvent.addHit(firstHit);
              eEvent.addHit(secondHit);
              eEvent.addHit(thirdHit); 
          }
        }
      }
    }
  if(eEvent.getHits().size() > 0) {
	cmbEvent = eEvent;    
	return true;
  }
  return false;
}

bool EventCategorizerTools::checkForScatter(const JPetEvent& event, JPetEvent& selEvent, double cutValue){

if(event.getHits().size() < 3) return false;
int c = 0;

int size = event.getHits().size();

  for(int i = 0; i < size; i++){
    for(int j = i+1; j < size; j++){
      auto hit1 = event.getHits().at(i);
      auto hit2	= event.getHits().at(j);

      double fTOF = EventCategorizerTools::calculateTOFByConvention(hit1, hit2);
      double fDistance = EventCategorizerTools::calculateDistance(hit1, hit2);
      double scatterAngle = EventCategorizerTools::calculateScatteringAngle(hit1, hit2);
      double fTimeDiff = fTOF - fDistance/kLightVelocity_cm_ps;

      if(fTimeDiff > cutValue){
        c++;
      }
   }
}
  if(c == 0) {    //means no scattering 
    return true;
  }
  return false;
}

bool EventCategorizerTools::thetaCheck0(const JPetEvent& event, JPetEvent& ttEvent, double &fTT0){
if(event.getHits().size() < 3) return false;
ttEvent = JPetEvent();

for (uint i = 0; i < event.getHits().size(); i++)
    {
      for (uint j = i + 1; j < event.getHits().size(); j++)
      {
        for (uint k = j + 1; k < event.getHits().size(); k++)
        {
          auto firstHit = event.getHits().at(i);
	        auto secondHit = event.getHits().at(j);
	        auto thirdHit = event.getHits().at(k);

          std::vector<double> relativeAngles = {calcAngles(firstHit, secondHit), calcAngles(secondHit, thirdHit), calcAngles(thirdHit, firstHit)};  
          std::sort(relativeAngles.begin(), relativeAngles.end());
          double transformedX = relativeAngles.at(1) + relativeAngles.at(0);
          if( transformedX > fTT0){
            ttEvent.addHit(firstHit);
            ttEvent.addHit(secondHit);
            ttEvent.addHit(thirdHit);
            return true;
          }
        }
      }
    }
  return false;
}

bool EventCategorizerTools::thetaCheck1(const JPetEvent& event, JPetEvent& ttEvent, double &fTT1){
if(event.getHits().size() < 3) return false;
ttEvent = JPetEvent();

for (uint i = 0; i < event.getHits().size(); i++)
    {
      for (uint j = i + 1; j < event.getHits().size(); j++)
      {
        for (uint k = j + 1; k < event.getHits().size(); k++)
        {
          auto firstHit = event.getHits().at(i);
	        auto secondHit = event.getHits().at(j);
	        auto thirdHit = event.getHits().at(k);

          std::vector<double> relativeAngles = {calcAngles(firstHit, secondHit), calcAngles(secondHit, thirdHit), calcAngles(thirdHit, firstHit)};  
          std::sort(relativeAngles.begin(), relativeAngles.end());

          double transformedY = relativeAngles.at(1) - relativeAngles.at(0);
          if( transformedY < fTT1){

            ttEvent.addHit(firstHit);
            ttEvent.addHit(secondHit);
            ttEvent.addHit(thirdHit);
            return true;
          }
        }
      }
    }
  return false;
}

bool EventCategorizerTools::checkForCut(const JPetEvent& event, JPetEvent& selEvent, int cutNr, int mp, double cut, 
const std::string& plotPfix, bool (*funcCheck)(const JPetEvent&, JPetEvent&, double&), JPetStatistics& stats )
{ 
  bool checkForCondition = false;
  JPetEvent tempEvent;
  checkForCondition = (*funcCheck)(event, tempEvent, cut);
  
  if (checkForCondition)
  {
    selEvent = tempEvent;
    stats.fillHistogram("cutStatus", cutNr);
    Plots1hit(tempEvent, stats, plotPfix);
    Plots2hits(tempEvent, stats, plotPfix);
    Plots3hits(tempEvent, stats, plotPfix);
    return true;
  }
  return false;
}



double EventCategorizerTools::calcET(const JPetHit& hit){
  TVector3 SourcePos(0.0, 0.0, 0.0);
  double regTime = ((hit.getPos() - SourcePos).Mag())/ kLightVelocity_cm_ps; 
  return (hit.getTime() - regTime); 
}
/*************plotting with the conditions cut */


bool EventCategorizerTools::checkForCut3(const JPetEvent& event, JPetStatistics& stats, JPetEvent& sel3Event, int &m)
{
  bool totCheckDeexHits = false;
  bool totCheckHits = false;

  bool furtherCheck = false;
  int mcity = 0;
  std::pair<bool, int> separation = std::make_pair(false, 0);

  JPetEvent totEventAnn;
  JPetEvent totEventDeex;

  JPetEvent Event1;
  JPetEvent Event2;
  JPetEvent Event3;

 // totCheckDeexHits = totCheck(event, totEventDeex, fToTCutDeexMin, fToTCutDeexMax);
  totCheckHits = totCheck(event, totEventAnn, fToTCutAnniMin, fToTCutAnniMax);
  
  if(totCheckHits && totEventAnn.getHits().size() != 0)
  {
    separation = ETSTEST(totEventAnn, Event1, Event2, Event3);
    furtherCheck = separation.first;
    mcity = separation.second;

    //sel3Event = totEventAnn;
    stats.fillHistogram("cutStatus", 3);
    Plots3hits(totEventAnn, stats, "totCut");
    Plots2hits(totEventAnn, stats, "totCut");
    Plots1hit(totEventAnn, stats, "totCut");
    
     /*if(Event1.getHits().size() != 0){
        
      }
      if(Event2.getHits().size() != 0){
        
      }
*/
    if(Event3.getHits().size() != 0){
        Plots1hit(Event3, stats, "CombSel");
        Plots2hits(Event3, stats,   "CombSel");
        Plots3hits(Event3, stats,  "CombSel");
        sel3Event = Event3;
         m = mcity;
    }
    return true;
  }
return false;
}

bool EventCategorizerTools::checkForCut2(const JPetEvent& event, JPetStatistics& stats, JPetEvent& sel2Event, int &m)
{ 
  bool checkForSC = false;
  JPetEvent scEvent;
  int minHits = 1;
  checkForSC = ScIDCheck(event, scEvent);
  
  if(checkForSC && scEvent.getHits().size() != 0)
  {
    sel2Event = scEvent;
    Plots1hit(scEvent, stats, "SCCut");
    Plots2hits(scEvent, stats, "SCCut");
    Plots3hits(scEvent, stats, "SCCut");
    return true;
  }
  return false;
}

bool EventCategorizerTools::checkForScat(const JPetEvent& event,  JPetStatistics& stats, JPetEvent& selScEvent , int m)
{ 
  bool checkForSc = false;
  JPetEvent ScEvent;
  double scatValue = -1800.0;
  checkForSc = checkForScatter(event, ScEvent, scatValue);
  
  if(checkForSc)
  {
    selScEvent = ScEvent;
    stats.fillHistogram("cutStatus", 7);
    Plots1hit(event, stats, "STCut");
    Plots2hits(event, stats, "STCut");
    Plots3hits(event, stats, "STCut");
    return true;
  }

    if(!checkForSc && event.getHits().size()!= 0)
  {
    selScEvent = ScEvent;
    Plots1hit(event, stats, "STCutRej");
    Plots2hits(event, stats, "STCutRej");
    Plots3hits(event, stats, "STCutRej");
    return true;
  }
  return false;
}

/*bool EventCategorizerTools::checkForCut4(const JPetEvent& event, JPetStatistics& stats, JPetEvent& sel4Event, int m, double val)
{ 
  if(event.getHits().size() < 3) return false;

  bool checkForCoP = false;
  JPetEvent CoPEvent;
  checkForCoP = coplanarCheck(event, stats, CoPEvent, val);
  
  if(checkForCoP && CoPEvent.getHits().size()!= 0)
  {
    sel4Event = CoPEvent;
    stats.fillHistogram("cutStatus", 4);
    Plots1hit(CoPEvent,  stats, "copCut");
    Plots2hits(CoPEvent, stats,"copCut");
    Plots3hits(CoPEvent,  stats, "copCut");
    return true;
  }
  return false;
}

bool EventCategorizerTools::checkForCut5(const JPetEvent& event, JPetStatistics& stats, JPetEvent& sel5Event, int m, double val)
{ 
   if(event.getHits().size() < 3) return false;
  bool checkForETS = false;
  JPetEvent ETSEvent;
  checkForETS = emissionTimeSpreadCheck(event, stats, ETSEvent, val);
  
  if(checkForETS && ETSEvent.getHits().size()!= 0)
  {
    sel5Event = ETSEvent;
    stats.fillHistogram("cutStatus", 5);
    Plots1hit(ETSEvent, stats, "ETSCut");
    Plots2hits(ETSEvent, stats, "ETSCut");
    Plots3hits(ETSEvent, stats, "ETSCut");
    return true;
  }
  return false;
}

bool EventCategorizerTools::checkForCut5A(const JPetEvent& event, JPetEvent& sel5AEvent)
{ 
  bool checkForETS = false;  ///for combined cut 
  JPetEvent ETSEvent;
  checkForETS = combinedCheck(event, ETSEvent);
  
  if(checkForETS && ETSEvent.getHits().size()!= 0)
  {
    sel5AEvent = ETSEvent;
    return true;
  }
  return false;
}

bool EventCategorizerTools::checkForCut6(const JPetEvent& event, JPetStatistics& stats, JPetEvent& sel6Event, int m, double val)
{ 
  bool checkFortt = false;
  JPetEvent tt0Event;
  checkFortt = thetaCheck0(event, stats, tt0Event, val);
  
  if(checkFortt)
  {
    sel6Event = tt0Event;
    stats.fillHistogram("cutStatus", 6);
    Plots1hit(tt0Event, stats, "tt0Cut");
    Plots2hits(tt0Event, stats, "tt0Cut");
    Plots3hits(tt0Event, stats, "tt0Cut");
    return true;
  }
  return false;
}

bool EventCategorizerTools::checkForCut7(const JPetEvent& event, JPetStatistics& stats, JPetEvent& sel7Event, int m, double val)
{ 
  bool checkFortt = false;
  JPetEvent ttEvent;
  checkFortt = thetaCheck1(event, stats, ttEvent, val);
  
  if(checkFortt)
  {
    sel7Event = ttEvent;
    stats.fillHistogram("cutStatus", 7);
    Plots1hit(ttEvent, stats, "tt1Cut");
    Plots2hits(ttEvent, stats, "tt1Cut");
    Plots3hits(ttEvent, stats, "tt1Cut");
    return true;
  }
  return false;
}
bool EventCategorizerTools::checkForCut1(const JPetEvent& event, JPetStatistics& stats, JPetEvent& sel1Event)
{ 
  bool checkForZ = false;
  JPetEvent zEvent;
  int minHits = 1;
  checkForZ = zTest(event, zEvent, minHits);
  
  if(checkForZ)
  {
    stats.fillHistogram("cutStatus", 2);
    sel1Event = zEvent;
    Plots1hit(zEvent, stats, "zCut");
    Plots2hits(zEvent, stats, "zCut");
    Plots3hits(zEvent, stats, "zCut");
    return true;
  }
  return false;
}
*/
/**
* Method for determining type of event - back to back 2 gamma
*/
bool EventCategorizerTools::checkFor2Gamma(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos,
  double b2bSlotThetaDiff, double b2bTimeDiff
)
{
  if (event.getHits().size() < 2) {
    return false;
  }
  for (uint i = 0; i < event.getHits().size(); i++) {
    for (uint j = i + 1; j < event.getHits().size(); j++) {
      JPetHit firstHit, secondHit;
      if (event.getHits().at(i).getTime() < event.getHits().at(j).getTime()) {
        firstHit = event.getHits().at(i);
        secondHit = event.getHits().at(j);
      } else {
        firstHit = event.getHits().at(j);
        secondHit = event.getHits().at(i);
      }
      // Checking for back to back
      double timeDiff = fabs(firstHit.getTime() - secondHit.getTime());
      double deltaLor = (secondHit.getTime() - firstHit.getTime()) * kLightVelocity_cm_ps / 2.;
      double theta1 = min(firstHit.getBarrelSlot().getTheta(), secondHit.getBarrelSlot().getTheta());
      double theta2 = max(firstHit.getBarrelSlot().getTheta(), secondHit.getBarrelSlot().getTheta());
      double thetaDiff = min(theta2 - theta1, 360.0 - theta2 + theta1);
      /*if (saveHistos) {
        stats.fillHistogram("2Gamma_Zpos", firstHit.getPosZ());
        stats.fillHistogram("2Gamma_Zpos", secondHit.getPosZ());
        stats.fillHistogram("2Gamma_TimeDiff", timeDiff / 1000.0);
        stats.fillHistogram("2Gamma_DLOR", deltaLor);
        stats.fillHistogram("2Gamma_ThetaDiff", thetaDiff);
        stats.fillHistogram("2Gamma_Dist", calculateDistance(firstHit, secondHit));
      }*/
      if (fabs(thetaDiff - 180.0) < b2bSlotThetaDiff && timeDiff < b2bTimeDiff) {
        /*if (saveHistos) {
          TVector3 annhilationPoint = calculateAnnihilationPoint(firstHit, secondHit);
          stats.fillHistogram("Annih_TOF", calculateTOFByConvention(firstHit, secondHit));
          stats.fillHistogram("AnnihPoint_XY", annhilationPoint.X(), annhilationPoint.Y());
          stats.fillHistogram("AnnihPoint_ZX", annhilationPoint.Z(), annhilationPoint.X());
          stats.fillHistogram("AnnihPoint_ZY", annhilationPoint.Z(), annhilationPoint.Y());
          stats.fillHistogram("Annih_DLOR", deltaLor);
        }*/
        return true;
      }
    }
  }
  return false;
}

/**
* Method for determining type of event - 3Gamma
*/
bool EventCategorizerTools::checkFor3Gamma(const JPetEvent& event, JPetStatistics& stats, bool saveHistos)
{
  if (event.getHits().size() < 3) return false;
  for (uint i = 0; i < event.getHits().size(); i++) {
    for (uint j = i + 1; j < event.getHits().size(); j++) {
      for (uint k = j + 1; k < event.getHits().size(); k++) {
        JPetHit firstHit = event.getHits().at(i);
        JPetHit secondHit = event.getHits().at(j);
        JPetHit thirdHit = event.getHits().at(k);

        vector<double> thetaAngles;
        thetaAngles.push_back(firstHit.getBarrelSlot().getTheta());
        thetaAngles.push_back(secondHit.getBarrelSlot().getTheta());
        thetaAngles.push_back(thirdHit.getBarrelSlot().getTheta());
        sort(thetaAngles.begin(), thetaAngles.end());

        vector<double> relativeAngles;
        relativeAngles.push_back(thetaAngles.at(1) - thetaAngles.at(0));
        relativeAngles.push_back(thetaAngles.at(2) - thetaAngles.at(1));
        relativeAngles.push_back(360.0 - thetaAngles.at(2) + thetaAngles.at(0));
        sort(relativeAngles.begin(), relativeAngles.end());
        double transformedX = relativeAngles.at(1) + relativeAngles.at(0);
        double transformedY = relativeAngles.at(1) - relativeAngles.at(0);

        if (saveHistos) {
         // stats.fillHistogram("3Gamma_Angles", transformedX, transformedY);
        }
      }
    }
  }
  return true;
}

/**
* Method for determining type of event - prompt
*/
bool EventCategorizerTools::checkForPrompt(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos,
  double deexTOTCutMin, double deexTOTCutMax, std::string fTOTCalculationType)
{
  for (unsigned i = 0; i < event.getHits().size(); i++) {
    double tot = HitFinderTools::calculateTOT(event.getHits().at(i), 
                                              HitFinderTools::getTOTCalculationType(fTOTCalculationType));
    if (tot > deexTOTCutMin && tot < deexTOTCutMax) {
      if (saveHistos) {
        //stats.fillHistogram("Deex_TOT_cut", tot);
      }
      return true;
    }
  }
  return false;
}

/**
* Method for determining type of event - scatter
*/
bool EventCategorizerTools::checkForScatter(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos, double scatterTOFTimeDiff, 
  std::string fTOTCalculationType)
{
  if (event.getHits().size() < 2) {
    return false;
  }
  for (uint i = 0; i < event.getHits().size(); i++) {
    for (uint j = i + 1; j < event.getHits().size(); j++) {
      JPetHit primaryHit, scatterHit;
      if (event.getHits().at(i).getTime() < event.getHits().at(j).getTime()) {
        primaryHit = event.getHits().at(i);
        scatterHit = event.getHits().at(j);
      } else {
        primaryHit = event.getHits().at(j);
        scatterHit = event.getHits().at(i);
      }

      double scattAngle = calculateScatteringAngle(primaryHit, scatterHit);
      double scattTOF = calculateScatteringTime(primaryHit, scatterHit);
      double timeDiff = scatterHit.getTime() - primaryHit.getTime();

      if (saveHistos) {
        stats.fillHistogram("ScatterTOF_TimeDiff", fabs(scattTOF - timeDiff));
      }

      if (fabs(scattTOF - timeDiff) < scatterTOFTimeDiff) {
        if (saveHistos) {
          stats.fillHistogram("ScatterAngle_PrimaryTOT", scattAngle, HitFinderTools::calculateTOT(primaryHit, 
                                                        HitFinderTools::getTOTCalculationType(fTOTCalculationType)));
          stats.fillHistogram("ScatterAngle_ScatterTOT", scattAngle, HitFinderTools::calculateTOT(scatterHit, 
                                                        HitFinderTools::getTOTCalculationType(fTOTCalculationType)));
        }
        return true;
      }
    }
  }
  return false;
}

/**
* Calculation of distance between two hits
*/
double EventCategorizerTools::calculateDistance(const JPetHit& hit1, const JPetHit& hit2)
{
  return (hit1.getPos() - hit2.getPos()).Mag();
}

/**
* Calculation of time that light needs to travel the distance between primary gamma
* and scattered gamma. Return value in picoseconds.
*/
double EventCategorizerTools::calculateScatteringTime(const JPetHit& hit1, const JPetHit& hit2)
{
  return calculateDistance(hit1, hit2) / kLightVelocity_cm_ps;
}
double EventCategorizerTools::calcAngles(const JPetHit& hit1, const JPetHit& hit2)
{
   return (TMath::RadToDeg() * hit1.getPos().Angle(hit2.getPos()));
}
/**
* Calculation of scatter angle between primary hit and scattered hit.
* This function assumes that source of first gamma was in (0,0,0).
* Angle is calculated from scalar product, return value in degrees.
*/
double EventCategorizerTools::calculateScatteringAngle(const JPetHit& hit1, const JPetHit& hit2)
{
  return TMath::RadToDeg() * hit1.getPos().Angle(hit2.getPos() - hit1.getPos());
}
/**
 * calculate 2D angle between hits
 * */
double EventCategorizerTools::calculate2DOpenAngles(const JPetHit& hit1, const JPetHit& hit2)
{
double x1 = hit1.getPosX();
double x2 = hit2.getPosX();
double y1 = hit1.getPosY();
double y2 = hit2.getPosY();
double dot = x1*x2 + y1*y2;
double det = x1*y2 - y1*x2 ;

       return TMath::RadToDeg() * TMath::ATan2(det, dot);
}


/**
* Calculation point in 3D, where annihilation occured
*/
TVector3 EventCategorizerTools::calculateAnnihilationPoint(const JPetHit& hitA, const JPetHit& hitB)
{
  double tof = EventCategorizerTools::calculateTOF(hitA, hitB);
  return calculateAnnihilationPoint(hitA.getPos(), hitB.getPos(), tof);
}

TVector3 EventCategorizerTools::calculateAnnihilationPoint(const TVector3& hitA, const TVector3& hitB, double tof)
{
  TVector3 middleOfLOR = 0.5 * (hitA + hitB);
  TVector3 versorOnLOR = (hitB - hitA).Unit()  ;

  double shift = 0.5 * tof  * kLightVelocity_cm_ps;
  TVector3 annihilationPoint(middleOfLOR.X() + shift * versorOnLOR.X(),
                             middleOfLOR.Y() + shift * versorOnLOR.Y(),
                             middleOfLOR.Z() + shift * versorOnLOR.Z());
  return annihilationPoint;
}

double EventCategorizerTools::calculateTOFByConvention(const JPetHit& hitA, const JPetHit& hitB)
{
  if (hitA.getBarrelSlot().getTheta() < hitB.getBarrelSlot().getTheta()) {
    return calculateTOF(hitA, hitB);
  } else {
    return calculateTOF(hitB, hitA);
  }
}

double EventCategorizerTools::calculateTOF(const JPetHit& hitA, const JPetHit& hitB)
{
  return EventCategorizerTools::calculateTOF(hitA.getTime(), hitB.getTime());
}

double EventCategorizerTools::calculateTOF(double time1, double time2)
{
  return (time1 - time2);
}

/**
* Calculating distance from the center of the decay plane
*/
double EventCategorizerTools::calculatePlaneCenterDistance(
  const JPetHit& firstHit, const JPetHit& secondHit, const JPetHit& thirdHit)
{
  TVector3 crossProd = (secondHit.getPos() - firstHit.getPos()).Cross(thirdHit.getPos() - secondHit.getPos());
  double distCoef = -crossProd.X() * secondHit.getPosX() - crossProd.Y() * secondHit.getPosY() - crossProd.Z() * secondHit.getPosZ();
  if (crossProd.Mag() != 0) {
    return (distCoef) / crossProd.Mag();
  } else {
    ERROR("One of the hit has zero position vector - unable to calculate distance from the center of the surface");
    return -1.;
  }
}
