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

          stats.fillHistogram(("relAngles12_23 " + postfix).c_str(), fabs(calcAngles(firstHit, secondHit)), fabs(calcAngles(secondHit, thirdHit)));  
          std::sort(relativeAngles.begin(), relativeAngles.end());

          double transformedX = relativeAngles.at(1) + relativeAngles.at(0);
          double transformedY = relativeAngles.at(1) - relativeAngles.at(0);
          stats.fillHistogram(("3g_rel_angles " + postfix).c_str(), transformedX, transformedY);
        }
      }
    }
  return true;
}

bool EventCategorizerTools::Plots2hits(const JPetEvent& event, JPetStatistics& stats, const std::string& postfix){
if(event.getHits().size() < 3) return false;
//std::vector<double> M(100);
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
      //double fscatDist = abs(fDistance - fTOF);
     // M.push_back(fTimeDiff);

      stats.fillHistogram(("timeDiff " + postfix).c_str(), fHitTimeDiff);
      stats.fillHistogram(("TOF " + postfix).c_str(), fTOF); 
      stats.fillHistogram(("Distance " + postfix).c_str(), fDistance);
      stats.fillHistogram(("openAngles " + postfix).c_str(), openAng);

	    stats.fillHistogram(("timeDiff_vs_angle " + postfix).c_str(), fHitTimeDiff, openAng);
      stats.fillHistogram(("Distance_vs_angle " + postfix).c_str(), fDistance, openAng);
      stats.fillHistogram(("timeDiff_vs_Distance " + postfix).c_str(), fHitTimeDiff, fDistance);
      
      //for checking scattering
      stats.fillHistogram(("M_ij " + postfix).c_str(), fTimeDiff);
      //stats.fillHistogram(("MDist_vs_angle " + postfix).c_str(), fscatDist, openAng);
      stats.fillHistogram(("Mdt_vs_angle " + postfix).c_str(), fTimeDiff, openAng);

      stats.fillHistogram(("ap_XY " + postfix).c_str(), ap.X(), ap.Y() );
      //stats.fillHistogram(("ap_XZ " + postfix).c_str(), ap.X(), ap.Z() );
      //stats.fillHistogram(("ap_YZ " + postfix).c_str(), ap.Y(), ap.Z() );
     // stats.fillHistogram(("TOF_vs_distance " + postfix).c_str(), fTOF, fDistance);
     // stats.fillHistogram(("TOF_vs_timeDifference " + postfix).c_str(), fTOF,fTimeDiff);   
      
    }   
  }
   /*sort(M.begin(), M.end());
   stats.fillHistogram(("M_ij_minimum " + postfix).c_str(), M.front());
   stats.fillHistogram(("M_ij_maximum " + postfix).c_str(), M.back());
   M.clear(); */
return true;
}

bool EventCategorizerTools::Plots1hit(const JPetEvent& event, JPetStatistics& stats, const std::string& postfix)
{
  if(event.getHits().size() < 3) return false;
  double fSumEnergy = 0.0;
  double smallestEnergy = std::numeric_limits<double>::max();
  double largestEnergy = std::numeric_limits<double>::lowest();

  for(int i = 0; i < event.getHits().size(); i++)
  {
    auto hit = event.getHits().at(i);
    auto energy =  hit.getEnergy()/1000.0;

    int scinID = hit.getScintillator().getID();

    if(energy < smallestEnergy){smallestEnergy = energy;}
	  if(energy > largestEnergy){largestEnergy = energy; }
    
    fSumEnergy += energy;

    stats.fillHistogram( ("HM " + postfix).c_str(), event.getHits().size());
    stats.fillHistogram(("TOT " + postfix).c_str(), energy);
    stats.fillHistogram(("Z " + postfix).c_str(), hit.getPos().Z());
   // stats.fillHistogram(("Phi " + postfix).c_str(),TMath::RadToDeg() * hit->getPos().Phi());
   // stats.fillHistogram(("Theta " + postfix).c_str(), TMath::RadToDeg() * hit->getPos().Theta());
    //stats.fillHistogram(("hit_XY " + postfix).c_str(), hit->getPos().X(), hit->getPos().Y() );
    stats.fillHistogram(("scinID_VS_TOT " + postfix).c_str(), scinID, energy);
  }
   // stats.fillHistogram(("SmallestEnergy " + postfix).c_str(), smallestEnergy/1000);
   // stats.fillHistogram(("LargestEnergy " + postfix).c_str(), largestEnergy/1000);
  
    stats.fillHistogram(("smallest_vs_largest " + postfix).c_str(), smallestEnergy, largestEnergy);
    stats.fillHistogram(("EnergySum " + postfix).c_str(), fSumEnergy);
    //stats.fillHistogram(("EnergySum-largest " + postfix).c_str(), (fSumEnergy - largestEnergy)/1000);
    //stats.fillHistogram(("EnergySum-Smallest " + postfix).c_str(), (fSumEnergy - smallestEnergy)/1000);
    stats.fillHistogram(("EnergySum-Smallest_vs_EnergySum-largest "+ postfix).c_str(), (fSumEnergy - smallestEnergy), (fSumEnergy - largestEnergy));

    fSumEnergy = 0.0;
    smallestEnergy = std::numeric_limits<double>::max();
    largestEnergy = std::numeric_limits<double>::lowest();
  return true;
}

/****************************************************** */
bool EventCategorizerTools::zTest(const JPetEvent& event, JPetEvent& zEvent, int minC){
 if(event.getHits().size() < 3) return false;
 bool isZValid = false;
  zEvent =  JPetEvent();
  int count = 0;

  for(auto & hit: event.getHits())
  {
    if(abs(hit.getPos().Z()) <= 23)
    {
      zEvent.addHit(hit);
      count++;		
    }
    else
      continue;
  }
  
  if(count >= minC)
  {
    isZValid = true;
  }
  return isZValid;
}

bool EventCategorizerTools::totCheck(const JPetEvent& event, JPetEvent& totEvent, double totCutMin, double totCutMax, int atleastNHit){
 if(event.getHits().size() < 3) return false;
 int Count = 0;
 totEvent = JPetEvent();

 for(auto & hit: event.getHits()){
     auto tot =  hit.getEnergy()/1000.0;
     if(tot >= totCutMin && tot <= totCutMax ){ 
      totEvent.addHit(hit);
      Count++;
      }
   }
 if (Count >= atleastNHit){ return true; }

 return false;
}

bool EventCategorizerTools::ScIDCheck(const JPetEvent& event, JPetEvent& scEvent)
{
  if(event.getHits().size() < 3) return false;
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



bool EventCategorizerTools::checkForCut1(const JPetEvent& event, JPetStatistics& stats, JPetEvent& sel1Event)
{ 
  bool checkForZ = false;
  JPetEvent zEvent;
  int minHits = 1;
  checkForZ = zTest(event, zEvent, minHits);
  
  if(checkForZ)
  {
    sel1Event = zEvent;
    Plots1hit(zEvent, stats, "zCut");
    Plots2hits(zEvent, stats, "zCut");
    Plots3hits(zEvent, stats, "zCut");
    return true;
  }
  return false;
}


bool EventCategorizerTools::checkForCut2(const JPetEvent& event, JPetStatistics& stats, JPetEvent& sel3Event, int n)
{
  bool totCheckDeexHits = false;
  bool totCheckHits = false;
  bool hitTypeSel = false;

  JPetEvent totEventAnn;
  JPetEvent totEventDeex;

  totCheckDeexHits = totCheck(event, totEventDeex, fToTCutDeexMin, fToTCutDeexMax, 1);
  totCheckHits = totCheck(event, totEventAnn, fToTCutAnniMin, fToTCutAnniMax, n);
  
  if (!totCheckDeexHits || !totCheckHits) 
  {
    return false;
  }
  
  if(totCheckDeexHits && totCheckHits)
  {
    sel3Event = totEventAnn;

    Plots3hits(totEventAnn, stats, "totCut");
    Plots2hits(totEventAnn, stats, "totCut");
    Plots1hit(totEventAnn, stats, "totCut");
    return true;
  }
return false;
}

bool EventCategorizerTools::checkForCut3(const JPetEvent& event, JPetStatistics& stats, JPetEvent& sel3Event)
{ 
  bool checkForSC = false;
  JPetEvent scEvent;
  int minHits = 1;
  checkForSC = ScIDCheck(event, scEvent);
  
  if(checkForSC)
  {
    sel3Event = scEvent;
    Plots1hit(scEvent, stats, "SCCut");
    Plots2hits(scEvent, stats, "SCCut");
    Plots3hits(scEvent, stats, "SCCut");
    return true;
  }
  return false;
}

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
      if (saveHistos) {
        stats.fillHistogram("2Gamma_Zpos", firstHit.getPosZ());
        stats.fillHistogram("2Gamma_Zpos", secondHit.getPosZ());
        stats.fillHistogram("2Gamma_TimeDiff", timeDiff / 1000.0);
        stats.fillHistogram("2Gamma_DLOR", deltaLor);
        stats.fillHistogram("2Gamma_ThetaDiff", thetaDiff);
        stats.fillHistogram("2Gamma_Dist", calculateDistance(firstHit, secondHit));
      }
      if (fabs(thetaDiff - 180.0) < b2bSlotThetaDiff && timeDiff < b2bTimeDiff) {
        if (saveHistos) {
          TVector3 annhilationPoint = calculateAnnihilationPoint(firstHit, secondHit);
          stats.fillHistogram("Annih_TOF", calculateTOFByConvention(firstHit, secondHit));
          stats.fillHistogram("AnnihPoint_XY", annhilationPoint.X(), annhilationPoint.Y());
          stats.fillHistogram("AnnihPoint_ZX", annhilationPoint.Z(), annhilationPoint.X());
          stats.fillHistogram("AnnihPoint_ZY", annhilationPoint.Z(), annhilationPoint.Y());
          stats.fillHistogram("Annih_DLOR", deltaLor);
        }
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
          stats.fillHistogram("3Gamma_Angles", transformedX, transformedY);
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
    return fabs(distCoef) / crossProd.Mag();
  } else {
    ERROR("One of the hit has zero position vector - unable to calculate distance from the center of the surface");
    return -1.;
  }
}
