
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
#include <tuple>

using namespace std;
int ha,hb,hcounta,hcountb=0;
/**

* Method for determining type of event - back to back 2 gamma
*/
bool EventCategorizerTools::checkFor2Gamma(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos,
  double b2bSlotThetaDiff, double b2bTimeDiff)
{
  if (event.getHits().size() < 2)
    {
      return false;
    }

  stats.fillHistogram("Hit_multiplicity_2g", event.getHits().size());
  for (auto i = 0; i < event.getHits().size(); i++)
    {
    double tot = event.getHits().at(i).getEnergy();
    if (tot > 65000)
      {
	return false;
      }
    // stats.fillHistogram("simple_tot", tot);
     if (fabs(event.getHits().at(i).getPosZ()) > 23)
       {
         return false;
       }
    }

  for (auto i = 0; i < event.getHits().size(); i++){
    for (auto j = i + 1; j < event.getHits().size(); j++) {
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
      double absthetadiff = fabs(thetaDiff - 180);
      
      if (saveHistos) {
        stats.fillHistogram("2Gamma_Zpos", firstHit.getPosZ());
        stats.fillHistogram("2Gamma_Zpos", secondHit.getPosZ());
        stats.fillHistogram("2Gamma_TimeDiff", timeDiff / 1000.0);
        stats.fillHistogram("2Gamma_DLOR", deltaLor);
        stats.fillHistogram("2Gamma_ThetaDiff", thetaDiff);
        stats.fillHistogram("2Gamma_Dist", calculateDistance(firstHit, secondHit));
	stats.fillHistogram("2Gamma_absThetaDiff", absthetadiff);
      }
      if (fabs(thetaDiff - 180.0) < b2bSlotThetaDiff && timeDiff < b2bTimeDiff) {
        if (saveHistos) {
          TVector3 annhilationPoint = calculateAnnihilationPoint(firstHit, secondHit);
          stats.fillHistogram("Annih_TOF", calculateTOFByConvention(firstHit, secondHit));
          stats.fillHistogram("AnnihPoint_XY", annhilationPoint.X(), annhilationPoint.Y());
          stats.fillHistogram("AnnihPoint_ZX", annhilationPoint.Z(), annhilationPoint.X());
          stats.fillHistogram("AnnihPoint_ZY", annhilationPoint.Z(), annhilationPoint.Y());
          stats.fillHistogram("Annih_DLOR", deltaLor);
	  stats.fillHistogram("2Gamma_TimeDiff_after", timeDiff / 1000.0);
          stats.fillHistogram("2Gamma_ThetaDiff_after", thetaDiff);	  
        }
        return true;
      }
     
    }
  }
  return true;
}

/**
* Method for determining type of event - 3Gamma
*/
bool EventCategorizerTools::checkFor3Gamma(const JPetEvent& event, JPetStatistics& stats, bool saveHistos)
{
  if (event.getHits().size() < 3) return false;
  stats.fillHistogram("Hit_multiplicity_3g", event.getHits().size());
  
  for (auto i = 0; i < event.getHits().size(); i++)
    {      
    double tot = event.getHits().at(i).getEnergy();
    if (tot > 65000)
      {
	return false;
      }
    
    if (fabs(event.getHits().at(i).getPosZ()) > 23)
       {
         return false;
       }
    }
  
    for (auto i = 0; i < event.getHits().size(); i++){
      for (auto j = i + 1; j < event.getHits().size(); j++) {
	for (auto k = j + 1; k < event.getHits().size(); k++) {
        JPetHit firstHit = event.getHits().at(i);
        JPetHit secondHit = event.getHits().at(j);
        JPetHit thirdHit = event.getHits().at(k);

        vector<double> thetaAngles;
        thetaAngles.push_back(firstHit.getBarrelSlot().getTheta());
        thetaAngles.push_back(secondHit.getBarrelSlot().getTheta()) ;
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

  for (auto i = 0; i < event.getHits().size(); i++) {
    if(fabs(event.getHits().at(i).getPosZ()) > 23)
      {
        return false;
      }
 
    double tot = event.getHits().at(i).getEnergy();
    stats.fillHistogram("SYNC_TOT", tot);
    
       
    if (tot > deexTOTCutMin && tot < deexTOTCutMax){
	if (saveHistos){
	 stats.fillHistogram("Deex_TOT_cut", tot);
	 stats.fillHistogram("Hit_multiplicity_prompt", event.getHits().size());
	}
	return true;
    }
  }
  return true;
}


/** Method for determining type of event - Annihilation*/
std::tuple<int, int, bool> EventCategorizerTools::checkForAnnihilation(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos)
{
  if (event.getHits().size() < 2)
    {
      return std::make_tuple(hcountb, hcounta, false);
    }
 
  for (auto i = 0; i < event.getHits().size(); i++) {
             
       double tot = event.getHits().at(i).getEnergy();
       stats.fillHistogram("Ann_TOT_before_cut", tot);
       if (tot > 65000)
	 {
	   return std::make_tuple(hcountb, hcounta, false);
	 }
       
       double z = fabs( event.getHits().at(i).getPosZ());
       if (z > 23)
	 {
	   return std::make_tuple(hcountb, hcounta, false);
	 }
       
     if (saveHistos)
       {
	 stats.fillHistogram("Ann_TOT", tot); 
         stats.fillHistogram("Hit_multiplicity_ann", event.getHits().size());
       }
          
  /* hcountb=hcountb+hb;
    hcounta=hcounta+ha;
    hb=0;
    ha=0;*/
    }
  return std::make_tuple(hcountb, hcounta, true);
}



/** Method for removing neighbouring hits*/
bool EventCategorizerTools::removeNeighbourhits(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos,
  std::string fTOTCalculationType)
{
  if (event.getHits().size() < 2)
    {
    return false;
    }

  else if (event.getHits().size() == 2)
   {
     
  stats.fillHistogram("Hit_multiplicity_n2", event.getHits().size());
  for (auto i = 0; i < 2; i++)
    {
       double tot = event.getHits().at(i).getEnergy();
       if (tot > 65000)
         return false;
     }

     vector<JPetHit> hits = event.getHits();
     if (fabs(hits[0].getPosZ()) > 23 && fabs(hits[1].getPosZ()) > 23)
       {
         return false;
       }

  JPetHit firstHit = event.getHits().at(0);
  JPetHit secondHit = event.getHits().at(1);
  int scinID1 = firstHit.getScintillator().getID();
  int scinID2 = secondHit.getScintillator().getID();
  TVector3 v1 = firstHit.getPos();
  TVector3 v2 = secondHit.getPos();
  int sc_diff = fabs(scinID1-scinID2);
  double distance = fabs((v2-v1).Mag());
  double open_angles = TMath::RadToDeg() * v1.Angle(v2);
  double time = fabs((firstHit.getTime()-secondHit.getTime()));
  double dx = v2.X()-v1.X();
  double dy = v2.Y()-v1.Y();
  double dz = v2.Z()-v1.Z();
  double del_time = fabs((distance/kLightVelocity_cm_ps) -time);

  if (saveHistos) {
            stats.fillHistogram("Opening_angle_before", open_angles);
            stats.fillHistogram("Scintillator_before", sc_diff);
            stats.fillHistogram("opening_angle_Vs_distance_before", open_angles, distance);
            stats.fillHistogram("opening_angle_Vs_z_before", open_angles, dz);
            stats.fillHistogram("opening_angle_Vs_scinID1", open_angles, scinID1);
            stats.fillHistogram("opening_angle_Vs_scinID2", open_angles, scinID2);
            stats.fillHistogram("delta_time_before", del_time);
            stats.fillHistogram("scID Distribution", scinID1,scinID2,open_angles);
            stats.fillHistogram("scinID1_Vs_scinID2",scinID1, scinID2);
            stats.fillHistogram("time_diff_Vs_distance_before", distance, del_time);
          }

  if (scinID1 == scinID2)  {return false;}
           else  {
        if (saveHistos) {
            stats.fillHistogram("Opening_angle", open_angles);
            stats.fillHistogram("Scintillator", sc_diff);
            stats.fillHistogram("opening_angle_Vs_distance", open_angles, distance);
            stats.fillHistogram("opening_angle_Vs_z", open_angles, dz);
            stats.fillHistogram("delta_time", del_time);
            stats.fillHistogram("scID Distribution_after", scinID1,scinID2,open_angles);
            stats.fillHistogram("time_diff_Vs_distance", distance, del_time);
                }
           }
  return true;

  
   }


 else if (event.getHits().size() == 3)
   {
     stats.fillHistogram("Hit_multiplicity_n3", event.getHits().size());
     for (auto i = 0; i < 3; i++){
       double tot = event.getHits().at(i).getEnergy();
       if (tot > 65000)
	 return false;
     }

     vector<JPetHit> hits = event.getHits();
     if (fabs(hits[0].getPosZ()) > 23 && fabs(hits[1].getPosZ()) > 23 && fabs(hits[2].getPosZ()) > 23)
       {
	 return false;
       }
     
     TVector3 v1 =  hits[0].getPos();
     TVector3 v2 =  hits[1].getPos();
     TVector3 v3 =  hits[2].getPos();
     double t1 = hits[0].getTime();
     double t2 = hits[1].getTime();
     double t3 = hits[2].getTime();

     vector<double> Angles;
     Angles.push_back(TMath::RadToDeg() * v1.Angle(v2));
     Angles.push_back(TMath::RadToDeg() * v2.Angle(v3));
     Angles.push_back(TMath::RadToDeg() * v1.Angle(v3));
     sort(Angles.begin(), Angles.end());
     double theta_sum = Angles.at(1)+Angles.at(0);
     double theta_diff = Angles.at(1)-Angles.at(0);

     vector<double> dist;
     vector<double> dt;

     for(auto i = 0; i<3; i++)
       {
	 for(auto j = i+1; j<3; j++)
	   {

	     double d =fabs((hits[i].getPos()-hits[j].getPos()).Mag());
	     double t =fabs(hits[i].getTime()-hits[j].getTime());
	     dist.push_back(d);
	     dt.push_back(fabs((d/kLightVelocity_cm_ps) -t));
	   }
       }

     if (saveHistos)
       {
	 stats.fillHistogram("qqq1", t3-t2);
	 stats.fillHistogram("qqq2", t2-t1);
	 stats.fillHistogram("qqq3", t3-t1);
	 stats.fillHistogram("sum_diff_angle_dist", theta_sum,theta_diff);
	 stats.fillHistogram("xyz1",dt.at(0),dist.at(0));
	 stats.fillHistogram("xyz2",dt.at(1),dist.at(1));
	 stats.fillHistogram("xyz3",dt.at(2),dist.at(2));
	 stats.fillHistogram("t1",dt.at(0));
	 stats.fillHistogram("t2",dt.at(1));
	 stats.fillHistogram("t3",dt.at(2));
	 stats.fillHistogram("d1",dist.at(0));

       }
     return true;
   }

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
  for (auto i = 0; i < event.getHits().size(); i++) {
    for (auto j = i + 1; j < event.getHits().size(); j++) {
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
