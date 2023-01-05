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
    double tot = event.getHits().at(i).getEnergy();					    
    if (saveHistos){
      stats.fillHistogram("SYNC_TOT", tot);
       if (tot > deexTOTCutMin && tot < deexTOTCutMax)
	 stats.fillHistogram("Deex_TOT_cut", tot);
    }
    return true;
     }                                                                                 
  return false;
}

/** Method for determining type of event - Annihilation*/
std::tuple<int, int, bool> EventCategorizerTools::checkForAnnihilation(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos)
{

  if (event.getHits().size() < 2) {
       return std::make_tuple(hcountb, hcounta, false);
   }
  for (unsigned i = 0; i < event.getHits().size(); i++) {
       hb++;
       double tot = event.getHits().at(i).getEnergy();
     if (saveHistos){stats.fillHistogram("Ann_TOT_before_cut", tot);}

     if (tot > 65000) {
      return std::make_tuple(hcountb, hcounta, false);
    }

    if (saveHistos){ 
      stats.fillHistogram("Ann_TOT", tot);
      ha++;}
    hcountb=hcountb+hb;
    hcounta=hcounta+ha;
    hb=0;
    ha=0;
    }
  return std::make_tuple(hcountb, hcounta, true);
}



/** Method for removing neighbouring hits*/
bool EventCategorizerTools::removeNeighbourhits(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos,
  std::string fTOTCalculationType)
{
  if (event.getHits().size() < 2){
    return false;
  }

  else if (event.getHits().size() == 2)
    {
  int nhit = event.getHits().size();
  for (uint i = 0; i < nhit; i++) {
    for (uint j = i + 1; j < nhit; j++) {
      
        JPetHit firstHit = event.getHits().at(i);
        JPetHit secondHit = event.getHits().at(j);
	int scinID1 = firstHit.getScintillator().getID();
        int scinID2 = secondHit.getScintillator().getID();

	int sc_diff = fabs(scinID1-scinID2);
	//	double distance = fabs((firstHit.getPos() - secondHit.getPos()).Mag());
	double time = fabs((firstHit.getTime()-secondHit.getTime()));
   
	TVector3 v1 = firstHit.getPos();
	TVector3 v2 = secondHit.getPos();
	double open_angles = TMath::RadToDeg() * v1.Angle(v2);
	 
       double dx = v2.X()-v1.X();
       double dy = v2.Y()-v1.Y();
       double dz = v2.Z()-v1.Z();
              double distance = sqrt(dx * dx + dy * dy + dz * dz);
	      // double distance = fabs((v2-v1).Mag()); 
	double del_time = fabs((distance/kLightVelocity_cm_ps) -time);
	
       	  if (saveHistos) {
            stats.fillHistogram("Opening_angle_before", open_angles);
            stats.fillHistogram("Scintillator_before", sc_diff);  
	    stats.fillHistogram("opening_angle_Vs_distance_before", open_angles, distance);
	    stats.fillHistogram("opening_angle_Vs_z_before", open_angles, dz);
	    stats.fillHistogram("delta_time_before", del_time);


	  }

           if (scinID1 == scinID2)  {return false;}
	   else  {		 
        if (saveHistos) {
	    stats.fillHistogram("Opening_angle", open_angles);
	    stats.fillHistogram("Scintillator", sc_diff);
	    stats.fillHistogram("opening_angle_Vs_distance", open_angles, distance);
	    stats.fillHistogram("opening_angle_Vs_z", open_angles, dz);
	    stats.fillHistogram("delta_time", del_time);

	}	
	     }
    }
  }//for loop for 2 hits ends here
         return true;
    }//2 hits check ends here
  
 else if (event.getHits().size() == 3)
    {
        int nhit = event.getHits().size();
        for (unsigned i = 0; i < nhit; i++) {
             double tot = event.getHits().at(i).getEnergy();
              if (tot < 90000){

                  for (uint i = 0; i < nhit; i++) {
                      for (uint j = i + 1; j < nhit; j++) {
                         for (uint k = j + 1; k < nhit; k++) {
                               JPetHit hit1 = event.getHits().at(i);
                               JPetHit hit2 = event.getHits().at(j);
                               JPetHit hit3 = event.getHits().at(k);

			       

                               vector<double> Angles;
                               Angles.push_back(hit1.getBarrelSlot().getTheta());
                               Angles.push_back(hit2.getBarrelSlot().getTheta());
                               Angles.push_back(hit3.getBarrelSlot().getTheta());
                               sort(Angles.begin(), Angles.end());

                               vector<double> relAngles;
                               relAngles.push_back(Angles.at(1) - Angles.at(0));
                               relAngles.push_back(Angles.at(2) - Angles.at(1));
                               relAngles.push_back(360.0 - Angles.at(2) + Angles.at(0));
                               sort(relAngles.begin(), relAngles.end());

                               TVector3 vec1 = hit2.getPos() - hit1.getPos();
	                       TVector3 vec2 = hit3.getPos() - hit1.getPos();
                               TVector3 cross = vec1.Cross(vec2);
        if (cross.Mag() == 0) {
	  if (saveHistos) {
            for (auto angle : relAngles) {
        stats.fillHistogram("Opening_angle_3", angle); }
	    }
	  }
	 else {return false;}//cross product is non zero  
	}
      }
		      //     else {return false;}//cross product is non zero
    }
  }
	      }// 3 for loop ends here
	      return true; } // tot check ends here
	return true;} // 3 hit checks end here




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
