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
#include "TEfficiency.h"
#include <iostream>
#include <TFile.h>


using namespace std;
int hcounta;
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
	  // stats.fillHistogram("Hit_multiplicity_3g", event.getHits().size());
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
  double deexTOTCutMin, double deexTOTCutMax, std::string fTOTCalculationType, int atleastNprompt)
{
  int NumberOfPrompts = 0;
  // stats.fillHistogram("Hit_multiplicity_prompt", event.getHits().size());
  double tot = 0.0;

  if (event.getHits().size() < atleastNprompt)
    {
      return false;
    }
  
  for (auto i = 0; i < event.getHits().size(); i++) {
    
    tot = event.getHits().at(i).getEnergy();
    stats.fillHistogram("SYNC_TOT", tot);
    
    if (tot > deexTOTCutMin && tot < deexTOTCutMax){
      NumberOfPrompts++;
      stats.fillHistogram("Deex_TOT_cut", tot);
    }
    
    if (NumberOfPrompts >= atleastNprompt)
      {	
	return true;
      }
  }
  
  stats.fillHistogram("Hit_multiplicity_prompt", event.getHits().size());
    return false;
}


/** Method for determining type of event - Annihilation*/

bool EventCategorizerTools::checkForAnnihilation(
const JPetEvent& event, JPetStatistics& stats, bool saveHistos, int atLeastNAnihilationHits, double TOT_Cut)
{
  int numberOfAnnihilation = 0;
  double tot = 0.0;
  double tota = 0.0;
  
  // stats.fillHistogram("Hit_multiplicity_ann0", event.getHits().size());
  if (event.getHits().size() < atLeastNAnihilationHits)
    {
      return false;
    }
  
  for (auto i = 0; i < event.getHits().size(); i++) {
  tot = event.getHits().at(i).getEnergy();     
  stats.fillHistogram("Ann_TOT_before_cut", tot);
  
  if (tot < TOT_Cut)
	 {
	   numberOfAnnihilation++;
	 }
       tota = event.getHits().at(i).getEnergy();
       stats.fillHistogram("Hit_multiplicity_ann0", event.getHits().size());
      
       if (numberOfAnnihilation>=atLeastNAnihilationHits)
	 {
	   stats.fillHistogram("Ann_TOT", tota);
	   return true;
	 }
  }
 
  stats.fillHistogram("Hit_multiplicity_ann1", event.getHits().size());
  return false;
}

/**Method for initial cuts**/


bool EventCategorizerTools::initialCut(                                   
const JPetEvent& event, JPetStatistics& stats, bool saveHistos, Counter& hitCounter)
{
 
  stats.fillHistogram("Hit_multiplicity_0", event.getHits().size());
  if (event.getHits().size() < 2 )
    {
      hitCounter.totalNumber++;
      return false;
    }                                                                                                                                                                                                   
  stats.fillHistogram("Hit_multiplicity_cut1", event.getHits().size()); 
  bool isZLessThan23 = true;
  for (auto i = 0; i < event.getHits().size(); i++) {
    hitCounter.totalNumber++;
    if (fabs( event.getHits().at(i).getPosZ()) > 23)
      {
	isZLessThan23 = false;
      }
    else{
      hitCounter.totalAccepted++;
    }
  }
  if (!isZLessThan23)	
    {
      return false;
    }

  stats.fillHistogram("Hit_multiplicity_cut2", event.getHits().size());
  return true;
}
  



/** Method for removing neighbouring hits*/
bool EventCategorizerTools::removeNeighbourhits(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos,
  std::string fTOTCalculationType)
{
  vector<JPetHit> hits = event.getHits();
 if (event.getHits().size() == 2)
   {
     
  stats.fillHistogram("Hit_multiplicity_n2", event.getHits().size());
  JPetHit firstHit = hits[0];
  JPetHit secondHit = hits[1];
  int scinID1 = firstHit.getScintillator().getID();
  int scinID2 = secondHit.getScintillator().getID();
  TVector3 v1 = firstHit.getPos();
  TVector3 v2 = secondHit.getPos();
  int sc_diff = fabs(scinID1-scinID2);
  double distance = fabs((v2-v1).Mag());
  double open_angles = TMath::RadToDeg() * v1.Angle(v2);
  double dt = fabs((firstHit.getTime()-secondHit.getTime()));
  double dx = v2.X()-v1.X();
  double dy = v2.Y()-v1.Y();
  double dz = v2.Z()-v1.Z();
  double del_time = fabs((distance/kLightVelocity_cm_ps) - dt);
  double del_time2 = fabs((distance/kLightVelocity_cm_ps));
  double del_time3 = fabs(dt);

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
            stats.fillHistogram("distance_vs_time_diff_before", distance, del_time);
	    stats.fillHistogram("distance_vs_time_diff_before2(dis)", distance, del_time2);
	    stats.fillHistogram("distance_vs_time_diff_before3(dt)", distance, dt);
	    stats.fillHistogram("time_difference", dt);
	    stats.fillHistogram("distance", distance);
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
            stats.fillHistogram("distance_vs_time_diff", distance, del_time);
                }
           }
  return true;
  
  
   }


 else if (event.getHits().size() == 3)
   {
     stats.fillHistogram("Hit_multiplicity_n3", event.getHits().size());
     
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
     double tot3 = 0.0;
     

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
     for (auto i =0; i < 3; i++)
       {
	 tot3 = hits[i].getEnergy();
       }

     if (saveHistos)
       {
	 stats.fillHistogram("hit_order1", t3-t2);
	 stats.fillHistogram("hit_order2", t2-t1);
	 stats.fillHistogram("hit_order3", t3-t1);
	 stats.fillHistogram("sum_diff_angle_dist", theta_sum,theta_diff);
	 stats.fillHistogram("xyz1",dist.at(0),dt.at(0));
	 stats.fillHistogram("xyz2",dist.at(1),dt.at(1));
	 stats.fillHistogram("xyz3",dist.at(2),dt.at(2));
	 stats.fillHistogram("t1",dt.at(0));
	 stats.fillHistogram("t2",dt.at(1));
	 stats.fillHistogram("t3",dt.at(2));
	 stats.fillHistogram("3hit_tot",tot3);
	 //	 stats.fillHistogram("xy",dt.at(0),dt.at(2));

       }
     // Hit_Eff->Fill(true,2);
     return true;
     
   }
     
 if (event.getHits().size()==5)
    {
      vector<TVector3> pos;
      vector<double> time;
      vector<double> dist;
      vector<double> dt;
      double tot5 = 0.0;
      stats.fillHistogram("Hit_multiplicity_n5", event.getHits().size());
      
      for(auto i = 0; i< 5; i++)
        {
          pos.push_back(hits[i].getPos());
          time.push_back(hits[i].getTime());
	}

      for(auto i = 0; i<5; i++)
        {
          for(auto j = i+1; j<5; j++)
            {
              double d =fabs((hits[i].getPos()-hits[j].getPos()).Mag());
              double t =fabs(hits[i].getTime()-hits[j].getTime());
             dist.push_back(d);
              dt.push_back(fabs((d/kLightVelocity_cm_ps) -t));
            }
        }

      for (auto i =0; i < 3; i++)
       {
         tot5 = hits[i].getEnergy();
       }

      if (saveHistos)
	{
            for(auto k =1; k<=10; k++)
              {
		// stats.fillHistogram(Form("D_vs_t%d", k), dist[k-1],dt[k-1]);
		//	stats.fillHistogram(Form("dt_vs_dt%d", k), dt[k-1],dt[k]);
		stats.fillHistogram("5hit_tot", tot5 );
		stats.fillHistogram("D_vs_dt", dist[k-1], dt[k-1]);
		stats.fillHistogram("dt_vs_dt", dt[k-1],dt[k]);
		  
              }
          }
      return true;
    }
 return true;
}

/*
std::pair<int, bool> EventCategorizerTools::checkFor5gamma(
const JPetEvent& event, JPetStatistics& stats, bool saveHistos)

{
  if (event.getHits().size() < 3)
    {
      return std::make_pair(false, n);
    }

  n = event.getHits().size();
  int c = (n*(n-1))/2;

  
  for(auto i = 0; i<event.getHits().size(); i++)
    {
      auto tot = event.getHits().at(i).getEnergy();
      if (tot > 65000)
	{
	  return std::make_pair(false, n);
	}
    }
  
  vector<JPetHit> hits = event.getHits();
  
  for(auto i = 0; i< event.getHits().size(); i++)
    {
     if (fabs(hits[i].getPosZ()) > 23)
       {
         return std::make_pair(false, c);
       }
    }
  
  if (event.getHits().size()==5)
    {
      vector<TVector3> pos;
      vector<double> time;
      vector<double> dist;
      vector<double> dt;
      for(auto i = 0; i< event.getHits().size(); i++)
	{
	  pos.push_back(hits[i].getPos());
	  time.push_back(hits[i].getTime());
	}
      

      for(auto i = 0; i<event.getHits().size(); i++)
	{
	  for(auto j = i+1; j<event.getHits().size(); j++)
	    {
	      double d =fabs((hits[i].getPos()-hits[j].getPos()).Mag());
	      double t =fabs(hits[i].getTime()-hits[j].getTime());
             dist.push_back(d);
              dt.push_back(fabs((d/kLightVelocity_cm_ps) -t));
	    }
	}
      
      if (saveHistos)
	  {
	    for(auto k =1; k<=10; k++)
	      {
		stats.fillHistogram(Form("D_vs_t%d", k), dist[k-1],dt[k-1]);
	      }
	  }
      return std::make_pair(true, n);
    }

}
*/

/**
* Method for determining type of event - scatter
*/
bool EventCategorizerTools::checkForScatter(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos, double scatterTOFTimeDiff, 
  std::string fTOTCalculationType)
{  
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

      if (abs(primaryHit.getScintillator().getID() - scatterHit.getScintillator().getID()) <= 2)
        {
          return true;
        }


      double scattAngle = calculateScatteringAngle(primaryHit, scatterHit);
      double scattTOF = calculateScatteringTime(primaryHit, scatterHit);
      double timeDiff = scatterHit.getTime() - primaryHit.getTime();

 

      if (saveHistos) {
        stats.fillHistogram("ScatterTOF_TimeDiff", fabs(scattTOF - timeDiff));
	stats.fillHistogram("ScatterAngle_PrimaryTOT_before", scattAngle, primaryHit.getEnergy());
	stats.fillHistogram("ScatterAngle_ScatterTOT_before", scattAngle, scatterHit.getEnergy());

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
