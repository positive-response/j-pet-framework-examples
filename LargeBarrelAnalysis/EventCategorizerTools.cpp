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
bool EventCategorizerTools::calculateMLParamsBefore(const JPetEvent& event, JPetStatistics& stats, bool saveHistos)
{
int eventSize = event.getHits().size();
	double fTOF = 0.0;
	double fDistance = 0.0;
        double fTimeDiff = 0.0;
	double f3DOpenAngle = 0.0;
	double f2DOpenAngle = 0.0;
	double fSmallestEnergy = 0.0;
	double fLargestEnergy = 0.0;
	double fSumEnergy = 0.0;
        double fDiffSumEnergyLarge = 0.0;
	double fDiffSumEnergySmall = 0.0;
	double fTheta = 0.0;
	double fPhi = 0.0;
	vector <double> hitEnergy{};
	vector <double> M{};

for(int i = 0; i < eventSize; i++)
{
	auto hitFirst = event.getHits().at(i);
	TVector3 v1(hitFirst.getPos());
	
	for(int j = i+1 ; j < eventSize; j++)
	{
		auto& hitSecond = event.getHits().at(j);
		TVector3 v2(hitSecond.getPos());
		fTOF = EventCategorizerTools::calculateTOFByConvention(hitFirst, hitSecond);
		fDistance = EventCategorizerTools::calculateDistance(hitFirst, hitSecond);
		fTimeDiff = fTOF - fDistance/kLightVelocity_cm_ps;
		f3DOpenAngle = EventCategorizerTools::calculateScatteringAngle(hitFirst, hitSecond);
		f2DOpenAngle = EventCategorizerTools::calculate2DOpenAngles(v1, v2);

		M.push_back(fTOF);
		stats.fillHistogram("TOFB", fTOF);
		stats.fillHistogram("TOF_vs_distanceB", fTOF, fDistance);
		stats.fillHistogram("TOF_vs_timeDifferenceB", fTOF,fTimeDiff);
		stats.fillHistogram("3DOpenAngleB", f3DOpenAngle);
		stats.fillHistogram("2DopenAngleB", f2DOpenAngle);

		if(eventSize == 5)
		{

			         stats.fillHistogram("TOF", fTOF);
		 		 stats.fillHistogram("TOF_vs_distance", fTOF, fDistance);
			 	 stats.fillHistogram("TOF_vs_timeDifference", fTOF,fTimeDiff);
				 stats.fillHistogram("3DOpenAngle", f3DOpenAngle);
                                 stats.fillHistogram("2DopenAngle", f2DOpenAngle);
		}
		
	}
	fSumEnergy += hitFirst.getEnergy();
	hitEnergy.push_back(hitFirst.getEnergy());
	fTheta = TMath::RadToDeg() * hitFirst.getPos().Theta();
	fPhi = TMath::RadToDeg() * hitFirst.getPos().Phi();
	stats.fillHistogram("ThetaB", fTheta);
	stats.fillHistogram("PhiB", fPhi);

}
for(int hi = 0; hi < eventSize; hi++)
{
	JPetHit hit = event.getHits().at(hi);
	hitEnergy.push_back(hit.getEnergy());
	fSumEnergy += hit.getEnergy();
}
sort(M.begin(), M.end());
sort(hitEnergy.begin(), hitEnergy.end());

fSmallestEnergy = hitEnergy.front();
fLargestEnergy = hitEnergy.back();
fDiffSumEnergyLarge = fSumEnergy - fSmallestEnergy;
fDiffSumEnergySmall = fSumEnergy - fLargestEnergy;

stats.fillHistogram("M_ijB", M.front());
stats.fillHistogram("SmallestEnergyB", fSmallestEnergy);
stats.fillHistogram("LargestEnergyB", fLargestEnergy);
stats.fillHistogram("EnergySumB", fSumEnergy);
stats.fillHistogram("EnergySum-largestB", fDiffSumEnergyLarge);
stats.fillHistogram("EnergySum-SmallestB", fDiffSumEnergySmall);
stats.fillHistogram("EnergySum-Smallest_vs_EnergySum-largestB", fDiffSumEnergySmall, fDiffSumEnergyLarge);
hitEnergy.clear();
fSumEnergy = 0.0;
return true;
}



bool EventCategorizerTools::calculateMLParamsAfter(const JPetEvent& event, JPetStatistics& stats, bool saveHistos)
{
int eventSize = event.getHits().size();
	double fTOF = 0.0;
	double fDistance = 0.0;
        double fTimeDiff = 0.0;
	double f3DOpenAngle = 0.0;
	double f2DOpenAngle = 0.0;
	double fSmallestEnergy = 0.0;
	double fLargestEnergy = 0.0;
	double fSumEnergy = 0.0;
        double fDiffSumEnergyLarge = 0.0;
	double fDiffSumEnergySmall = 0.0;
	double fTheta = 0.0;
	double fPhi = 0.0;
	vector <double> hitEnergy{};
	vector <double> M{};

for(int i = 0; i < eventSize; i++)
{
	auto hitFirst = event.getHits().at(i);
	TVector3 v1(hitFirst.getPos());
	
	for(int j = i+1 ; j < eventSize; j++)
	{
		auto& hitSecond = event.getHits().at(j);
		TVector3 v2(hitSecond.getPos());
		fTOF = EventCategorizerTools::calculateTOFByConvention(hitFirst, hitSecond);
		fDistance = EventCategorizerTools::calculateDistance(hitFirst, hitSecond);
		fTimeDiff = fTOF - fDistance/kLightVelocity_cm_ps;
		f3DOpenAngle = EventCategorizerTools::calculateScatteringAngle(hitFirst, hitSecond);
		f2DOpenAngle = EventCategorizerTools::calculate2DOpenAngles(v1, v2);

		M.push_back(fTimeDiff);
		stats.fillHistogram("TOF", fTOF);
		stats.fillHistogram("TOF_vs_distance", fTOF, fDistance);
		stats.fillHistogram("TOF_vs_timeDifference", fTOF,fTimeDiff);
		stats.fillHistogram("3DOpenAngle", f3DOpenAngle);
		stats.fillHistogram("2DopenAngle", f2DOpenAngle);
		
	}
	fSumEnergy += hitFirst.getEnergy();
	hitEnergy.push_back(hitFirst.getEnergy());
	fTheta = TMath::RadToDeg() * hitFirst.getPos().Theta();
	fPhi = TMath::RadToDeg() * hitFirst.getPos().Phi();
	stats.fillHistogram("Theta", fTheta);
	stats.fillHistogram("Phi", fPhi);

}
for(int hi = 0; hi < eventSize; hi++)
{
	JPetHit hit = event.getHits().at(hi);
	hitEnergy.push_back(hit.getEnergy());
	fSumEnergy += hit.getEnergy();
}
sort(M.begin(), M.end());
sort(hitEnergy.begin(), hitEnergy.end());

fSmallestEnergy = hitEnergy.front();
fLargestEnergy = hitEnergy.back();
fDiffSumEnergyLarge = fSumEnergy - fSmallestEnergy;
fDiffSumEnergySmall = fSumEnergy - fLargestEnergy;

stats.fillHistogram("M_ij", M.front());
stats.fillHistogram("SmallestEnergy", fSmallestEnergy);
stats.fillHistogram("LargestEnergy", fLargestEnergy);
stats.fillHistogram("EnergySum", fSumEnergy);
stats.fillHistogram("EnergySum-largest", fDiffSumEnergyLarge);
stats.fillHistogram("EnergySum-Smallest", fDiffSumEnergySmall);
stats.fillHistogram("EnergySum-Smallest_vs_EnergySum-largest", fDiffSumEnergySmall, fDiffSumEnergyLarge);
hitEnergy.clear();
fSumEnergy = 0.0;
return true;
}


bool EventCategorizerTools::checkForInitialCuts(const JPetEvent& event, JPetStatistics& stats, bool saveHistos, double fLowEnergyCut, double fAnnihilationEnergyCut)
{
	int eventSize = event.getHits().size();
	bool isInitialCut = true;
	int n = 0;
	for(auto & hit: event.getHits())
	{
		stats.fillHistogram("Hit_Z_POS(before)", hit.getPosZ());
		stats.fillHistogram("Energy(before)", hit.getEnergy());
		if(abs(hit.getPosZ()) < 23 && hit.getEnergy() < fAnnihilationEnergyCut && hit.getEnergy() > fLowEnergyCut)
		{
			n++;
			isInitialCut = true;
			stats.fillHistogram("Hit_Z_POS(after)", hit.getPosZ());
			stats.fillHistogram("Energy(after)", hit.getEnergy());
		}
		else	
		{
			isInitialCut = false;
		
		}
	}
	if (n != eventSize)
		return false;
	else
		return true;
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
        stats.fillHistogram("Deex_TOT_cut", tot);
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
double EventCategorizerTools::calculate2DOpenAngles(const TVector3& hit1, const TVector3& hit2)
{
     TVector3 hitFirst (hit1.X(), hit1.Y(), 0.0); 
     TVector3 hitSecond (hit2.X(), hit2.Y(), 0.0);
       return TMath::RadToDeg() * hitFirst.Angle(hitSecond);
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
