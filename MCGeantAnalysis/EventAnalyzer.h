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
 *  @file EventAnalyzer.h
 */

#ifndef EVENTANALYZER_H
#define EVENTANALYZER_H

#include <JPetEvent/JPetEvent.h>
#include <JPetHit/JPetHit.h>
#include <JPetUserTask/JPetUserTask.h>
#include "../LargeBarrelAnalysis/EventCategorizerTools.h"
#include "TRandom.h"


class EventAnalyzer : public JPetUserTask {
public:
  EventAnalyzer(const char *name);
  virtual ~EventAnalyzer();
  virtual bool init() override;
  virtual bool exec() override;
  virtual bool terminate() override;

  static const double kUnknownEventType;

protected:
  void fillResolutionHistograms(const JPetEvent &event,const JPetTimeWindowMC *tw);
  bool hitParser(const JPetEvent &event,const JPetTimeWindowMC *tw, JPetEvent& sel1Event);
  double calcAngles(const JPetHit& hit1, const JPetHit& hit2);
  bool Twohit(const JPetEvent &event);
  bool Fourhit(const JPetEvent &event);
  //bool saveEvt(const JPetEvent &event, const JPetTimeWindowMC *tw);

	TRandom r;
	int random = r.Rndm();
	//std::string Name = "mltest"+std::to_string(random)+".root";
  std::string Name = "mltest"+std::to_string(10)+".root";
  TFile *file = TFile::Open(Name.c_str(), "RECREATE");	
	TTree* pTree22 = nullptr;

  /*std::vector<float> fPosX {};
	std::vector<float> fPosY {};
	std::vector<float> fPosZ {};
	std::vector<float> fEnergy {};
	std::vector<float> fTime {};
	std::vector<unsigned int> fHitType {};
	std::vector<unsigned int> fVtxIndex {};
  */
  float fPosX;
	float fPosY;
	float fPosZ;
	float fEnergy;
	float fTime;
	int fHitType;
	//int fVtxIndex;

                        
  bool fIsMC = false;
};
#endif /* !EVENTANALYZER_H */
