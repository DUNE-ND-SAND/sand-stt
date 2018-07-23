#include <map>
#include <iostream>
#include <vector>

#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4Event.h"
#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4HitSegment.h"

#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoTrd2.h>
#include <TChain.h>

#include "/mnt/nas01/users/mtenti/wd/analysis/KLOEcal/loader/loader.C"

void STTDigi()
{
  std::cout << "DigitizeSTT(const char* finname, const char* foutname)" << std::endl;
  std::cout << "finname could contain wild card" << std::endl;
}

void Cluster(TG4Event* ev, TGeoManager* geo, std::map<std::string,std::vector<hit> >& cluster_map)
{
  cluster_map.clear();
    
  for(unsigned int j = 0; j < ev->SegmentDetectors["StrawTracker"].size(); j++)
  {
    const TG4HitSegment& hseg = ev->SegmentDetectors["StrawTracker"].at(j);
    
    double x = 0.5 * (hseg.Start.X() + hseg.Stop.X());
    double y = 0.5 * (hseg.Start.Y() + hseg.Stop.Y());
    double z = 0.5 * (hseg.Start.Z() + hseg.Stop.Z());
    
    std:string sttname = geo->FindNode(x,y,z)->GetName();
    
    hit h;
    h.det = sttname;
    h.x1 = hseg.Start.X();
    h.y1 = hseg.Start.Y();
    h.z1 = hseg.Start.Z();
    h.t1 = hseg.Start.T();
    h.x2 = hseg.Stop.X();
    h.y2 = hseg.Stop.Y();
    h.z2 = hseg.Stop.Z();
    h.t2 = hseg.Stop.T();
    h.de = hseg.EnergyDeposit;
    h.pid = hseg.PrimaryId;
    h.index = j;
    
    std::string cluster_name(sttname);
    cluster_name += "_" + std::to_string(hseg.PrimaryId);
    
    cluster_map[cluster_name].push_back(h);
  }
}

void Cluster2Digit(std::map<std::string,std::vector<hit> >& cluster_map, std::vector<digit>& digit_vec)
{    
  for(std::map<std::string,std::vector<hit> >::iterator it = cluster_map.begin(); it != cluster_map.end(); ++it)
  {      
    digit d;
    d.de = 0;    
    d.det = it->second[0].det;  
  
    for(unsigned int k = 0; k < it->second.size(); k++)
    {
      d.hindex.push_back(it->second[k].index);
      d.de += it->second[k].de;
    }
    
    d.hor = (d.det.find("STTPlane1FULL") != std::string::npos) ? false : true;
    
    std::sort(it->second.begin(), it->second.end(), isHitBefore);
    
    d.t = it->second.at(0).t1;
    d.x = 0.5 * (it->second.front().x1 + it->second.back().x2);
    d.y = 0.5 * (it->second.front().y1 + it->second.back().y2);
    d.z = 0.5 * (it->second.front().z1 + it->second.back().z2);
    
    digit_vec.push_back(d);
  }
}

void Digitize(TG4Event* ev, TGeoManager* geo, std::vector<digit>& digit_vec)
{
  std::map<std::string,std::vector<hit> > cluster_map;
  digit_vec.clear();
  
  Cluster(ev, geo, cluster_map);
  Cluster2Digit(cluster_map, digit_vec);
}

void DigitizeSTT(const char* finname, const char* foutname)
{
  
  TChain* t = new TChain("EDepSimEvents","EDepSimEvents");
  t->Add(finname);
  TFile f(t->GetListOfFiles()->At(0)->GetTitle());
  TGeoManager* geo = (TGeoManager*) f.Get("EDepSimGeometry");
  
  TG4Event* ev = new TG4Event;
  
  t->SetBranchAddress("Event",&ev);
  
  std::vector<digit> digit_vec;
  
  TFile fout(foutname,"RECREATE");
  TTree tstt("tSttDigi","KLOE Stt Digitization");
  tstt.Branch("Stt","std::vector<digit>",&digit_vec);
    
  const int nev = t->GetEntries();
  
  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;

  for(int i = 0; i < nev; i++)
  {
    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
      
    t->GetEntry(i);
    
    Digitize(ev, geo, digit_vec);
    
    tstt.Fill();
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
    
  fout.cd();
  tstt.Write();
  fout.Close();
  
  delete t;
  f.Close();
}