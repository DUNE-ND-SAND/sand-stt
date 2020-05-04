#include <TGeoManager.h>
#include <TGeoBBox.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TParameter.h>

#include <iostream>
#include <vector>
#include <map>

#include "/wd/dune-it/ext_bkg/kloe-simu/src/display.cpp"

void erase_element(int val, std::vector<int>& vec)
{
  std::vector<int>::iterator ite = std::find(vec.begin(),vec.end(),val);
  if(ite == vec.end())
  {
    std::cout << "problem: " << val << " not found in vec (size: " << vec.size() << ")" << std::endl;
    for(std::vector<int>::iterator it = vec.begin(); it != vec.end(); ++it)
    {
      std::cout << distance(vec.begin(), it) << " " << *it << std::endl;
    }
    exit(1);
  } 
  else
  {
    vec.erase(ite);
  }
}

void processMtrVct(std::vector<int>& regMtrX, std::vector<int>& reg0trX, std::map<int,int>& regionsX, TH1I& hmultX, const int minreg)
{
    while(regMtrX.size() != 0)
    {
      std::map<int,int>::iterator this_el = regionsX.find(regMtrX.front());
      std::map<int,int>::iterator next_el = std::next(this_el);
      std::map<int,int>::iterator prev_el = std::prev(this_el);
      std::map<int,int>::iterator next2_el = std::next(this_el,2);
      std::map<int,int>::iterator prev2_el = std::prev(this_el,2);
      bool mantain = false;
      
      if(next_el != regionsX.end() && this_el != regionsX.begin())
      {
        if(this_el->second < minreg && hmultX.GetBinContent(prev_el->first) == 1 && hmultX.GetBinContent(next_el->first) == 1)
        {
          for(int k = 0; k < this_el->second; k++)
          {
            hmultX.SetBinContent(this_el->first+k,1);
          }
          prev_el->second += this_el->second + next_el->second;
          regionsX.erase(this_el);
          regionsX.erase(next_el);
          regMtrX.erase(regMtrX.begin());        
          
          continue;
        }
      }
      
      if(next_el != regionsX.end())
      {        
        if(next_el->second < minreg)
        {
          mantain = true;
          
          if(next2_el != regionsX.end())
          {
            if(hmultX.GetBinContent(next2_el->first) == 2)
            {
              if(hmultX.GetBinContent(next_el->first) == 0)
                erase_element(next_el->first, reg0trX);
              for(int k = 0; k < next_el->second; k++)
              {
                hmultX.SetBinContent(next_el->first+k,2);
              }
              this_el->second += next_el->second + next2_el->second;
              erase_element(next2_el->first, regMtrX);
              regionsX.erase(next_el);
              regionsX.erase(next2_el);
            }
            else
            {
              if(hmultX.GetBinContent(next_el->first) == 0)
                erase_element(next_el->first, reg0trX);
              for(int k = 0; k < next_el->second; k++)
              {
                hmultX.SetBinContent(next_el->first+k,2);
              }
              this_el->second += next_el->second;
              regionsX.erase(next_el);
            }
          }
          else
          {
            if(hmultX.GetBinContent(next_el->first) == 0)
              erase_element(next_el->first, reg0trX);
            for(int k = 0; k < next_el->second; k++)
            {
              hmultX.SetBinContent(next_el->first+k,2);
            }
            this_el->second += next_el->second;
            regionsX.erase(next_el);
          }
        }
      }
      
      if(this_el != regionsX.begin())
      {
        if(prev_el->second < minreg)
        {
          mantain = true;
          
          if(prev_el != regionsX.begin())
          {
            if(hmultX.GetBinContent(prev2_el->first) == 2)
            {
              if(hmultX.GetBinContent(prev_el->first) == 0)
                erase_element(prev_el->first, reg0trX);
              for(int k = 0; k < prev_el->second; k++)
              {
                hmultX.SetBinContent(prev_el->first+k,2);
              }
              prev2_el->second += this_el->second + prev_el->second;
              erase_element(prev2_el->first, regMtrX);
              regionsX.erase(prev_el);
              regionsX.erase(this_el);
              regMtrX.front() = prev2_el->first;
              
            }
            else
            {
              if(hmultX.GetBinContent(prev_el->first) == 0)
                erase_element(prev_el->first, reg0trX);
              for(int k = 0; k < prev_el->second; k++)
              {
                hmultX.SetBinContent(prev_el->first+k,2);
              }
              prev_el->second += this_el->second;
              regionsX.erase(this_el);
              regMtrX.front() = prev_el->first;
            }
          }
          else
          {
            if(hmultX.GetBinContent(prev_el->first) == 0)
              erase_element(prev_el->first, reg0trX);
            for(int k = 0; k < prev_el->second; k++)
            {
              hmultX.SetBinContent(prev_el->first+k,2);
            }
            prev_el->second += this_el->second;
            regionsX.erase(this_el);
            regMtrX.front() = prev_el->first;
          }
        }
      }
      if(mantain == false)
      {
        regMtrX.erase(regMtrX.begin());
      }
    }
}

void process0trVct(std::vector<int>& reg0trX, std::map<int,int>& regionsX, TH1I& hmultX, const int minreg)
{
    while(reg0trX.size() != 0)
    {
      std::map<int,int>::iterator this_el = regionsX.find(reg0trX.front());
      std::map<int,int>::iterator next_el = std::next(this_el);
      std::map<int,int>::iterator prev_el = std::prev(this_el);
      std::map<int,int>::iterator next2_el = std::next(this_el,2);
      std::map<int,int>::iterator prev2_el = std::prev(this_el,2);
      bool mantain = false;      
      
      if(next_el != regionsX.end() && this_el != regionsX.begin())
      {
        if(this_el->second < minreg && hmultX.GetBinContent(prev_el->first) == 1 && hmultX.GetBinContent(next_el->first) == 1)
        {
          for(int k = 0; k < this_el->second; k++)
          {
            hmultX.SetBinContent(this_el->first+k,1);
          }
          prev_el->second += this_el->second + next_el->second;
          regionsX.erase(this_el);
          regionsX.erase(next_el); 
          reg0trX.erase(reg0trX.begin());
          continue;
        }
      }
      
      if(next_el != regionsX.end())
      {        
        if(next_el->second < minreg)
        {
          mantain = true;
          
          if(next2_el != regionsX.end())
          {
            if(hmultX.GetBinContent(next2_el->first) == 0)
            {
              for(int k = 0; k < next_el->second; k++)
              {
                hmultX.SetBinContent(next_el->first+k,0);
              }
              this_el->second += next_el->second + next2_el->second;
              erase_element(next2_el->first, reg0trX);
              regionsX.erase(next_el);
              regionsX.erase(next2_el);
              
            }
            else
            {
              for(int k = 0; k < next_el->second; k++)
              {
                hmultX.SetBinContent(next_el->first+k,hmultX.GetBinContent(next2_el->first));
              }
              next_el->second += next2_el->second;
              regionsX.erase(next2_el);   
            }
          }
          else
          {
            for(int k = 0; k < next_el->second; k++)
            {
              hmultX.SetBinContent(next_el->first+k,0);
            }
            this_el->second += next_el->second;
            regionsX.erase(next_el);
          }
        }
      }
      
      if(this_el != regionsX.begin())
      {        
        if(prev_el->second < minreg)
        {
          mantain = true;
          
          if(prev_el != regionsX.begin())
          {
            if(hmultX.GetBinContent(prev2_el->first) == 0)
            {
              for(int k = 0; k < prev_el->second; k++)
              {
                hmultX.SetBinContent(prev_el->first+k,0);
              }
              prev2_el->second += prev_el->second + this_el->second;
              erase_element(prev2_el->first, reg0trX);
              regionsX.erase(prev_el);
              regionsX.erase(this_el);  
              reg0trX.front() = prev2_el->first;        
            }
            else
            {
              for(int k = 0; k < prev_el->second; k++)
              {
                hmultX.SetBinContent(prev_el->first+k,hmultX.GetBinContent(prev2_el->first));
              }
              prev_el->second += this_el->second;
              regionsX.erase(this_el);
              reg0trX.front() = prev_el->first;
            }
          }
          else
          {
            for(int k = 0; k < prev_el->second; k++)
            {
              hmultX.SetBinContent(prev_el->first+k,0);
            }
            prev_el->second += this_el->second;
            regionsX.erase(this_el);
            reg0trX.front() = prev_el->first;
          }
        }
      }
      if(mantain == false)
      {
        reg0trX.erase(reg0trX.begin());
      }
    }
}

void MeanAndRMS(std::vector<digit>* digits, TH1D& hmeanX, TH1D& hrmsX, TH1I& hnX, TH1I& hmultX, TH1D& hmeanY, TH1D& hrmsY, TH1I& hnY, TH1I& hmultY)
{
  
    double mean, mean2, var, rms;
    int n;
  
    for(unsigned int j = 0; j < digits->size(); j++)
    {
      if(digits->at(j).hor == 1)
      {
        hrmsY.Fill(digits->at(j).z,digits->at(j).y*digits->at(j).y);
        hmeanY.Fill(digits->at(j).z,digits->at(j).y);
        hnY.Fill(digits->at(j).z);
      }
      else
      {
        hrmsX.Fill(digits->at(j).z,digits->at(j).x*digits->at(j).x);
        hmeanX.Fill(digits->at(j).z,digits->at(j).x);
        hnX.Fill(digits->at(j).z);
      }
    }
    
    for(unsigned int j = 0; j < hmeanX.GetNbinsX(); j++)
    {
      mean = -1;
      rms = -1;
      n = hnX.GetBinContent(j+1);
      
      if(n  != 0)
      {
        mean = hmeanX.GetBinContent(j+1)/n;
        mean2 = hrmsX.GetBinContent(j+1)/n;
        var = mean2 - mean * mean;
        rms = sqrt(var);
      }
      if(n>1) n=2;
      
      hrmsX.SetBinContent(j+1,rms);
      hmeanX.SetBinContent(j+1,mean);
      hmultX.SetBinContent(j+1,n);
      
      mean = -1;
      rms = -1;
      n = hnY.GetBinContent(j+1);
      if(n  != 0)
      {
        mean = hmeanY.GetBinContent(j+1)/n;
        mean2 = hrmsY.GetBinContent(j+1)/n;
        var = mean2 - mean * mean;
        rms = sqrt(var);
      }
      if(n>1) n=2;
      
      hrmsY.SetBinContent(j+1,rms);
      hmeanY.SetBinContent(j+1,mean);
      hmultY.SetBinContent(j+1,n);
    }
}

void fillRegionXY(TH1I& hmult, TH1I& hmultX, TH1I& hmultY, std::map<int,int>& regions)
{
    int last_v = -1;
    int last_b = -1;
    
    double mult;
    
    for(unsigned int j = 0; j < hmult.GetNbinsX(); j++)
    {
      mult = std::max(hmultX.GetBinContent(j+1),hmultY.GetBinContent(j+1));
      hmult.SetBinContent(j+1, mult);
      
      if(mult == last_v)
      {
        regions[last_b]++;
      }
      else
      {
        regions[j+1] = 1;
        last_b = j+1;
        last_v = mult;
      }
    }
}

void fillRegion(TH1I& hmultX, std::map<int,int>& regionsX)
{
    int last_vx = -1;
    int last_bx = -1;
    
    for(unsigned int j = 0; j < hmultX.GetNbinsX(); j++)
    {
      if(hmultX.GetBinContent(j+1) == last_vx)
      {
        regionsX[last_bx]++;
      }
      else
      {
        regionsX[j+1] = 1;
        last_bx = j+1;
        last_vx = hmultX.GetBinContent(j+1);
      }
    }
}

void orderByWidth(TH1I& hmultX, std::map<int,int>& regionsX, std::vector<int>& regMtrX, std::vector<int>& reg0trX)
{
    for (std::map<int,int>::iterator it=regionsX.begin(); it!=regionsX.end(); ++it)
    {
      if(hmultX.GetBinContent(it->first) == 2)
      {
        std::vector<int>::iterator ite;
        for (ite=regMtrX.begin(); ite!=regMtrX.end(); ++ite)
        {
          if(regionsX[*ite] < it->second)
            break;
        }
        regMtrX.insert(ite,it->first);
      }
      if(hmultX.GetBinContent(it->first) == 0)
      {
        std::vector<int>::iterator ite;
        for (ite=reg0trX.begin(); ite!=reg0trX.end(); ++ite)
        {
          if(regionsX[*ite] < it->second)
            break;
        }
        reg0trX.insert(ite,it->first);
      }
    }
}

void filterDigitsModuleY(std::map<int, double>& yd, std::vector<int>& toremove, const double epsilon)
{
  double mean = 0.;
  double var = 0.;
  double rms = 0.;
  double rms_thr = 0.;
  
  for(std::map<int, double>::iterator it = yd.begin(); it != yd.end(); ++it)
  {
    mean += it->second;
    var += it->second*it->second;
  }
  
  mean /= yd.size();
  var /= yd.size();
  var -= mean *mean;
  rms = sqrt(var)* (1. + epsilon);
  
  double maxd = 0.;
  std::map<int, double>::iterator idx; 
  
  for(std::map<int, double>::iterator it = yd.begin(); it != yd.end(); ++it)
  {
    if(abs(it->second - mean) > maxd)
    {
      maxd = abs(it->second - mean);
      idx = it;
    }
  }
  
  if(maxd > rms && rms > rms_thr)
  {
    toremove.push_back(idx->first);
    yd.erase(idx);
    filterDigitsModuleY(yd, toremove, epsilon);
  }
}

void filterDigitsModuleX(std::map<int, double>& xd, std::vector<int>& toremove, const double epsilon)
{
  double mean = 0.;
  double var = 0.;
  double rms = 0.;
  double rms_thr = 0.;
  
  for(std::map<int, double>::iterator it = xd.begin(); it != xd.end(); ++it)
  {
    mean += it->second;
    var += it->second*it->second;
  }
  
  mean /= xd.size();
  var /= xd.size();
  var -= mean *mean;
  rms = sqrt(var)* (1. + epsilon);
  
  double maxd = 0.;
  std::map<int, double>::iterator idx; 
  
  for(std::map<int, double>::iterator it = xd.begin(); it != xd.end(); ++it)
  {
    if(abs(it->second - mean) > maxd)
    {
      maxd = abs(it->second - mean);
      idx = it;
    }
  }
  
  if(maxd > rms && rms > rms_thr)
  {
    toremove.push_back(idx->first);
    xd.erase(idx);
    filterDigitsModuleY(xd, toremove, epsilon);
  }
}

void filterDigits(std::vector<digit>* digits, TH1D& hdummy, const double epsilon)
{
  std::map<int, std::map<int, double> > dy;
  std::map<int, std::map<int, double> > dx;
  std::vector<int> toremove;
  for(unsigned int i = 0; i < digits->size(); i++)
  {
    if(digits->at(i).hor == 1)
    {
      dy[hdummy.FindBin(digits->at(i).z)][i] = digits->at(i).y;
    }
    else
    {
      dx[hdummy.FindBin(digits->at(i).z)][i] = digits->at(i).x;
    }
  }
  
  for(std::map<int, std::map<int, double> >::iterator it = dy.begin(); it != dy.end(); ++it)
  {
    filterDigitsModuleY(it->second, toremove, epsilon);
  }
  
  for(std::map<int, std::map<int, double> >::iterator it = dx.begin(); it != dx.end(); ++it)
  {
    filterDigitsModuleX(it->second, toremove, epsilon);
  }
  
  std::sort(toremove.begin(),toremove.end());
  
  for(int i = toremove.size() - 1; i >= 0; i--)
  {
    digits->erase(digits->begin()+toremove[i]);
  }
}

void VtxFinding(int minreg = 3, double epsilon = 0.1, bool wideMultReg = false)
{
  gROOT->SetBatch();

  bool display = false;
  bool save2root = true && display;
  bool save2pdf = true && display;

// root [0] TGeoManager::Import("../geo/nd_hall_kloe_empty.gdml")
// (TGeoManager *) 0x3037200
// root [12] gGeoManager->cd("volWorld/rockBox_lv_0/volDetEnclosure_0/volKLOE_0")
// (bool) true
// root [16] double vorigin[3]
// (double [3]) { 0.0000000, 0.0000000, 0.0000000 }
// root [17] double vmaster[3]
// (double [3]) { 0.0000000, 0.0000000, 0.0000000 }
// root [18] gGeoManager->LocalToMaster(vorigin,vmaster)
// root [19] vmaster
// (double [3]) { 0.0000000, -238.47300, 2391.0000 }
// root [7] v = gGeoManager->GetVolume("volSTTFULL")
// (TGeoVolume *) 0x17fc5b10
// root [9] TGeoTube* tb = (TGeoTube*) v->GetShape()
// root [11] tb->GetRmin()
// (double) 0.0000000
// root [12] tb->GetRmax()
// (double) 200.00000
// root [13] tb->GetDz()
// (double) 169.00000

  double kloe_center[] = { 0.0000000, -2384.7300, 23910.000 };
  double kloe_size[] = { 2 * 1690.00000, 2 * 2000.00000, 2 * 2000.00000};

  TFile f("../files/reco/numu_internal_10k.0.reco.root");
  
  TTree* tDigit = (TTree*) f.Get("tDigit");
  TTree* tMC = (TTree*) f.Get("EDepSimEvents");
  TGeoManager* g = (TGeoManager*) f.Get("EDepSimGeometry");
  
  TString path_prefix = "volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/MagIntVol_volume_PV_0/volSTTFULL_PV_0/";
  TGeoVolume* v = g->FindVolumeFast("volSTTFULL_PV");
  
  double origin[3];
  double master[3];
  double last_z = 0.;
  double dz;
  std::vector<double> binning;
  bool is_first = true;
  
  origin[0] = 0.;
  origin[1] = 0.;
  origin[2] = 0.;
  
  // assuming they are order by Z position
  for(int i = 0; i < v->GetNdaughters(); i++)
  {
    TString name = v->GetNode(i)->GetName();
    
    if(name.Contains("sttmod") || name.Contains("volfrontST"))
    {
      TString path = path_prefix + name;
      g->cd(path.Data());
      g->LocalToMaster(origin,master);
      TGeoBBox* b = (TGeoBBox*) v->GetNode(i)->GetVolume()->GetShape();
      dz = b->GetDX();
    
      if(is_first)
      {
        is_first = false;
        binning.push_back(master[2] - dz);
      }
      else
      {
        binning.push_back(0.5 * (last_z + master[2] - dz));
      }
      last_z = master[2] + dz;
    }
  }
  binning.push_back(last_z);
  
  const double z_sampl = 44.4;
  const double xy_sampl = 5;
  
  int nall = 0;
  int nok = 0;
  
  double vrmsX = 0.;
  double vrmsY = 0.;
  double vrmsZ = 0.;
  double vrms3D = 0.;
  double vmeanX = 0.;
  double vmeanY = 0.;
  double vmeanZ = 0.;
  double vmean3D = 0.;
  double norm_dist = 0.;
  
  std::cout << std::setw(20) << "Z sampling (mm)" << std::setw(20) << "XY sampling (mm)" << std::endl;
  std::cout << std::setw(20) << z_sampl << std::setw(20) << xy_sampl << std::endl;
  
  /*
  for(unsigned int i = 0; i < binning.size()-1; i++)
  {
    std::cout << i << " " << binning.data()[i] << " " << binning.data()[i+1] << " " << binning.data()[i+1] - binning.data()[i] << std::endl;
  }
  */
    
  std::vector<digit>* digits = new std::vector<digit>;
  TG4Event* ev = new TG4Event;
  
  tDigit->SetBranchAddress("Stt",&digits);
  tMC->SetBranchAddress("Event",&ev);
  
  TFile fout("vtx.root","RECREATE");
  TTree tv("tv","tv");
  
  double xvtx_true, yvtx_true, zvtx_true;
  int nint;
  double x_int[256];
  double y_int[256];
  double z_int[256];
  
  double xvtx_reco, yvtx_reco, zvtx_reco;
    
  unsigned int n0 = 0;
  unsigned int n1 = 0;
  unsigned int nM = 0;
  
  int isCC;
  int VtxMulti;
  int VtxMono;
    
  tv.Branch("xvtx_true",&xvtx_true,"xvtx_true/D");
  tv.Branch("yvtx_true",&yvtx_true,"yvtx_true/D");
  tv.Branch("zvtx_true",&zvtx_true,"zvtx_true/D");
  
  tv.Branch("xvtx_reco",&xvtx_reco,"xvtx_reco/D");
  tv.Branch("yvtx_reco",&yvtx_reco,"yvtx_reco/D");
  tv.Branch("zvtx_reco",&zvtx_reco,"zvtx_reco/D");
  
  tv.Branch("VtxMulti",&VtxMulti,"VtxMulti/I");
  tv.Branch("VtxMono",&VtxMono,"VtxMono/I");
  tv.Branch("isCC",&isCC,"isCC/I");
  //tv.Branch("nint",&nint,"nint/I"); 
  tv.Branch("n0",&n0,"n0/I"); 
  tv.Branch("n1",&n1,"n1/I"); 
  tv.Branch("nM",&nM,"nM/I"); 
  //tv.Branch("x_int",x_int,"x_int[nint]/D"); 
  //tv.Branch("y_int",y_int,"y_int[nint]/D"); 
  //tv.Branch("z_int",z_int,"z_int[nint]/D"); 
  //tv.Branch("digits","std::vector<digit>",&digits); 
  
  TH1D hrmsX("hrmsX","rmsX;Z (mm); rmsX (mm)",binning.size()-1,binning.data());
  TH1D hrmsY("hrmsY","rmsY;Z (mm); rmsY (mm)",binning.size()-1,binning.data());
  hrmsY.SetLineColor(kRed);
    
  TH1D hmeanX("hmeanX","meanX;Z (mm); meanX (mm)",binning.size()-1,binning.data());
  TH1D hmeanY("hmeanY","meanY;Z (mm); meanY (mm)",binning.size()-1,binning.data());
  hmeanY.SetLineColor(kRed);
    
  TH1I hnX("hnX","nX;Z (mm); nX",binning.size()-1,binning.data());
  TH1I hnY("hnY","nY;Z (mm); nY",binning.size()-1,binning.data());
  hnY.SetLineColor(kRed);
  
  TH1D h0trX("h0trX","h0trX;Z (mm); rmsX (mm)",binning.size()-1,binning.data());
  TH1D h1trX("h1trX","h1trX;Z (mm); rmsX (mm)",binning.size()-1,binning.data());
  TH1D hmtrX("hmtrX","hmtrX;Z (mm); rmsX (mm)",binning.size()-1,binning.data());
  
  TH1D h0trY("h0trY","h0trY;Z (mm); rmsY (mm)",binning.size()-1,binning.data());
  TH1D h1trY("h1trY","h1trY;Z (mm); rmsY (mm)",binning.size()-1,binning.data());
  TH1D hmtrY("hmtrY","hmtrY;Z (mm); rmsY (mm)",binning.size()-1,binning.data());
  
  TH1D h0tr("h0tr","h0tr;Z (mm); rmsY (mm)",binning.size()-1,binning.data());
  TH1D h1tr("h1tr","h1tr;Z (mm); rmsY (mm)",binning.size()-1,binning.data());
  TH1D hmtr("hmtr","hmtr;Z (mm); rmsY (mm)",binning.size()-1,binning.data());
  
  TH1I hmult("hmult","hmult;Z (mm); multipliciy",binning.size()-1,binning.data());
  TH1I hmultX("hmultX","hmultX;Z (mm); multipliciy",binning.size()-1,binning.data());
  TH1I hmultY("hmultY","hmultY;Z (mm); multipliciy",binning.size()-1,binning.data());
  
  TH1D hdummy("hdummy","hdummy;Z (mm); multipliciy",binning.size()-1,binning.data());
  
  h0trX.SetFillColorAlpha(kRed, 0.15);
  h1trX.SetFillColorAlpha(kBlue, 0.15);
  hmtrX.SetFillColorAlpha(kGreen, 0.15);
  
  h0trY.SetFillColorAlpha(kRed, 0.15);
  h1trY.SetFillColorAlpha(kBlue, 0.15);
  hmtrY.SetFillColorAlpha(kGreen, 0.15);
  
  h0tr.SetFillColorAlpha(kRed, 0.15);
  h1tr.SetFillColorAlpha(kBlue, 0.15);
  hmtr.SetFillColorAlpha(kGreen, 0.15);
  
  h0trX.SetLineColor(0);
  h1trX.SetLineColor(0);
  hmtrX.SetLineColor(0);
  
  h0trX.SetLineWidth(0);
  h1trX.SetLineWidth(0);
  hmtrX.SetLineWidth(0);
  
  h0trY.SetLineColor(0);
  h1trY.SetLineColor(0);
  hmtrY.SetLineColor(0);
  
  h0trY.SetLineWidth(0);
  h1trY.SetLineWidth(0);
  hmtrY.SetLineWidth(0);
  
  h0tr.SetLineColor(0);
  h1tr.SetLineColor(0);
  hmtr.SetLineColor(0);
  
  h0tr.SetLineWidth(0);
  h1tr.SetLineWidth(0);
  hmtr.SetLineWidth(0);
  
  hrmsX.SetStats(false);
  hrmsY.SetStats(false);
  hmeanX.SetStats(false);
  hmeanY.SetStats(false);
  hnX.SetStats(false);
  hnY.SetStats(false);
  h0trX.SetStats(false);
  h1trX.SetStats(false);
  hmtrX.SetStats(false);
  h0trY.SetStats(false);
  h1trY.SetStats(false);
  hmtrY.SetStats(false);
  h0tr.SetStats(false);
  h1tr.SetStats(false);
  hmtr.SetStats(false);
  hmult.SetStats(false);
  hmultX.SetStats(false);
  hmultY.SetStats(false);
  
  TParameter<double> xv("xv",0);
  TParameter<double> yv("yv",0);
  TParameter<double> zv("zv",0);
  
  init("../files/reco/numu_internal_10k.0.reco.root");
  
  int first = 0;
  int last = 9999;/*tDigit->GetEntries();*/
  int nev = last - first + 1;
    
  vector<double> zvX;
  vector<double> zvY; 
  vector<double> zvtx;
  
  vector<TH1D> vharctg; 
  vector<TGraph>  vguv;
  vector<TH2D> vhHT; 
  
  TCanvas c;
  c.SaveAs("rms.pdf(");

  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;
  
  for(int i = first; i < last+1; i++)
  {
    
    tDigit->GetEntry(i);
    tMC->GetEntry(i);

    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i) / nev * 100)
              << "%]" << std::flush;
    
    TString reaction = ev->Primaries.at(0).GetReaction();
    VtxMulti = 0;
    VtxMono = 0;
    
    if(reaction.Contains("CC") == false)
      isCC = 0;
    else
      isCC = 1;
    
    xv.SetVal(ev->Primaries.at(0).GetPosition().X());
    yv.SetVal(ev->Primaries.at(0).GetPosition().Y());
    zv.SetVal(ev->Primaries.at(0).GetPosition().Z());
    
    xvtx_true = xv.GetVal();
    yvtx_true = yv.GetVal();
    zvtx_true = zv.GetVal();
    
    TDirectoryFile* fd;
    
    if(save2root)
      fd = new TDirectoryFile(TString::Format("ev_%d",i).Data(),TString::Format("ev_%d",i).Data(),"",&fout);
    
    hrmsX.Reset("ICESM");
    hrmsY.Reset("ICESM");
    hmeanX.Reset("ICESM");
    hmeanY.Reset("ICESM");
    hnX.Reset("ICESM");
    hnY.Reset("ICESM");
    h0trX.Reset("ICESM");
    h1trX.Reset("ICESM");
    hmtrX.Reset("ICESM");
    h0trY.Reset("ICESM");
    h1trY.Reset("ICESM");
    hmtrY.Reset("ICESM");
    h0tr.Reset("ICESM");
    h1tr.Reset("ICESM");
    hmtr.Reset("ICESM");
    
    hrmsX.SetTitle(TString::Format("event: %d",i).Data());
    hrmsY.SetTitle(TString::Format("event: %d",i).Data());
    hmeanX.SetTitle(TString::Format("event: %d",i).Data());
    hmeanY.SetTitle(TString::Format("event: %d",i).Data());
    hnX.SetTitle(TString::Format("event: %d",i).Data());
    hnY.SetTitle(TString::Format("event: %d",i).Data());
    h0trX.SetTitle(TString::Format("event: %d",i).Data());
    h1trX.SetTitle(TString::Format("event: %d",i).Data());
    hmtrX.SetTitle(TString::Format("event: %d",i).Data());
    h0trY.SetTitle(TString::Format("event: %d",i).Data());
    h1trY.SetTitle(TString::Format("event: %d",i).Data());
    hmtrY.SetTitle(TString::Format("event: %d",i).Data());
    h0tr.SetTitle(TString::Format("event: %d",i).Data());
    h1tr.SetTitle(TString::Format("event: %d",i).Data());
    hmtr.SetTitle(TString::Format("event: %d",i).Data());
    
    // filter digits
    // digits distant more than rms from mean are removed
    filterDigits(digits, hdummy, epsilon);
    
    // eveluate rms, mean and multiplicity in Y and X as a function of the STT module
    MeanAndRMS(digits, hmeanX, hrmsX, hnX, hmultX, hmeanY, hrmsY, hnY, hmultY);
    
    // find regions of: 0 tracks, 1 track, multi tracks
    std::map<int,int> regions;
    std::map<int,int> regionsX;
    std::map<int,int> regionsY;
    std::vector<int> regMtr;
    std::vector<int> regMtrX;
    std::vector<int> regMtrY;
    std::vector<int> reg0tr;
    std::vector<int> reg0trX;
    std::vector<int> reg0trY;
    
    fillRegion(hmultX, regionsX);
    fillRegion(hmultY, regionsY);
    fillRegionXY(hmult, hmultX, hmultY, regions);
    
    // order regions by width
    orderByWidth(hmultX, regionsX, regMtrX, reg0trX);
    orderByWidth(hmultY, regionsY, regMtrY, reg0trY);
    orderByWidth(hmult, regions, regMtr, reg0tr);
    
    // merge regions less than threshold [in number of modules] (threshold to be optimized) 
    // to the adiacent ones starting from wider multi track regions
    // after multi tracks regions => wider 0 tracks regions
    // if regions of multi tracks and 0 tracks less than 
    // threshold are in the middle of 1 track region, the former
    // is merged to the latter
    
    processMtrVct(regMtr, reg0tr, regions, hmult, minreg);
    process0trVct(reg0tr, regions, hmult, minreg);
      
    n0 = 0;  
    n1 = 0;  
    nM = 0;  
    
    for(std::map<int,int>::iterator it = regions.begin(); it != regions.end(); ++it)
    {
      if(hmult.GetBinContent(it->first) == 0)
      {
        n0++;
      }
      else if(hmult.GetBinContent(it->first) == 1)
      {
        n1++;
      }
      else if(hmult.GetBinContent(it->first) == 2)
      {
        nM++;
      }      
    }
    
    if(n0 + n1 + nM != regions.size())
      std::cout << "warning line:" << __LINE__ << std::endl;
    
    // find vertex as modules with less spreas between tracks
    int minb = 0;
    int maxb = hmult.GetNbinsX();
    
    regMtr.clear();
    reg0tr.clear();
    
    if(wideMultReg)
    {
      orderByWidth(hmult, regions, regMtr, reg0tr);
      if(regMtr.size()>0)
      {
        std::map<int,int>::iterator it = regions.find(regMtr.front());
        minb = it->first;
        maxb = it->first + it->second;
      }
    }
    
    int idx = 0;
    double rms = 10000.; 
    double rmsX, rmsY;
    VtxMulti = 0;
    for(int j = minb; j < maxb; j++)
    {
      if(hmult.GetBinContent(j+1) == 2 && hrmsX.GetBinContent(j+1) > 0. && hrmsY.GetBinContent(j+1) > 0.)
      {
        rmsX = hrmsX.GetBinContent(j+1);
        rmsY = hrmsY.GetBinContent(j+1);
        
        if(rmsX*rmsX+rmsY*rmsY < rms)
        {
          rms = rmsX*rmsX+rmsY*rmsY;
          idx = j+1;
          VtxMulti = 1;
        }
      }
    }
    
    if(VtxMulti == 0)
    {
      for(int j = 0; j < hmult.GetNbinsX(); j++)
      {
        if(hmult.GetBinContent(j+1) == 1 && hrmsX.GetBinContent(j+1) > -1. && hrmsY.GetBinContent(j+1) > -1.)
        {
          idx = j+1;
          VtxMono = 1;
          break;
        }
      }
    }
    
    if(VtxMulti == 1 || VtxMono == 1)
    {
      xvtx_reco = hmeanX.GetBinContent(idx);
      yvtx_reco = hmeanY.GetBinContent(idx);
      zvtx_reco = hmult.GetXaxis()->GetBinLowEdge(idx);
    }
    else
    {
      xvtx_reco = -1E10;
      yvtx_reco = -1E10;
      zvtx_reco = -1E10;
    }
    
    const int tol = 3;
    
    if(VtxMulti == 1)
    {
      nall++;
      if(abs(xvtx_reco-xvtx_true) < tol*xy_sampl && abs(yvtx_reco-yvtx_true) < tol*xy_sampl && abs(zvtx_reco-zvtx_true) < tol*z_sampl)
        nok++;
      
      norm_dist = sqrt(pow((xvtx_reco-xvtx_true)/xy_sampl,2)+pow((yvtx_reco-yvtx_true)/xy_sampl,2)+pow((zvtx_reco-zvtx_true)/z_sampl,2));
      
      vmeanX += xvtx_reco-xvtx_true;
      vmeanY += yvtx_reco-yvtx_true;
      vmeanZ += zvtx_reco-zvtx_true;
      vmean3D += norm_dist;
      
      vrmsX += (xvtx_reco-xvtx_true)*(xvtx_reco-xvtx_true);
      vrmsY += (yvtx_reco-yvtx_true)*(yvtx_reco-yvtx_true);
      vrmsZ += (zvtx_reco-zvtx_true)*(zvtx_reco-zvtx_true);
      vrms3D += norm_dist * norm_dist;
    }
    
    /*
    // find surface between:
    // - 0 tracks => 1 tracks : mono prong vertex
    // - 0 tracks => multi tracks: multi prong vertex
    // - 1 track => multi tracks: re-interaction or vertex with back-scattering
    zvX.clear();
    zvY.clear();
    
    for(std::map<int,int>::iterator it = regionsX.begin(); std::next(it) != regionsX.end(); ++it)
    {
      if(hmultX.GetBinContent(it->first) == 0 && hmultX.GetBinContent(std::next(it)->first) != 0)
      {
        zvX.push_back(hmultX.GetXaxis()->GetBinLowEdge(std::next(it)->first));
      }
      else if(hmultX.GetBinContent(it->first) == 1 && hmultX.GetBinContent(std::next(it)->first) == 2)
      {
        zvX.push_back(hmultX.GetXaxis()->GetBinLowEdge(std::next(it)->first));
      }
    }
    
    for(std::map<int,int>::iterator it = regionsY.begin(); std::next(it) != regionsY.end(); ++it)
    {
      if(hmultY.GetBinContent(it->first) == 0 && hmultY.GetBinContent(std::next(it)->first) != 0)
      {
        zvY.push_back(hmultY.GetXaxis()->GetBinLowEdge(std::next(it)->first));
      }
      else if(hmultY.GetBinContent(it->first) == 1 && hmultY.GetBinContent(std::next(it)->first) == 2)
      {
        zvY.push_back(hmultY.GetXaxis()->GetBinLowEdge(std::next(it)->first));
      }
    }*/
    
    // merge vertex between X and Y plan
    // if they are in between 200 mm => to be optimized
    /*
    std::map<double,int> vv;
    for(unsigned int j = 0; j < zvX.size(); j++)
    {
      vv[zvX[j]] = 1;
    }
    for(unsigned int j = 0; j < zvY.size(); j++)
    {
      vv[zvY[j]] = -1;
    }
    
    zvtx.clear();
    const double max_dist = 0.; //mm
    while(vv.size() != 0)
    {
      std::map<double,int>::iterator it = vv.begin();
      std::map<double,int>::iterator itt = it;
      bool found = false;
      while(++itt != vv.end())
      {
        if(abs(itt->first-it->first) < max_dist)
        {
          if(it->second != itt->second)
          {
            if(std::next(itt) != vv.end())
            {
              if(abs(itt->first-it->first) < abs(itt->first-std::next(it)->first) || itt->second == std::next(itt)->second)
              {
                found = true;
              }
            }
            else
            {
              found = true;
            }
            break;
          }
        }
      }
      if(found)
      {
        zvtx.push_back(0.5*(itt->first+it->first));
        vv.erase(itt);
        vv.erase(it);
      }
      else
      {
        zvtx.push_back(it->first);
        vv.erase(it);
      }
    }
    */
    /*
    double u, v, d;
    double zc, yc;
    double zz, yy;
    
    vguv.clear();
    vharctg.clear();
    vhHT.clear();
    
    nint = zvtx.size();
    for(unsigned int k = 0; k < zvtx.size(); k++)
    {
      x_int[k] = hmeanX.GetBinContent(hmeanX.FindBin(zvtx[k]));
      y_int[k] = hmeanY.GetBinContent(hmeanY.FindBin(zvtx[k]));
      z_int[k] = zvtx[k];
      
      if(x_int[k] == -1. || y_int[k] == -1.)
      {
        std::cout << "event: " << i << " ; z_int[k]   = " << z_int[k] << std::endl;
        std::cout << "warning....: x_int[k]   = " << x_int[k] << " ; y_int[k]   = " << y_int[k] << std::endl;
        std::cout << "next   ....: x_int[k+1] = " << hmeanX.GetBinContent(hmeanX.FindBin(zvtx[k])+1) << " ; y_int[k+1] = " << hmeanY.GetBinContent(hmeanY.FindBin(zvtx[k])+1) << std::endl;
        continue;
      }
      TH2D hHT(TString::Format("hHT_%d",k).Data(),TString::Format("zv: %.1f - reco: %.1f; zc (mm); yc (mm)",zv.GetVal(),z_int[k]).Data(),100,-2000,2000,100,-2000,2000);
      
      TH1D harctg(TString::Format("harctg_%d",k).Data(),TString::Format("zv: %.1f - reco: %.1f; #phi (rad)",zv.GetVal(),z_int[k]).Data(), 200, -TMath::Pi(), TMath::Pi());
      
      std::vector<double> vu;
      std::vector<double> vv;
      
      for(unsigned int j = 0; j < digits->size(); j++)
      {
        if(digits->at(j).hor == 1)
        {
          
          u = (digits->at(j).z - z_int[k]);
          v = (y_int[k] - digits->at(j).y);
          d = (u*u + v*v);
          
          u /=d;
          v /=d;
          
          vu.push_back(u);
          vv.push_back(v);
          
          double y = digits->at(j).y - y_int[k];
          double z = digits->at(j).z - z_int[k];
          double R = z*z + y*y; 
          
          for(int kk = 0; kk < hHT.GetNbinsX(); kk++)
          {
            zz = hHT.GetXaxis()->GetBinCenter(kk+1);
            yy = (-z*zz + 0.5*R)/y;
            hHT.Fill(zz, yy);
          }
          
          harctg.Fill(TMath::ATan2(v,u));
        }
      } 
      TGraph guv(vu.size(), vu.data(),vv.data());
      guv.SetName(TString::Format("huv_%d",k).Data());
      guv.SetTitle(TString::Format("zv: %.1f - reco: %.1f; u; v",zv.GetVal(),z_int[k]).Data()); 
      
      vguv.push_back(std::move(guv));
      vharctg.push_back(std::move(harctg));
      vhHT.push_back(std::move(hHT));
    }
    */
    
    if(display)
    {
      
      for (std::map<int,int>::iterator it=regionsX.begin(); it!=regionsX.end(); ++it)
      {
        for(int k = 0; k < it->second; k++)
        {
          h0trX.SetBinContent(it->first + k, hmultX.GetBinContent(it->first) == 0 ? 1 : 0);
          h1trX.SetBinContent(it->first + k, hmultX.GetBinContent(it->first) == 1 ? 1 : 0);
          hmtrX.SetBinContent(it->first + k, hmultX.GetBinContent(it->first) == 2 ? 1 : 0);
        }
      }
      
      for (std::map<int,int>::iterator it=regionsY.begin(); it!=regionsY.end(); ++it)
      {
        for(int k = 0; k < it->second; k++)
        {
          h0trY.SetBinContent(it->first + k, hmultY.GetBinContent(it->first) == 0 ? 1 : 0);
          h1trY.SetBinContent(it->first + k, hmultY.GetBinContent(it->first) == 1 ? 1 : 0);
          hmtrY.SetBinContent(it->first + k, hmultY.GetBinContent(it->first) == 2 ? 1 : 0);
        }
      }
      
      for (std::map<int,int>::iterator it=regions.begin(); it!=regions.end(); ++it)
      {
        for(int k = 0; k < it->second; k++)
        {
          h0tr.SetBinContent(it->first + k, hmult.GetBinContent(it->first) == 0 ? 1 : 0);
          h1tr.SetBinContent(it->first + k, hmult.GetBinContent(it->first) == 1 ? 1 : 0);
          hmtr.SetBinContent(it->first + k, hmult.GetBinContent(it->first) == 2 ? 1 : 0);
        }
      }
      
      // display of the event
      show(i,true,false,true);
      
      TCanvas* cev = (TCanvas*) gROOT->FindObject("cev");
      
      // histograms with regions
      TLine lv1;
      lv1.SetLineStyle(2);
      lv1.SetLineColor(kBlack);
      lv1.SetLineWidth(2);
      lv1.SetX1(zv.GetVal());
      lv1.SetX2(zv.GetVal());
      
      TLine lv2;
      lv2.SetLineStyle(2);
      lv2.SetLineColor(kBlack);
      lv2.SetLineWidth(2);
      lv2.SetX1(zv.GetVal());
      lv2.SetX2(zv.GetVal());
      
      TLine* lvx;
      TLine* lvy;
      
      c.Clear();
      c.Divide(1,2);
      
      h0tr.Scale((hrmsX.GetMaximum() <= 0 && hrmsY.GetMaximum() <= 0) ? 1 : std::max(hrmsX.GetMaximum(),hrmsY.GetMaximum()));
      h1tr.Scale((hrmsX.GetMaximum() <= 0 && hrmsY.GetMaximum() <= 0) ? 1 : std::max(hrmsX.GetMaximum(),hrmsY.GetMaximum()));
      hmtr.Scale((hrmsX.GetMaximum() <= 0 && hrmsY.GetMaximum() <= 0) ? 1 : std::max(hrmsX.GetMaximum(),hrmsY.GetMaximum()));
      
      c.cd(1);
      h0tr.Draw("B");
      h1tr.Draw("BSAME");
      hmtr.Draw("BSAME");
      hrmsX.Draw("HISTSAME");
      lv1.SetY1(hrmsX.GetMinimum() == -1 ? 0 : hrmsX.GetMinimum());
      lv1.SetY2((hrmsX.GetMaximum() <= 0 && hrmsY.GetMaximum() <= 0) ? 1 : std::max(hrmsX.GetMaximum(),hrmsY.GetMaximum()));
      lv1.Draw();
      
      c.cd(2);
      h0tr.Draw("B");
      h1tr.Draw("BSAME");
      hmtr.Draw("BSAME");
      hrmsY.Draw("HISTSAME");
      lv2.SetY1(hrmsY.GetMinimum() == -1 ? 0 : hrmsY.GetMinimum());
      lv2.SetY2((hrmsX.GetMaximum() <= 0 && hrmsY.GetMaximum() <= 0) ? 1 : std::max(hrmsX.GetMaximum(),hrmsY.GetMaximum()));
      lv2.Draw();
      
      // save to pdf
      if(save2pdf)
      {
        if(cev)
        {
          cev->SaveAs("rms.pdf");
        }
        
        c.SaveAs("rms.pdf");
      
        c.Clear();
        c.cd();
        
        for(unsigned int k = 0; k < vharctg.size(); k++)
        {
          vharctg[k].Draw();
          c.SaveAs("rms.pdf");
          vguv[k].Draw("ap*");
          c.SaveAs("rms.pdf");
          vhHT[k].Draw("colz");
          c.SaveAs("rms.pdf");
        }
      }
      
      c.Clear();
      c.Divide(2,1);
      
      TMarker* m;
      TLine* l;
      TBox* b;
      
      c.cd(1)->DrawFrame(21500,-5000,26500,0);
      c.cd(1);
      m = new TMarker(zvtx_reco,yvtx_reco,8);
      m->Draw();
      m = new TMarker(zv.GetVal(),yv.GetVal(),21);
      m->SetMarkerColor(kBlue);
      m->Draw();
      for(unsigned int j = 0; j < digits->size(); j++)
      {
        if(digits->at(j).hor == 1)
        {
          m = new TMarker(digits->at(j).z,digits->at(j).y,1);
          m->Draw();
        }
      }
      
      for(int j = 0; j < hmeanY.GetNbinsX(); j++)
      {
        if(hmeanY.GetBinContent(j+1) != -1)
        {
          l = new TLine(hmeanY.GetXaxis()->GetBinLowEdge(j+1),hmeanY.GetBinContent(j+1),hmeanY.GetXaxis()->GetBinUpEdge(j+1),hmeanY.GetBinContent(j+1));
          l->SetLineColor(kRed);
          l->Draw();
        }
        if(hrmsY.GetBinContent(j+1) > 0.)
        {
          b = new TBox(hmeanY.GetXaxis()->GetBinLowEdge(j+1),hmeanY.GetBinContent(j+1)+hrmsY.GetBinContent(j+1),hmeanY.GetXaxis()->GetBinUpEdge(j+1),hmeanY.GetBinContent(j+1)-hrmsY.GetBinContent(j+1));
          b->SetFillColorAlpha(kBlue,0.15);
          b->SetLineColor(kBlue);
          b->Draw();
        }
      }
      
      c.cd(2)->DrawFrame(21500,-2500,26500,2500);
      c.cd(2);
      m = new TMarker(zvtx_reco,xvtx_reco,8);
      m->Draw();
      m = new TMarker(zv.GetVal(),xv.GetVal(),21);
      m->SetMarkerColor(kBlue);
      m->Draw();
      for(unsigned int j = 0; j < digits->size(); j++)
      {
        if(digits->at(j).hor != 1)
        {
          m = new TMarker(digits->at(j).z,digits->at(j).x,1);
          m->Draw();
        }
      }
      
      for(int j = 0; j < hmeanX.GetNbinsX(); j++)
      {
        if(hmeanX.GetBinContent(j+1) != -1)
        {
          l = new TLine(hmeanX.GetXaxis()->GetBinLowEdge(j+1),hmeanX.GetBinContent(j+1),hmeanX.GetXaxis()->GetBinUpEdge(j+1),hmeanX.GetBinContent(j+1));
          l->SetLineColor(kRed);
          l->Draw();
        }
        if(hrmsX.GetBinContent(j+1) > 0.)
        {
          b = new TBox(hmeanX.GetXaxis()->GetBinLowEdge(j+1),hmeanX.GetBinContent(j+1)+hrmsX.GetBinContent(j+1),hmeanX.GetXaxis()->GetBinUpEdge(j+1),hmeanX.GetBinContent(j+1)-hrmsX.GetBinContent(j+1));
          b->SetFillColorAlpha(kBlue,0.15);
          b->SetLineColor(kBlue);
          b->Draw();
        }
      }
      
      c.SaveAs("rms.pdf");
    }
    
    if(save2root)
    {
      fd->Add(&hmeanX);
      fd->Add(&hnX);
      fd->Add(&hrmsX);
      
      fd->Add(&hmeanY);
      fd->Add(&hnY);
      fd->Add(&hrmsY);
      
      fd->Add(&xv);
      fd->Add(&yv);
      fd->Add(&zv);
      /*
      for(unsigned int k = 0; k < vharctg.size(); k++)
      {
        fd->Add(&vharctg[k]);
        fd->Add(&vguv[k]);
        fd->Add(&vhHT[k]);
      }*/
      
      fout.cd();
      fd->Write();
    }
    
    tv.Fill();
  }
  
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
  
  vmeanX /= nall;
  vmeanY /= nall;
  vmeanZ /= nall;
  vmean3D /= nall;
  
  vrmsX /= nall;
  vrmsY /= nall;
  vrmsZ /= nall;
  vrms3D /= nall;
  
  vrmsX -= vmeanX*vmeanX;
  vrmsY -= vmeanY*vmeanY;
  vrmsZ -= vmeanZ*vmeanZ;
  vrms3D -= vmean3D*vmean3D;
  
  vrmsX = sqrt(vrmsX);
  vrmsY = sqrt(vrmsY);
  vrmsZ = sqrt(vrmsZ);
  vrms3D = sqrt(vrms3D);
  
  std::cout << std::setw(10) << "minreg" << std::setw(10) << "epsilon" << std::setw(15) << "wideMultReg" << std::setw(15) << "gefficiency" << std::setw(15) << "efficiency" << std::setw(10) << "rmsX" << std::setw(10) << "rmsY" << std::setw(10) << "rmsZ" << std::setw(10) << "rms3D" << std::endl;
  std::cout << std::setw(10) << minreg << std::setw(10) << epsilon << std::setw(15) << wideMultReg << std::setw(15) << double(nall)/nev << std::setw(15) << double(nok)/nall << std::setw(10) << vrmsX << std::setw(10) << vrmsY << std::setw(10) << vrmsZ << std::setw(10) << vrms3D << std::endl;  
  
  c.SaveAs("rms.pdf)");
  fout.cd();
  tv.Write();
  fout.Close();
}

void findvtx()
{
  /*
  VtxFinding(1, 0.0, false);
  VtxFinding(1, 0.00000001, false);
  VtxFinding(1, 0.0000001, false);
  VtxFinding(1, 0.000001, false);
  VtxFinding(1, 0.00001, false);
  VtxFinding(1, 0.0001, false);
  VtxFinding(1, 0.001, false);
  VtxFinding(1, 0.01, false);
  VtxFinding(1, 0.1, false);
  VtxFinding(1, 0.2, false);
  VtxFinding(1, 0.3, false);
  VtxFinding(1, 0.4, false);
  VtxFinding(1, 0.5, false);
  VtxFinding(1, 0.6, false);
  VtxFinding(1, 0.7, false);
  VtxFinding(1, 0.8, false);
  VtxFinding(1, 0.9, false);
  VtxFinding(1, 1., false);
  VtxFinding(1, 2., false);
  VtxFinding(1, 0.5, true);
  VtxFinding(2, 0.5, true);
  VtxFinding(3, 0.5, true);
  VtxFinding(4, 0.5, true);
  VtxFinding(5, 0.5, true);
  VtxFinding(6, 0.5, true);
  VtxFinding(7, 0.5, true);
  VtxFinding(8, 0.5, true);
  VtxFinding(9, 0.5, true);
  VtxFinding(10, 0.5, true);
  VtxFinding(2, 0.5, false);
  VtxFinding(3, 0.5, false);
  VtxFinding(4, 0.5, false);
  VtxFinding(5, 0.5, false);
  VtxFinding(6, 0.5, false);
  VtxFinding(7, 0.5, false);
  VtxFinding(8, 0.5, false);
  VtxFinding(9, 0.5, false);
  VtxFinding(10, 0.5, false);*/
  
  VtxFinding(1, 0.5, false);
  
  gStyle->SetPalette(53);
  
  TFile f("vtx.root");
  TTree* tv = (TTree*) f.Get("tv");
  TCanvas c;
  TH1D* h;
  
  TCut vtxReco = "VtxMono == 1 || VtxMulti == 1";
  TCut monoprong = "VtxMono == 1";
  TCut multiprong = "VtxMulti == 1";
  TCut notLAr = "zvtx_true > 22100";
  
  // x VS z
  tv->Draw("xvtx_true:zvtx_true>>hxVSz","","");
  TH2D* hxVSz = (TH2D*) gROOT->FindObject("hxVSz");
  hxVSz->SetTitle(";Z (mm); X (mm); #Deltar_{3D} (mm)");
  hxVSz->Draw();
  c.SaveAs("vtx.pdf(");
  
  // y VS z
  tv->Draw("yvtx_true:zvtx_true>>hyVSz","","");
  TH2D* hyVSz = (TH2D*) gROOT->FindObject("hyVSz");
  hyVSz->SetTitle(";Z (mm); Y (mm); #Deltar_{3D} (mm)");
  hyVSz->Draw();
  c.SaveAs("vtx.pdf");
  
  // x VS z reco
  tv->Draw("xvtx_reco:zvtx_reco>>hxVSz_reco",vtxReco,"");
  TH2D* hxVSz_reco = (TH2D*) gROOT->FindObject("hxVSz_reco");
  hxVSz_reco->SetTitle(";Z (mm); X (mm); #Deltar_{3D} (mm)");
  hxVSz_reco->Draw();
  c.SaveAs("vtx.pdf");
  
  // y VS z reco
  tv->Draw("yvtx_reco:zvtx_reco>>hyVSz_reco",vtxReco,"");
  TH2D* hyVSz_reco = (TH2D*) gROOT->FindObject("hyVSz_reco");
  hyVSz_reco->SetTitle(";Z (mm); Y (mm); #Deltar_{3D} (mm)");
  hyVSz_reco->Draw();
  c.SaveAs("vtx.pdf");
  
  // x VS z VS mean dr
  tv->Draw("xvtx_reco:zvtx_reco>>hxVSz_n(25,21500,26500,25,-2000,2000)","VtxMulti == 1","colz");
  TH2D* hxVSz_n = new TH2D(*((TH2D*) gROOT->FindObject("hxVSz_n")));
  tv->Draw("xvtx_reco:zvtx_reco>>hxVSz_dv(25,21500,26500,25,-2000,2000)","(VtxMulti == 1)*sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2)+pow(zvtx_reco-zvtx_true,2))","colz");
  TH2D* hxVSz_dv = new TH2D(*((TH2D*) gROOT->FindObject("hxVSz_dv")));
  TH2D hxVSz_meandv = (*hxVSz_dv)/(*hxVSz_n);
  hxVSz_meandv.SetStats(false);
  hxVSz_meandv.SetTitle(";Z (mm); Y (mm); #Deltar_{3D} (mm)");
  hxVSz_meandv.Draw("colz");
  c.SaveAs("vtx.pdf");
  
  // y VS z VS mean dr
  tv->Draw("yvtx_reco:zvtx_reco>>hyVSz_n(25,21500,26500,25,-5000,0)","VtxMulti == 1","colz");
  TH2D* hyVSz_n = new TH2D(*((TH2D*) gROOT->FindObject("hyVSz_n")));
  tv->Draw("yvtx_reco:zvtx_reco>>hyVSz_dv(25,21500,26500,25,-5000,0)","(VtxMulti == 1)*sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2)+pow(zvtx_reco-zvtx_true,2))","colz");
  TH2D* hyVSz_dv = new TH2D(*((TH2D*) gROOT->FindObject("hyVSz_dv")));
  TH2D hyVSz_meandv = (*hyVSz_dv)/(*hyVSz_n);
  hyVSz_meandv.SetStats(false);
  hyVSz_meandv.SetTitle(";Z (mm); Y (mm); #Deltar_{3D} (mm)");
  hyVSz_meandv.Draw("colz");
  c.SaveAs("vtx.pdf");
  
  
  c.SetLogy(true);
  // dx  
  tv->Draw("xvtx_reco-xvtx_true>>dvx_multi(100,-4000,4000)",multiprong,"");
  TH1D* dvx_multi = new TH1D(*((TH1D*) gROOT->FindObject("dvx_multi")));
  dvx_multi->SetFillColorAlpha(kBlue,0.3);
  
  tv->Draw("xvtx_reco-xvtx_true>>dvx_mono(100,-4000,4000)",monoprong,"");
  TH1D* dvx_mono = new TH1D(*((TH1D*) gROOT->FindObject("dvx_mono")));
  dvx_mono->SetFillColorAlpha(kRed,0.3);
  
  THStack dvx_all("dvx_all",";xreco-xtrue (mm)");
  dvx_all.Add(dvx_multi);
  dvx_all.Add(dvx_mono);
  
  TLegend leg(0.79,0.79,0.99,0.99);
  leg.AddEntry(dvx_multi,"mtr","fl");
  leg.AddEntry(dvx_mono,"1tr","fl");
  
  dvx_all.Draw("nostack");
  leg.Draw();
  c.SaveAs("vtx.pdf");
  
  // dy
  tv->Draw("yvtx_reco-yvtx_true>>dvy_multi(100,-5000,5000)",multiprong,"");
  TH1D* dvy_multi = new TH1D(*((TH1D*) gROOT->FindObject("dvy_multi")));
  dvy_multi->SetFillColorAlpha(kBlue,0.3);
  
  tv->Draw("yvtx_reco-yvtx_true>>dvy_mono(100,-5000,5000)",monoprong,"");
  TH1D* dvy_mono = new TH1D(*((TH1D*) gROOT->FindObject("dvy_mono")));
  dvy_mono->SetFillColorAlpha(kRed,0.3);
  
  THStack dvy_all("dvy_all",";yreco-ytrue (mm)");
  dvy_all.Add(dvy_multi);
  dvy_all.Add(dvy_mono);
  
  dvy_all.Draw("nostack");
  leg.Draw();
  c.SaveAs("vtx.pdf");
  
  //dz
  tv->Draw("zvtx_reco-zvtx_true>>dvz_multi(100,-4000,4000)",multiprong,"");
  TH1D* dvz_multi = new TH1D(*((TH1D*) gROOT->FindObject("dvz_multi")));
  dvz_multi->SetFillColorAlpha(kBlue,0.3);
  
  tv->Draw("zvtx_reco-zvtx_true>>dvz_mono(100,-4000,4000)",monoprong,"");
  TH1D* dvz_mono = new TH1D(*((TH1D*) gROOT->FindObject("dvz_mono")));
  dvz_mono->SetFillColorAlpha(kRed,0.3);
  
  THStack dvz_all("dvz_all",";zreco-ztrue (mm)");
  dvz_all.Add(dvz_multi);
  dvz_all.Add(dvz_mono);
  
  dvz_all.Draw("nostack");
  leg.Draw();
  c.SaveAs("vtx.pdf");
  
  // dxy
  tv->Draw("sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2))>>dvabsxy_multi(200,0,5000)",multiprong,"");
  TH1D* dvabsxy_multi = new TH1D(*((TH1D*) gROOT->FindObject("dvabsxy_multi")));
  dvabsxy_multi->SetFillColorAlpha(kBlue,0.3);
  
  tv->Draw("sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2))>>dvabsxy_mono(200,0,5000)",monoprong,"");
  TH1D* dvabsxy_mono = new TH1D(*((TH1D*) gROOT->FindObject("dvabsxy_mono")));
  dvabsxy_mono->SetFillColorAlpha(kRed,0.3);
  
  THStack dvabsxy_all("dvabsxy_all",";#Deltar_{xy} (mm)");
  dvabsxy_all.Add(dvabsxy_multi);
  dvabsxy_all.Add(dvabsxy_mono);
  
  dvabsxy_all.Draw("nostack");
  leg.Draw();
  c.SaveAs("vtx.pdf");
  
  c.SetLogy(false);
  
  // dxy cumulative  
  TH1D hintxy(*dvabsxy_multi);
  TH1D hintxy_mono(*dvabsxy_multi);
  TH1D hintxy_multi(*dvabsxy_multi);
  hintxy.SetFillColor(0);
  hintxy.SetStats(false);
  hintxy_multi.SetFillColor(0);
  hintxy_multi.SetLineColor(kBlue);
  hintxy_multi.SetLineStyle(2);
  hintxy_multi.SetStats(false);
  hintxy_mono.SetFillColor(0);
  hintxy_mono.SetLineColor(kRed);
  hintxy_mono.SetLineStyle(2);
  hintxy_mono.SetStats(false);
  for(int i = 0; i < hintxy.GetNbinsX(); i++)
  {
    hintxy.SetBinContent(i+1,(dvabsxy_multi->Integral(1,i+1)+dvabsxy_mono->Integral(1,i+1))/(dvabsxy_multi->Integral()+dvabsxy_mono->Integral()));
    hintxy_mono.SetBinContent(i+1,dvabsxy_mono->Integral(1,i+1)/dvabsxy_mono->Integral());
    hintxy_multi.SetBinContent(i+1,dvabsxy_multi->Integral(1,i+1)/dvabsxy_multi->Integral());
  }
  c.DrawFrame(0,0,1000,1,"fraction of events with #Deltar_{xy} < #Deltar^{thr}_{xy};#Deltar^{thr}_{xy} (mm)");
  hintxy.Draw("csame");
  hintxy_mono.Draw("csame");
  hintxy_multi.Draw("csame");
  
  TLegend leg2(0.69,0.11,0.89,0.31);
  leg2.AddEntry(&hintxy,"all vtx","l");
  leg2.AddEntry(&hintxy_multi,"mtr","l");
  leg2.AddEntry(&hintxy_mono,"1tr","l");
  leg2.Draw();
  
  
  c.SaveAs("vtx.pdf");
  
  c.SetLogy(true);
  
  // dz
  tv->Draw("abs(zvtx_reco-zvtx_true)>>dvabsz_multi(100,0,4000)",multiprong,"");
  TH1D* dvabsz_multi = new TH1D(*((TH1D*) gROOT->FindObject("dvabsz_multi")));
  dvabsz_multi->SetFillColorAlpha(kBlue,0.3);
  
  tv->Draw("abs(zvtx_reco-zvtx_true)>>dvabsz_mono(100,0,4000)",monoprong,"");
  TH1D* dvabsz_mono = new TH1D(*((TH1D*) gROOT->FindObject("dvabsz_mono")));
  dvabsz_mono->SetFillColorAlpha(kRed,0.3);
  
  THStack dvabsz_all("dvabsz_all",";|zreco-ztrue| (mm)");
  dvabsz_all.Add(dvabsz_multi);
  dvabsz_all.Add(dvabsz_mono);
  
  dvabsz_all.Draw("nostack");
  leg.Draw();
  c.SaveAs("vtx.pdf");
  
  c.SetLogy(false);
  // dz cumulative
  TH1D hintz(*dvabsz_multi);
  TH1D hintz_mono(*dvabsz_multi);
  TH1D hintz_multi(*dvabsz_multi);
  hintz.SetFillColor(0);
  hintz.SetStats(false);
  hintz_mono.SetFillColor(0);
  hintz_mono.SetLineColor(kRed);
  hintz_mono.SetLineStyle(2);
  hintz_mono.SetStats(false);
  hintz_multi.SetFillColor(0);
  hintz_multi.SetLineColor(kBlue);
  hintz_multi.SetLineStyle(2);
  hintz_multi.SetStats(false);
  for(int i = 0; i < hintz.GetNbinsX(); i++)
  {
    hintz.SetBinContent(i+1,(dvabsz_multi->Integral(1,i+1)+dvabsz_mono->Integral(1,i+1))/(dvabsz_multi->Integral()+dvabsz_mono->Integral()));
    hintz_mono.SetBinContent(i+1,dvabsz_mono->Integral(1,i+1)/dvabsz_mono->Integral());
    hintz_multi.SetBinContent(i+1,dvabsz_multi->Integral(1,i+1)/dvabsz_multi->Integral());
  }
  c.DrawFrame(0,0,1000,1,"fraction of events with |zreco-ztrue| < |zreco-ztrue|_{thr};|zreco-ztrue|_{thr} (mm)");
  hintz.Draw("csame");
  hintz_mono.Draw("csame");
  hintz_multi.Draw("csame");
  leg2.Draw();
  c.SaveAs("vtx.pdf");
  c.SetLogy(true);
  
  // 3D
  tv->Draw("sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2)+pow(zvtx_reco-zvtx_true,2))>>dv3D_multi(250,0,5500)",multiprong,"");
  TH1D* dv3D_multi = new TH1D(*((TH1D*) gROOT->FindObject("dv3D_multi")));
  dv3D_multi->SetFillColorAlpha(kBlue,0.3);
  
  tv->Draw("sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2)+pow(zvtx_reco-zvtx_true,2))>>dv3D_mono(250,0,5500)",monoprong,"");
  TH1D* dv3D_mono = new TH1D(*((TH1D*) gROOT->FindObject("dv3D_mono")));
  dv3D_mono->SetFillColorAlpha(kRed,0.3);
  
  THStack dv3D_all("dv3D_all",";#Deltar_{3D} (mm)");
  dv3D_all.Add(dv3D_multi);
  dv3D_all.Add(dv3D_mono);
  
  dv3D_all.Draw("nostack");
  leg.Draw();
  c.SaveAs("vtx.pdf");
  
  c.SetLogy(false);
  
  // 3D cumulative
  TH1D hint3D(*dv3D_multi);
  TH1D hint3D_mono(*dv3D_multi);
  TH1D hint3D_multi(*dv3D_multi);
  hint3D.SetFillColor(0);
  hint3D.SetStats(false);
  hint3D_mono.SetFillColor(0);
  hint3D_mono.SetLineColor(kRed);
  hint3D_mono.SetLineStyle(2);
  hint3D_mono.SetStats(false);
  hint3D_multi.SetFillColor(0);
  hint3D_multi.SetLineColor(kBlue);
  hint3D_multi.SetLineStyle(2);
  hint3D_multi.SetStats(false);
  for(int i = 0; i < hint3D.GetNbinsX(); i++)
  {
    hint3D.SetBinContent(i+1,(dv3D_multi->Integral(1,i+1)+dv3D_mono->Integral(1,i+1))/(dv3D_multi->Integral()+dv3D_mono->Integral()));
    hint3D_mono.SetBinContent(i+1,dv3D_mono->Integral(1,i+1)/dv3D_mono->Integral());
    hint3D_multi.SetBinContent(i+1,dv3D_multi->Integral(1,i+1)/dv3D_multi->Integral());
  }
  c.DrawFrame(0,0,1000,1,"fraction of events with #Deltar_{3D} < #Deltar^{thr}_{3D};#Deltar^{thr}_{3D} (mm)");
  hint3D.Draw("csame");
  hint3D_mono.Draw("csame");
  hint3D_multi.Draw("csame");
  leg2.Draw();
  c.SaveAs("vtx.pdf");
  
  // xy
  tv->Draw("yvtx_reco-yvtx_true:xvtx_reco-xvtx_true>>dvxy(500,-1000,1000,500,-1000,1000)",vtxReco,"colz");
  TH2D* h2D = (TH2D*) gROOT->FindObject("dvxy");
  h2D->SetTitle(";xreco-xtrue (mm);yreco-ytrue (mm)");
  h2D->Draw("colz");
  TPad p("p","",0.12,0.55,0.45,0.88);
  p.cd();
  p.DrawFrame(-100,-100,100,100);
  h2D->SetStats(false);
  h2D->Draw("colsame");
  c.cd();
  p.Draw();
  c.SaveAs("vtx.pdf)");
  
  // dr 2D (dr < 50 cm)
  int nn = tv->Draw("sqrt(pow(yvtx_reco-yvtx_true,2)+pow(zvtx_reco-zvtx_true,2))>>h","sqrt(pow(yvtx_reco-yvtx_true,2)+pow(zvtx_reco-zvtx_true,2))<500",vtxReco);
  h = (TH1D*) gROOT->FindObject("h");
  h->SetTitle(";#Deltar_{zy} (mm)");
  h->SetFillColor(38);
  h->Draw();
  std::cout << "dr_zy < 500 mm => " << double(nn)/tv->GetEntries() << " rms: " << h->GetRMS() << std::endl;
  
  // dr_zy (dr_zy < 5 cm)
  nn = tv->Draw("sqrt(pow(yvtx_reco-yvtx_true,2)+pow(zvtx_reco-zvtx_true,2))>>h","sqrt(pow(yvtx_reco-yvtx_true,2)+pow(zvtx_reco-zvtx_true,2))<50",vtxReco);
  h = (TH1D*) gROOT->FindObject("h");
  h->SetTitle(";#Deltar_{zy} (mm)");
  h->SetFillColor(38);
  h->Draw();
  std::cout << "dr_zy < 50 mm => " << double(nn)/tv->GetEntries() << " rms: " << h->GetRMS() << std::endl;
  
  // dr_3D (dr_3D < 50 cm)
  nn = tv->Draw("sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2)+pow(zvtx_reco-zvtx_true,2))>>h","sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2)+pow(zvtx_reco-zvtx_true,2))<500",vtxReco);
  h = (TH1D*) gROOT->FindObject("h");
  h->SetTitle(";#Deltar_{3D} (mm)");
  h->SetFillColor(38);
  h->Draw();
  std::cout << "dr_3D < 500 mm => " << double(nn)/tv->GetEntries() << " rms: " << h->GetRMS() << std::endl;
  
  nn = tv->Draw("sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2)+pow(zvtx_reco-zvtx_true,2))>>h","sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2)+pow(zvtx_reco-zvtx_true,2))<50",vtxReco);
  h = (TH1D*) gROOT->FindObject("h");
  h->SetTitle(";#Deltar_{3D} (mm)");
  h->SetFillColor(38);
  h->Draw();
  std::cout << "dr_3D < 50 mm => " << double(nn)/tv->GetEntries() << " rms: " << h->GetRMS() << std::endl;
  c.SetLogy(false);
}