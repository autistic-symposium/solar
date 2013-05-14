/********************************************************
 * 
 * For informatiob about the code, please read README 
 * 
 * steinkirch, spring/2012
 * 
 *******************************************************/

#include "TVirtualFFT.h"
#include <iostream>
using namespace std;


/**********************************************
 * 
 * Main Function starts here
 * 
 * *******************************************/

void DataRedux( double  db = -25., char * outfile = "output/sat_power.root", char * name_file = "SAT",  double  number_files = 5,   double xnbins = 1000,  double xbin_min = 0,  double xbin_max = 15,  double ynbins = 1000,  double ybin_min = -6,  double ybin_max = -8 )
//void DataRedux( double  db = -25., char * outfile = "output/sun_power.root", char * name_file = "SUN",  double  number_files = 5,  int xnbins = 1000,  double xbin_min = 0,  double xbin_max = 20.,  int ynbins = 50,  double ybin_min = -6.8,  double ybin_max = -7.4)
{

/**********************************************
* Global Variables
* ********************************************/

  const int numberfiles = number_files;
  const char namefile[1000] = name_file;
  const char output[1000] = outfile;
  
  char name[1000];
  char title[10000];
  
  TH2D* h[numberfiles];
  TH1D* h1D[numberfiles];
  TH1D* h1Dcut[numberfiles];
  
  ifstream f[numberfiles]; 
  
  TFile *output_root = new TFile(output,"RECREATE");


/***********************************************
 Opening ascii files and creating the root file 
************************************************/
 for(int i = 0; i < number_files; i++)
 {
   sprintf(name,"%s/%s%i.txt",namefile,name_file,i+1);
   f[i].open(name);
 }
 
  
/************************************************
  Filling histograms
*************************************************/

  for(int i = 0; i < numberfiles; i++)
  {
    double t, v, p, nlines=0, first;
    
    sprintf(name,"%s_%i", namefile, i+1);
    sprintf(title,"%s - Baseline %i", namefile, i+1);
    h[i] = new TH2D(name,title,xnbins, xbin_min, xbin_max, ynbins, ybin_min, ybin_max);

    sprintf(name,"%s_1D_%i", namefile, i+1);
    h1D[i] = new TH1D(name,title,xnbins, xbin_min, xbin_max);    
    
    sprintf(name,"%s_1D_cut_%i", namefile, i+1);
    h1Dcut[i] = new TH1D(name,title,xnbins, xbin_min, xbin_max);
    
    while (1) {
      f[i] >> t >> v ; 
      if (!f[i].good()) break;
      if(nlines==0) first = t; 
      p = VToP(v,db);
      h[i]->Fill(t-first,p);
      h1D[i]->Fill(t-first,p);
      Cut(h1Dcut[i],i,t-first,p);
      nlines++; 
    }
  }
  
 

/****************************************
 * Drawing and Saving 
 *****************************************/
  
  sprintf(name,"%s - P(t) vs t",namefile);
  TCanvas *c = new TCanvas("c",name,900,600);  
  
  TPad *pad1 = new TPad("pad1","",0.03,0.62,0.50,0.92,21);
  pad1->Draw();
  pad1->SetGridy();
  c->Clear();
  c->Divide(2,5);
 
  Style_Cosmetics();
  
  int canvas_number = 1;
  for(i=0; i<numberfiles; i++){
    c->SetGridy();
    // Plot Power vc t
    gStyle->SetOptLogy(1);
    c->cd(canvas_number);
      
      Style_Cosmetics();
    Power_Cosmetics(h1D[i], i);
    c->Update();

    // Plot FFT Magnitude
    gStyle->SetOptLogy(0);
    c->cd(canvas_number+1);
    TH1 *hm =0;
    FFT_Mag(h1Dcut[i],hm, i, namefile);
    c->Update();
    
    // Increment number camvas
    canvas_number = canvas_number + 2;
  }

  sprintf(name,"plots/%s_power.ps",namefile);
  c->SaveAs(name);
   
  output_root->Write();  
  f->close();  
}



/*********************************************
 * 
 *                Aux Functions
 * 
 **********************************************/

double VToP(const double v, const double db) {
  double pow = v*db;
  double p = 10**((pow-30)/10);
  return p;
}

void Cut(TH1D * h,int  i, double  t, double  p){
  if(i==0){
    if((p<15) ||(p>4)){
          h->Fill(t,p);
    }
  }
  else if (i==1){
    if((p<15) || (p>4)){
          h->Fill(t,p);
    }
  }
  else if (i==2){
    if((p<15) || (p>4)){
          h->Fill(t,p);
    }
  }
  else if (i==3){
    if((p<15) && (p>4)){
          h->Fill(t,p);
    }
  }   
  else if (i==4){
    if((p<15) || (p>4)){
          h->Fill(t,p);
    }
  }
}

void Power_Cosmetics(TH1D *h, int  i){
    h->SetMarkerColor(i+1);  
    h->SetMarkerSize(0.5);
    h->SetLineStyle(1);
    h->SetLineColor(1);
    h->SetTitleSize(0.03, "x");
    h->SetXTitle("Time (s)"); 
    h->SetYTitle("log(Intensity)[digit]");
    h->GetXaxis()->CenterTitle(1);
    h->GetYaxis()->CenterTitle(1);
    h->GetXaxis()->SetLabelSize(0.05); 
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleOffset(0.9);
    float xmax = h->GetXaxis()->GetXmax();
    float xmin = h->GetXaxis()->GetXmin();
    h->GetXaxis()->SetLimits(xmin,xmax);
    h->GetYaxis()->SetMoreLogLabels();
    h->Draw("P");
  //  h->GetXaxis()->SetLabelOffset(0.15);
      TLine *l1 = new TLine();
  l1->SetLineWidth(0.5);
  l1->SetLineStyle(2);
  l1->SetLineColor(2);
  l1->DrawLine(3.05, 0, 3.05, 10);
  l1->DrawLine(5, 0, 5, 14);
}


double fpeaks(double *x, double *par) {
   double result = par[0] + par[1]*x[0];
   npeaks = 5;
   for (int p=0;p<npeaks;p++) {
      double norm  = par[3*p+2];
      double mean  = par[3*p+3];
      double sigma = par[3*p+4];
      result += norm*TMath::Gaus(x[0],mean,sigma);
   }
   return result;
}

void FindPeak(TH1 *hm, int * i, char * namefile){
  int np =5, p, max = 0;
  int npeaks = TMath::Abs(np);
  double par[3000];
  par[0] = 0.8;
  par[1] = -0.6/1000;

  for (p=0;p<npeaks;p++) {
    par[3*p+2] = 1;
    par[3*p+3] = 10+gRandom->Rndm()*980;
    par[3*p+4] = 3+2*gRandom->Rndm();
  }
 
  TSpectrum *s = new TSpectrum(2*npeaks,1);
  int nfound = s->Search(hm,2,"",0.10);
  printf("Found %d candidate peaks to fit\n",nfound);

  TH1 *hb = s->Background(hm,20,"same");
  if (np <0) return;

  // loope over peaks
  TF1 *fline = new TF1("fline","pol1",0,1000);
  hm->Fit("fline","qn");
  par[0] = fline->GetParameter(0);
  par[1] = fline->GetParameter(1);
  npeaks = 0;
  float *xpeaks = s->GetPositionX();
  for (p=0;p<nfound;p++) {
    float xp = xpeaks[p];
    int bin = hm->GetXaxis()->FindBin(xp);
    float yp = hm->GetBinContent(bin);
    if (yp-TMath::Sqrt(yp) < fline->Eval(xp)) continue;
    par[3*npeaks+2] = yp;
    par[3*npeaks+3] = xp;
    par[3*npeaks+4] = 3;
    npeaks++;
  }
  printf("Found %d useful peaks to fit\n",npeaks); 
  printf("Now fitting: Be patient\n");
  if (max < npeaks) max = npeaks;
  TF1 *fit = new TF1("fit",fpeaks,0,1000,2+3*npeaks);
  TVirtualFitter::Fitter(hm,10+3*npeaks);
  fit->SetParameters(par);
  fit->SetNpx(1000);
  hm->Fit("fit"); 
}


void FFT_Mag(TH1D *h, TH1 *hm, int *  i,char * namefile){
  TVirtualFFT::SetTransform(0);
  char name[10000];
  sprintf(name,"%s_MAG_%i", namefile,i);
  TH1F *h2 = (TH1F*)h->Clone("h2");
   hm= h2->FFT(hm,  name);
   hm->SetTitle("Magnitude of the 1st transform");
   double scale = 1/h->Integral();
 // h->Scale(scale);
   hm->Draw();
   FindPeak(hm,  i,  name);
//  hm->SetLineColor(i+1);
   // for "real" frequencies you have to divide the x-axes range with the range of your function 
   //(in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n)
//  hm->SetStats(kFALSE); 
}



void Style_Cosmetics(){    
  // stat box
  gStyle->SetStatX(.997);  //  mid x-point for a stat box
  gStyle->SetStatW(0.16);  //  width of a stat box
  gStyle->SetOptStat(1111);  //values displayed in a stat box
  gStyle->SetStatColor(91);   //  set the stat box fill color
  gStyle->SetStatBorderSize(2);   //  set the stat box border size
  gStyle->SetStatFontSize(0.044);   // font size for a stat box

  // pad
  gStyle->SetPadTopMargin(0.1);  //top margin for a graph
  gStyle->SetPadRightMargin(0.03); //  right margin for a graph
  gStyle->SetPadBottomMargin(0.11);  //   bottom margin for a graph
  gStyle->SetPadLeftMargin(0.12);  //  left margin for a graph

  // title
  gStyle->SetTitleColor(2,"xy");  //  graph title colors 
  gStyle->SetLabelSize(0.09,"xy");   //size of graph labels
  gStyle->SetTitleOffset(.2,"x");  //  position of title,
  gStyle->SetTitleOffset(1.6,"y");   //  position of  title, y-axis

 // gStyle->UseCurrentStyle();  
}


