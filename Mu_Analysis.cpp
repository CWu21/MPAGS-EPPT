#include <iomanip>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>

#include "TH1F.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"


std::vector<std::map<int, std::vector<double>>> HitsDC1L(5);
std::vector<std::map<int, std::vector<double>>> HitsDC2L(5);
vector<int> flag(1000);
vector<TH1F*> H2;

auto pi = TMath::Pi();

void DefineHistos(){

    H2.push_back(new TH1F("Start_Point_x","Start Point (z = -8m);x(m);Entries",100,-0.0005,0.0005));
    H2.push_back(new TH1F("Start_Point_y","Start Point (z = -8m);y(m);Entries",100,-0.0005,0.0005));
    H2.push_back(new TH1F("Start_Angle_xz","Start Angle in xz plane;/phi (degree);Entries",100,-0.005,0.005));
    H2.push_back(new TH1F("Start_Angle_yz","Start Angle in yz plane;/theta (degree);Entries",100,-0.005,0.005));
    H2.push_back(new TH1F("Momentum","Momentum;P (GeV);Entries",100,95,105));


}

void DrawHitsLocation(){

    std::ifstream file1("./Hits.txt");

    std::cout << "Parsing HitsInfo" << std::endl;
    int row, instance, Dc1, Dc2;
    int I=0;
    vector<double> Info(6);

    vector<TH1F*> H;

    H.push_back(new TH1F("DC1X0","DC1X0;x;Hits",200,-1,1));
    H.push_back(new TH1F("DC1Y0","DC1Y0;y;Hits",200,-1,1));
    H.push_back(new TH1F("DC1Z0","DC1Z0;z;Hits",5,0,5));
    H.push_back(new TH1F("DC1X1","DC1X1;x;Hits",200,-1,1));
    H.push_back(new TH1F("DC1Y1","DC1Y1;y;Hits",200,-1,1));
    H.push_back(new TH1F("DC1Z1","DC1Z1;z;Hits",5,0,5));
    H.push_back(new TH1F("DC1X2","DC1X2;x;Hits",200,-1,1));
    H.push_back(new TH1F("DC1Y2","DC1Y2;y;Hits",200,-1,1));
    H.push_back(new TH1F("DC1Z2","DC1Z2;z;Hits",5,0,5));
    H.push_back(new TH1F("DC1X3","DC1X3;x;Hits",200,-1,1));
    H.push_back(new TH1F("DC1Y3","DC1Y3;y;Hits",200,-1,1));
    H.push_back(new TH1F("DC1Z3","DC1Z3;z;Hits",5,0,5));
    H.push_back(new TH1F("DC1X4","DC1X4;x;Hits",200,-1,1));
    H.push_back(new TH1F("DC1Y4","DC1Y4;y;Hits",200,-1,1));
    H.push_back(new TH1F("DC1Z4","DC1Z4;z;Hits",5,0,5));

    H.push_back(new TH1F("DC2X0","DC2X0;x;Hits",200,-15,-7));
    H.push_back(new TH1F("DC2Y0","DC2Y0;y;Hits",200,-1,1));
    H.push_back(new TH1F("DC2Z0","DC2Z0;z;Hits",5,0,5));
    H.push_back(new TH1F("DC2X1","DC2X1;x;Hits",200,-15,-7));
    H.push_back(new TH1F("DC2Y1","DC2Y1;y;Hits",200,-1,1));
    H.push_back(new TH1F("DC2Z1","DC2Z1;y;Hits",5,0,5));
    H.push_back(new TH1F("DC2X2","DC2X2;x;Hits",200,-15,-7));
    H.push_back(new TH1F("DC2Y2","DC2Y2;y;Hits",200,-1,1));
    H.push_back(new TH1F("DC2Z2","DC2Z2;z;Hits",5,0,5));
    H.push_back(new TH1F("DC2X3","DC2X3;x;Hits",200,-15,-7));
    H.push_back(new TH1F("DC2Y3","DC2Y3;y;Hits",200,-1,1));
    H.push_back(new TH1F("DC2Z3","DC2Z3;z;Hits",5,0,5));
    H.push_back(new TH1F("DC2X4","DC2X4;x;Hits",200,-15,-7));
    H.push_back(new TH1F("DC2Y4","DC2Y4;y;Hits",200,-1,1));
    H.push_back(new TH1F("DC2Z4","DC2Z4;z;Hits",5,0,5));

    for(std::string line; getline(file1, line ); ){
        std::istringstream iss(line);

        if ((iss>>row>>instance>>Info[0]>>Info[1]>>Info[2]>>Info[3]>>Info[4]>>Info[5]>>Dc1>>Dc2)){
        flag[row]=1;
       // std::cout << "row= " << row << " instance= " << instance <<" xyz= " <<Info[0]<<" "<< Info[1]<<" "<< Info[2]<<" "<<Info[3]<<" "<<Info[4]<<" "<<Info[5]<<std::endl;
        if(Dc1==5&&Dc2==5){  // get rid of non-muon event
        for(int j=0; j<3; j++){
                switch(instance){
                    case 0: 
                        HitsDC1L[0][row].push_back(Info[j]);
                        HitsDC2L[0][row].push_back(Info[j+3]);
                        H[j]->Fill(Info[j]);
                        H[j+15]->Fill(Info[j+3]);
                        break;
                    case 1: 
                        HitsDC1L[1][row].push_back(Info[j]);
                        HitsDC2L[1][row].push_back(Info[j+3]);
                        H[j+3]->Fill(Info[j]);
                        H[j+18]->Fill(Info[j+3]);
                        break;
                    case 2: 
                        HitsDC1L[2][row].push_back(Info[j]);
                        HitsDC2L[2][row].push_back(Info[j+3]);
                        H[j+6]->Fill(Info[j]);
                        H[j+21]->Fill(Info[j+3]);
                        break;
                    case 3: 
                        HitsDC1L[3][row].push_back(Info[j]);
                        HitsDC2L[3][row].push_back(Info[j+3]);
                        H[j+9]->Fill(Info[j]);
                        H[j+24]->Fill(Info[j+3]);
                        break;
                    case 4: 
                        HitsDC1L[4][row].push_back(Info[j]);
                        HitsDC2L[4][row].push_back(Info[j+3]);
                        H[j+12]->Fill(Info[j]);
                        H[j+27]->Fill(Info[j+3]);
                        break;
                }
                
            }
        } else flag[row]=0;  
        
        }
    }

    file1.close();

//----------------Drawing-----------------
 
/*   for(int j=0; j<30; j++) {
            string Save(H[j]->GetName());
            TString SaveAs="Hits/"+Save+".png";
            TCanvas* C=new TCanvas();
            H[j]->Draw();
            C->SaveAs(SaveAs);
    }

*/
}

double GetMomentum(double x1, double y1, double x2, double y2, double phi1, double phi2){

    double d = sqrt(pow((x2 - x1),2)+2*2);
    double phi = fabs(phi1-phi2);
    double R = d*0.5/sin(phi/2);
    double Pt = 0.3*0.5*R; // Transverse Momentum Pt = 0.3*z*B*R

    double y = fabs(y1-y2);
    double P3 = 0.3*0.5*y/phi; // Momentum along the magnetic field

    double P = sqrt(pow(Pt,2)+pow(P3,2));

    //std::cout <<"------------------------"<<std::endl;
    //std::cout <<" R = "<<R<<std::endl;
    //std::cout <<" phi = "<<phi<<std::endl;
    //std::cout <<" d = "<<d<<std::endl;
    //std::cout <<" PT = "<<Pt<<std::endl;
    //std::cout <<" P3 = "<<P3<<std::endl;

    return P;
}   

void ReadParameter(TGraph *gr1, TGraph *gr2, TGraph *gr3, TGraph *gr4){
    
    TF1 *func1 = gr1->GetFunction("pol1");
    TF1 *func2 = gr2->GetFunction("pol1");
    double a1 = func1->GetParameter(0);
    double b1 = func1->GetParameter(1);
    double a2 = func2->GetParameter(0);
    double b2 = func2->GetParameter(1);
    
    TF1 *func3 = gr3->GetFunction("pol1");
    TF1 *func4 = gr4->GetFunction("pol1");
    double c1 = func3->GetParameter(0);
    double d1 = func3->GetParameter(1);
    double c2 = func4->GetParameter(0);
    double d2 = func4->GetParameter(1);

    double x0 = a1+b1*(-8);   // start point5.7
    double y0 = c1+d1*(-8);
    double phi1 = atan(b1); // angles in DC1
    double theta1 = atan(d1);
    double phi2 = atan(b2); // angles in DC2
    double theta2 = atan(d2);

    double x1 = a1+b1*(-1);
    double y1 = c1+d1*(-1);
    double x2 = a2+b2*(1);
    double y2 = c2+d2*(1);  // in and out point of magnetic field

    double P = GetMomentum(x1, y1, x2, y2, phi1, phi2);

     //std::cout <<" x0 = "<< x0 <<" angle = "<<theta1-theta2<<std::endl;
     //std::cout <<" in = "<< x1 <<" out = "<<x2<<std::endl;

    H2[0]->Fill(x0);
    H2[1]->Fill(y0);
    H2[2]->Fill(phi1*180/pi);
    H2[3]->Fill(theta1*180/pi);
    H2[4]->Fill(P);

    //std::cout <<" P = "<<P<<std::endl;

}

void  DrawTracks(){
   
   vector<TGraph*> gr1, gr2, gr3, gr4;
   int row;
   Int_t j=0;
   TCanvas *c1 = new TCanvas("Tracks","Tracks");
   TMultiGraph *mg1 = new TMultiGraph();
   TMultiGraph *mg2 = new TMultiGraph();
   TMultiGraph *mg3 = new TMultiGraph();
   TMultiGraph *mg4 = new TMultiGraph();
    mg1->SetTitle("DC1 Tracks x;z(m);x(m)");
    mg2->SetTitle("DC2 Tracks x;z(m);x(m)");
    mg3->SetTitle("DC1 Tracks y;z(m);y(m)");
    mg4->SetTitle("DC2 Tracks y;z(m);y(m)");

   for(row=0;row<1000;row++){

    if(flag[row]==1){
        Double_t x[5], y[5], z[5], x1[5], y1[5], z1[5];
        Int_t n = 5;
        for (Int_t i=0;i<n;i++) {

            x[i] = HitsDC1L[i][row][0]*0.001;
            y[i] = HitsDC1L[i][row][1]*0.001;
            z[i] = HitsDC1L[i][row][2]*0.5-6;
            x1[i]= HitsDC2L[i][row][0]*0.001;
            y1[i]= HitsDC2L[i][row][1]*0.001;
            z1[i]= HitsDC1L[i][row][2]*0.5+2.5;
      }

  //----------for debugging-----
 /* std::cout << "x=["; 
    std::copy(std::begin(x),std::end(x),std::ostream_iterator<double>(std::cout,","));
    std::cout << "\b]" << std::endl;
    std::cout << "x1=["; 
    std::copy(std::begin(x1),std::end(x1),std::ostream_iterator<double>(std::cout,","));
    std::cout << "\b]" << std::endl;
*/

    gr1.push_back(new TGraph(n,z,x));
    gr2.push_back(new TGraph(n,z1,x1));
    gr3.push_back(new TGraph(n,z,y));
    gr4.push_back(new TGraph(n,z1,y1));

    gr1[j]->Fit("pol1","Q");
    gr2[j]->Fit("pol1","Q");
    gr3[j]->Fit("pol1","Q");
    gr4[j]->Fit("pol1","Q");

    gr1[j]->SetLineColor(4);
    gr2[j]->SetLineColor(4);
    gr3[j]->SetLineColor(4);
    gr4[j]->SetLineColor(4);

    mg1->Add(gr1[j]);
    mg2->Add(gr2[j]);
    mg3->Add(gr3[j]);
    mg4->Add(gr4[j]);

    ReadParameter(gr1[j], gr2[j], gr3[j], gr4[j]);
    
    j++;
   }
   }

    mg1->Draw("AC*");
    c1->SaveAs("Hits/Dc1Tracks_x.png");
    mg2->Draw("AC*");
    c1->SaveAs("Hits/Dc2Tracks_x.png");
    mg3->Draw("AC*");
    c1->SaveAs("Hits/Dc1Tracks_y.png");
    mg4->Draw("AC*");
    c1->SaveAs("Hits/Dc2Tracks_y.png");
}



void Mu_Analysis(){

    DefineHistos();
    DrawHitsLocation();
    DrawTracks();

   for(int j=0; j<H2.size(); j++) {
            string Save(H2[j]->GetName());
            TString SaveAs="Hits/"+Save+".png";
            TCanvas* C=new TCanvas();

            std::cout <<" Fit for "<< Save <<" : "<<std::endl;
            H2[j]->Fit("gaus");
            H2[j]->Draw();
            C->SaveAs(SaveAs);
   }

}