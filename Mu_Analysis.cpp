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
 
    for(int j=0; j<30; j++) {
            string Save(H[j]->GetName());
            TString SaveAs="Hits/"+Save+".png";
            TCanvas* C=new TCanvas();
            H[j]->Draw();
            C->SaveAs(SaveAs);
    }


}

void  DrawTracks(){
   
   vector<TGraph*> gr1,gr2;
   int row;
   Int_t j=0;
   TCanvas *c1 = new TCanvas("Tracks","Tracks",200,10,600,400);
   TMultiGraph *mg1 = new TMultiGraph();
   TMultiGraph *mg2 = new TMultiGraph();
    mg1->SetTitle("DC1 Tracks;z(m);x");
    mg2->SetTitle("DC2 Tracks;z(m);x");

   for(row=0;row<1000;row++){

    if(flag[row]==1){
        Double_t x[5], z[5], x1[5], z1[5];
        Int_t n = 5;
        for (Int_t i=0;i<n;i++) {

            x[i] = HitsDC1L[i][row][0];
            z[i] = HitsDC1L[i][row][2]*0.5;
            x1[i]= HitsDC2L[i][row][0];
            z1[i]= HitsDC1L[i][row][2]*0.5;
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
    gr1[j]->SetLineColor(4);
    mg1->Add(gr1[j]);

    gr2.push_back(new TGraph(n,z1,x1));
    gr2[j]->SetLineColor(4);
    mg2->Add(gr2[j]);
    
    j++;
   }
   }

    mg1->Draw("AC*");
    c1->SaveAs("Hits/Dc1Tracks.png");
    mg2->Draw("AC*");
    c1->SaveAs("Hits/Dc2Tracks.png");
}




void Mu_Analysis(){

    DrawHitsLocation();
    DrawTracks();

}