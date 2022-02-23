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
#include <tuple>

map <int, vector<double>> EC, HC, DC1X, DC1Y, DC1Z, DC2X, DC2Y, DC2Z;
vector<double> ECtot, HCtot, a, b, c, d;


vector<TH1F*> H2;
TRandom3*rndm=new TRandom3(0);


void DefineHistos(){

    H2.push_back(new TH1F("Momentum","Momentum;p_{T} (GeV);Entries",100,50,150));// 
    H2.push_back(new TH1F("ECEnergy_x","EC;x (cm);Entries",100,-30,30));
    H2.push_back(new TH1F("ECEnergy_y","EC;y (cm);Entries",100,-30,30));
    H2.push_back(new TH1F("HCEnergy_x","HC;x (cm);Entries",100,-150,150));
    H2.push_back(new TH1F("HCEnergy_y","HC;y (cm);Entries",50,-30,30));

    H2.push_back(new TH1F("ECEnergy","EC;Energy(GeV);Entries",100,90,100));
    H2.push_back(new TH1F("HCEnergy","HC;Energy(GeV);Entries",100,0,10));
    H2.push_back(new TH1F("TotalEnergy","Energy;Energy(GeV);Entries",100,90,100));

}

void DefineBranches(){
    
    TFile *myFile = TFile::Open("./Positron_100GeV_0p5T.root");

    TTreeReader R("B5", myFile);

    TTreeReaderValue<int>  Dc1Hits(R,"Dc1Hits");
    TTreeReaderValue<int>  Dc2Hits(R,"Dc2Hits");
    TTreeReaderValue<double>  ECEnergy(R,"ECEnergy");
    TTreeReaderValue<double>  HCEnergy(R,"HCEnergy");
    TTreeReaderValue<double>  Time1(R,"Time1");
    TTreeReaderValue<double>  Time2(R,"Time2");
    TTreeReaderArray<double>  ECEnergyVector(R,"ECEnergyVector");
    TTreeReaderArray<double>  HCEnergyVector(R,"HCEnergyVector");
    TTreeReaderArray<double>  Dc1HitsVector_x(R,"Dc1HitsVector_x");
    TTreeReaderArray<double>  Dc1HitsVector_y(R,"Dc1HitsVector_y");
    TTreeReaderArray<double>  Dc1HitsVector_z(R,"Dc1HitsVector_z");
    TTreeReaderArray<double>  Dc2HitsVector_x(R,"Dc2HitsVector_x");
    TTreeReaderArray<double>  Dc2HitsVector_y(R,"Dc2HitsVector_y");
    TTreeReaderArray<double>  Dc2HitsVector_z(R,"Dc2HitsVector_z");


//---------------Store Branches---------------------

    int I = 0;
    while(R.Next()){
        double E1=0, E2=0;

        for(int j = 0; j < 80; j++){
            EC[I].push_back(ECEnergyVector.At(j));
            E1 += ECEnergyVector.At(j);
        }
        for(int k = 0; k < 20; k++){
            HC[I].push_back(HCEnergyVector.At(k));
            E2 += HCEnergyVector.At(k);
        }
        for(int x1 = 0; x1 < Dc1HitsVector_x.GetSize(); x1++){
            DC1X[I].push_back(Dc1HitsVector_x.At(x1));
        }       
        for(int y1 = 0; y1 < Dc1HitsVector_y.GetSize(); y1++){
            DC1Y[I].push_back(Dc1HitsVector_y.At(y1));
        } 
        for(int z1 = 0; z1 < Dc1HitsVector_z.GetSize(); z1++){
            DC1Z[I].push_back(Dc1HitsVector_z.At(z1));
        }
        for(int x2 = 0; x2 < Dc2HitsVector_x.GetSize(); x2++){
            DC2X[I].push_back(Dc2HitsVector_x.At(x2));
        }       
        for(int y2 = 0; y2 < Dc2HitsVector_y.GetSize(); y2++){
            DC2Y[I].push_back(Dc2HitsVector_y.At(y2));
        } 
        for(int z2 = 0; z2 < Dc2HitsVector_z.GetSize(); z2++){
            DC2Z[I].push_back(Dc2HitsVector_z.At(z2));
        } 

        ECtot.push_back(E1);
        HCtot.push_back(E2);
        //ECtot.push_back(*ECEnergy.Get());
        //HCtot.push_back(*HCEnergy.Get());

        I++;

        //if(E1>4000||*ECEnergy.Get()>4000)std::cout <<"I= "<<I<<" ECSum = "<<E1<<" ECtot = "<<*ECEnergy.Get()<<std::endl;

    }

}

double GetMomentum(double x1, double y1, double x2, double y2, double phi1, double phi2){

    double d = sqrt(pow((x2 - x1),2)+2*2);
    double phi = fabs(phi1-phi2);
    double R = d*0.5/sin(phi/2);
    double Pt = 0.3*0.5*R; // Transverse Momentum Pt = 0.3*z*B*R

    //double y = fabs(y1-y2);
    //double P3 = 0.3*0.5*y/phi; // Momentum along the magnetic field

    //double P = sqrt(pow(Pt,2)+pow(P3,2));

    //std::cout <<"------------------------"<<std::endl;
    //std::cout <<" R = "<<R<<std::endl;
    //std::cout <<" phi = "<<phi<<std::endl;
    //std::cout <<" d = "<<d<<std::endl;
    //std::cout <<" PT = "<<Pt<<std::endl;
    //std::cout <<" P3 = "<<P3<<std::endl;

    return Pt;
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

    double x0 = a1+b1*(-8);   // start point
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
    H2[0]->Fill(P);

    a.push_back(a2);
    b.push_back(b2);
    c.push_back(c2);
    d.push_back(d2);

}

double GetFloatPrecision(double value, double precision)
{
    return (floor((value * pow(10, precision) + 0.5))*1.0 / pow(10, precision)); 
}


void DrawTracks(){

    vector<TGraph*> gr1, gr2, gr3, gr4;
    int j=0;
    TCanvas *c1 = new TCanvas("Tracks","Tracks");
    TMultiGraph *mg1 = new TMultiGraph();
    TMultiGraph *mg2 = new TMultiGraph();
    TMultiGraph *mg3 = new TMultiGraph();
    TMultiGraph *mg4 = new TMultiGraph();
    mg1->SetTitle("DC1 Tracks x;z(m);x(m)");
    mg2->SetTitle("DC2 Tracks x;z(m);x(m)");
    mg3->SetTitle("DC1 Tracks y;z(m);y(m)");
    mg4->SetTitle("DC2 Tracks y;z(m);y(m)");

    for(j = 0; j < 1000; j++){

        Double_t x[5], y[5], z[5], x1[5], y1[5], z1[5];
        int n = 5;
        int k = 0, m = 0;

            for(int i=0; i<DC1Z[j].size() && k < 5;i++){
                if(DC1Z[j][i] != k) {
                    continue;
                }
                x[k] = DC1X[j][i]*0.001;
                //x[k] = GetFloatPrecision(x[k], 4);
                x[k] += rndm->Gaus(0.,0.0001);
                y[k] = DC1Y[j][i]*0.001;
                z[k] = DC1Z[j][i]*0.5-6;
                k++;
            }
/*
        std::cout << "Event : "<<j<<"===============For Debug================="<<std::endl;
        std::cout << "x=["; 
            std::copy(std::begin(x),std::end(x),std::ostream_iterator<double>(std::cout,","));
        std::cout << "\b]" << std::endl;
        std::cout << "z=["; 
        std::copy(std::begin(z),std::end(z),std::ostream_iterator<double>(std::cout,","));
        std::cout << "\b]" << std::endl;
        std::cout << "DC1X=["; 
            std::copy(std::begin(DC1X[j]),std::end(DC1X[j]),std::ostream_iterator<double>(std::cout,","));
        std::cout << "\b]" << std::endl;
        std::cout << "DC1Z=["; 
            std::copy(std::begin(DC1Z[j]),std::end(DC1Z[j]),std::ostream_iterator<double>(std::cout,","));
        std::cout << "\b]" << std::endl;
*/        

            for(int i=0; i<DC2Z[j].size() && m < 5;i++){
                if(DC2Z[j][i] != m) {
                    continue;
                }
                x1[m] = DC2X[j][i]*0.001;
                //x1[m] = GetFloatPrecision(x1[m], 4);
                x1[m] += rndm->Gaus(0.,0.0001);
                y1[m] = DC2Y[j][i]*0.001;
                z1[m] = DC2Z[j][i]*0.5+2.5;
                m++;
            }

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

void DrawEMC(){
    TCanvas *c2 = new TCanvas("ECHits","ECHits");
    TGraph *grEMC1 = new TGraph();
    TGraph2D *grEMC2 = new TGraph2D();

    int n = 1000;
    double x[n], y[n];
    for(int i = 0; i < n; i++){
        x[i] = a[i] + b[i]*6.72;
        y[i] = c[i] + d[i]*6.72;

        grEMC2->SetPoint(i, x[i], y[i], ECtot[i]);
    } 
    *grEMC1 = TGraph(n,x,y);
    grEMC1->SetTitle("Hits on EMC;x(m);y(m)");
    grEMC1->SetMarkerStyle(7);
    grEMC1->SetMarkerColor(4);
    grEMC1->Draw("AP");
    c2->SaveAs("Hits/EMCHits.png");

    grEMC2->SetTitle("Hits on EMC;x(m);y(m);Energy(Mev)");
    grEMC2->SetMarkerStyle(7);
    grEMC2->Draw("pcol"); //pcol
    c2->SaveAs("Hits/EMCHits2.png");

}

void DrawHC(){
    TCanvas *c3 = new TCanvas("HCHits","HCHits");
    TGraph *grHC1 = new TGraph();
    TGraph2D *grHC2 = new TGraph2D();

    int n = 1000;
    double x[n], y[n];
    for(int i = 0; i < n; i++){
        x[i] = a[i] + b[i]*7.5;
        y[i] = c[i] + d[i]*7.5;

        grHC2->SetPoint(i, x[i], y[i], HCtot[i]);
    } 
    *grHC1 = TGraph(n,x,y);
    grHC1->SetTitle("Hits on HC;x(m);y(m)");
    grHC1->SetMarkerStyle(7);
    grHC1->SetMarkerColor(4);
    grHC1->Draw("AP");
    c3->SaveAs("Hits/HCHits.png");

    grHC2->SetTitle("Hits on HC;x(m);y(m);Energy(Mev)");
    grHC2->SetMarkerStyle(7);
    grHC2->Draw("pcol"); //pcol
    c3->SaveAs("Hits/HCHits2.png");

}

void DrawEnergyDistribution(){
    TCanvas *c = new TCanvas();
    TH2* h1 = new TH2D("ECEnergy","ECEnergy", 50, -150, 150, 50, -30, 30);
    TH2* h2 = new TH2D("HCEnergy","HCEnergy", 50, -150, 150, 50, -30, 30);

    for(int j = 0; j < 1000; j++){
        int ECx, ECy, HCx, HCy;
        double ECxE{0}, ECyE{0}, ECEx{0}, ECEy{0};
        double HCxE{0}, HCyE{0}, HCEx{0}, HCEy{0};
        for(int i = 0; i < 80; i++){
            ECx = i/4;
            ECy = i%4;
            ECxE += ECx*EC[j][i];
            ECyE += ECy*EC[j][i];
        }
        //std::cout << "Event : "<<j<<" x = " <<ECxE/ECE<<std::endl;
        ECEx = (ECxE/ECtot[j]-9.5)*15;
        ECEy = (ECyE/ECtot[j]-1.5)*15;

        H2[1]->Fill(ECEx,1);
        H2[2]->Fill(ECEy,1);
        h1->Fill(ECEx,ECEy,1);

        for(int k = 0; k < 20; k++){
            HCx = k/2;
            HCy = k%2;
            HCxE += HCx*HC[j][k];
            HCyE += HCy*HC[j][k];            
        }
        HCEx = (HCxE/HCtot[j]-4.5)*30;
        HCEy = (HCyE/HCtot[j]-0.5)*30;
        
        H2[3]->Fill(HCEx,1);
        H2[4]->Fill(HCEy,1);
        h2->Fill(HCEx,HCEy,1);

        H2[5]->Fill(ECtot[j]/1000,1);
        H2[6]->Fill(HCtot[j]/1000,1);
        H2[7]->Fill((ECtot[j]+HCtot[j])/1000,1);

    }
    h1->GetXaxis()->SetTitle("x(cm)");
    h1->GetYaxis()->SetTitle("y(cm)");
    h1->GetZaxis()->SetTitle("Entries");
    h1->Draw("surf3");
    c->SaveAs("Hits/ECEnergy_xy.png");

    h2->GetXaxis()->SetTitle("x(cm)");
    h2->GetYaxis()->SetTitle("y(cm)");
    h2->GetZaxis()->SetTitle("Entries");
    h2->Draw("surf3");
    c->SaveAs("Hits/HCEnergy_xy.png");



}



void Analysis(){

    //TFile * File = new TFile("test.root","recreate");
    DefineHistos();
    DefineBranches();
    DrawTracks();
    //DrawEMC();
    //DrawHC();
    DrawEnergyDistribution();
    
    for(int j=0; j<H2.size(); j++) {
            string Save(H2[j]->GetName());
            TString SaveAs="Hits/"+Save+".png";
            TCanvas* C=new TCanvas();

            //std::cout <<" Fit for "<< Save <<" : "<<std::endl;
            //H2[j]->Fit("gaus");
            H2[j]->Draw();
            C->SaveAs(SaveAs);
   }


}