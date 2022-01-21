#include <iomanip>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>

std::map<int, std::vector<int>> HitsDC1L0;
std::map<int, std::vector<int>> HitsDC1L1;
std::map<int, std::vector<int>> HitsDC1L2;
std::map<int, std::vector<int>> HitsDC1L3;
std::map<int, std::vector<int>> HitsDC1L4;

std::map<int, std::vector<int>> HitsDC2L0;
std::map<int, std::vector<int>> HitsDC2L1;
std::map<int, std::vector<int>> HitsDC2L2;
std::map<int, std::vector<int>> HitsDC2L3;
std::map<int, std::vector<int>> HitsDC2L4;


void ReadHitsLocation(){

    std::ifstream file1("./Hits.txt");

    std::cout << "Parsing HitsInfo" << std::endl;
    int row, instance;
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

        if ((iss>>row>>instance>>Info[0]>>Info[1]>>Info[2]>>Info[3]>>Info[4]>>Info[5])){
       // std::cout << "row= " << row << " instance= " << instance <<" xyz= " <<Info[0]<<" "<< Info[1]<<" "<< Info[2]<<" "<<Info[3]<<" "<<Info[4]<<" "<<Info[5]<<std::endl;
        for(int j=0; j<3; j++){
                switch(instance){
                    case 0: 
                        HitsDC1L0[row].push_back(Info[j]);
                        HitsDC2L0[row].push_back(Info[j+3]);
                        H[j]->Fill(Info[j]);
                        H[j+15]->Fill(Info[j+3]);
                        break;
                    case 1: 
                        HitsDC1L1[row].push_back(Info[j]);
                        HitsDC2L1[row].push_back(Info[j+3]);
                        H[j+3]->Fill(Info[j]);
                        H[j+18]->Fill(Info[j+3]);
                        break;
                    case 2: 
                        HitsDC1L2[row].push_back(Info[j]);
                        HitsDC2L2[row].push_back(Info[j+3]);
                        H[j+6]->Fill(Info[j]);
                        H[j+21]->Fill(Info[j+3]);
                        break;
                    case 3: 
                        HitsDC1L3[row].push_back(Info[j]);
                        HitsDC2L3[row].push_back(Info[j+3]);
                        H[j+9]->Fill(Info[j]);
                        H[j+24]->Fill(Info[j+3]);
                        break;
                    case 4: 
                        HitsDC1L4[row].push_back(Info[j]);
                        HitsDC2L4[row].push_back(Info[j+3]);
                        H[j+12]->Fill(Info[j]);
                        H[j+27]->Fill(Info[j+3]);
                        break;
                }
                
            }
        }
    }

    file1.close();

    for(int j=0; j<30; j++) {
            string Save(H[j]->GetName());
            TString SaveAs="Hits/"+Save+".png";
            TCanvas* C=new TCanvas();
            H[j]->Draw();
            C->SaveAs(SaveAs);}

    

}

void Analysis(){

     ReadHitsLocation();

}