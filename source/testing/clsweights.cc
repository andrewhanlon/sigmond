#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdint>
#include <vector>
#include <cmath>

using namespace std;


void read_openqcd_binary(const string file_name, bool v12)
{
 bool first_file = !nrw;

 fstream file(file_name, ios::binary | ios::in);
 if (!file.is_open()){
    cout << "could not open file " << file_name << endl;
    exit(1);}
 
 int32_t temp_nrw;
 file.read(reinterpret_cast<char *>(&temp_nrw), sizeof(temp_nrw));
 if (!first_file && temp_nrw != nrw){
    cout << "Mismatch in number of reweighting factors among files" << endl;
    exit(1);}
 nrw = temp_nrw;

 vector<int32_t> temp_nfct(nrw,1);
 if (!v12)
    file.read(reinterpret_cast<char *>(temp_nfct.data()), temp_nfct.size()*sizeof(int32_t));
 if (!first_file && temp_nfct != nfct){
    cout << "Mismatch in number of Hasenbusch factors among reweighting files" << endl;
    exit(1);}
 nfct = temp_nfct;

 vector<int32_t> temp_nsrc(nrw);
 file.read(reinterpret_cast<char *>(temp_nsrc.data()), temp_nsrc.size()*sizeof(int32_t));
 if (!first_file && temp_nsrc != nsrc){
    cout << "Mismatch in number of source fields among reweighting files" << endl;
    exit(1);}
 nsrc = temp_nsrc;

 vector<vector<vector<double> > > sqn_nc(nrw);
 vector<vector<vector<double> > > lnr_nc(nrw);

 int32_t nc;
 file.read(reinterpret_cast<char *>(&nc), sizeof(nc));
 while (!file.eof()){
    vector<double> rw_nc(nrw, 1.);
    for (int32_t irw=0; irw < nrw; irw++){
      sqn_nc[irw].resize(nfct[irw]);
      lnr_nc[irw].resize(nfct[irw]);
      for (int32_t ifct=0; ifct < nfct[irw]; ifct++){
        sqn_nc[irw][ifct].resize(nsrc[irw]);
        lnr_nc[irw][ifct].resize(nsrc[irw]);
        file.read(reinterpret_cast<char *>(sqn_nc[irw][ifct].data()), sqn_nc[irw][ifct].size()*sizeof(double));
        file.read(reinterpret_cast<char *>(lnr_nc[irw][ifct].data()), lnr_nc[irw][ifct].size()*sizeof(double));

        double sum = 0.;
        for (int32_t isrc=0; isrc < nsrc[irw]; isrc++)
          sum += exp(-lnr_nc[irw][ifct][isrc]);
        rw_nc[irw] *= sum/(nsrc[irw]);}}
    sqn.push_back(sqn_nc);
    lnr.push_back(lnr_nc);
    rw.push_back(rw_nc);
    file.read(reinterpret_cast<char *>(&nc), sizeof(nc));}

 file.close();
}

