#include "reweightings_handler.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdint>
#include <vector>
#include <cmath>
#include <Endian.h>

using namespace std;

void ReweightingsHandler::getData(RVector& result)
{
 read_data();
 if (number_of_meas != rw_prod.size())
    throw(std::runtime_error("Discrepancy between number of measurements and number of reweighting factors."));
}

bool ReweightingsHandler::queryData()
{
 read_data();
 if (number_of_meas != rw_prod.size()) return false;
 return true;
}

//*******************************************************************

void ReweightingsHandler::read_data()
{
 if (nrw) return; // data has already been read

 void (ReweightingsHandler::*read_file_ptr)(const string file_name);
 if (file_format == "OPENQCD")
    read_file_ptr = &ReweightingsHandler::read_openqcd;
 else if (file_format == "OPENQCD_12")
    read_file_ptr = &ReweightingsHandler::read_openqcd12;
 else if (file_format == "ASCII")
    read_file_ptr = &ReweightingsHandler::read_ascii;
 else{
    cout << "unrecognized reweightings file format " << file_format << endl;
    exit(1);}

 for (vector<string>::const_iterator fn=file_names.begin(); fn!=file_names.end(); ++fn){
    (this->*read_file_ptr)(*fn);}

 rw_prod.resize((uint)rw.size());
 for (uint n = 0; n < rw.size(); n++){
   rw_prod[n] = 1.;
   for (vector<double>::const_iterator it=rw[n].begin(); it!=rw[n].end(); it++)
     rw_prod[n] *= *it;}
}

void ReweightingsHandler::read_openqcd(const string file_name)
{
 read_openqcd_binary(file_name, false);
}

void ReweightingsHandler::read_openqcd12(const string file_name)
{
 read_openqcd_binary(file_name, true);
}

void ReweightingsHandler::read_ascii(const string file_name)
{
 fstream file(file_name, ios::in);
 if (!file.is_open()){
    cout << "could not open file " << file_name << endl;
    exit(1);}

 string line;
 while (getline(file, line)){
    if (line.front() == '#') continue;
    if (!nrw){
      istringstream ss(line);
      string token;
      int tokens = 0;
      while (ss >> token) tokens++;
      nrw = tokens - 2;
    }
    istringstream ss(line);
    string nc;
    ss >> nc;
    vector<double> rw_nc(nrw);
    for (vector<double>::iterator it=rw_nc.begin(); it!=rw_nc.end(); it++)
      ss >> *it;

    rw.push_back(rw_nc);
 }
 file.close();
}

void ReweightingsHandler::read_openqcd_binary(const string file_name, bool v12)
{
 if (endianness() == BigEndian){
    cout << "Big endian currently not supported...sorry" << endl;
    exit(1);}

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

