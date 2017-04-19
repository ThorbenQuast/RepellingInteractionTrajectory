#include "simResults.h"


simResult::simResult() {
  _data["px_e1_in"] = std::vector<double>();
  _data["py_e1_in"] = std::vector<double>();
  _data["pz_e1_in"] = std::vector<double>();
  _data["E_e1_in"] = std::vector<double>();
  _data["px_e1_out"] = std::vector<double>();
  _data["py_e1_out"] = std::vector<double>();
  _data["pz_e1_out"] = std::vector<double>();
  _data["E_e1_out"] = std::vector<double>();
  _data["mass_q1"] = std::vector<double>();
  _data["IP1_x"] = std::vector<double>();
  _data["IP1_y"] = std::vector<double>();
  _data["IP1_z"] = std::vector<double>();
  _data["px_q1_in"] = std::vector<double>();
  _data["py_q1_in"] = std::vector<double>();
  _data["pz_q1_in"] = std::vector<double>();
  _data["E_q1_in"] = std::vector<double>();
  _data["px_q1_out"] = std::vector<double>();
  _data["py_q1_out"] = std::vector<double>();
  _data["pz_q1_out"] = std::vector<double>();
  _data["E_q1_out"] = std::vector<double>();
  _data["Q2_1"] = std::vector<double>();

  _data["px_e2_in"] = std::vector<double>();
  _data["py_e2_in"] = std::vector<double>();
  _data["pz_e2_in"] = std::vector<double>();
  _data["E_e2_in"] = std::vector<double>();
  _data["mass_q2"] = std::vector<double>();
  _data["IP2_x"] = std::vector<double>();
  _data["IP2_y"] = std::vector<double>();
  _data["IP2_z"] = std::vector<double>();
  _data["px_e2_out"] = std::vector<double>();
  _data["py_e2_out"] = std::vector<double>();
  _data["pz_e2_out"] = std::vector<double>();
  _data["E_e2_out"] = std::vector<double>();
  _data["px_q2_in"] = std::vector<double>();
  _data["py_q2_in"] = std::vector<double>();
  _data["pz_q2_in"] = std::vector<double>();
  _data["E_q2_in"] = std::vector<double>();
  _data["px_q2_out"] = std::vector<double>();
  _data["py_q2_out"] = std::vector<double>();
  _data["pz_q2_out"] = std::vector<double>();
  _data["E_q2_out"] = std::vector<double>();
  _data["Q2_2"] = std::vector<double>();

  _data["DeltaR_initial"] = std::vector<double>();;
  _data["relPt_initial"] = std::vector<double>();;

  _data["px_q1_final"] = std::vector<double>();
  _data["py_q1_final"] = std::vector<double>();
  _data["pz_q1_final"] = std::vector<double>();
  _data["E_q1_final"] = std::vector<double>();
  _data["px_q2_final"] = std::vector<double>();
  _data["py_q2_final"] = std::vector<double>();
  _data["pz_q2_final"] = std::vector<double>();
  _data["E_q2_final"] = std::vector<double>();
  _data["relPt_final"] = std::vector<double>();
  _data["event_weight"] = std::vector<double>();
};

void simResult::addEntry(std::string key, double value) {
  _data[key].push_back(value);
};

void simResult::printHeader(std::ofstream &outfile) {
  std::map<std::string, std::vector<double> >::iterator it = _data.begin();
    size_t column_counter = 0;
    for (std::map<std::string, std::vector<double> >::iterator it = _data.begin()=_data.begin(); it != _data.end(); it++) {
      if (column_counter>0) outfile<<",";
      outfile<<it->first;
      column_counter++;
    }
    outfile<<std::endl;
};

void simResult::printEntries(std::ofstream &outfile) {
  std::map<std::string, std::vector<double> >::iterator it = _data.begin();
  int N_entries = it->second.size();
  for (size_t i = 0; i<N_entries; i++) {
    size_t column_counter = 0;
    for (it=_data.begin(); it != _data.end(); it++) {
      if (N_entries != it->second.size()) {
        std::cout<<"INVALID output table format-->writing aborted!"<<std::endl;
        break;
      }
      if (column_counter>0) outfile<<",";
      outfile<<it->second[i];
      column_counter++;
    }
    outfile<<std::endl;
  }
};