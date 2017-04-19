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
};

void simResult::addEntry(std::string key, double value) {
  _data[key].push_back(value);
};

void simResult::printHeader() {

};
void simResult::printEntries() {
  std::map<std::string, std::vector<double> >::iterator it = _data.begin();
  int N_entries = it->second.size();
  std::cout<<"Number of simulated events: "<<N_entries<<std::endl;
  for (size_t i = 0; i<N_entries; i++) {
    std::cout<<"******"<<std::endl;
    std::cout<<"EVENT NUMBER "<<i+1<<std::endl;
    for (it=_data.begin(); it != _data.end(); it++) {
      if (N_entries != it->second.size()) {
        std::cout<<"INVALID output table format!"<<std::endl;
        break;
      }
      std::cout<<it->first<<": "<<it->second[i]<<std::endl;
    }
  }
  /*
  for(size_t i=0; i<_data["relPt_final"].size(); i++) {
    std::cout<<std::endl<<std::endl;
    std::cout<<"Event: "<<i<<std::endl;
    std::cout<<"px_e1_in: "<<_data["px_e1_in"][i]<<std::endl;
    std::cout<<"py_e1_in: "<<_data["py_e1_in"][i]<<std::endl;
    std::cout<<"pz_e1_in: "<<_data["pz_e1_in"][i]<<std::endl;
    std::cout<<"E_e1_in: "<<_data["E_e1_in"][i]<<std::endl;
    std::cout<<"px_e1_out: "<<_data["px_e1_out"][i]<<std::endl;
    std::cout<<"py_e1_out: "<<_data["py_e1_out"][i]<<std::endl;
    std::cout<<"pz_e1_out: "<<_data["pz_e1_out"][i]<<std::endl;
    std::cout<<"E_e1_out: "<<_data["E_e1_out"][i]<<std::endl;
    std::cout<<"mass_q1: "<<_data["mass_q1"][i]<<std::endl;
    std::cout<<"px_q1_in: "<<_data["px_q1_in"][i]<<std::endl;
    std::cout<<"py_q1_in: "<<_data["py_q1_in"][i]<<std::endl;
    std::cout<<"pz_q1_in: "<<_data["pz_q1_in"][i]<<std::endl;
    std::cout<<"E_q1_in: "<<_data["E_q1_in"][i]<<std::endl;
    std::cout<<"px_q1_out: "<<_data["px_q1_out"][i]<<std::endl;
    std::cout<<"py_q1_out: "<<_data["py_q1_out"][i]<<std::endl;
    std::cout<<"pz_q1_out: "<<_data["pz_q1_out"][i]<<std::endl;
    std::cout<<"E_q1_out: "<<_data["E_q1_out"][i]<<std::endl;
    std::cout<<"Q2_1: "<<_data["Q2_1"][i]<<std::endl;
    std::cout<<"px_e2_in: "<<_data["px_e2_in"][i]<<std::endl;
    std::cout<<"py_e2_in: "<<_data["py_e2_in"][i]<<std::endl;
    std::cout<<"pz_e2_in: "<<_data["pz_e2_in"][i]<<std::endl;
    std::cout<<"E_e2_in: "<<_data["E_e2_in"][i]<<std::endl;
    std::cout<<"mass_q2: "<<_data["mass_q2"][i]<<std::endl;
    std::cout<<"px_e2_out: "<<_data["px_e2_out"][i]<<std::endl;
    std::cout<<"py_e2_out: "<<_data["py_e2_out"][i]<<std::endl;
    std::cout<<"pz_e2_out: "<<_data["pz_e2_out"][i]<<std::endl;
    std::cout<<"E_e2_out: "<<_data["E_e2_out"][i]<<std::endl;
    std::cout<<"px_q2_in: "<<_data["px_q2_in"][i]<<std::endl;
    std::cout<<"py_q2_in: "<<_data["py_q2_in"][i]<<std::endl;
    std::cout<<"pz_q2_in: "<<_data["pz_q2_in"][i]<<std::endl;
    std::cout<<"E_q2_in: "<<_data["E_q2_in"][i]<<std::endl;
    std::cout<<"px_q2_out: "<<_data["px_q2_out"][i]<<std::endl;
    std::cout<<"py_q2_out: "<<_data["py_q2_out"][i]<<std::endl;
    std::cout<<"pz_q2_out: "<<_data["pz_q2_out"][i]<<std::endl;
    std::cout<<"E_q2_out: "<<_data["E_q2_out"][i]<<std::endl;
    std::cout<<"Q2_2: "<<_data["Q2_2"][i]<<std::endl;
    std::cout<<"DeltaR_initial: "<<_data["DeltaR_initial"][i]<<std::endl;
    std::cout<<"relPt_initial: "<<_data["relPt_initial"][i]<<std::endl;
    std::cout<<"px_q1_final: "<<_data["px_q1_final"][i]<<std::endl;
    std::cout<<"py_q1_final: "<<_data["py_q1_final"][i]<<std::endl;
    std::cout<<"pz_q1_final: "<<_data["pz_q1_final"][i]<<std::endl;
    std::cout<<"E_q1_final: "<<_data["E_q1_final"][i]<<std::endl;
    std::cout<<"px_q2_final: "<<_data["px_q2_final"][i]<<std::endl;
    std::cout<<"py_q2_final: "<<_data["py_q2_final"][i]<<std::endl;
    std::cout<<"pz_q2_final: "<<_data["pz_q2_final"][i]<<std::endl;
    std::cout<<"E_q2_final: "<<_data["E_q2_final"][i]<<std::endl;
    std::cout<<"relPt_final: "<<_data["relPt_final"][i]<<std::endl;
  }
  */
};