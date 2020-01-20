#include "DataSet.hpp"

std::set<std::string> DataSet::getUnique(const Rcpp::CharacterVector& col) {
  std::set<std::string> unique;
  std::string val;
  for (int i = 0; i < n; i++) {
    val = (std::string) col[i];
    if (val == "*" || val=="-99")
      continue;
    unique.insert(val);
  }
  return unique;
}

DataSet::DataSet(const Rcpp::DataFrame& df, const int& maxDiscrete) {
  this->maxDiscrete=maxDiscrete;
  const Rcpp::CharacterVector names = df.names();
  this->m = names.length();
  this->n = df.nrows();
  this->data.set_size(this->n,this->m);
  int numUnique;
  std::string val, curName;
  
  for (int i = 0; i < m; i++) {

    Rcpp::CharacterVector col = df[i];

    curName = (std::string) names[i];

    std::set<std::string> unique = getUnique(col);
    
    numUnique = unique.size();
      
    if (numUnique > maxDiscrete) {
      
      variables.push_back(new ContinuousVariable(curName));

    } else {
      
      std::vector<std::string> categories;
      
      for (std::set<std::string>::iterator it=unique.begin(); it!=unique.end(); it++)
	categories.push_back(*it);

      std::sort(categories.begin(), categories.end());
      
      variables.push_back(new DiscreteVariable(curName, categories));
    }

    variableNames.push_back(curName);
    name2idx.insert(std::pair<std::string, int>(curName, i));


    for (int j = 0; j < n; j++) {

      val = (std::string) col[j];

      if (variables.at(i)->isMissingValue(val)) {

	throw std::runtime_error("Missing values not yet implemented");
      
      } else if (variables.at(i)->isContinuous()) {
	
	if (((ContinuousVariable*) variables.at(i))->checkValue(val)) {
	  this->data(j,i) = std::stod(val);
	  
	} else {
	  
	  throw std::runtime_error("invalid value encountered in row " + std::to_string(j) + " of variable " + curName);
	  
	}
      } else if (variables.at(i)->isDiscrete()) {
	
	if (((DiscreteVariable*) variables.at(i))->checkValue(val)) {
	  this->data(j,i) = (double) ((DiscreteVariable*) variables.at(i))->getIndex(val);

	} else {
	  
	  throw std::runtime_error("invalid value encountered in row " + std::to_string(j) + " of variable " + curName);
	}
      } else {
	
	throw std::runtime_error("invalid variable type");
      }
    }
  }
}

DataSet::~DataSet() {
  for (int i = 0; i < variables.size(); i++)
    delete variables.at(i);
}

// [[Rcpp::export]]
Rcpp::RObject rCausalMGMData(const Rcpp::DataFrame& df, const int& maxDiscrete=5) {
  DataSet *ds = new DataSet(df, maxDiscrete);

  Rcpp::XPtr<DataSet> data = Rcpp::XPtr<DataSet>(ds);
  
  data.attr("data") = df;
  
  return Rcpp::as<Rcpp::RObject>(data);
}
