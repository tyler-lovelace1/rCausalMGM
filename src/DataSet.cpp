#include "DataSet.hpp"

std::set<std::string> DataSet::getUnique(const Rcpp::CharacterVector &col)
{
    std::set<std::string> unique;
    std::string val;
    for (int i = 0; i < n; i++)
    {
        val = (std::string)col[i];
        if (val == "*" || val == "-99")
            continue;
        unique.insert(val);
    }
    return unique;
}

bool DataSet::checkCensoring(const Rcpp::CharacterVector& col,
			     arma::vec& values,
			     arma::uvec& censor) {
    bool censored = false;
    std::string val;
    for (int i = 0; i < n; i++) {
	val = (std::string)col[i];
	censor[i] = 1;
	if (val.back()=='+') {
	    censor[i] = 0;
	    censored = true;
	    // values[i] = std::stod(val.substr(0,val.length()-1));
	} 
	values[i] = std::stod(val);
    }
    return censored;
}


DataSet::DataSet(const Rcpp::DataFrame &df) {
    this->maxDiscrete = maxDiscrete;
    const Rcpp::CharacterVector names = df.names();
    this->m = names.length();
    this->n = df.nrows();
    // this->data.set_size(this->n, this->m);
    this->data = arma::mat(n,m,arma::fill::zeros);
    this->missing = arma::umat(n,m,arma::fill::zeros);
    this->name2idx = std::unordered_map<std::string, int>();
    this->var2idx = std::unordered_map<Node, int>();
    int numUnique;
    std::string val, curName;

    for (int i = 0; i < m; i++) {
	
	curName = (std::string) names[i];

	// Rcpp::CharacterVector xVec = Rcpp::clone(df[i]);
	Rcpp::RObject x = df[i];

	// Rcpp::LogicalVector isna = Rcpp::is_na(xVec);
	
	// missing.col(i) = Rcpp::as<arma::uvec>(isna);

	if (Rcpp::is<Rcpp::NumericVector>(x)){
	    if (Rf_isMatrix(x))  {
		// Rcpp::Rcout << "NumericMatrix\n";
		if (x.inherits("Surv")) {
		    // Rcpp::Rcout << "Surv Object\n";
		    arma::mat surv = Rcpp::as<arma::mat>(x);
		    variables.push_back(Node(new CensoredVariable(curName)));
		    arma::vec values(surv.col(0));
		    arma::uvec censor(arma::conv_to<arma::uvec>::from(surv.col(1)));
		    arma::uvec nonmissing = arma::intersect(arma::find_finite(values),
							    arma::find_finite(censor));
		    variables[i].setCensor(values(nonmissing), censor(nonmissing));
		    data.col(i) = values;
		    // missing.col(i) = Rcpp::as<arma::uvec>(Rcpp::is_na(Rcpp::NumericVector(values.begin(), values.end())));
		} else {
		    Rcpp::stop(curName + " is a numeric matrix but is not a Surv object");
		}
	    } else {
		// Rcpp::Rcout << "NumericVector\n";
		arma::vec values = Rcpp::as<arma::vec>(x);
		// Rcpp::Rcout << values.t();
		arma::vec uniqVals = arma::unique(values(arma::find_finite(values)));
		// Rcpp::Rcout << uniqVals.t();
		if (uniqVals.n_elem < 10) {
		    Rcpp::Rcout << "Warning: Variable " + curName + " has fewer than 10 unique values and is numeric. " + curName + " is being treated as a continuous variable. If intended to be categorical, convert " + curName + " to a factor.\n";
		}
		variables.push_back(Node(new ContinuousVariable(curName)));
		data.col(i) = values;
		// missing.col(i) = Rcpp::as<arma::uvec>(Rcpp::is_na(Rcpp::as<Rcpp::NumericVector>(x)));
	    }
	} else if (Rcpp::is<Rcpp::IntegerVector>(x)) {
	    if(Rf_isFactor(x)) {
		// Rcpp::Rcout << "Factor " << curName << "\n";
		std::vector<std::string> tempLevels = Rcpp::as<std::vector<std::string>>(x.attr("levels"));
		arma::vec values = Rcpp::as<arma::vec>(x);

		arma::vec uniqVals = arma::sort(arma::unique(values(arma::find_finite(values))));

		// Rcpp::Rcout << values.t() << std::endl;

		// Rcpp::Rcout << uniqVals.t() << std::endl;


		std::vector<std::string> levels;

		for (double val : uniqVals) {
		    levels.push_back(tempLevels[(int) (val-1)]);
		}

		arma::vec mappedValues(n, arma::fill::zeros);

		for (int idx = 0; idx < n; idx++) {
		    for (int cat = 0; cat < levels.size(); cat++) {
			if (tempLevels[(int) (values[idx]-1)] == levels[cat]) {
			    mappedValues[idx] = cat;
			}
		    }
		}

		// for ( std::string val : levels) {
		//     Rcpp::Rcout << val << "\t";
		// }
		// Rcpp::Rcout << std::endl;
		
		// Rcpp::Rcout << mappedValues.t() << std::endl;

		if (levels.size() >= 10) {
		    Rcpp::warning("Categorical variable " + curName + " has 10 or more categories. Fitting models with large numbers of categories is not recommended. If " + curName + " is intended to be a continuous variable, convert it to numeric.");
		}	    
	    
		variables.push_back(Node(new DiscreteVariable(curName, levels)));
		data.col(i) = mappedValues;
		// missing.col(i) = Rcpp::as<arma::uvec>(Rcpp::is_na(Rcpp::as<Rcpp::IntegerVector>(x)));
	    }
	    else {
		// Rcpp::Rcout << "IntegerVector\n";
		arma::vec values = Rcpp::as<arma::vec>(x);
		arma::vec uniqVals = arma::unique(values(arma::find_finite(values)));
		if (uniqVals.n_elem < 10) {
		    Rcpp::warning("Variable " + curName + " has fewer than 10 unique values and is numeric. " + curName + " is being treated as a continuous variable. If intended to be categorical, convert " + curName + " to a factor.");
		}
		variables.push_back(Node(new ContinuousVariable(curName)));
		data.col(i) = values;	
		// missing.col(i) = Rcpp::as<arma::uvec>(Rcpp::is_na(Rcpp::as<Rcpp::IntegerVector>(x)));
	    }
	}
	else if (Rcpp::is<Rcpp::CharacterVector>(x)) {
	    // Rcpp::Rcout << "CharacterVector\n";
	    Rcpp::CharacterVector levels = Rcpp::sort_unique(Rcpp::as<Rcpp::CharacterVector>(x));
	    Rcpp::CharacterVector values = Rcpp::as<Rcpp::CharacterVector>(x);

	    arma::vec mappedValues(n, arma::fill::zeros);

	    for (int idx = 0; idx < n; idx++) {
		for (int cat = 0; cat < levels.size(); cat++) {
		    if (values[idx] == levels[cat]) {
			mappedValues[idx] = cat;
		    }
		}
	    }

	    if (levels.size() >= 10) {
	        Rcpp::warning("Categorical variable " + curName + " has 10 or more categories. Fitting models with large numbers of categories is not recommended. If " + curName + " is intended to be a continuous variable, convert it to numeric.");
	    }
	    
	    variables.push_back(Node(new DiscreteVariable(curName, Rcpp::as<std::vector<std::string>>(levels))));
	    data.col(i) = mappedValues;
	    
	    // missing.col(i) = Rcpp::as<arma::uvec>(Rcpp::is_na(Rcpp::as<Rcpp::CharacterVector>(x)));
	    
	}
	variableNames.push_back(curName);
        name2idx.insert(std::pair<std::string, int>(curName, i));
        var2idx.insert(std::pair<Node, int>(variables[i], i));
    }

    // Rcpp::Rcout << data << std::endl;

    missing.elem(arma::find_nonfinite(data)).ones();

    // if (arma::accu(missing) > 0) {
    // 	Rcpp::warning("Missing values detected. Samples with missing values will be dropped.");
    // }
}


DataSet::DataSet(const Rcpp::DataFrame &df, const int maxDiscrete) {
    this->maxDiscrete = maxDiscrete;
    const Rcpp::CharacterVector names = df.names();
    this->m = names.length();
    this->n = df.nrows();
    this->data.set_size(this->n, this->m);
    this->name2idx = std::unordered_map<std::string, int>();
    this->var2idx = std::unordered_map<Node, int>();
    int numUnique;
    std::string val, curName;

    for (int i = 0; i < m; i++)
    {

	// Rcpp::Rcout << (std::string)names[i] << std::endl;
	// Rcpp::Rcout << TYPEOF(df[i]) << std::endl;

        Rcpp::CharacterVector col = df[i];

        curName = (std::string)names[i];

        std::set<std::string> unique = getUnique(col);

        numUnique = unique.size();

        if (numUnique > maxDiscrete) {
	    arma::vec values(n);
	    arma::uvec censor(n);
	    if (checkCensoring(col, values, censor)) {
		variables.push_back(Node(new CensoredVariable(curName)));
	        variables[i].setCensor(values, censor);
	    } else {
		variables.push_back(Node(new ContinuousVariable(curName)));
	    }
        }
        else {

            std::vector<std::string> categories;

            for (auto it = unique.begin(); it != unique.end(); it++)
                categories.push_back(*it);

            std::sort(categories.begin(), categories.end());

            variables.push_back(Node(new DiscreteVariable(curName, categories)));
        }

        variableNames.push_back(curName);
        name2idx.insert(std::pair<std::string, int>(curName, i));
        var2idx.insert(std::pair<Node, int>(variables[i], i));

        for (int j = 0; j < n; j++)
        {

            val = (std::string)col[j];

            if (variables[i].isMissingValue(val))
            {
                throw std::runtime_error("Missing values not yet implemented");
            }
            else if (variables[i].isContinuous())
            {
                if (variables[i].checkValue(val))
                {
                    this->data(j, i) = std::stod(val);
                }
                else
                {
                    throw std::runtime_error("invalid value encountered in row " + std::to_string(j) + " of continuous variable " + curName + 
                        " (HINT: If this variable is intended to be discrete, consider raising the maxDiscrete parameter)");
                }
            }
            else if (variables[i].isDiscrete())
            {
                if (variables[i].checkValue(val))
                {
                    this->data(j, i) = (double) variables[i].getIndex(val);
                }
                else
                {
                    throw std::runtime_error("invalid value encountered in row " + std::to_string(j) + " of discrete variable " + curName);
                }
            }
	    else if (variables[i].isCensored()) {
		if (variables[i].checkValue(val)) {
		    this->data(j, i) = std::stod(val);
		    if (this->data(j, i) <= 0) {
			throw std::runtime_error("invalid value encountered in row " + std::to_string(j) + " of censored variable " + curName + " (HINT: Time values of censored variables must be greater than 0)");
		    }
		  // this->data(j, i) = std::log(this->data(j, i));
		}
		else {
		    throw std::runtime_error("invalid value encountered in row " + std::to_string(j) + " of censored variable " + curName + " (HINT: If this variable is intended to be discrete, consider raising the maxDiscrete parameter)");
		}
	    }
            else
            {
                throw std::runtime_error("invalid variable type");
            }
        }
    }
}

void DataSet::dropMissing() {
    if (arma::accu(missing) > 0) {
	Rcpp::warning("Missing values detected. Samples with missing values will be dropped.");

	arma::uvec missingByRow = arma::sum(missing, 1);
	arma::uvec completeSamples = arma::find(missingByRow == 0);

	n = completeSamples.n_elem;
        data = data.rows(completeSamples);
	
	for (int j = 0; j < m; j++) {
	    if (variables[j].isCensored()) {
		arma::uvec censor = variables[j].getCensor();
		variables[j].setCensor(data.col(j), censor.elem(completeSamples));
	    }
	}
    }
}

void DataSet::npnTransform() {
    for (int j = 0; j < m; j++) {
	if (variables[j].isContinuous()) {
	    arma::vec vals(data.col(j));
	    double mu = arma::mean(vals);
	    double sd = arma::stddev(vals);
	    
	    arma::uvec indices = arma::sort_index(vals);
	    int marker = 0;

	    for (int i = 0; i < n; i++) {
		if (vals(indices(i)) > vals(indices(marker))) {
		    vals(indices.subvec(marker, i-1)).fill(i);
		    marker = i;
		}
	    }
	    vals(indices.subvec(marker, n-1)).fill(n);
	    vals /= ((double) (n+1));
	    // Rcpp::Rcout << vals(indices).t() << std::endl;
	    // Rcpp::Rcout << vals.t() << std::endl;
	    // Rcpp::Rcout << data.col(j).t() << std::endl;

	    vals = Rcpp::qnorm(Rcpp::NumericVector(vals.begin(), vals.end()), mu, sd);
	    data.col(j) = vals;
	    
	    // Rcpp::Rcout << vals.t() << std::endl;
	}
    }
}

void DataSet::addVariable(Node v) {
    // Rcpp::Rcout << "adding variable " << v->getName() << "..." << std::endl;
    variables.push_back(v);
    arma::vec col = arma::vec(n);
    for (int i = 0; i < n; i++)
        col[i] = arma::datum::nan;

    variableNames.push_back(v.getName());
    name2idx.insert(std::pair<std::string, int>(v.getName(), m));
    var2idx.insert(std::pair<Node, int>(v, m));

    data.insert_cols(m, col);
    m++;

    // Rcpp::Rcout << "variable " << v->getName() << " inserted" << std::endl;
}

void DataSet::addVariable(int i, Node v)
{
    if (i > m || i < 0)
        throw std::invalid_argument("index " + std::to_string(i) + " out of bounds");

    variables.insert(variables.begin() + i, v);
    arma::vec col = arma::vec(n);
    for (int i = 0; i < n; i++)
        col[i] = arma::datum::nan;

    variableNames.push_back(v.getName());
    name2idx.insert(std::pair<std::string, int>(v.getName(), i));
    var2idx.insert(std::pair<Node, int>(v, i));

    data.insert_cols(i, col);
    m++;
}

// DataSet::~DataSet() {
//     for (int i = 0; i < variables.size(); i++)
// 	delete variables[i];
// }

// // WARNING: only call when done with an R-exposed function
// void DataSet::deleteVariables() {
//     for (int i = 0; i < variables.size(); i++)
//         delete variables[i];
// }


// // Deep copy
// std::vector<Variable *> DataSet::copyVariables() {
//     std::vector<Variable *> result = std::vector<Variable *>(m);

//     for (int i = 0; i < m; i++) {
//         if (variables[i]->isContinuous())
//             result[i] = new ContinuousVariable(*((ContinuousVariable*) variables[i]));
//         else if (variables[i]->isDiscrete())
//             result[i] = new DiscreteVariable(*((DiscreteVariable*) variables[i]));
//     }

//     return result;
// }

// DataSet::DataSet(const DataSet& ds) {
//     // this->maxDiscrete=ds.maxDiscrete;
//     // this->m = ds.m;
//     // this->n = ds.n;
//     // // this->variables = ds.variables;
//     // this->variableNames = ds.variableNames;
//     // this->name2idx = ds.name2idx;
//     // // this->var2idx = ds.var2idx;
//     // this->data = ds.data;
  
//     maxDiscrete=ds.maxDiscrete;
//     m = ds.m;
//     n = ds.n;
//     variables = std::vector<Variable*>();
//     variableNames = ds.variableNames;
//     name2idx = ds.name2idx;
//     // this->variableNames = std::vector<std::string>();
//     // this->name2idx = std::unordered_map<std::string, int>();
//     var2idx = std::unordered_map<Variable*, int>();
//     data = ds.data;

//     for (int j = 0; j < m; j++) {
// 	if (ds.variables.at(j)->isDiscrete())
// 	    variables.push_back(new DiscreteVariable(*((DiscreteVariable*) ds.variables[j])));
// 	else if (ds.variables.at(j)->isContinuous())
// 	    variables.push_back(new ContinuousVariable(*((ContinuousVariable*) ds.variables[j])));

// 	var2idx.insert(std::pair<Variable*, int>(variables[j], j));
//     }
// }

DataSet::DataSet(const DataSet& ds, const arma::urowvec& rows) {
    maxDiscrete = ds.maxDiscrete;
    m = ds.m;
    n = rows.n_elem;
    variables = ds.variables;
    variableNames = ds.variableNames;
    name2idx = ds.name2idx;
    var2idx = ds.var2idx;
    data = ds.data.rows(rows);

    for (int j = 0; j < m; j++) {
        if (variables[j].isCensored()) {
	    arma::uvec censor = variables[j].getCensor();
	    variables[j].setCensor(data.col(j), censor.elem(rows));
	}
    }
}

// DataSet& DataSet::operator=(const DataSet& ds) {
//     // this->maxDiscrete=ds.maxDiscrete;
//     // this->m = ds.m;
//     // this->n = ds.n;
//     // this->variables = ds.variables;
//     // this->variableNames = ds.variableNames;
//     // this->name2idx = ds.name2idx;
//     // this->var2idx = ds.var2idx;
//     // this->data = ds.data;

//     maxDiscrete=ds.maxDiscrete;
//     m = ds.m;
//     n = ds.n;
//     variables = std::vector<Variable*>();
//     variableNames = ds.variableNames;
//     name2idx = ds.name2idx;
//     // this->variableNames = std::vector<std::string>();
//     // this->name2idx = std::unordered_map<std::string, int>();
//     var2idx = std::unordered_map<Variable*, int>();
//     data = ds.data;

//     for (int j = 0; j < m; j++) {
// 	if (ds.variables.at(j)->isDiscrete())
// 	    variables.push_back(new DiscreteVariable(*((DiscreteVariable*) ds.variables[j])));
// 	else if (ds.variables.at(j)->isContinuous())
// 	    variables.push_back(new ContinuousVariable(*((ContinuousVariable*) ds.variables[j])));

// 	var2idx.insert(std::pair<Variable*, int>(variables[j], j));
//     }
  
//     return *this;
// }

// DataSet::DataSet(DataSet&& ds) {
//     // this->maxDiscrete=ds.maxDiscrete;
//     // this->m = ds.m;
//     // this->n = ds.n;
//     // this->variables = ds.variables;
//     // this->variableNames = ds.variableNames;
//     // this->name2idx = ds.name2idx;
//     // this->var2idx = ds.var2idx;
//     // this->data = ds.data;
  
//     // this->variables.clear();
//     // this->var2idx.clear();

//     maxDiscrete=ds.maxDiscrete;
//     m = ds.m;
//     n = ds.n;
//     variables = std::vector<Variable*>();
//     variableNames = ds.variableNames;
//     name2idx = ds.name2idx;
//     // this->variableNames = std::vector<std::string>();
//     // this->name2idx = std::unordered_map<std::string, int>();
//     var2idx = std::unordered_map<Variable*, int>();
//     data = ds.data;

//     for (int j = 0; j < m; j++) {
// 	if (ds.variables.at(j)->isDiscrete())
// 	    variables.push_back(new DiscreteVariable(*((DiscreteVariable*) ds.variables[j])));
// 	else if (ds.variables.at(j)->isContinuous())
// 	    variables.push_back(new ContinuousVariable(*((ContinuousVariable*) ds.variables[j])));

// 	var2idx.insert(std::pair<Variable*, int>(variables[j], j));
//     }
  
//     // for (int j = 0; j < this->m; j++) {
//     //   if (ds.variables.at(j)->isDiscrete())
//     //     this->variables.push_back(new DiscreteVariable(*((DiscreteVariable*)ds.variables.at(j))));
//     //   else if (ds.variables.at(j)->isContinuous())
//     //     this->variables.push_back(new ContinuousVariable(*((ContinuousVariable*)ds.variables.at(j))));

//     //   var2idx.insert(std::pair<Variable*, int>(this->variables.at(j), j));
//     // }
  
// }

// DataSet& DataSet::operator=(DataSet&& ds) {
//     // this->maxDiscrete=ds.maxDiscrete;
//     // this->m = ds.m;
//     // this->n = ds.n;
//     // this->variables = ds.variables;
//     // this->variableNames = ds.variableNames;
//     // this->name2idx = ds.name2idx;
//     // this->var2idx = ds.var2idx;
//     // this->data = ds.data;

//     maxDiscrete=ds.maxDiscrete;
//     m = ds.m;
//     n = ds.n;
//     variables = std::vector<Variable*>();
//     variableNames = ds.variableNames;
//     name2idx = ds.name2idx;
//     // this->variableNames = std::vector<std::string>();
//     // this->name2idx = std::unordered_map<std::string, int>();
//     var2idx = std::unordered_map<Variable*, int>();
//     data = ds.data;

//     for (int j = 0; j < m; j++) {
// 	if (ds.variables.at(j)->isDiscrete())
// 	    variables.push_back(new DiscreteVariable(*((DiscreteVariable*) ds.variables[j])));
// 	else if (ds.variables.at(j)->isContinuous())
// 	    variables.push_back(new ContinuousVariable(*((ContinuousVariable*) ds.variables[j])));

// 	var2idx.insert(std::pair<Variable*, int>(variables[j], j));
//     }

//     // this->variables.clear();
//     // this->var2idx.clear();
  
//     // for (int j = 0; j < this->m; j++) {
//     //   if (ds.variables.at(j)->isDiscrete())
//     //     this->variables.push_back(new DiscreteVariable(*((DiscreteVariable*)ds.variables.at(j))));
//     //   else if (ds.variables.at(j)->isContinuous())
//     //     this->variables.push_back(new ContinuousVariable(*((ContinuousVariable*)ds.variables.at(j))));

//     //   var2idx.insert(std::pair<Variable*, int>(this->variables.at(j), j));
//     // }
//     return *this;
// }

std::vector<Node> DataSet::getContinuousVariables() {
    std::vector<Node> result = std::vector<Node>();

    for (int i = 0; i < m; i++) {
        if (variables[i].isContinuous())
            result.push_back(variables[i]);
    }  

    return result;
}

std::vector<Node> DataSet::getDiscreteVariables() {
    std::vector<Node> result = std::vector<Node>();

    for (int i = 0; i < m; i++) {
        if (variables[i].isDiscrete())
            result.push_back(variables[i]);
    } 

    return result;
}

std::vector<Node> DataSet::getCensoredVariables() {
    std::vector<Node> result = std::vector<Node>();

    for (int i = 0; i < m; i++) {
        if (variables[i].isCensored())
            result.push_back(variables[i]);
    } 

    return result;
}


bool DataSet::isMixed() {
    bool hasCont = false;
    bool hasDisc = false;

    for (const Node& var : variables) {
	if (!hasCont)
	    hasCont = var.isContinuous();
	
	if (!hasDisc)
	    hasDisc = var.isDiscrete();
	
	if (hasCont && hasDisc)
	    break;
    }
    
    return hasCont && hasDisc;
}

bool DataSet::isContinuous() {
    bool hasCont = false;
    bool hasDisc = false;

    for (const Node& var : variables) {
	if (!hasCont)
	    hasCont = var.isContinuous();
	
	if (!hasDisc)
	    hasDisc = var.isDiscrete();
	
	if (hasCont && hasDisc)
	    break;
    }
    
    return hasCont && !hasDisc;
}

bool DataSet::isDiscrete() {
    bool hasCont = false;
    bool hasDisc = false;

    for (const Node& var : variables) {
	if (!hasCont)
	    hasCont = var.isContinuous();
	
	if (!hasDisc)
	    hasDisc = var.isDiscrete();
	
	if (hasCont && hasDisc)
	    break;
    }
    
    return !hasCont && hasDisc;
}

bool DataSet::isCensored() {
    bool hasCens = false;

    for (const Node& var : variables) {
	if (!hasCens)
	    hasCens = var.isCensored();
	else
	    break;
    }
    
    return hasCens;
}



int DataSet::getInt(int row, int col) {

    // Node var = getVariable(col);
    
    if (!variables[col].isDiscrete()) {
	throw std::invalid_argument("Column indicated is not of type DISCRETE");
    }
    return (int) data(row, col);
}

// std::vector<Variable *> DataSet::copyContinuousVariables() {
//     std::vector<Variable *> result = std::vector<Variable *>();

//     for (int i = 0; i < m; i++) {
//         if (variables[i]->isContinuous())
//             result.push_back(new ContinuousVariable(*((ContinuousVariable*) variables[i])));
//     }  

//     return result;
// }

// std::vector<Variable *> DataSet::copyDiscreteVariables() {
//     std::vector<Variable *> result = std::vector<Variable *>();

//     for (int i = 0; i < m; i++) {
//         if (variables[i]->isDiscrete())
//             result.push_back(new DiscreteVariable(*((DiscreteVariable*) variables[i])));
//     } 

//     return result;
// }

arma::mat DataSet::getContinuousData() {
    std::vector<arma::uword> continuousColumns = std::vector<arma::uword>();

    for (arma::uword i = 0; i < m; i++) {
        if (variables[i].isContinuous())
            continuousColumns.push_back(i);
    }

    return data.cols(arma::uvec(continuousColumns));
}

arma::mat DataSet::getDiscreteData() {
    std::vector<arma::uword> discreteColumns = std::vector<arma::uword>();

    for (arma::uword i = 0; i < m; i++) {
        if (variables[i].isDiscrete())
            discreteColumns.push_back(i);
    }

    return data.cols(arma::uvec(discreteColumns));
}

arma::mat DataSet::getCensoredData() {
    std::vector<arma::uword> censoredColumns = std::vector<arma::uword>();

    for (arma::uword i = 0; i < m; i++) {
        if (variables[i].isCensored())
            censoredColumns.push_back(i);
    }

    return data.cols(arma::uvec(censoredColumns));
}


std::vector<int> DataSet::getDiscLevels() {
    std::vector<int> result;

    for (arma::uword i = 0; i < m; i++) {
        if (variables[i].isDiscrete())
            result.push_back(variables[i].getNumCategories());
    }

    return result;
}

// [[Rcpp::export]]
void DataSetTest(const Rcpp::DataFrame &df, const int maxDiscrete = 5) {
    DataSet ds = DataSet(df);
    // Rcpp::Rcout << ds;

    ds.dropMissing();

    std::vector<DataSet> dsVec;
    for (int i = 0; i < 5; i++)
	dsVec.push_back(ds);

    ds.addVariable(Node(new ContinuousVariable("Dummy1")));

    // Rcpp::Rcout << ds;

    int col = ds.getColumn(ds.getVariable("Dummy1"));
    for (int i = 0; i < ds.getNumRows(); i++)
        ds.set(i, col, ((double)i) / 100.0);

    // Rcpp::Rcout << ds;

    col = 3;
    ds.addVariable(col, Node(new DiscreteVariable("Dummy2", 4)));

    // Rcpp::Rcout << ds;

    col = ds.getColumn(ds.getVariable("Dummy2"));
    for (int i = 0; i < ds.getNumRows(); i++)
        ds.set(i, col, i % 4);

    // Rcpp::Rcout << ds;

    std::vector<Node> nodes;
    for (const Node& var : ds.getVariables()) {
	nodes.push_back(var);
    }

    
    for (const Node& n : nodes) {
	if (n.isContinuous())
            Rcpp::Rcout << "C:";
        else if (n.isDiscrete())
            Rcpp::Rcout << "D:";
	else if (n.isCensored())
            Rcpp::Rcout << "Cens:";
        Rcpp::Rcout << n.getName();
        Rcpp::Rcout << "\t";
    }
    Rcpp::Rcout << std::endl;

    Rcpp::Rcout << "Copying all nodes\n";
    std::vector<Node> nodes2(nodes);
    
    for (const Node& n : nodes2) {
	if (n.isContinuous())
            Rcpp::Rcout << "C:";
        else if (n.isDiscrete())
            Rcpp::Rcout << "D:";
	else if (n.isCensored())
            Rcpp::Rcout << "Cens:";
        Rcpp::Rcout << n.getName();
        Rcpp::Rcout << "\t";
    }
    Rcpp::Rcout << std::endl;

    Rcpp::Rcout << "Moving all nodes\n";
    std::vector<Node> nodes3 = std::move(nodes);
    
    for (const Node& n : nodes3) {
	if (n.isContinuous())
            Rcpp::Rcout << "C:";
        else if (n.isDiscrete())
            Rcpp::Rcout << "D:";
	else if (n.isCensored())
            Rcpp::Rcout << "Cens:";
        Rcpp::Rcout << n.getName();
        Rcpp::Rcout << "\t";
    }
    Rcpp::Rcout << std::endl << std::endl;

    // Rcpp::Rcout << ds << std::endl;

    // Rcpp::Rcout << "std::vector use test\n";

    ds.npnTransform();

    for (int i = 0; i < 5; i++)
	dsVec.push_back(ds);

    for (int i = 0; i < 10; i++) 
	Rcpp::Rcout << "DataSet " << i << ":\n" << dsVec[i] << std::endl;

}

std::ostream &operator<<(std::ostream &os, DataSet &ds) {
    for (int i = 0; i < ds.getNumColumns(); i++) {
        if (ds.variables[i].isContinuous())
            os << "C:";
        else if (ds.variables[i].isDiscrete())
            os << "D:";
	else if (ds.variables[i].isCensored())
		os << "Cens:";
        os << ds.variables[i].getName();
        os << "\t";
    }

    os << "\n";

    for (int i = 0; i < std::min(ds.getNumRows(), 10); i++) {
        for (int j = 0; j < ds.getNumColumns(); j++) {
            os << ds.data(i, j);
	    if (ds.variables[j].isCensored()) 
		if (!ds.variables[j].getCensor(i))
		    os << "+";
	    
            os << "\t";
        }
        os << "\n";
    }
    os << "(" << ds.getNumRows() << ", " << ds.getNumColumns() << ")\n";

    return os;
}

