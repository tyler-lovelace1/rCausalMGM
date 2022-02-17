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

DataSet::DataSet(const Rcpp::DataFrame &df, const int maxDiscrete)
{
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

        Rcpp::CharacterVector col = df[i];

        curName = (std::string)names[i];

        std::set<std::string> unique = getUnique(col);

        numUnique = unique.size();

        if (numUnique > maxDiscrete)
        {

            variables.push_back(Node(new ContinuousVariable(curName)));
        }
        else
        {

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
            else
            {
                throw std::runtime_error("invalid variable type");
            }
        }
    }
}

void DataSet::addVariable(Node v)
{
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

std::vector<int> DataSet::getDiscLevels() {
    std::vector<int> result;

    for (arma::uword i = 0; i < m; i++) {
        if (variables[i].isDiscrete())
            result.push_back(variables[i].getNumCategories());
    }

    return result;
}

//[[Rcpp::export]]
void DataSetTest(const Rcpp::DataFrame &df, const int maxDiscrete = 5) {
    DataSet ds = DataSet(df, maxDiscrete);
    // Rcpp::Rcout << ds;

    std::vector<DataSet> dsVec;
    for (int i = 0; i < 5; i++)
	dsVec.push_back(ds);

    ds.addVariable(Node(new ContinuousVariable("Y1")));

    // Rcpp::Rcout << ds;

    int col = ds.getColumn(ds.getVariable("Y1"));
    for (int i = 0; i < ds.getNumRows(); i++)
        ds.set(i, col, ((double)i) / 100.0);

    // Rcpp::Rcout << ds;

    col = 3;
    ds.addVariable(col, Node(new DiscreteVariable("Y2", 4)));

    // Rcpp::Rcout << ds;

    col = ds.getColumn(ds.getVariable("Y2"));
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
        Rcpp::Rcout << n.getName();
        Rcpp::Rcout << "\t";
    }
    Rcpp::Rcout << std::endl << std::endl;

    // Rcpp::Rcout << ds << std::endl;

    // Rcpp::Rcout << "std::vector use test\n";

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
        os << ds.variables[i].getName();
        os << "\t";
    }

    os << "\n";

    for (int i = 0; i < std::min(ds.getNumRows(), 10); i++)
    {
        for (int j = 0; j < ds.getNumColumns(); j++)
        {
            os << ds.data(i, j);
            os << "\t";
        }
        os << "\n";
    }
    os << "(" << ds.getNumRows() << ", " << ds.getNumColumns() << ")\n";

    return os;
}

