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
    this->var2idx = std::unordered_map<Variable *, int>();
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

            variables.push_back(new ContinuousVariable(curName));
        }
        else
        {

            std::vector<std::string> categories;

            for (std::set<std::string>::iterator it = unique.begin(); it != unique.end(); it++)
                categories.push_back(*it);

            std::sort(categories.begin(), categories.end());

            variables.push_back(new DiscreteVariable(curName, categories));
        }

        variableNames.push_back(curName);
        name2idx.insert(std::pair<std::string, int>(curName, i));
        var2idx.insert(std::pair<Variable *, int>(variables[i], i));

        for (int j = 0; j < n; j++)
        {

            val = (std::string)col[j];

            if (variables[i]->isMissingValue(val))
            {
                throw std::runtime_error("Missing values not yet implemented");
            }
            else if (variables[i]->isContinuous())
            {
                if (((ContinuousVariable *)variables[i])->checkValue(val))
                {
                    this->data(j, i) = std::stod(val);
                }
                else
                {
                    throw std::runtime_error("invalid value encountered in row " + std::to_string(j) + " of variable " + curName);
                }
            }
            else if (variables[i]->isDiscrete())
            {
                if (((DiscreteVariable *)variables[i])->checkValue(val))
                {
                    this->data(j, i) = (double)((DiscreteVariable *)variables[i])->getIndex(val);
                }
                else
                {
                    throw std::runtime_error("invalid value encountered in row " + std::to_string(j) + " of variable " + curName);
                }
            }
            else
            {
                throw std::runtime_error("invalid variable type");
            }
        }
    }
}

void DataSet::addVariable(Variable *v)
{
    variables.push_back(v);
    arma::vec col = arma::vec(n);
    for (int i = 0; i < n; i++)
        col[i] = arma::datum::nan;

    variableNames.push_back(v->getName());
    name2idx.insert(std::pair<std::string, int>(v->getName(), m));
    var2idx.insert(std::pair<Variable *, int>(v, m));

    data.insert_cols(m, col);
    m++;
}

void DataSet::addVariable(int i, Variable *v)
{
    if (i > m || i < 0)
        throw std::invalid_argument("index " + std::to_string(i) + " out of bounds");

    variables.insert(variables.begin() + i, v);
    arma::vec col = arma::vec(n);
    for (int i = 0; i < n; i++)
        col[i] = arma::datum::nan;

    variableNames.push_back(v->getName());
    name2idx.insert(std::pair<std::string, int>(v->getName(), i));
    var2idx.insert(std::pair<Variable *, int>(v, i));

    data.insert_cols(i, col);
    m++;
}

DataSet::~DataSet()
{
    for (int i = 0; i < variables.size(); i++)
        delete variables[i];
}


// Deep copy
std::vector<Variable *> DataSet::copyVariables() {
    std::vector<Variable *> result = std::vector<Variable *>(m);

    for (int i = 0; i < m; i++) {
        if (variables[i]->isContinuous())
            result[i] = new ContinuousVariable(*((ContinuousVariable*) variables[i]));
        else if (variables[i]->isDiscrete())
            result[i] = new DiscreteVariable(*((DiscreteVariable*) variables[i]));
    }

    return result;
}

std::vector<Variable *> DataSet::getContinuousVariables() {
    std::vector<Variable *> result = std::vector<Variable *>();

    for (int i = 0; i < m; i++) {
        if (variables[i]->isContinuous())
            result.push_back(variables[i]);
    }  

    return result;
}

std::vector<Variable *> DataSet::getDiscreteVariables() {
    std::vector<Variable *> result = std::vector<Variable *>();

    for (int i = 0; i < m; i++) {
        if (variables[i]->isDiscrete())
            result.push_back(variables[i]);
    } 

    return result;
}

std::vector<Variable *> DataSet::copyContinuousVariables() {
    std::vector<Variable *> result = std::vector<Variable *>();

    for (int i = 0; i < m; i++) {
        if (variables[i]->isContinuous())
            result.push_back(new ContinuousVariable(*((ContinuousVariable*) variables[i])));
    }  

    return result;
}

std::vector<Variable *> DataSet::copyDiscreteVariables() {
    std::vector<Variable *> result = std::vector<Variable *>();

    for (int i = 0; i < m; i++) {
        if (variables[i]->isDiscrete())
            result.push_back(new DiscreteVariable(*((DiscreteVariable*) variables[i])));
    } 

    return result;
}

arma::mat DataSet::getContinuousData() {
    std::vector<arma::uword> continuousColumns = std::vector<arma::uword>();

    for (arma::uword i = 0; i < m; i++) {
        if (variables[i]->isContinuous())
            continuousColumns.push_back(i);
    }

    return data.cols(arma::uvec(continuousColumns));
}

arma::mat DataSet::getDiscreteData() {
    std::vector<arma::uword> discreteColumns = std::vector<arma::uword>();

    for (arma::uword i = 0; i < m; i++) {
        if (variables[i]->isDiscrete())
            discreteColumns.push_back(i);
    }

    return data.cols(arma::uvec(discreteColumns));
}

std::vector<int> DataSet::getDiscLevels() {
    std::vector<int> result;

    for (arma::uword i = 0; i < m; i++) {
        if (variables[i]->isDiscrete())
            result.push_back(((DiscreteVariable*) variables[i])->getNumCategories());
    }

    return result;
}

// [[Rcpp::export]]
void DataSetTest(const Rcpp::DataFrame &df, const int maxDiscrete = 5)
{
    DataSet ds = DataSet(df, maxDiscrete);
    Rcpp::Rcout << ds;

    ds.addVariable(new ContinuousVariable("Y1"));

    Rcpp::Rcout << ds;

    int col = ds.getColumn(ds.getVariable("Y1"));
    for (int i = 0; i < ds.getNumRows(); i++)
        ds.set(i, col, ((double)i) / 100.0);

    Rcpp::Rcout << ds;

    col = 3;
    ds.addVariable(col, new DiscreteVariable("Y2", 4));

    Rcpp::Rcout << ds;

    col = ds.getColumn(ds.getVariable("Y2"));
    for (int i = 0; i < ds.getNumRows(); i++)
        ds.set(i, col, i % 4);

    Rcpp::Rcout << ds;
}

std::ostream &operator<<(std::ostream &os, DataSet &ds)
{
    for (int i = 0; i < ds.getNumColumns(); i++)
    {
        if (ds.variables[i]->isContinuous())
            os << "C:";
        else if (ds.variables[i]->isDiscrete())
            os << "D:";
        os << ds.variables[i]->getName();
        os << "\t";
    }

    os << "\n";

    for (int i = 0; i < ds.getNumRows(); i++)
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
