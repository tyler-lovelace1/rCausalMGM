#include "Bootstrap.hpp"

arma::umat Bootstrap::getBootstrapSamples() {

    if (subSize < 1)
        throw std::invalid_argument("Sample size must be > 0");
    
    arma::umat sampMat(B, subSize);

    for(arma::uword i = 0; i < B; i++) {
        arma::urowvec curSamp;// (
	    // arma::conv_to<arma::urowvec>::from(
	    // 	arma::floor(N * arma::rowvec(subSize, arma::fill::randu))
	    // 	)
	    // );
	Rcpp::Rcout << "Bootstrap sample " << i << std::endl;
	bool done = false;
	int attempts = 5000;
        while(!done) {
	    // Rcpp::Rcout << "Attempt " << 5000 - attempts << std::endl;
	    curSamp = arma::conv_to<arma::urowvec>::from(
		arma::floor(N * arma::rowvec(subSize, arma::fill::randu))
		);
            for (arma::uword j = 0; j < i; j++) {
                if (arma::all(curSamp == sampMat.row(j))) {
		    continue;
                }
            }
	    // Rcpp::Rcout << "Sample indices: " << curSamp << std::endl;
	    done = true;
	    DataSet subset(d, curSamp);
	    
	    if (StabilityUtils::checkForVariance(subset, d) != -1) {
		done = false;
	    }
	    attempts--;
	    if (attempts == 0) {
		// Rcpp::Rcout << "ERROR: Unable to find a subsampled dataset of size " << b << " where there are at least one category of every discrete variable" << std::endl;
		throw std::invalid_argument("Unable to find a bootstrapped dataset of size " + std::to_string(subSize) + " where there are at least two samples for each category of every discrete variable. The number of samples per subsampled dataset can be increased with the subSize parameter to address this problem.");
	    }
        }
        sampMat.row(i) = curSamp;
    }
    
    return sampMat;
}

EdgeListGraph Bootstrap::runBootstrap() {

    arma::umat subs = getBootstrapSamples();

    // Rcpp::Rcout << subs << std::endl;

    EdgeListGraph g;

    std::vector<EdgeListGraph> graphVec;

    for (arma::uword i = 0; i < B; i++) {

	DataSet subset(d, subs.row(i));

	if (alg == "mgm") {
	    
	    MGM mgm(subset, lambda);
	    g = mgm.search();
	    Rcpp::Rcout << g << std::endl;
	    
	}

	graphVec.push_back(g);

    }

    Rcpp::Rcout << graphVec.size() << std::endl;

    return makeEnsembleGraph(graphVec);
    
}


EdgeListGraph Bootstrap::makeEnsembleGraph(std::vector<EdgeListGraph>& graphVec) {

    Rcpp::Rcout << "makeEnsembleGraph\n";

    // string key: var1 var2
    // vec indices:  0    1    2    3    4    5    6
    //              ---  -->  <--  <->  o->  <-o  o-o
    std::unordered_map<std::string, arma::vec> edgeFreq;

    EdgeListGraph ensembleGraph(d.getVariables());

    for (EdgeListGraph g : graphVec) {
	Rcpp::Rcout << g << std::endl;
	
	for (Edge e : g.getEdges()) {

	    Rcpp::Rcout << e << std::endl;
	    Rcpp::Rcout << e.toString() << std::endl;

	    bool node1first = e.getNode1()->getName() < e.getNode2()->getName();
	    
	    std::string key = (node1first) ?
		e.getNode1()->getName() + " " + e.getNode2()->getName() :
		e.getNode2()->getName() + " " + e.getNode1()->getName();

	    int offset = (node1first) ? 0 : 1;

	    int typeIdx;

	    if (Edge::isUndirectedEdge(e))                typeIdx = 0;
	    else if (Edge::isDirectedEdge(e))             typeIdx = 1 + offset;
	    else if (Edge::isBidirectionalEdge(e))        typeIdx = 3;
	    else if (Edge::isPartiallyOrientedEdge(e))    typeIdx = 4 + offset;
	    else if (Edge::isNondirectedEdge(e))          typeIdx = 6;
	    else throw std::runtime_error("Invalid edge type in " + e.toString());
	    
	    Rcpp::Rcout << "key: " << key << std::endl;
	    if (!edgeFreq.count(key)) {
		edgeFreq[key] = arma::vec(7, arma::fill::zeros);
	    }
	    edgeFreq[key](typeIdx)++;
	}
    }
}
