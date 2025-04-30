#include "Bootstrap.hpp"

// arma::umat Bootstrap::getBootstrapSamples() {

//     if (N < 1)
//         throw std::invalid_argument("Sample size must be > 0");
    
//     arma::umat sampMat(B, N);

//     for(arma::uword i = 0; i < B; i++) {
//         arma::urowvec curSamp;// (
// 	    // arma::conv_to<arma::urowvec>::from(
// 	    // 	arma::floor(N * arma::rowvec(subSize, arma::fill::randu))
// 	    // 	)
// 	    // );
// 	// Rcpp::Rcout << "Bootstrap sample " << i << std::endl;
// 	bool done = false;
// 	int attempts = 5000;
//         while(!done) {
// 	    // Rcpp::Rcout << "Attempt " << 5000 - attempts << std::endl;
// 	    curSamp = arma::conv_to<arma::urowvec>::from(
// 		arma::floor(N * arma::rowvec(N, arma::fill::randu))
// 		);
//             for (arma::uword j = 0; j < i; j++) {
//                 if (arma::all(curSamp == sampMat.row(j))) {
// 		    continue;
//                 }
//             }
// 	    // Rcpp::Rcout << "Sample indices: " << curSamp << std::endl;
// 	    done = true;
// 	    DataSet subset(d, curSamp);
	    
// 	    if (StabilityUtils::checkForVariance(subset, d) != -1) {
// 		done = false;
// 	    }
// 	    attempts--;
// 	    if (attempts == 0) {
// 		// Rcpp::Rcout << "ERROR: Unable to find a subsampled dataset of size " << b << " where there are at least one category of every discrete variable" << std::endl;
// 		throw std::invalid_argument("Unable to find a bootstrapped dataset of size " + std::to_string(N) + " where there are at least two samples for each category of every discrete variable. The number of samples per subsampled dataset can be increased with the subSize parameter to address this problem.");
// 	    }
//         }
//         sampMat.row(i) = curSamp;
//     }
    
//     return sampMat;
// }

EdgeListGraph Bootstrap::runBootstrap() {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // arma::umat subs;
    
    if (replace)   samps = StabilityUtils::subSampleWithReplacement(d, N, B);
    else           samps = StabilityUtils::subSampleNoReplacement(d, N, B);

    if (verbose) Rcpp::Rcout << "  Bootstrapping..." << std::endl;

    // EdgeListGraph g;

    std::vector<EdgeListGraph> graphVec;

    bool mgmInit = alg.find("mgm") != std::string::npos;
    bool censFlag = d.isCensored();

    for (arma::uword i = 0; i < B; i++) {

	if (verbose) Rcpp::Rcout << "\r    Running bootstrap " << i+1 << "...";

	DataSet subset(d, samps.row(i));

	// Rcpp::Rcout << subset << std::endl;

	if (alg == "mgm") {
	    if (!censFlag) {
		MGM mgm(subset, lambda);
		graphVec.push_back(mgm.search());
	    } else {
		CoxMGM coxmgm(subset, lambda);
		graphVec.push_back(coxmgm.search());
	    }
	    // Rcpp::Rcout << g << std::endl;
	    
	} else {
	    EdgeListGraph ig;
	    IndTestMultiCox itm(subset, alpha);
	    
	    if (mgmInit) {
		// MGM mgm(subset, lambda);
	        // ig = mgm.search();
		if (!censFlag) {
		    MGM mgm(subset, lambda);
		    ig = mgm.search();	       
		} else {
		    // Rcpp::Rcout << "Starting CoxMGM...\n";
		    CoxMGM coxmgm(subset, lambda);
		    ig = coxmgm.search();
		    // Rcpp::Rcout << "CoxMGM Complete\n";
		}
	    }

	    if (alg.find("pc") != std::string::npos) {
		
		PcStable causalAlg((IndependenceTest*) &itm);
		if (threads > 0) causalAlg.setThreads(threads);

		causalAlg.setKnowledge(this->knowledge);
		causalAlg.setVerbose(false);
		causalAlg.setFDR(false);
		causalAlg.setOrientRule(orientRule);
	    
		if (mgmInit) {
		    causalAlg.setInitialGraph(&ig);
		}
	    
		graphVec.push_back(causalAlg.search());
		
	    } else if (alg.find("fci") != std::string::npos) {
		
		Fci causalAlg((IndependenceTest*) &itm);
		if (threads > 0) causalAlg.setThreads(threads);

		causalAlg.setKnowledge(this->knowledge);
		causalAlg.setVerbose(false);
		causalAlg.setFDR(false);
		causalAlg.setPossibleDsepSearchDone(false);
		causalAlg.setOrientRule(orientRule);
	    
		if (mgmInit) {
		    causalAlg.setInitialGraph(&ig);
		}
	    
		graphVec.push_back(causalAlg.search());

		// Rcpp::Rcout << "FCI Complete\n";
	    } else {
		throw std::runtime_error("Invalid algorithm selected");
	    }   
	}

	// graphVec.push_back(g);

    }

    // Rcpp::Rcout << graphVec.size() << std::endl;

    EdgeListGraph ensGraph = makeEnsembleGraph(graphVec);

    if (verbose) Rcpp::Rcout << std::endl << "  Finished" << std::endl;

    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    std::string fullAlg = graphVec.at(0).getAlgorithm() + " Bootstrap Ensemble";

    ensGraph.setAlgorithm(fullAlg);

    std::string graphType = graphVec.at(0).getGraphType();
    
    ensGraph.setGraphType(graphType);

    if (mgmInit) {
	ensGraph.setHyperParam("lambda", lambda);
    }

    if (alg != "mgm") {
	ensGraph.setHyperParam("alpha", { alpha });
    }

    return ensGraph;
    
}


EdgeListGraph Bootstrap::makeEnsembleGraph(std::vector<EdgeListGraph>& graphVec) {

    // Rcpp::Rcout << "makeEnsembleGraph\n";

    // vec indices:  0    1    2    3    4    5    6    7
    //              ---  -->  <--  <->  o->  <-o  o-o  None
    std::map<NodePair, arma::rowvec> edgeFreq;

    EdgeListGraph ensembleGraph(d.getVariables());

    for (auto&& g : graphVec) {
	//Rcpp::Rcout << g << std::endl;
	
	for (auto&& e : g.getEdges()) {

	    Node n1 = e.getNode1();
	    Node n2 = e.getNode2();

	    int offset = (n1 < n2) ? 0 : 1;

	    int typeIdx;

	    if (Edge::isUndirectedEdge(e))                typeIdx = 0;
	    else if (Edge::isDirectedEdge(e))             typeIdx = 1 + offset;
	    else if (Edge::isBidirectionalEdge(e))        typeIdx = 3;
	    else if (Edge::isPartiallyOrientedEdge(e))    typeIdx = 4 + offset;
	    else if (Edge::isNondirectedEdge(e))          typeIdx = 6;
	    else {
		std::ostringstream os;
		os << "Invalid edge type in " << e;
		throw std::runtime_error(os.str());
	    }

	    // Rcpp::Rcout << e << ":   " << typeIdx << std::endl;

	    NodePair pair = std::minmax(n1, n2);
	    
	    // Rcpp::Rcout << "key: " << key << std::endl;
	    if (!edgeFreq.count(pair)) {
		edgeFreq[pair] = arma::rowvec(8, arma::fill::zeros);
	    }
	    edgeFreq[pair](typeIdx)++;
	}
    }
    
    for (auto it = edgeFreq.begin(); it!=edgeFreq.end(); it++) {
	Node n1 = it->first.first;
	Node n2 = it->first.second;
	
	// if (n1 > n2)
	//     std::swap(n1, n2);

	it->second /= B;
	
	// Rcpp::Rcout << it->first.first << "  " << it->first.second << "," << it->second;
	
	double adjFreq = arma::sum(it->second);
	it->second(7) = 1 - adjFreq;
	
	// int maxIdx = -1;
	// if (adjFreq >= thresh) {
	int maxIdx = it->second.index_max();
	
	// if (maxIdx == 7) {
	//     if (adjFreq > 0.5) {
	// 	if (useNondirected)
	// 	    ensembleGraph.addNondirectedEdge(n1, n2);
	// 	else
	// 	    ensembleGraph.addUndirectedEdge(n1, n2);
	//     }
	//     continue;
	// }
	
	switch (maxIdx) {
	case 0:
	    ensembleGraph.addUndirectedEdge(n1, n2);
	    break;
	case 1:
	    // if (ensemble == "majority" && it->second(maxIdx) < 0.5) {
	    // 	if (useNondirected)
	    // 	    ensembleGraph.addNondirectedEdge(n1, n2);
	    // 	else
	    // 	    ensembleGraph.addUndirectedEdge(n1, n2);
	    // } else
	    ensembleGraph.addDirectedEdge(n1, n2);
	    break;
	case 2:
	    // if (ensemble == "majority" && it->second(maxIdx) < 0.5) {
	    // 	if (useNondirected)
	    // 	    ensembleGraph.addNondirectedEdge(n1, n2);
	    // 	else
	    // 	    ensembleGraph.addUndirectedEdge(n1, n2);
	    // } else
	    ensembleGraph.addDirectedEdge(n2, n1);
	    break;
	case 3:
	    // if (ensemble == "majority" && it->second(maxIdx) < 0.5) {
	    // 	if (useNondirected)
	    // 	    ensembleGraph.addNondirectedEdge(n1, n2);
	    // 	else
	    // 	    ensembleGraph.addUndirectedEdge(n1, n2);
	    // } else
	    ensembleGraph.addBidirectedEdge(n1, n2);
	    break;
	case 4:
	    // if (ensemble == "majority" && it->second(maxIdx) < 0.5) {
	    // 	if (useNondirected)
	    // 	    ensembleGraph.addNondirectedEdge(n1, n2);
	    // 	else
	    // 	    ensembleGraph.addUndirectedEdge(n1, n2);
	    // } else
	    ensembleGraph.addPartiallyOrientedEdge(n1, n2);
	    break;
	case 5:
	    // if (ensemble == "majority" && it->second(maxIdx) < 0.5) {
	    // 	if (useNondirected)
	    // 	    ensembleGraph.addNondirectedEdge(n1, n2);
	    // 	else
	    // 	    ensembleGraph.addUndirectedEdge(n1, n2);
	    // } else
	    ensembleGraph.addPartiallyOrientedEdge(n2, n1);
	    break;
	case 6:
	    // if (ensemble == "majority" && it->second(maxIdx) < 0.5) {
	    //     if (useNondirected)   ensembleGraph.addNondirectedEdge(n1, n2);
	    //     else                  ensembleGraph.addUndirectedEdge(n1, n2);
	    // }
	    ensembleGraph.addNondirectedEdge(n1, n2);
	    break;
	}
	// }
    }

    // std::vector<Edge> ensEdgeList = ensembleGraph.getEdgeList();

    // // Edge::sortEdges(ensEdgeList);

    // std::sort(ensEdgeList.begin(), ensEdgeList.end(),
    // 	      [](Edge e1, Edge e2) {
    // 		  // Node e1n1 = e1.getNode1();
    // 		  // Node e1n2 = e1.getNode2();
    // 		  // Node e2n1 = e2.getNode1();
    // 		  // Node e2n2 = e2.getNode2();
		  
    // 		  NodePair pair1 = std::minmax(e1.getNode1(), e1.getNode2());
    // 		  NodePair pair2 = std::minmax(e2.getNode1(), e2.getNode2());

    // 		  // if (e1n1 > e1n2) std::swap(e1n1, e1n2);
    // 		  // if (e2n1 > e2n2) std::swap(e2n1, e2n2);

    // 		  // if (e1n1 == e2n1)
    // 		  //     return e1n2 < e2n2;
    // 		  return pair1 < pair2;
    // 	      });

    std::vector<NodePair> ensEdgeList;
    ensEdgeList.reserve(edgeFreq.size());

    for (auto&& kv : edgeFreq) {
	ensEdgeList.push_back(kv.first);
    }

    std::sort(ensEdgeList.begin(), ensEdgeList.end(),
	      [&] (const NodePair& pair1, const NodePair& pair2) {
		  if (edgeFreq.at(pair1)(7) != edgeFreq.at(pair2)(7))
		      return edgeFreq.at(pair1)(7) < edgeFreq.at(pair2)(7);
		  return pair1 < pair2;			  
	      });

    Rcpp::StringVector var1 = {};
    Rcpp::StringVector var2 = {};
    Rcpp::StringVector interaction = {};
    Rcpp::NumericVector undir = {};
    Rcpp::NumericVector rdir = {};
    Rcpp::NumericVector ldir = {};
    Rcpp::NumericVector bidir = {};
    Rcpp::NumericVector rpdir = {};
    Rcpp::NumericVector lpdir = {};
    Rcpp::NumericVector nondir = {};
    Rcpp::NumericVector none = {};

    for (NodePair pair : ensEdgeList) {
	Node n1 = pair.first;
	Node n2 = pair.second;

	// bool flipDir = (n1 < n2) ? false : true;

	// NodePair pair = std::minmax(n1, n2);

	// if (n1 > n2)
	//     std::swap(n1, n2);
	
        var1.push_back(pair.first.getName());
	var2.push_back(pair.second.getName());

	if (ensembleGraph.isAdjacentTo(n1, n2)) {
	    Edge e = ensembleGraph.getEdge(n1, n2);
	    bool flipDir = (ensembleGraph.isDirectedFromTo(n1, n2)) ? false : true;
	    if (Edge::isUndirectedEdge(e)) {
		interaction.push_back("---");
	    } else if (Edge::isDirectedEdge(e)) {
		if (flipDir)
		    interaction.push_back("<--");
		else
		    interaction.push_back("-->");
	    }
	    else if (Edge::isBidirectionalEdge(e)) {
		interaction.push_back("<->");
	    } else if (Edge::isPartiallyOrientedEdge(e)) {
		if (flipDir)
		    interaction.push_back("<-o");
		else
		    interaction.push_back("o->");
	    } else if (Edge::isNondirectedEdge(e)) {
		interaction.push_back("o-o");
	    } else {
		throw std::runtime_error("Invalid edge type in " + e.toString());
	    }
	} else {
	    interaction.push_back("");
	}

        undir.push_back(edgeFreq[pair](0));
	rdir.push_back(edgeFreq[pair](1));
	ldir.push_back(edgeFreq[pair](2));
	bidir.push_back(edgeFreq[pair](3));
	rpdir.push_back(edgeFreq[pair](4));
	lpdir.push_back(edgeFreq[pair](5));
	nondir.push_back(edgeFreq[pair](6));
	none.push_back(edgeFreq[pair](7));
	    
    }

    stabs = Rcpp::DataFrame::create(Rcpp::_["var1"]=var1,
				    Rcpp::_["interaction"]=interaction,
				    Rcpp::_["var2"]=var2,
				    Rcpp::_["undir"]=undir,
				    Rcpp::_["right.dir"]=rdir,
				    Rcpp::_["left.dir"]=ldir,
				    Rcpp::_["bidir"]=bidir,
				    Rcpp::_["right.partdir"]=rpdir,
				    Rcpp::_["left.partdir"]=lpdir,
				    Rcpp::_["nondir"]=nondir,
				    Rcpp::_["none"]=none);
    
    return ensembleGraph;
}
