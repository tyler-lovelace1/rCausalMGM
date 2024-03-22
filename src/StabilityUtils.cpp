#include "StabilityUtils.hpp"

arma::umat StabilityUtils::subSampleNoReplacement(DataSet& data, int subSize, int numSub) {

    int sampSize = data.getNumRows();
    
    if (subSize <= 1)
        throw std::invalid_argument("Subsample size must be > 1");

    if (subSize >= sampSize) 
	throw std::invalid_argument("Subample size must be < " + std::to_string(sampSize));

    arma::urowvec indices(sampSize);
    for (arma::uword i = 0; i < sampSize; i++) {
        indices(i) = i;
    }

    int attempts = 5000;
    
    arma::umat sampMat(numSub, subSize);

    for(arma::uword i = 0; i < numSub; i++) {
        arma::urowvec curSamp;
        while(true) {
            SAMP:   
	    if (attempts == 0) {
		// Rcpp::Rcout << "ERROR: Unable to find a subsampled dataset of size " << b << " where there are at least one category of every discrete variable" << std::endl;
		throw std::invalid_argument("Unable to find a subsampled dataset of size " + std::to_string(subSize) + " where there are at least two samples for each category of every discrete variable and each continuous variable has a variance > 0. The number of samples per subsampled dataset can be increased with the subSize parameter to address this problem.");
	    }
	    indices = arma::shuffle(indices);
            curSamp = indices.head(subSize);
            for (arma::uword j = 0; j < i; j++) {
                if (arma::all(curSamp == sampMat.row(j))) {
                    goto SAMP;
                }
		DataSet subset(data, curSamp);
		if (checkForVariance(subset, data) != -1) {
		    attempts--;
		    goto SAMP;
		}
            }
            break;
        }
        sampMat.row(i) = curSamp;
    }
    
    return sampMat;
}

arma::umat StabilityUtils::subSampleLOO(DataSet& data) {
    arma::umat sampMat(data.getNumRows(), data.getNumRows()-1);

    for(int i = 0; i < sampMat.n_rows; i++) {
        int count = 0;
        for(int j = 0; j < sampMat.n_rows; j++) {
            if(i==j) goto A;
            sampMat(i, count) = j;
            count++;
            A:;
        }
    }
    
    // Print a warning if there's a variance issue
    for (int s = 0; s < sampMat.n_rows; s++) {
        DataSet dataSubSamp(data, sampMat.row(s));
        if (checkForVariance(dataSubSamp, data) != -1) {
            Rcpp::Rcout << "Variance issue with dataset: " << s << std::endl;
        }
    }

    return sampMat;
}

// arma::urowvec StabilityUtils::subSampleIndices(int N, int subSize) {

//     arma::urowvec indices(N);
//     for (arma::uword i = 0; i < N; i++) {
//         indices(i) = i;
//     }

//     indices = arma::shuffle(indices);
//     arma::urowvec samp(subSize);
//     for (arma::uword i = 0; i < subSize; i++) {
//         samp(i) = indices(i);
//     }

//     return samp;
// }


// arma::umat StabilityUtils::subSampleWithReplacement(int sampSize, int subSize, int numSub) {

//     if (subSize < 1)
//         throw std::invalid_argument("Sample size must be > 0");
  
//     arma::umat sampMat(numSub, subSize);

//     for(arma::uword i = 0; i < numSub; i++) {
//         arma::urowvec curSamp(arma::floor(sampSize * arma::vec(subSize, fill::randu)));
//         while(true) {
//             SAMP:
//             for (arma::uword j = 0; j < i; j++) {
//                 if (arma::all(curSamp == sampMat.row(j))) {
// 		    curSamp = arma::floor(sampSize * arma::vec(subSize, fill::randu))
//                     goto SAMP;
//                 }
//             }
//             break;
//         }
//         sampMat.row(i) = curSamp;
//     }
    
//     return sampMat;
// }


arma::umat StabilityUtils::subSampleWithReplacement(DataSet& data, int subSize, int numSub) {

    int sampSize = data.getNumRows();
    
    if (subSize <= 1)
        throw std::invalid_argument("Subsample size must be > 1");

    int attempts = 5000;
    
    arma::umat sampMat(numSub, subSize);

    for(arma::uword i = 0; i < numSub; i++) {
        arma::urowvec curSamp;
        while(true) {
            SAMP:   
	    if (attempts == 0) {
		// Rcpp::Rcout << "ERROR: Unable to find a subsampled dataset of size " << b << " where there are at least one category of every discrete variable" << std::endl;
		throw std::invalid_argument("Unable to find a subsampled dataset of size " + std::to_string(subSize) + " where there are at least two samples for each category of every discrete variable and each continuous variable has a variance > 0. The number of samples per subsampled dataset can be increased with the subSize parameter to address this problem.");
	    }
            curSamp = arma::conv_to<arma::urowvec>::from(arma::floor(sampSize * arma::vec(subSize, arma::fill::randu)));
            for (arma::uword j = 0; j < i; j++) {
                if (arma::all(curSamp == sampMat.row(j))) {
                    goto SAMP;
                }
		DataSet subset(data, curSamp);
		if (checkForVariance(subset, data) != -1) {
		    attempts--;
		    goto SAMP;
		}
            }
            break;
        }
        sampMat.row(i) = curSamp;
    }
    
    return sampMat;
}

//Get subsample size given sample size
//Uses the heuristic given in Learning Mixed Graphical Models with Separate Sparsity Parameters and Stability-Based Model Selection
int StabilityUtils::getSubSize(int sampleSize) {
    int b = (int) std::floor(10*std::sqrt(sampleSize));
    if (b >= sampleSize) {
        b = 3*sampleSize/4;
    }

    return b;
}

//returns undirected skeleton matrix (symmetric)
arma::mat StabilityUtils::skeletonToMatrix(EdgeListGraph& graph, DataSet& d) {
    int n = graph.getNumNodes();
    arma::mat matrix(n, n, arma::fill::zeros);

    // map nodes in order of appearance
    std::unordered_map<Node, int> map;
    for (Node node : graph.getNodes()) {
        map[node] = d.getColumn(node);
    }

    // mark edges
    for (Edge edge : graph.getEdges()) {
        // if directed find which is parent/child
        Node node1 = edge.getNode1();
        Node node2 = edge.getNode2();

        matrix(map[node1], map[node2]) = 1.0;
        matrix(map[node2], map[node1]) = 1.0;
    }

    return matrix;
}


//returns cube encoding of graph for stars
arma::cube StabilityUtils::graphToCube(EdgeListGraph& graph, DataSet& d) {
    int n = graph.getNumNodes();
    int m = 1;
    
    if (graph.getGraphType() == "markov equivalence class")
	m = 2;
    else if (graph.getGraphType() == "partial ancestral graph")
	m = 3;
    
    arma::cube cube(n, n, m, arma::fill::zeros);

    // map nodes in order of appearance
    std::unordered_map<Node, int> map;
    for (Node node : graph.getNodes()) {
        map[node] = d.getColumn(node);
    }

    
    if (m==1) { // undirected graph or unknown graph type: considers only skeleton
	// mark edges
	for (Edge edge : graph.getEdges()) {
	    // if directed find which is parent/child
	    Node node1 = edge.getNode1();
	    Node node2 = edge.getNode2();

	    int i = (map[node1] < map[node2]) ? map[node1] : map[node2];
	    int j = (map[node1] > map[node2]) ? map[node1] : map[node2];

	    Rcpp::Rcout << "  ( " << i << ", " << j << " )\n";

	    cube(i, j, 0) = 1.0;
	}
    } else if (m==2) { // cpdag graph: considers skeleton and orientations
	// mark edges
	for (Edge edge : graph.getEdges()) {
	    // if directed find which is parent/child
	    Node node1 = edge.getNode1();
	    Node node2 = edge.getNode2();
	    int i, j, k;

	    if (edge.getEndpoint1() == ENDPOINT_TAIL
		&& edge.getEndpoint2() == ENDPOINT_ARROW) {
		i = map[node2];
		j = map[node1];
		k = 0;
	    } else if (edge.getEndpoint1() == ENDPOINT_ARROW
		       && edge.getEndpoint2() == ENDPOINT_TAIL) {
		i = map[node1];
		j = map[node2];
		k = 0;
	    } else if (edge.getEndpoint1() == ENDPOINT_TAIL
		       && edge.getEndpoint2() == ENDPOINT_TAIL) {
		i = (map[node1] < map[node2]) ? map[node1] : map[node2];
		j = (map[node1] > map[node2]) ? map[node1] : map[node2];
		k = 1;
	    } else if (edge.getEndpoint1() == ENDPOINT_ARROW
		       && edge.getEndpoint2() == ENDPOINT_ARROW) {
		i = (map[node1] < map[node2]) ? map[node1] : map[node2];
		j = (map[node1] > map[node2]) ? map[node1] : map[node2];
		k = 1;
	    } else {
		throw std::runtime_error("Invalid edge type for CPDAGs encountered in graphToCube");
	    }

	    cube(i, j, k) = 1.0;
	}
	
    } else if (m==3) { // pag graph: considers skeleton and orientations
	// mark edges
	for (Edge edge : graph.getEdges()) {
	    // if directed find which is parent/child
	    Node node1 = edge.getNode1();
	    Node node2 = edge.getNode2();
	    int i, j, k;

	    if (edge.getEndpoint1() == ENDPOINT_TAIL
		&& edge.getEndpoint2() == ENDPOINT_ARROW) {
		i = map[node2];
		j = map[node1];
		k = 0;
	    } else if (edge.getEndpoint1() == ENDPOINT_ARROW
		       && edge.getEndpoint2() == ENDPOINT_TAIL) {
		i = map[node1];
		j = map[node2];
		k = 0;
	    } else if (edge.getEndpoint1() == ENDPOINT_CIRCLE
		&& edge.getEndpoint2() == ENDPOINT_ARROW) {
		i = map[node2];
		j = map[node1];
		k = 1;
	    } else if (edge.getEndpoint1() == ENDPOINT_ARROW
		       && edge.getEndpoint2() == ENDPOINT_CIRCLE) {
		i = map[node1];
		j = map[node2];
		k = 1;
	    } else if (edge.getEndpoint1() == ENDPOINT_CIRCLE
		       && edge.getEndpoint2() == ENDPOINT_CIRCLE) {
		i = (map[node1] < map[node2]) ? map[node1] : map[node2];
		j = (map[node1] > map[node2]) ? map[node1] : map[node2];
		k = 2;
	    } else if (edge.getEndpoint1() == ENDPOINT_ARROW
		       && edge.getEndpoint2() == ENDPOINT_ARROW) {
		i = (map[node1] > map[node2]) ? map[node1] : map[node2];
		j = (map[node1] < map[node2]) ? map[node1] : map[node2];
		k = 2;
	    } else {
		throw std::runtime_error("Invalid edge type for PAGs encountered in graphToCube");
	    }

	    cube(i, j, k) = 1.0;
	}
    }

    return cube;
}

int StabilityUtils::checkForVariance(DataSet& d, DataSet& full) {
    arma::mat t = d.getData();
    for (arma::uword i = 0; i < d.getNumColumns(); i++) {
        if (d.getVariable(i).isContinuous()) {
            arma::vec curr = standardizeData(t.col(i));

            double var = arma::var(curr);

            if (var <= 0.00001) {
		Rcpp::Rcout << "   " << d.getVariable(i).getName() << ":  " << var << std::endl;
                return i;
            }

        } else if (d.getVariable(i).isDiscrete()) {

	    if (d.getVariable(i).getNumCategories() < 2) {
		Rcpp::Rcout << "    " << d.getVariable(i).getName() << " has only one category\n";
		return i;
	    }
	    
            std::unordered_map<int, int> cats;
            for (arma::uword j = 0; j < full.getNumRows(); j++) {
                cats[full.getInt(j, i)] = 0;
            }

            for (arma::uword j = 0; j < d.getNumRows(); j++) {
                int currCat = d.getInt(j, i);
                if (cats.find(currCat) == cats.end()) {
                    throw std::invalid_argument("Found a category not in the full dataset");
                } else {
                    cats[currCat] += 1;
                }
            }

            for (std::pair<int, int> element : cats) {
                if (element.second < 3) {
		    Rcpp::Rcout << "    " << d.getVariable(i).getName() << ":  {";
		    for (std::pair<int, int> item : cats) {
			Rcpp::Rcout << item.first << " : " << item.second << ", ";
		    }
		    Rcpp::Rcout << "}\n";
                    return i;
                }
            }
        } else if (d.getVariable(i).isCensored()) {
	    if (d.getVariable(i).getNEvents() < 5) {
		Rcpp::Rcout << "   " << d.getVariable(i).getName() << ":  "
			    << d.getVariable(i).getNEvents()
			    << " events" << std::endl;
		return i;
	    }
	} else {
	    throw std::runtime_error("Invalid variable type for node " + d.getVariable(i).getName());
	}
    }
    return -1;
}

int StabilityUtils::checkForVariance(DataSet& d) {
    arma::mat t = d.getData();
    for (arma::uword i = 0; i < d.getNumColumns(); i++) {
        if (d.getVariable(i).isContinuous()) {
            arma::vec curr = standardizeData(t.col(i));

            double var = arma::var(curr);

            if (var <= 0.00001) {
		Rcpp::Rcout << "   " << d.getVariable(i).getName() << ":  " << var << std::endl;
                return i;
            }

        } else if (d.getVariable(i).isDiscrete()) {

	    if (d.getVariable(i).getNumCategories() < 2) {
		Rcpp::Rcout << "    " << d.getVariable(i).getName() << " has only one category\n";
		return i;
	    }
	  
            std::unordered_map<int, int> cats;
            for (arma::uword j = 0; j < d.getNumRows(); j++) {
                cats[d.getInt(j, i)] = 0;
            }

            for (arma::uword j = 0; j < d.getNumRows(); j++) {
                int currCat = d.getInt(j, i);
                if (cats.find(currCat) == cats.end()) {
                    throw std::invalid_argument("Found an invalid category");
                } else {
                    cats[currCat] += 1;
                }
            }

            for (std::pair<int, int> element : cats) {
                if (element.second < 5) {
		    Rcpp::Rcout << "    " << d.getVariable(i).getName() << ":  {";
		    for (std::pair<int, int> item : cats) {
			Rcpp::Rcout << item.first << " : " << item.second << ", ";
		    }
		    Rcpp::Rcout << "}\n";
                    return i;
                }
            }
        } else if (d.getVariable(i).isCensored()) {
	    if (d.getVariable(i).getNEvents() < 10) {
		Rcpp::Rcout << "   " << d.getVariable(i).getName() << ":  "
			    << d.getVariable(i).getNEvents()
			    << " events" << std::endl;
		return i;
	    }
	} else {
	    throw std::runtime_error("Invalid variable type for node " + d.getVariable(i).getName());
	}
    }
    return -1;
}

arma::vec StabilityUtils::standardizeData(const arma::vec& data) {
    arma::vec data2 = data;
    
    double mean = arma::mean(data);
    data2 -= mean;

    double norm = arma::stddev(data2);
    norm = (norm == 0) ? 1 : norm;
    data2 /= norm;

    return data2;
}

// arma::umat StabilityUtils::subsampData(DataSet& data, int N, int b) {
//     arma::umat samp = subSampleNoReplacement(data.getNumRows(), b, N);
//     int attempts = 5000;
//     bool done = false;
//     while(!done) {
//         // Rcpp::Rcout << "Attempt " << 5000 - attempts;
//         samp = subSampleNoReplacement(data.getNumRows(), b, N);
//         done = true;
//         for (int i = 0; i < samp.n_rows; i++) {
//             DataSet subset(data, samp.row(i));
//             if (checkForVariance(subset, data) != -1) {
//                 done = false;
//                 break;
//             }
//         }
//         attempts--;
//         if (attempts == 0) {
//             // Rcpp::Rcout << "ERROR: Unable to find a subsampled dataset of size " << b << " where there are at least one category of every discrete variable" << std::endl;
//             throw std::invalid_argument("Unable to find a subsampled dataset of size " + std::to_string(b) + " where there are at least two samples for each category of every discrete variable. The number of samples per subsampled dataset can be increased with the subSize parameter to address this problem.");
//         }
//     }
//     return samp;
// }

arma::mat StabilityUtils::stabilitySearchPar(DataSet& data, std::vector<double>& lambda, int num_threads, int N, int b) {
    arma::umat samp = subSampleNoReplacement(data, b, N);
    // int attempts = 5000;
    // bool done = false;
    // while(!done) {
    //     // Rcpp::Rcout << "Attempt " << 5000 - attempts;
    //     samp = subSampleNoReplacement(data.getNumRows(), b, N);
    //     done = true;
    //     for (int i = 0; i < samp.n_rows; i++) {
    //         DataSet subset(data, samp.row(i));
    //         if (checkForVariance(subset, data) != -1) {
    //             done = false;
    //             break;
    //         }
    //     }
    //     attempts--;
    //     if (attempts == 0) {
    //         // Rcpp::Rcout << "ERROR: Unable to find a subsampled dataset of size " << b << " where there are at least one category of every discrete variable" << std::endl;
    //         throw std::invalid_argument("Unable to find a subsampled dataset of size " + std::to_string(b) + " where there are at least two samples for each category of every discrete variable. The number of samples per subsampled dataset can be increased with the subSize parameter to address this problem.");
    //     }
    // }
    // Rcpp::Rcout << std::endl;

    return stabilitySearchPar(data, lambda, num_threads, samp);
}

// LOO
arma::mat StabilityUtils::stabilitySearchPar(DataSet& data, std::vector<double>& lambda, int num_threads) {
    arma::umat samps(data.getNumRows(), data.getNumRows()-1);

    for(int i = 0; i < samps.n_rows; i++) {
        int count = 0;
        for(int j = 0; j < samps.n_rows; j++) {
            if(i==j) goto A;
            samps(i, count) = j;
            count++;
            A:;
        }
    }
    
    // Print a warning if there's a variance issue
    for (int s = 0; s < samps.n_rows; s++) {
        DataSet dataSubSamp(data, samps.row(s));
        if (checkForVariance(dataSubSamp, data) != -1) {
            Rcpp::Rcout << "Variance issue with dataset: " << s << std::endl;
        }
    }
    
    return stabilitySearchPar(data, lambda, num_threads, samps);
}

arma::mat StabilityUtils::stabilitySearchPar(DataSet& data, std::vector<double>& lambda, int num_threads, arma::umat& subs) {
    BlockingQueue<int> taskQueue(subs.n_rows);
    if (num_threads <= 0) num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) {
        num_threads = 4;
        Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
    }

    int numVars = data.getNumColumns();
    arma::mat thetaMat(numVars, numVars, arma::fill::zeros);
    std::mutex matMutex; // For protecting thetaMat
    RcppThread::ThreadPool pool(num_threads);
    std::vector<MGM> mgmList;

    for (int i = 0; i < subs.n_rows; i++) {
    	DataSet dataSubSamp(data, subs.row(i));
    	mgmList.push_back(MGM(dataSubSamp));
    }


    pool.parallelForEach(mgmList,
			 [&] (MGM& m) {
			     m.setLambda(lambda);
			     EdgeListGraph g = m.search();
			     arma::mat curAdj = skeletonToMatrix(g, data);
				 
			     {
				 std::lock_guard<std::mutex> matLock(matMutex);
				 thetaMat += curAdj;
			     }
			 });

    pool.join();

    
    // auto producer = [&]() {
    //     for (int i = 0; i < subs.n_rows; i++) {
    //         taskQueue.push(i);

    // 	    if (RcppThread::isInterrupted()) {
    // 	      break;
    // 	    }
    //     }

    //     // Poison pills
    //     for (int i = 0; i < num_threads; i++) {
    //         taskQueue.push(-1);
    //     }
    // };

    // auto consumer = [&]() {
    //     while (true) {
    //         int sample = taskQueue.pop();

    //         if (sample < 0) break; // Poison pill

    // 	    if (RcppThread::isInterrupted()) {
    // 	      break;
    // 	    }

    //         DataSet dataSubSamp(data, subs.row(sample));
    //         MGM mgm(dataSubSamp, lambda);
    //         EdgeListGraph g = mgm.search();
    //         arma::mat curAdj = skeletonToMatrix(g, dataSubSamp);

    // 	    {
    //         std::lock_guard<std::mutex> matLock(matMutex);
    //         thetaMat += curAdj;
    // 	    }
    //     }
    // };


    // RcppThread::ThreadPool pool(num_threads);

    // pool.push(producer);

    // for (int i = 0; i < num_threads; i++) {
    //     pool.push(consumer);
    // }

    // pool.join();
    
    // std::vector<RcppThread::Thread> threads;

    // threads.push_back(RcppThread::Thread( producer ));

    // for (int i = 0; i < num_threads; i++) {
    //     threads.push_back(RcppThread::Thread( consumer ));
    // }

    // for (int i = 0; i < threads.size(); i++) {
    //     threads[i].join();
    // }

    return thetaMat / subs.n_rows;
}

double StabilityUtils::stabilitySearchStars(DataSet& data,
					    std::string& alg,
					    double param,
					    EdgeListGraph* initialGraph,
					    int num_threads,
					    bool adjacency,
					    int N,
					    int b) {
    arma::umat samp = subSampleNoReplacement(data, b, N);
    // if (alg != "mgm") {
    // 	b = std::min(4*data.getNumRows()/5, (int) std::floor(20*std::sqrt(data.getNumRows())));
    // 	Rcpp::Rcout << "Subsample size = " << b << std::endl;
    // }
    // int attempts = 5000;
    // bool done = false;
    // while(!done) {
    //     // Rcpp::Rcout << "Attempt " << 5000 - attempts;
    //     samp = subSampleNoReplacement(data.getNumRows(), b, N);
    //     done = true;
    //     for (int i = 0; i < samp.n_rows; i++) {
    //         DataSet subset(data, samp.row(i));
    //         if (checkForVariance(subset, data) != -1) {
    //             done = false;
    //             break;
    //         }
    //     }
    //     attempts--;
    //     if (attempts == 0) {
    //         Rcpp::Rcout << "ERROR: Unable to find a subsampled dataset of size " << b << " where there are at least one category of every discrete variable" << std::endl;
    //         throw std::invalid_argument("Unable to find a subsampled dataset of size " + std::to_string(b) + " where there are at least one category of every discrete variable");
    //     }
    // }
    // Rcpp::Rcout << std::endl;

    return stabilitySearchStars(data, alg, param, initialGraph, num_threads, adjacency, samp);
}

double StabilityUtils::stabilitySearchStars(DataSet& data,
					    std::string& alg,
					    double param,
					    EdgeListGraph* initialGraph,
					    int num_threads,
					    bool adjacency) {
    // LOO
    arma::umat samps(data.getNumRows(), data.getNumRows()-1);

    for(int i = 0; i < samps.n_rows; i++) {
        int count = 0;
        for(int j = 0; j < samps.n_rows; j++) {
            if(i==j) goto A;
            samps(i, count) = j;
            count++;
            A:;
        }
    }
    
    // Print a warning if there's a variance issue
    for (int s = 0; s < samps.n_rows; s++) {
        DataSet dataSubSamp(data, samps.row(s));
        if (checkForVariance(dataSubSamp, data) != -1) {
            Rcpp::Rcout << "Variance issue with dataset: " << s << std::endl;
        }
    }
    
    return stabilitySearchStars(data, alg, param, initialGraph, num_threads, adjacency, samps);
}

double StabilityUtils::stabilitySearchStars(DataSet& data,
					    std::string& alg,
					    double param,
					    EdgeListGraph* initialGraph,
					    int num_threads,
					    bool adjacency,
					    arma::umat& subs) {
    int numVars = data.getNumColumns();
    int layers;
    double numPossEdges=0;
    std::vector<double> lambda;
    double alpha;

    if (alg == "mgm") {
	adjacency = true;
	numPossEdges = 0.5 * numVars * (numVars - 1);

	int lamLength;

	if (data.isMixed()) {
	    lamLength = 3;
	    // if (data.isCensored()) {
	    // 	lamLength = 5;
	    // } else {
	    // 	lamLength = 3;
	    // }
	} else {
	    throw std::runtime_error("MGM is not implemented for purely continuous or purely discrete datasets.");
	}

	for (int i = 1; i < lamLength; i++) {
	    lambda.push_back(param);
	}
	
    } else if (adjacency) {
	if (initialGraph == NULL) {
	    numPossEdges = numVars * (numVars - 1) / 2;
	} else {
	    std::string ig_alg = initialGraph->getAlgorithm();
	    if (ig_alg.find("MGM")!=std::string::npos) {
		std::size_t start = ig_alg.find("[");
		std::size_t end = ig_alg.find("]");
		if (start != end && end != std::string::npos) {
		    std::stringstream ss(ig_alg.substr(start+1, end-start-1));
		    for (double lam; ss >> lam;) {
			lambda.push_back(lam);
			while (ss.peek()==',' || ss.peek()==' ')
			    ss.ignore();
		    }
		    
		    Rcpp::Rcout << lambda[0] << ", " << lambda[1] << ", " << lambda[2] << std::endl;
		    numPossEdges = numVars * (numVars - 1) / 2;
		}
	    }
	    if (numPossEdges == 0) {
		// std::unordered_set<Edge> edgeSet = initialGraph->getEdges();
		numPossEdges = initialGraph->getEdges().size();
	    }
	}
	alpha = param;
    } else if (alg == "pc" || alg == "cpc" || alg == "pcm") {
	layers = 2;
	if (initialGraph == NULL) {
	    numPossEdges = 1.5 * numVars * (numVars - 1);
	} else {
	    std::string ig_alg = initialGraph->getAlgorithm();
	    if (ig_alg.find("MGM")!=std::string::npos) {
		std::size_t start = ig_alg.find("[");
		std::size_t end = ig_alg.find("]");
		if (start != end && end != std::string::npos) {
		    std::stringstream ss(ig_alg.substr(start+1, end-start-1));
		    for (double lam; ss >> lam;) {
			lambda.push_back(lam);
			while (ss.peek()==',' || ss.peek()==' ')
			    ss.ignore();
		    }
		    
		    Rcpp::Rcout << lambda[0] << ", " << lambda[1] << ", " << lambda[2] << std::endl;
		    numPossEdges = 1.5 * numVars * (numVars - 1);
		}
	    }
	    if (numPossEdges == 0) {
		// std::unordered_set<Edge> edgeSet = initialGraph->getEdges();
		numPossEdges = 3.0 * initialGraph->getEdges().size();
	    }
	}
	alpha = param;
    } else if (alg == "fci" || alg == "cfci" || alg == "fcim") {
	layers = 3;
	if (initialGraph == NULL) {
	    numPossEdges = 3.0 * numVars * (numVars - 1);
	} else {
	    std::string ig_alg = initialGraph->getAlgorithm();
	    if (ig_alg.find("MGM")!=std::string::npos) {
		std::size_t start = ig_alg.find("[");
		std::size_t end = ig_alg.find("]");
		if (start != end && end != std::string::npos) {
		    std::stringstream ss(ig_alg.substr(start+1, end-start-1));
		    for (double lam; ss >> lam;) {
			lambda.push_back(lam);
			while (ss.peek()==',' || ss.peek()==' ')
			    ss.ignore();
		    }
		    
		    Rcpp::Rcout << lambda[0] << ", " << lambda[1] << ", " << lambda[2] << std::endl;
		    numPossEdges = 3.0 * numVars * (numVars - 1);
		}
	    }
	    if (numPossEdges == 0) {
		// std::unordered_set<Edge> edgeSet = initialGraph->getEdges();
		numPossEdges = 6.0 * initialGraph->getEdges().size();
	    }
	}
	alpha = param;
    } else {
	throw std::invalid_argument("Invalid algorithm: " + alg
				    + "\n   Algorithm must be in the list: { mgm, pc, cpc, pcm, fci, cfci, fcim }");
    }

    EdgeListGraph g;
    EdgeListGraph ig;
    arma::mat thetaMat;
    arma::cube thetaCube;

    if (adjacency)
	thetaMat = arma::mat(numVars, numVars, arma::fill::zeros);
    else
	thetaCube = arma::cube(numVars, numVars, layers, arma::fill::zeros);

    // for (int i = 0; i < subs.n_rows; i++) {
    // 	DataSet dataSubSamp(data, subs.row(i));
	
    // 	if (alg == "mgm") {
	    
    // 	    MGM mgm(dataSubSamp, lambda);
    // 	    g = mgm.search();
	    
    // 	} else if (alg == "pc") {
	    
    // 	    IndTestMulti itm(dataSubSamp, alpha);
	    
    // 	    PcStable pcs((IndependenceTest*) &itm);
    // 	    if (num_threads > 0) pcs.setThreads(num_threads);
	    
    // 	    pcs.setVerbose(false);

    // 	    if (!lambda.empty()) {
    // 		MGM mgm(dataSubSamp, lambda);
    // 		ig = mgm.search();
    // 		pcs.setInitialGraph(&ig);
    // 	    } else {
    // 		pcs.setInitialGraph(initialGraph);
    // 	    }
	    
    // 	    g = pcs.search();
	    
    // 	} else if (alg == "cpc") {
	    
    // 	    IndTestMulti itm(dataSubSamp, alpha);
	    
    // 	    CpcStable cpc((IndependenceTest*) &itm);
    // 	    if (num_threads > 0) cpc.setThreads(num_threads);
	    
    // 	    cpc.setVerbose(false);

    // 	    if (!lambda.empty()) {
    // 		MGM mgm(dataSubSamp, lambda);
    // 		ig = mgm.search();
    // 		cpc.setInitialGraph(&ig);
    // 	    } else {
    // 		cpc.setInitialGraph(initialGraph);
    // 	    }
	    
    // 	    g = cpc.search();
	    
    // 	} else if (alg == "pcm") {
	    
    // 	    IndTestMulti itm(dataSubSamp, alpha);
	    
    // 	    PcMax pcm((IndependenceTest*) &itm);
    // 	    if (num_threads > 0) pcm.setThreads(num_threads);
	    
    // 	    pcm.setVerbose(false);
	    
    // 	    if (!lambda.empty()) {
    // 		MGM mgm(dataSubSamp, lambda);
    // 		ig = mgm.search();
    // 		pcm.setInitialGraph(&ig);
    // 	    } else {
    // 		pcm.setInitialGraph(initialGraph);
    // 	    }
	    
    // 	    g = pcm.search();
	    
    // 	} else if (alg == "fci") {
	    
    // 	    IndTestMulti itm(dataSubSamp, alpha);
	    
    // 	    Fci fci((IndependenceTest*) &itm);
    // 	    if (num_threads > 0) fci.setThreads(num_threads);
	    
    // 	    fci.setVerbose(false);
	    
    // 	    if (!lambda.empty()) {
    // 		MGM mgm(dataSubSamp, lambda);
    // 		ig = mgm.search();
    // 		fci.setInitialGraph(&ig);
    // 	    } else {
    // 		fci.setInitialGraph(initialGraph);
    // 	    }
	    
    // 	    g = fci.search();
	    
    // 	} else if (alg == "cfci") {
	    
    // 	    IndTestMulti itm(dataSubSamp, alpha);
	    
    // 	    Cfci cfci((IndependenceTest*) &itm);
    // 	    if (num_threads > 0) cfci.setThreads(num_threads);
	    
    // 	    cfci.setVerbose(false);
	    
    // 	    if (!lambda.empty()) {
    // 		MGM mgm(dataSubSamp, lambda);
    // 		ig = mgm.search();
    // 		cfci.setInitialGraph(&ig);
    // 	    } else {
    // 		cfci.setInitialGraph(initialGraph);
    // 	    }
	    
    // 	    g = cfci.search();
	    
    // 	} else if (alg == "fcim") {
	    
    // 	    IndTestMulti itm(dataSubSamp, alpha);
	    
    // 	    FciMax fcim((IndependenceTest*) &itm);
    // 	    if (num_threads > 0) fcim.setThreads(num_threads);
	    
    // 	    fcim.setVerbose(false);

    // 	    if (!lambda.empty()) {
    // 		MGM mgm(dataSubSamp, lambda);
    // 		ig = mgm.search();
    // 		fcim.setInitialGraph(&ig);
    // 	    } else {
    // 		fcim.setInitialGraph(initialGraph);
    // 	    }
	    
    // 	    g = fcim.search();
	    
    // 	} else {
    // 	    throw std::invalid_argument("Unrecognized search algorithm");
    // 	}

    // 	// arma::cube curAdj = graphToCube(g, dataSubSamp);
    // 	if (adjacency)
    // 	    thetaMat += skeletonToMatrix(g, dataSubSamp);
    // 	else
    // 	    thetaCube += graphToCube(g, dataSubSamp);
    // }
    
    // if (adjacency) {
    // 	thetaMat = thetaMat / subs.n_rows;
    // 	Rcpp::Rcout << thetaMat << std::endl;
    // 	arma::mat destable = 2 * thetaMat % (1-thetaMat);
    // 	Rcpp::Rcout << destable << std::endl;
    // 	Rcpp::Rcout << arma::accu(destable) / 2 << " / " << numPossEdges << std::endl;
    // 	return (arma::accu(destable) / 2) / numPossEdges;
    // }
    
    thetaCube = thetaCube / subs.n_rows;
    arma::cube destable = 2 * thetaCube % (1-thetaCube);
    return arma::accu(destable) / numPossEdges;
}
