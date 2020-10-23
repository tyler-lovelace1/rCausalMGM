#include "StabilityUtils.hpp"

arma::umat StabilityUtils::subSampleNoReplacement(int sampSize, int subSize, int numSub) {

    if (subSize < 1)
        throw std::invalid_argument("Sample size must be > 0");

    arma::urowvec indices(sampSize);
    for (arma::uword i = 0; i < sampSize; i++) {
        indices(i) = i;
    }

    arma::umat sampMat(numSub, subSize);

    for(arma::uword i = 0; i < numSub; i++) {
        indices = arma::shuffle(indices);
        arma::urowvec curSamp;
        while(true) {
            SAMP:
            curSamp = subSampleIndices(sampSize, subSize);
            for (arma::uword j = 0; j < i; j++) {
                if (arma::all(curSamp == sampMat.row(j))) {
                    goto SAMP;
                }
            }
            break;
        }
        sampMat.row(i) = curSamp;
    }
    
    return sampMat;
}

arma::urowvec StabilityUtils::subSampleIndices(int N, int subSize) {

    arma::urowvec indices(N);
    for (arma::uword i = 0; i < N; i++) {
        indices(i) = i;
    }

    indices = arma::shuffle(indices);
    arma::urowvec samp(subSize);
    for (arma::uword i = 0; i < subSize; i++) {
        samp(i) = indices(i);
    }

    return samp;
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
    std::unordered_map<Variable*, int> map;
    for (Variable* node : graph.getNodes()) {
        map[node] = d.getColumn(node);
    }

    // mark edges
    for (Edge edge : graph.getEdges()) {
        // if directed find which is parent/child
        Variable* node1 = edge.getNode1();
        Variable* node2 = edge.getNode2();

        matrix(map[node1], map[node2]) = 1.0;
        matrix(map[node2], map[node1]) = 1.0;
    }

    return matrix;
}

int StabilityUtils::checkForVariance(DataSet& d, DataSet& full) {
    arma::mat t = d.getData();
    for (arma::uword i = 0; i < d.getNumColumns(); i++) {
        if (d.getVariable(i)->isContinuous()) {
            arma::vec curr = standardizeData(t.col(i));

            double var = arma::var(curr);

            if (var <= 0.0001) {
                return i;
            }

        } else {
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
                if (element.second < 2) {
                    return i;
                }
            }
        }
    }
    return -1;
}

arma::vec StabilityUtils::standardizeData(const arma::vec& data) {
    arma::vec data2 = data;
    
    double mean = arma::mean(data);
    data2 -= mean;

    double norm = arma::stddev(data2);
    data2 /= norm;

    return data2;
}

arma::mat StabilityUtils::stabilitySearchPar(DataSet& data, std::vector<double>& lambda, int N, int b) {
    
    arma::umat samp = subSampleNoReplacement(data.getNumRows(), b, N);
    int attempts = 5000;
    bool done = false;
    while(!done) {
        //TODO the ordering here is weird
        done = true;
        for (int i = 0; i < samp.n_rows; i++) {
            DataSet subset(data, samp.row(i));
            if (checkForVariance(subset, data) != -1) done = false;
        }
        samp = subSampleNoReplacement(data.getNumRows(), b, N);
        attempts--;
        if (attempts == 0) {
            Rcpp::Rcout << "Unable to find a subsampled dataset of size " << b << " where there are at least one category of every discrete variable" << std::endl;
            exit(-1);
        }
    }

    return stabilitySearchPar(data, lambda, samp);
}

// LOO
arma::mat StabilityUtils::stabilitySearchPar(DataSet& data, std::vector<double>& lambda) {
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
    
    return stabilitySearchPar(data, lambda, samps);
}

arma::mat StabilityUtils::stabilitySearchPar(DataSet& data, std::vector<double>& lambda, arma::umat& subs) {
    BlockingQueue<int> taskQueue(10000);
    int parallelism = std::thread::hardware_concurrency();
    if (parallelism == 0) {
        parallelism = 4;
        Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
    }

    int numVars = data.getNumColumns();
    arma::mat thetaMat(numVars, numVars, arma::fill::zeros);
    std::mutex matMutex; // For protecting thetaMat

    auto producer = [&]() {
        for (int i = 0; i < subs.n_rows; i++) {
            taskQueue.push(i);
        }

        // Poison pills
        for (int i = 0; i < parallelism; i++) {
            taskQueue.push(-1);
        }
    };

    auto consumer = [&]() {
        while (true) {
            int sample = taskQueue.pop();

            if (sample < 0) break; // Poison pill

            DataSet dataSubSamp(data, subs.row(sample));
            MGM mgm(dataSubSamp, lambda);
            EdgeListGraph g = mgm.search();
            arma::mat curAdj = skeletonToMatrix(g, dataSubSamp);
            
            std::lock_guard<std::mutex> matLock(matMutex);
            thetaMat += curAdj;
        }
    };

    std::vector<std::thread> threads;

    threads.push_back(std::thread( producer ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(std::thread( consumer ));
    }

    for (int i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

    return thetaMat / subs.n_rows;
}
