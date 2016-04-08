#include <fstream>
#include <istream>
#include <string>
#include <vector>
#include <numeric>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <limits> // std::numeric_limits<double>::max();
#include <algorithm> // min()
#include <errno.h>
#include <cstring> // strerror()
#include <sys/time.h>
#include <ios>

std::ostream& operator<<(std::ostream& os, std::vector<double>& d) {
    os << "[";
    for (int i; i < (int)d.size(); i++) {
        os << d[i] << " ";
    }
    os << "]" << std::endl;
    return os;
}


// FOR DTW TAKE BIN SIZE AND RUN PAA ON QUERY FILE
void usage() {
    std::cout << "use this program in one of the following ways" << std::endl;
    std::cout << "./PDTW -PAA file blockSize" << std::endl;
    std::cout << "./PDTW -DTW dataFile queryFile" << std::endl;
    std::cout << "./PDTW -OUTLIERS dataFile length" << std::endl;
}

struct DTWData {
    double distance;
    int timeSeries;
    int locStart;
    int blkSize;
    DTWData(double in_distance, int in_timeSeries, int in_locStart, int in_blkSize) {
        distance = in_distance;
        timeSeries = in_timeSeries;
        locStart = in_locStart;
        blkSize = in_blkSize;
    }

    DTWData() {
        distance = 0;
        timeSeries = 0;
        locStart = 0;
        blkSize = 0;
    }

    bool empty() {
        return (distance == 0) && (timeSeries == 0) && (locStart == 0) && (blkSize == 0);
    }
};

std::ostream& operator<<(std::ostream& os, DTWData& d) {
    os << "Time Series: " << d.timeSeries << std::endl;
    os << "Location: " << d.locStart << " to " << d.locStart+d.blkSize << std::endl;
    os << "Distance: " << d.distance << std::endl;
    return os;
}

// takes the given file and adds '_PAA' to the end
// returns the new filepath
std::string createPAAFilepath(std::string in) {
    std::size_t point = in.find(".");
    if (point == std::string::npos)
        return in + "_PAA";
    return in.substr(0, point) + "_PAA" + in.substr(point);
}

std::string getPAAFilename(std::string incoming) {
    if (incoming.find("_PAA") != std::string::npos) { // found "_PAA"
        return incoming;
    }
    return createPAAFilepath(incoming);
}

std::vector<double> timeSeriesToVector(std::string series) {
    std::vector<double> vec;
    std::size_t prev = 0, pos;
    while ((pos = series.find_first_of(" ,", prev)) != std::string::npos) {
        if (pos > prev)
            vec.push_back(stod(series.substr(prev, pos-prev)));
        prev = pos+1;
    }
    if (prev < series.length())
        vec.push_back(stod(series.substr(prev, std::string::npos)));

    return vec;
}

//simple Euclidean distance between two points
double distFunc(double x, double y) {
    return sqrt(pow((x - y), 2));
}

//DTW calculation between two vectors. one is query and other Data sequence
double simpleDTW(const std::vector<double>& t1, const std::vector<double>& t2) {
//    std::cout << "t1" << std::endl;
//    for (int k=0; k< t1.size(); k++) {
//        std::cout << t1[k] << " ";
//    }
//    std::cout << std::endl << "t2" << std::endl;
//    for (int l=0; l< t2.size(); l++ ) {
//        std::cout << t2[l] << " ";
//    }
//    std::cout << std::endl << std::endl;

    int m = t1.size();
    int n = t2.size();
    // create cost matrix
    double cost[m][n];
    cost[0][0] = distFunc(t1[0], t2[0]);

    // calculate first row
    for(int i = 1; i < m; i++) {
        cost[i][0] = cost[i-1][0] + distFunc(t1[i], t2[0]);
    }
    // calculate first column
    for(int j = 1; j < n; j++) {
        cost[0][j] = cost[0][j-1] + distFunc(t1[0], t2[j]);
    }
    // fill matrix
    for(int i = 1; i < m; i++) {
        for(int j = 1; j < n; j++) {
            cost[i][j] = distFunc(t1[i],t2[j])+ std::min(cost[i-1][j],std::min(cost[i][j-1], cost[i-1][j-1]));
        }
    }

//    std::cout << "dist: " << cost[m-1][n-1] << std::endl;
    return cost[m-1][n-1];
}

// both incoming files should be PAA'd already
DTWData DTWaFile(std::string dataFile, std::string queryFile) {
    std::ifstream data(getPAAFilename(dataFile).c_str()); // bigger file
    std::ifstream query(getPAAFilename(queryFile).c_str()); // smaller file

    if (!data) {
            std::cout << "ERROR: ifstream failed on " << dataFile << ": " << strerror(errno) << std::endl;
            return DTWData();
    }
    if (!query) {
            std::cout << "ERROR: ifstream failed on " << queryFile << ": " << strerror(errno) << std::endl;
            return DTWData();
    }

    std::string dataTimePoints;
    std::string queryData;

    std::vector<double> dataVector;
    std::vector<double> queryVector;
    int curTimeSeries = -1;
    int bestMatchTimeSeries = 1;
    int bestMatchIdx = 1;
    int bestMatchBlkSz = 2;
    double bestMatchDistance = std::numeric_limits<double>::max();

    // turn query into vector
    std::getline(query, queryData);
    std::size_t prev = 0, pos;
    while ((pos = queryData.find_first_of(" ,", prev)) != std::string::npos) {
        if (pos > prev)
            queryVector.push_back(stod(queryData.substr(prev, pos-prev)));
        prev = pos+1;
    }
    if (prev < queryData.length())
        queryVector.push_back(stod(queryData.substr(prev, std::string::npos)));

    while (dataTimePoints.empty() && data.good())
        std::getline(data, dataTimePoints);

    while (data.good()) {
        // get the next data timeseries
        curTimeSeries++;
        // split the timeseries numbers on space or comma
        std::size_t prev = 0, pos;
        while ((pos = dataTimePoints.find_first_of(" ,", prev)) != std::string::npos) {
            if (pos > prev)
                dataVector.push_back(stod(dataTimePoints.substr(prev, pos-prev)));
            prev = pos+1;
        }
        if (prev < dataTimePoints.length())
            dataVector.push_back(stod(dataTimePoints.substr(prev, std::string::npos)));

        // run through all combinations from query
        for (int blkSz = 2; blkSz <= (int)dataVector.size(); blkSz++) {
            for (int startIdx = 0; startIdx+blkSz <= (int)dataVector.size(); startIdx++) {
                std::vector<double> subVec(dataVector.begin()+startIdx, dataVector.begin()+startIdx+blkSz);
                double newBest = std::min(simpleDTW(queryVector, subVec), bestMatchDistance);
		if (newBest != bestMatchDistance) {
                    bestMatchDistance = newBest;
                    bestMatchIdx = startIdx;
                    bestMatchBlkSz = blkSz;
                    bestMatchTimeSeries = curTimeSeries;
		}
            }
        }

        // get the next data timeseries
        dataVector.clear();
        dataTimePoints = "";
        while (dataTimePoints.empty() && data.good())
            std::getline(data, dataTimePoints);
    }
    return DTWData(bestMatchDistance, bestMatchTimeSeries, bestMatchIdx, bestMatchBlkSz);
}

// just take length, compute all sections of a timeseries with that length
// check against all of the sections of all other timeseries.
// sum the max.

// the incoming file should be PAA'd already
DTWData detectOutliers(std::string dataFile, int length) {
    std::ifstream dataStream(getPAAFilename(dataFile).c_str());

    if (!dataStream) {
        std::cout << "ERROR: ifstream failed on " << dataFile << ": " << strerror(errno) << std::endl;
        return DTWData();
    }

    std::string dataTimePoints;
    std::string queryData;
    std::vector<double> dataVector;
    std::vector<double> queryVector;
    int dataPos;
    int curQTimeSeries = -1;
    int curDTimeSeries = 0;
    std::vector<double> diffs;
    int queryStart = 0, queryEnd = length;
    int dataStart = 0, dataEnd = length;

    // get first time series
    while (queryData.empty())
        std::getline(dataStream, queryData);

    // run through each time series as query
    while (dataStream.good()) {
        curQTimeSeries++;
        if (curQTimeSeries == 0) // first one
            diffs.push_back(0.0); // for each time series, add an entry
        // turn time series into vector
        queryVector = timeSeriesToVector(queryData);

        // grab location in file
        dataPos = dataStream.tellg();

        // grab next data time series
        while (dataTimePoints.empty() && dataStream.good())
            std::getline(dataStream, dataTimePoints);

        // reset the data time series count
        curDTimeSeries = curQTimeSeries;

        // run through the rest of the time series
        while (dataStream.good()) {
            double worstDist = 0;
            int worstDataStart = 0;
            int worstQueryStart = 0;
            // reset the query start and end
            queryStart = 0;
            queryEnd = length;
            curDTimeSeries++;

            if (curQTimeSeries == 0) // first one
                diffs.push_back(0.0); // for each time series, add an entry
            // turn time series into vector
            dataVector = timeSeriesToVector(dataTimePoints);

            // run through all chunks of the query
            while (queryEnd <= (int) queryVector.size()) {
                // reset the data start and end for the next query chunk
                dataStart = 0;
                dataEnd = length;

                // run through all chunks of the data
                while (dataEnd <= (int) dataVector.size()) {
                    std::vector<double> subQVec(queryVector.begin()+queryStart, queryVector.begin()+queryEnd);
                    std::vector<double> subDVec(dataVector.begin()+dataStart, dataVector.begin()+dataEnd);
                    double oldWorst = worstDist;
                    worstDist = std::max(simpleDTW(subQVec, subDVec), worstDist);
                    if (oldWorst != worstDist) {
                        worstDataStart = dataStart;
                        worstQueryStart = queryStart;
                    }

                    dataStart++;
                    dataEnd++;
                }
                queryStart++;
                queryEnd++;
            }

            diffs[curQTimeSeries] += worstDist; // add new worst to
            diffs[curDTimeSeries] += worstDist;
            dataVector.clear();
            dataTimePoints = "";

            // get next time series
            while (dataTimePoints.empty() && dataStream.good())
                std::getline(dataStream, dataTimePoints);
        }

        dataStream.clear();
        dataStream.seekg(dataPos, dataStream.beg);

        queryData = "";
        dataTimePoints = "";
        dataVector.clear();
        queryVector.clear();
        // turn "next" data time series into query vector
        while (queryData.empty() && dataStream.good())
            std::getline(dataStream, queryData);
    }

    double worst =0;
    int loc = 0;
    for (int k=0; k < (int)diffs.size(); k++) {
        worst = std::max(diffs[k], worst);
        loc = k;
    }

    std::cout << "size: " << diffs.size() << std::endl;

    return DTWData(worst, loc, 0, 0); // FIX THIS
}

// how big will the incoming numbers be?
double PAA(std::vector<double> inData) {
    double ave = 0;
    for (int i=0; i< (int)inData.size(); i++) {
            ave += inData[i];
    }
    ave /= (int)inData.size();
    return ave;
}

// runs through the given file and calculates PAA with a N points in each block
std::string PAAaFile(std::string inData, int N) {
    if (N <= 0) {
        std::cout << "ERROR: block size is " << N << ": should be greater than 0" << std::endl;
        return "";
    }
    std::ifstream inFile(inData.c_str());
    std::string PAAFile = createPAAFilepath(inData);
    std::ofstream outFile(PAAFile.c_str());
    std::string timePoint;
    std::size_t newLine = std::string::npos;
    std::vector<double> data;
    char delim = ',';

    if (!inFile) {
        std::cout << "ERROR: ifstream failed on " << inData << ": " << strerror(errno) << std::endl;
        return "";
    }

    int len = inFile.tellg();
    std::getline(inFile, timePoint);
    if (timePoint.find(",") == std::string::npos) { // no commas
        delim = ' ';
    }
    inFile.seekg(len, std::ios_base::beg);

    while (inFile.good()) {
        // run until N time points are grabbed OR
        // the end of the time series is hit
        for (int i=0; (i < N) && (newLine == std::string::npos); i++ ) {
            // get one time value
            std::getline(inFile, timePoint, delim);
            // check if we hit the end of the time series
            // ex) "number1\nnumber2"
            newLine = timePoint.find("\n");
            data.push_back(stod(timePoint));
        }

        // send one chunk of data to PAA and write to outFile
        outFile << PAA(data) << " ";
        // clear data to grab the next chunk
        data.clear();

        // if you hit the end of a time series "number1\nnumber2"
        if (newLine != std::string::npos) {
            // timer series break
            outFile << std::endl << std::endl;
            // add the second number to data
            if (timePoint.substr(newLine).size() > 1) {
                    data.push_back(stod(timePoint.substr(newLine)));
            }
            // reset newLine
            newLine = std::string::npos;
        }
    }
    return PAAFile;
}


int main(int argc, char** argv) {
    if (argc != 4) {
        usage();
        return -1;
    }

    clock_t start, stop;

    if (!std::string(argv[1]).compare("-PAA")) {
        start = clock();
        std::string PAAFile = PAAaFile(argv[2], atoi(argv[3]));
        stop = clock();
        if (PAAFile.empty()) {
            std::cout << "PAAaFile failed" << std::endl;
            return -1;
        }
        std::cout << "PAA computed for " << argv[2] << ". Written to " << PAAFile << std::endl;
        int s = (stop - start) / CLOCKS_PER_SEC;
        long t = ((stop - start) % CLOCKS_PER_SEC)/1000;
        std::cout << "PAA took: " << s << " seconds and " << t << " milliseconds to compute" << std::endl;
    } else if (!std::string(argv[1]).compare("-DTW")) {
        start = clock();
        DTWData d = DTWaFile(argv[2], argv[3]);
        stop = clock();
        if (d.empty()) {
            std::cout << "DTWaFile failed" << std::endl;
            return -1;
        }
        std::cout << "DTW computed for data " << argv[2] << " and query " << argv[3];
        std::cout << ". Best distance found was: " << std::endl << d << std::endl;
        int s = (stop - start)/ CLOCKS_PER_SEC;
        long t = ((stop - start) % CLOCKS_PER_SEC)/1000;
        std::cout << "DTW took: " << s << " seconds and " << t << " milliseconds to compute" << std::endl;
    } else if (!std::string(argv[1]).compare("-OUTLIERS")) {
        start = clock();
        DTWData d = detectOutliers(argv[2], atoi(argv[3]));
        stop = clock();
        if (d.empty()) {
            std::cout << "detectOutliers failed" << std::endl;
            return -1;
        }
        std::cout << "Outlier detection computed for data " << argv[2];
        std::cout << ". Worst distance found was: " << std::endl << d << std::endl;
        int s = (stop - start)/ CLOCKS_PER_SEC;
        long t = ((stop - start) % CLOCKS_PER_SEC)/1000;
        std::cout << "DTW took: " << s << " seconds and " << t << " milliseconds to compute" << std::endl;
    } else {
        usage();
        return -1;
    }

    return 1;
}
