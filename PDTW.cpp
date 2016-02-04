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

// Take big set and multiply it by 4-5 to make huge dataset

// DTW doesn't need bin number, instead take every possible subsequence for each timeseries

void usage() {
	std::cout << "./PDTW OPTION blockSize" << std::endl;
	std::cout << "OPTIONS" << std::endl;
	std::cout << "\t-PAA file" << std::endl;
	std::cout << "\t-DTW dataFile queryFile" << std::endl;
}

struct DTWData {
	double distance;
	int timeSeries;
	int location;
	DTWData(double in_distance, int in_timeSeries, int in_location) {
		distance = in_distance;
		timeSeries = in_timeSeries;
		location = in_location;
	}

	DTWData() {
		distance = 0;
		timeSeries = 0;
		location = 0;
	}

	bool empty() {
		return (distance == 0) && (timeSeries == 0) && (location == 0);
	}
};

std::ostream& operator<<(std::ostream& os, DTWData& d) {
		os << "Time Series: " << d.timeSeries << std::endl;
		os << "Location: " << d.location << std::endl;
		os << "Distance: " << d.distance << std::endl;
		return os;
}

//simple Euclidean distance between two points
double distFunc(double x, double y) {
    return sqrt(pow((x - y), 2));
}

//DTW calculation between two vectors. one is query and other Data sequence
double simpleDTW(const std::vector<double>& t1, const std::vector<double>& t2) {
	std::cout << "t1" << std::endl;
	for (int k=0; k < t1.size(); k++) {
		std::cout << t1[k] << " ";
	}
	std::cout << std::endl << "t2" << std::endl;
	for (int m=0; m < t2.size(); m++) {
		std::cout << t2[m] << " ";
	}
	std::cout << std::endl << std::endl;

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
    return cost[m-1][n-1];
}

// both incoming files should be PAA'd already
DTWData DTWaFile(std::string dataFile, std::string queryFile, int N) {
	if (N <= 0) {
		std::cout << "ERROR: block size is " << N << ": should be greater than 0" << std::endl;
	}

	std::ifstream data(dataFile.c_str());
	std::ifstream query(queryFile.c_str());

	if (!data) {
		std::cout << "ERROR: ifstream failed on " << dataFile << ": " << strerror(errno) << std::endl;
		return DTWData();
	}
	if (!query) {
		std::cout << "ERROR: ifstream failed on " << queryFile << ": " << strerror(errno) << std::endl;
		return DTWData();
	}

	std::string dataTimePoints;
	std::string queryTimePoint;

	std::vector<double> dataVector;
	std::vector<double> queryVector;
	int curTimeSeries = 0;
	int curLocation = 1;

	int bestMatchTimeSeries = 1;
	int bestMatchLocation = 1;
	double bestMatchDistance = std::numeric_limits<double>::max();

	// DTW takes the whole data series and moves the query
	// DTW takes 1 timeseries from the data file and moves the query of length N up one every time.
	// TODO: grab N values from data and N values from query
	while (data.good() && query.good()) {
		dataVector.clear();
		std::size_t queryNewLine = std::string::npos;
		curTimeSeries++;
		curLocation = 1;
		dataTimePoints = "";
		// ------------------------------------------ FULL DATA TIMESERIES
		// get the next timeseries
		while (dataTimePoints.empty())
			std::getline(data, dataTimePoints);

		// split the string into numbers on space or comma
		std::size_t prev = 0, pos;
	    while ((pos = dataTimePoints.find_first_of(" ,", prev)) != std::string::npos) {
	        if (pos > prev)
	            dataVector.push_back(stod(dataTimePoints.substr(prev, pos-prev)));
	        prev = pos+1;
	    }
	    if (prev < dataTimePoints.length())
	        dataVector.push_back(stod(dataTimePoints.substr(prev, std::string::npos)));

	    // ------------------------------------------- INITIAL N FROM QUERY
		// get initial N from query
		for (int j=(int)queryVector.size(); j < N; j++ ) {
			// get one time value
			std::getline(query, queryTimePoint, ' ');
			if (queryTimePoint.find("\n") != std::string::npos) {	// hit the end of the time series
				std::cout << "Hit end of timeseries. Invalid block size" << std::endl;
				return DTWData(); // failed
			}
			queryVector.push_back(stod(queryTimePoint));
		}

		// check best distance for first set of data
		// std::cout << "at timeseries: " << curTimeSeries << " location: " << curLocation << std::endl;
		double newBest = std::min(simpleDTW(dataVector, queryVector), bestMatchDistance);
		if (newBest != bestMatchDistance) {
			bestMatchDistance = newBest;
			bestMatchLocation = curLocation;
			bestMatchTimeSeries = curTimeSeries;
		}

		std::getline(query, queryTimePoint, ' ');
		// run through the rest of the query timeseries
		while (query.good()) {
			curLocation++;
			// check if we hit the end of the time series
			// ex) "number1\nnumber2"
			queryNewLine = queryTimePoint.find("\n");
			if (queryNewLine != std::string::npos) {	// hit the end of the timeseries
				queryVector.clear();
				queryVector.push_back(stod(queryTimePoint));
				break;
			}
			// remove the first data point
			queryVector.erase(queryVector.begin());
			// add the new data point
			queryVector.push_back(stod(queryTimePoint));

			// check new best distance
			// std::cout << "at timeseries: " << curTimeSeries << " location: " << curLocation << std::endl;
			double newBest = std::min(simpleDTW(dataVector, queryVector), bestMatchDistance);
			if (newBest != bestMatchDistance) {
				bestMatchDistance = newBest;
				bestMatchLocation = curLocation;
				bestMatchTimeSeries = curTimeSeries;
			}

			// get next time value
			std::getline(query, queryTimePoint, ' ');
		}
	}
	return DTWData(bestMatchDistance, bestMatchTimeSeries, bestMatchLocation);
}

// takes the given file and adds '_PAA' to the end
// returns the new filepath
std::string createPAAFilepath(std::string in) {
	std::size_t point = in.find(".");
	if (point == std::string::npos)
		return in + "_PAA";
	return in.substr(0, point) + "_PAA" + in.substr(point);
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


// CREATE NEW FLAG FOR OUTLIER.
// OUTLIER WILL CHECK A FILE AGAINST ITSELF TO FIND THE MAX DTW DISTANCE
int main(int argc, char** argv) {
	if (argc != 4 && argc != 5) {
		usage();
		return -1;
	}

	struct timeval start, stop;
	start.tv_usec = 0;
	stop.tv_usec = 0;

	if ((!std::string(argv[1]).compare("-PAA")) && (argc == 4)) {
		gettimeofday(&start, NULL);
		std::string PAAFile = PAAaFile(argv[2], atoi(argv[3]));
		gettimeofday(&stop, NULL);
		if (PAAFile.empty())
			return -1;
		std::cout << "PAA computed for " << argv[2] << ". Written to " << PAAFile << std::endl;
		long s = (stop.tv_sec - start.tv_sec);
		long t = (stop.tv_usec - start.tv_usec)/1000;
		std::cout << "PAA took: " << s << " seconds and " << t << " milliseconds to compute" << std::endl;
	} else if ((!std::string(argv[1]).compare("-DTW")) && (argc == 5)) {
		gettimeofday(&start, NULL);
		DTWData d = DTWaFile(argv[2], argv[3], atoi(argv[4]));
		gettimeofday(&stop, NULL);
		if (d.empty())
			return -1;
		std::cout << "DTW computed for data " << argv[2] << " and query " << argv[3];
		std::cout << ". Best distance found was: " << std::endl << d << std::endl;
		long s = (stop.tv_sec - start.tv_sec);
		long t = (stop.tv_usec - start.tv_usec)/1000;
		std::cout << "DTW took: " << s << " seconds and " << t << " microseconds to compute" << std::endl;
	} else {
		usage();
		return -1;
	}

	return 1;
}
