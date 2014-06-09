/*
	Filename:		ProxAn_main.cpp
	Author:			Andrew Kienapple / Bradley Walker
	Date:			November 18, 2013
	Description:	An application that generates statistics about the proximity of 
					residential addresses in the City of Toronto to the available 
					Food Banks.
	Changelog:		November 23, 2013
					Updated with parallel capabilities, changed Msg_t percentage variable
					from double to float to fix overload error in passing floating point
					through MPI message passing
*/


#include <iostream>
#include <string>
#include <utility> //std::pair
#include <math.h>
#include <mpi.h>
#include <vector>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <ctime>
using namespace std;

//input file names
const string FOODBANKS_FN = "foodbanks.dat";
const string RESIDENCES_FN = "residences.dat";

//process tags and msg buffer
const int TAG_DATA = 0, TAG_QUIT = 1;
const int MSG_SIZE = 100;

typedef struct {
	int c1, c2, c3, c4;
	float p1, p2, p3, p4;
} Msg_t;

//c1 counts the distances that are 0-1 km away, c2 counts the distances that are 1-2 km and so on
//p1 stores the percentage of 0-1 km foodbanks and so on.

bool getCoordinates( ifstream &in, vector<pair<double,double>> &vCoords) {
	//Fill a vector with the coordinate values
	istream_iterator<double> start(in), end;
	vector<double> numbers(start, end);	
	double xCoord, yCoord;

	//Pair the values together and place back into vCoords
	for( vector<double>::iterator it = numbers.begin(); it != numbers.end(); ++it ) {
		xCoord = *it;
		++it;
		if( it != numbers.end() )
			yCoord = *it;
		else {
			cerr << "Error reading input file " << FOODBANKS_FN << endl;
			cerr << "There may be an odd number of values. Check file." << endl;
			return false;
		}

		vCoords.push_back( make_pair( xCoord, yCoord ) );
	}
	return true;
}

double calcDistance( double xRes, double yRes, double xFood, double yFood ) {
	return abs( sqrt( pow( xFood - xRes, 2 ) + pow( yFood - yRes, 2 ) ) ) / 1000;
}

bool countClosestProximities( ifstream &in, 
							 vector<pair<double, double>> &vFoodBanks, 
							 Msg_t &counts, int numProcs, int rank  ){
	vector<double> vDists;
	int nextIt = rank;
	double xRes, yRes;
	double shortestDist = 50000.0; //I chose 50000 because I need a number guaranteed 
								   //to be longer than any calculated distance in order to
								   //find the shortest distance. 50000 km is about as long as the Earth's circumference.
	int i = 0;
	// Start reading the files and doing calculations
	while( in >> xRes >> yRes ) {

		if(i == nextIt) // 
		{
			nextIt += numProcs;

			for( vector<pair<double,double>>::iterator it = vFoodBanks.begin(); it != vFoodBanks.end(); ++it ) {
				vDists.push_back( calcDistance( xRes, yRes, it->first, it->second ) );

			}
		
			//find the shortest distance that was calculated for this residence
			for( vector<double>::iterator it = vDists.begin(); it != vDists.end(); ++it ) {
				if( *it < shortestDist )
					shortestDist = *it;
			}
		
			//Increment the appropriate counter value
			if( shortestDist < 1 )
				++counts.c1;
			else if( shortestDist >= 1 && shortestDist < 2 )
				++counts.c2;
			else if( shortestDist >= 2 && shortestDist < 5 )
				++counts.c3;
			else if( shortestDist >= 5 )
				++counts.c4;

			vDists.clear();
		}
		
		shortestDist = 50000.0;
	}

	return true;
}

// Set the percentage values
void calculatePercentages( Msg_t &counts ) {
	int total = counts.c1 + counts.c2 + counts.c3 + counts.c4;

	counts.p1 = (counts.c1 / (float)total) * 100.0f;
	counts.p2 = (counts.c2 / (float)total) * 100.0f;
	counts.p3 = (counts.c3 / (float)total) * 100.0f;
	counts.p4 = (counts.c4 / (float)total) * 100.0f;
}

// Master Process Only
void printResults( Msg_t counts, double elapsedTime, int rank ) {

	if (rank == 1)
	{
		cout << "Proximity of Residential Addresses to Foodbanks in Toronto" << endl <<
				"----------------------------------------------------------" << endl;

		cout <<  "Elapsed Time in Seconds:\t" << elapsedTime << endl << endl;
	}

	// if rank is 0 (to be passed at end of master process) then use final statement
	if (rank != 0)
		cout << "Process #" << (rank) << " results for " << (counts.c1 + counts.c2 + counts.c3 + counts.c4) << " addresses..." << endl;
	else
		cout << "Aggregated results for all " << (counts.c1 + counts.c2 + counts.c3 + counts.c4) << " addresses..." << endl;

	cout << "Nearest Foodbank (km)" 
		<< right << setfill(' ') << setw(19) << "# of Addresses" 
		<< right << setfill(' ') << setw(17) << "% of Addresses" << endl
		<< "---------------------"
		<< right << setfill(' ') << setw(19) << "--------------"
		<< right << setfill(' ') << setw(17) << "--------------" << endl;

	cout << left << setfill(' ') << setw(21) << "        0 - 1" 
		<< right << setfill(' ') << setw(18) << counts.c1
		<< right << setfill(' ') << setw(16) << fixed << setprecision(4) << counts.p1 << endl;

	cout << left << setfill(' ') << setw(21) << "        1 - 2" 
		<< right << setfill(' ') << setw(18) << counts.c2
		<< right << setfill(' ') << setw(16) << fixed << setprecision(4) << counts.p2 << endl;

	cout << left << setfill(' ') << setw(21) << "        2 - 5" 
		<< right << setfill(' ') << setw(18) << counts.c3
		<< right << setfill(' ') << setw(16) << fixed << setprecision(4) << counts.p3 << endl;

	cout << left << setfill(' ') << setw(21) << "      n    > 5" 
		<< right << setfill(' ') << setw(18) << counts.c4
		<< right << setfill(' ') << setw(16) << fixed << setprecision(2) << counts.p4 << endl << endl;

}

// Creates the derived type
MPI_Datatype createMsgType() {
	MPI_Datatype new_type;
	int count = 2;
	int blocklens[] = { 4, 4 };
	
	MPI_Aint indices[2];
	indices[0] = 0;
	MPI_Type_extent( MPI_FLOAT, &indices[1] );
	indices[1] *= 4;

	MPI_Datatype old_types[] = { MPI_INT, MPI_FLOAT };

	//Call datatype constructor
	MPI_Type_struct( count, blocklens, indices, old_types, &new_type);
	MPI_Type_commit(&new_type);

	return new_type;
}

void processSlave(int numProcs, int rank)
{
	MPI_Datatype msgType = createMsgType();
	
	vector<pair<double, double>> vFoodBanks;

	double elapsedTime = (double)time(0);

	Msg_t counts;
	counts.c1 = 0; counts.c2 = 0; counts.c3 = 0; counts.c4 = 0;

	ifstream rin( RESIDENCES_FN );
	ifstream fin( FOODBANKS_FN );

	// Check stream integrities
	if( !rin.is_open() || rin.fail() ) {
		cerr << "Error opening input file " << RESIDENCES_FN << endl;
	}
	else if (!fin.is_open() || fin.fail() ) {
		cerr << "Error opening input file " << FOODBANKS_FN << endl;
	}

	if( !getCoordinates( fin, vFoodBanks) )
		cerr << "Error opening input file at Process #" << rank << endl;
	fin.close();

	countClosestProximities( rin, vFoodBanks, counts, numProcs, rank );

	calculatePercentages( counts );

	rin.close();

	MPI_Send(&counts, 2, msgType, 0, TAG_DATA,
						MPI_COMM_WORLD);

	//MPI_Send();
}

void processMaster(int numProcs)
{
	// MPI Datatypes
	MPI_Status status;
	MPI_Datatype msgType = createMsgType();

	vector<pair<double, double>> vFoodBanks;

	double elapsedTime = (double)time(0);

	// Create and initialize the message struct and container
	Msg_t counts;
	Msg_t msgCounts; // counts from message passing
	Msg_t finalCounts;
	vector<pair<Msg_t, int>> vecCount;

	counts.c1 = 0; counts.c2 = 0; counts.c3 = 0; counts.c4 = 0;
	finalCounts.c1 = 0; finalCounts.c2 = 0; finalCounts.c3 = 0; finalCounts.c4 = 0;

	ifstream rin( RESIDENCES_FN );
	ifstream fin( FOODBANKS_FN );

	// Check stream integrities
	if( !rin.is_open() || rin.fail() ) {
		cerr << "Error opening input file " << RESIDENCES_FN << endl;
	}
	else if (!fin.is_open() || fin.fail() ) {
		cerr << "Error opening input file " << FOODBANKS_FN << endl;
	}

	if( !getCoordinates( fin, vFoodBanks) )
		cerr << "Error opening input file at Master Process" << endl;
	fin.close();

	countClosestProximities( rin, vFoodBanks, counts, numProcs, 0);

	calculatePercentages( counts );

	// push master process calculations to vector
	vecCount.push_back(pair<Msg_t, int>(counts, 1));

	rin.close();

	// message count
	int i = 1;

	// enter only if slave processes are being utilized in instance
	while(numProcs != 1)
	{
		MPI_Recv(&msgCounts, 2, msgType, MPI_ANY_SOURCE, MPI_ANY_TAG, 
			MPI_COMM_WORLD, &status);

		if( status.MPI_TAG == TAG_DATA )
		{
			// push recieved slave calculations and process number to vector
			vecCount.push_back(pair<Msg_t, int>(msgCounts, (status.MPI_SOURCE+1)));
			++i;
		}

		// if messages recieved equal the total of slave processes, then break. OR if number of slace processes is just one, break.
		if (i == numProcs && numProcs != 2)
			break;
	}

	elapsedTime = (double)time(0) - elapsedTime;

	// Print results for each process
	for(auto & p : vecCount)
	{
		finalCounts.c1 += p.first.c1; finalCounts.c2 += p.first.c2; finalCounts.c3 += p.first.c3; finalCounts.c4 += p.first.c4;

		printResults(p.first, elapsedTime, p.second);
	}

	// Calculate precentages for total addresses
	calculatePercentages( finalCounts );

	// Final Print
	printResults( finalCounts, elapsedTime, 0);
}

int main( int argc, char* argv[] ) {	
	vector<pair<double, double>> vFoodBanks;

	// Create and initialize the message struct
	Msg_t counts;
	counts.c1 = 0; counts.c2 = 0; counts.c3 = 0; counts.c4 = 0;

	ifstream rin( RESIDENCES_FN );
	ifstream fin( FOODBANKS_FN );

	// Check stream integrities
	if( !rin.is_open() || rin.fail() ) {
		cerr << "Error opening input file " << RESIDENCES_FN << endl;
		return EXIT_FAILURE;
	}
	else if (!fin.is_open() || fin.fail() ) {
		cerr << "Error opening input file " << FOODBANKS_FN << endl;
		return EXIT_FAILURE;
	}

	if( MPI_Init(&argc, &argv) == MPI_SUCCESS )
	{
		// Obtain the rank and the # processes
		int numProcs, rank; 
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		if( rank == 0 )
			processMaster(numProcs);
		else
			processSlave(numProcs, rank);

		MPI_Finalize();
	}

	return 0;
}
