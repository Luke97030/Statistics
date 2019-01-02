// Program: Haoyu Kong
// Start Date: 2018/06/10
// End Date: 2018/06/14


#include <iostream>
#include <vector>
#include <string>
#include <exception>
#include <iterator>
#include <algorithm>
#include <iomanip>  // for manipulating the input and output
#include <tuple> 

using namespace	std;



/*
*	@ Purpose: to check the subStr of a string is a numberial data
*	@ Parameter: taking the substring as parameter
*	@ Return: return true or false depends on the result of converting
*/
bool subStrIsANum(string s)
{

	std::string::size_type sz;
	bool is_valid = false;
	double parsed_number = 0.0;

	try
	{
		parsed_number = std::stod(s, &sz);
		is_valid = true;
	}
	catch (...)
	{
		is_valid = false;
	}

	return is_valid;
}

/*
*	@ Purpose: to find the minimum number in two different vectors
*	@ Parameter 1: takes the first vector data
*	@ Parameter 2: takes the second vector data
*	@ Return: return two minimum number for two vectors respectively
*/
vector<double> findMinNum(vector<double> v1, vector<double> v2)
{
	vector<double> minV;
	double v1Min = *min_element(v1.begin(), v1.end());
	double v2Min = *min_element(v2.begin(), v2.end());
	minV.push_back(v1Min);
	minV.push_back(v2Min);
	return minV;
}

/*
*	@ Purpose: to find the max number in two different vectors
*	@ Parameter 1: takes the first vector data
*	@ Parameter 2: takes the second vector data
*	@ Return: return two minimum number for two vectors respectively
*/
vector<double> findMaxNum(vector<double> v1, vector<double> v2)
{
	vector<double> maxV;
	double v1Max = *max_element(v1.begin(), v1.end());
	double v2Max = *max_element(v2.begin(), v2.end());
	maxV.push_back(v1Max);
	maxV.push_back(v2Max);
	return maxV;
}


/*
*	@ Purpose: partition algorithm returns the high index of the lower partition
*	@ Parameter 1: vector
*	@ Parameter 2: 0
*	@ Parameter 3: vector.end() - 1
*	@ Return: int
*/
double partition(vector<double> &vi, int low, int up)
{
	auto pivot = vi[up];
	int i = low - 1;
	for (int j = low; j < up; j++)
	{
		if (vi[j] <= pivot)
		{
			i++;
			swap(vi[i], vi[j]);
		}
	}
	swap(vi[i + 1], vi[up]);
	return i + 1;
}

/*
*	@ Purpose: call the partition function and call it-self
*	@ Parameter 1: vector
*	@ Parameter 1: 0
*	@ Parameter 1: vi.size() - 1
*	@ Return: none
*/
void quickSort(vector<double> &vi, int low, int up)
{
	if (low < up)
	{
		int mid = int(partition(vi, low, up));
		// The mid position is on the place, so we don't need to consider it again.  
		//That's why below is mid-1, not mid! Otherwise it will occur overflow error!!!  
		quickSort(vi, low, mid - 1);
		quickSort(vi, mid + 1, up);
	}
}

/*
*	@ Purpose: call the qucikSort function
*	@ Parameter 1: vector
*	@ Return: none
*/
void qSort(vector<double> &vi)
{
	quickSort(vi, 0, vi.size() - 1);
}

/*
*	@ Purpose: find out the mod between vector elements
*	@ Parameter 1: vector
*	@ Parameter 2: size of vector
*	@ Return: a turple, bacause i want to return the frequence as an int and return the mode number as double
*/
tuple<int, double> searchVectorForMod(vector<double> v, int size)
{
	// create a new vector with same size as parameter 1
	vector<int> modCounter;
	modCounter.resize(size);
	// push 1 to the vector as frequence
	for (int i = 0; i < size; i++)
	{
		modCounter[i] = 1;
	}

	int maxAccureCount = 0;
	double mode = 0.0;
	int printCount = 0;


	// got through vector
	for (int i = 1; i < size; i++)
	{
		// v[1] = v[0]
		if (v[i] == v[i - 1])
		{
			// modeCounter + 1 if find same elements in the vector
			modCounter[i] = modCounter[i - 1] + 1;
		}

		if (maxAccureCount < modCounter[i])
			maxAccureCount = modCounter[i];
	}

	// there are duplicate numner in vector
	if (maxAccureCount > 1)
	{
		for (int i = 0; i < size; i++)
		{
			if (modCounter[i] == maxAccureCount)
			{
				mode = v[i];
				printCount++;
			}
		}
	}
	if (printCount > 1) {
		mode = 0;
		maxAccureCount = 0;

	}

	return tuple<int, double>(maxAccureCount, mode);
}


void checkMode(tuple<int, double> v1Mode, tuple<int, double> v2Mode)
{
	// for vector 1 mode 
	if (std::get<0>(v1Mode) == 1 && std::get<0>(v2Mode) == 1)
	{
		cout << "mode \t\t\t\t" << "no mode" << "\t\t" << "no mode" << endl;
	}
	else if (std::get<0>(v1Mode) > 1 && std::get<0>(v2Mode) == 1)
	{
		cout << "mode\t\t\t\tfreq.= " << std::get<0>(v1Mode) << "\t" << "no mode" << endl;
		cout << "\t\t\t\t" << std::get<1>(v1Mode) << endl;
	}
	else if (std::get<0>(v1Mode) > 1 && std::get<0>(v2Mode) == 0)
	{
		cout << "mode\t\t\t\tfreq.= " << std::get<0>(v1Mode) << "\t" << "no mode" << endl;
		cout << "\t\t\t\t" << std::get<1>(v1Mode) << endl;
	}
	else if (std::get<0>(v1Mode) == 1 && std::get<0>(v2Mode) > 1)
	{
		cout << "mode\t\t\t\t" << "no mode" << "\t\tfreq.= " << std::get<0>(v2Mode) << endl;
		cout << "\t\t\t\t\t\t" << std::get<1>(v2Mode) << endl;
	}
	else if (std::get<0>(v1Mode) == 0 && std::get<0>(v2Mode) > 1)
	{
		cout << "mode\t\t\t\t" << "no mode" << "\t\tfreq.= " << std::get<0>(v2Mode) << endl;
		cout << "\t\t\t\t\t\t" << std::get<1>(v2Mode) << endl;
	}
	else if (std::get<0>(v1Mode) > 1 && std::get<0>(v2Mode) > 1)
	{
		cout << "mode\t\t\t\tfreq.= " << std::get<0>(v1Mode) << "\tfreq.= " << std::get<0>(v2Mode) << endl;
		cout << "\t\t\t\t" << std::get<1>(v1Mode) << "\t\t" << std::get<1>(v2Mode) << endl;
	}

}





/*
*	@ Purpose: return absolute value
*	@ Parameter 1: the result of an element in vector - mean
*	@ Return: double
*/
double absoluteValue(double d)
{
	if (d < 0)
		return 0 - d;
	else
		return d;
}



/*
*	@ Purpose: find out the MeanA bsoulte Deviation About Mode for two data sets
*	@ Parameter 1: vector
*	@ Parameter 2: can be mode
*	@ Parameter 3: size of vector
*	@ Return: double
*/
void checkModeForDeviation(tuple<int, double> v1Mode, tuple<int, double> v2Mode, vector<double> vec1, vector<double> vec2)
{
	double result;
	double temp = 0.0;
	double temp2 = 0.0;
	double result2;
	// no mode for both data sets
	if (std::get<0>(v1Mode) == 1 && std::get<0>(v2Mode) == 1)
	{
		cout << "mode \t\t\t\t" << "no mode" << "\t\t" << "no mode" << endl;
	}
	// vector 1 has mode, vector 2 not
	else if (std::get<0>(v1Mode) > 1 && std::get<0>(v2Mode) == 1)
	{
		for (unsigned i = 0; i < vec1.size(); i++)
		{
			temp += absoluteValue(vec1[i] - std::get<1>(v1Mode));
		}
		result = temp / vec1.size();
		cout << "->about the mode\t\t" << result << "\t\t" << "no mode" << endl;
	}
	// vector 1 has mode, vector 2 not
	else if (std::get<0>(v1Mode) > 1 && std::get<0>(v2Mode) == 0)
	{
		for (unsigned i = 0; i < vec1.size(); i++)
		{
			temp += absoluteValue(vec1[i] - std::get<1>(v1Mode));

		}
		result = temp / vec1.size();
		cout << "->about the mode\t\t" << result << "\t\t" << "no mode" << endl;
	}

	// vector 2 has mode, vector 1 not
	else if (std::get<0>(v1Mode) == 1 && std::get<0>(v2Mode) > 1)
	{
		for (unsigned i = 0; i < vec2.size(); i++)
		{
			temp += absoluteValue(vec2[i] - std::get<1>(v2Mode));

		}
		result = temp / vec2.size();
		cout << "->about the mode\t\t" << "no mode" << "\t\t" << result << endl;

	}
	// vector 2 has mode, vector 1 not
	else if (std::get<0>(v1Mode) == 0 && std::get<0>(v2Mode) > 1)
	{
		for (unsigned i = 0; i < vec2.size(); i++)
		{
			temp += absoluteValue(vec2[i] - std::get<1>(v2Mode));
		}
		result = temp / vec2.size();
		cout << "->about the mode\t\t" << "no mode" << "\t\t" << result << endl;
	}

	// both data sets have mode 
	else if (std::get<0>(v1Mode) > 1 && std::get<0>(v2Mode) > 1)
	{
		for (unsigned i = 0; i < vec1.size(); i++)
		{
			temp += absoluteValue(vec1[i] - std::get<1>(v1Mode));
		}
		result = temp / vec1.size();

		for (unsigned i = 0; i < vec2.size(); i++)
		{
			temp2 += absoluteValue(vec2[i] - std::get<1>(v2Mode));
		}
		result2 = temp2 / vec2.size();
		cout << "->about the mode\t\t" << result << "\t\t" << result2 << endl;

	}

}


/*
*	@ Purpose: find out the Variance between vector elements
*	@ Parameter 1: vector
*	@ Parameter 2: mean of the vector
*	@ Parameter 3: size of vector
*	@ Return: double
*/
double calcVariance(vector<double> vec, double mean, int size)
{
	double difference = 0;
	double variance = 0;

	for (int i = 0; i < size; ++i)
	{
		difference = absoluteValue(vec[i] - mean);
		variance += difference * difference;
	}

	variance = variance / size;

	return variance;
}


/*
*	@ Purpose: find out the Deviation by sqrt variance
*	@ Parameter 1: variance
*	@ Return: double
*/
double calculateDeviation(double variance)
{
	return sqrt(variance);
}

/*
*	@ Purpose: find out the MeanA bsoulte Deviation About Mean, Median for two data sets
*	@ Parameter 1: vector
*	@ Parameter 2: can be mean, median
*	@ Parameter 3: size of vector
*	@ Return: double
*/
double calculateMeanAbsoulteDeviation(vector<double> vec, double mmm, int size)
{
	double result;
	double temp = 0.0;
	for (int i = 0; i < size; i++)
	{
		temp += absoluteValue(vec[i] - mmm);
		
	}
	result = temp / size;
	return result;
}



/*
*	@ Purpose: find out regression line a
*	@ Parameter 1: vector1,
*	@ Parameter 2: vector2,
*	@ Parameter 3: vector1 size,
*	@ Parameter 4: vector2 size,
*	@ Return: double
*/
double findRegressionLineB(vector<double> v1, vector<double> v2)
{	
	double regressionB = 0.0;
	// only can calculate the regression line a and b when user enter two same size data sets
	if (v1.size() == v2.size())
	{
		int n = v1.size();
		double sumX = 0.0;
		double sumY = 0.0;
		double sumDoubleX = 0.0;
		double sumDoubleY = 0.0;
		double sumXY = 0.0;
		
		for (unsigned i = 0; i < v1.size(); i++)
		{
			sumX += v1[i];
		}

		for (unsigned i = 0; i < v1.size(); i++)
		{
			sumDoubleX += v1[i] * v1[i];
		}
		for (unsigned i = 0; i < v2.size(); i++)
		{
			sumY += v2[i];
		}

		for (unsigned i = 0; i < v2.size(); i++)
		{
			sumDoubleY += (v2[i] * v2[i]);
		}

		for (unsigned i = 0; i < v2.size(); i++)
		{
			sumXY += v1[i] * v2[i];
		}
		regressionB = (n * sumXY - sumX * sumY) / (n * sumDoubleX - sumX * sumX);
	}
	else
	{
		cout << "The number of X and Y are not equal, so cannot calculate regression line a and b" << endl;
		return 0;
	}

	return regressionB;
}


/*
*	@ Purpose: find out regression line a
*	@ Parameter 1: vector1,
*	@ Parameter 2: vector2,
*	@ Parameter 3: vector1 size,
*	@ Parameter 4: vector2 size,
*	@ Return: double
*/
double findRegressionLineA(double regressionB, double xMean, double yMean) 
{
	double regressionA = 0.0;
	regressionA = yMean - regressionB * xMean;
	return regressionA;
}


/*
*	@ Purpose: find out outliers by comparing with 2 * deviation
*	@ Parameter 1: vector,
*	@ Parameter 2: mean,
*	@ Parameter 3: deviation,
*	@ Return: none
*/
void find2DOutliers(vector<double> vec, double mean, double deviation)
{
	vector<double> newVec;
	vector<double> largeThan2DVec;
	double finalDeviation = 0.0;
	for (unsigned i = 0; i < vec.size(); i++)
	{
		newVec.push_back(absoluteValue(vec[i] - mean));
	}
	finalDeviation = 2 * deviation;

	for (unsigned i = 0; i < newVec.size(); i++)
	{
		if (newVec[i] > finalDeviation)
		{
			largeThan2DVec.push_back(newVec[i]);
		}
	}

	// this means there are outliers
	if (largeThan2DVec.size() > 0)
	{
		cout << "# outliers = " << largeThan2DVec.size() << endl;
		cout << "\t\t\t\t\t\t";
		for (unsigned i = 0; i < largeThan2DVec.size(); i++)
		{
			cout << vec[i] << ", ";
		}
		cout << endl;
	}
	else
	{
		cout << "no outliers\t";
	}

}


/*
*	@ Purpose: find out outliers by comparing with 3 * deviation
*	@ Parameter 1: vector,
*	@ Parameter 2: mean,
*	@ Parameter 3: deviation,
*	@ Return: none
*/
void find3DOutliers(vector<double> vec, double mean, double deviation)
{
	vector<double> newVec;
	vector<double> largeThan3DVec;
	double finalDeviation = 0.0;
	for (unsigned i = 0; i < vec.size(); i++)
	{
		newVec.push_back(absoluteValue(vec[i] - mean));
	}
	finalDeviation = 3 * deviation;

	for (unsigned i = 0; i < newVec.size(); i++)
	{
		if (newVec[i] > finalDeviation)
		{
			largeThan3DVec.push_back(newVec[i]);
		}
	}

	// this means there are outliers
	if (largeThan3DVec.size() > 0)
	{
		cout << "# outliers = " << largeThan3DVec.size() << endl;
		cout << "\t\t\t\t\t\t";
		for (unsigned i = 0; i < largeThan3DVec.size(); i++)
		{
			cout << vec[i] << ", ";
		}
		cout << endl;
	}
	else
	{
		cout << "no outliers\t";
	}

}

int main()
{
	// all the variables 
	vector<double> firstV;
	vector<double> secondV;
	string wholeInput = "";
	double firstVInput = 0.0;
	double seoncdVInput = 0.0;
	double v1Min = 0.0;
	double v2Min = 0.0;
	double v1Max = 0.0;
	double v2Max = 0.0;
	double v1Total = 0.0;
	double v2Total = 0.0;
	double v1Mean = 0.0;
	double v2Mean = 0.0;
	double v1Median = 0.0;
	double v2Median = 0.0;
	tuple<int, double> v1Mode;
	tuple<int, double> v2Mode;
	double v1Variance = 0.0;
	double v2Variance = 0.0;
	double v1Deviation = 0.0;
	double v2Deviation = 0.0;
	double v1MeanAbsoulteDeviationAboutMean = 0.0;
	double v2MeanAbsoulteDeviationAboutMean = 0.0;
	double v1MeanAbsoulteDeviationAboutMedian = 0.0;
	double v2MeanAbsoulteDeviationAboutMedian = 0.0;
	double v1MeanAbsoulteDeviationAboutMode= 0.0;
	double v2MeanAbsoulteDeviationAboutMode = 0.0;
	double v1Outliers2D = 0.0;
	double v2Outliers2D = 0.0;
	double v1Outliers3D = 0.0;
	double v2Outliers3D = 0.0;
	double a = 0.0;
	double b = 0.0;



	cout << "States, (c) 2018 Haoyu Kong" << endl;
	cout << "!!!!Explaination: for the statistics program, "<< endl;
	cout << "!!!!you have to press ENTER after you enter two number which separate by comma" << endl;
	cout << "!!!!Therefore, Please follow the rules :)" << endl;
	cout << endl << endl <<endl;
	cout << "Enter a list of comma-separated real number paris terminated by ^Z" << endl;
	cout << endl;


	while (!cin.eof()) {
		bool lIsNum = false;
		bool rIsNum = false;

		// get the user input, using getline to allow user to enter any sapce as they want before they enter ENTER key
		std::getline(std::cin, wholeInput);
		// the program allows user enter white space at anywhere when they are inputing 
		std::string::iterator end_pos = remove(wholeInput.begin(), wholeInput.end(), ' ');
		// remove all the white space between a string
		// get the wholeInput without any space.
		wholeInput.erase(end_pos, wholeInput.end());
		// find the index of "," in a string
		std::size_t found = wholeInput.find(',');

		// from the beginning of a string to one digit left hand side of common 
		string subStringFirstV = wholeInput.substr(0, found);
		lIsNum = subStrIsANum(subStringFirstV);
		// from one digit right hand side of common to the end
		string subStringSecondV = wholeInput.substr(found + 1, wholeInput.length());
		rIsNum = subStrIsANum(subStringSecondV);

		//convert the two substr to double and assign them to local variable 
		if (lIsNum == true && rIsNum == true)
		{
			firstVInput = atof(subStringFirstV.c_str());
			seoncdVInput = atof(subStringSecondV.c_str());
			firstV.push_back(firstVInput);
			secondV.push_back(seoncdVInput);
		}
		else if (firstV.size() == 0 && secondV.size() == 0)
		{
			cout << "Empty data set!!!" << endl;
			break;
		}
		else
		{
			break;
		}

	}


	cout << "Results: " << endl;
	cout << "----------------------------------------------------------------" << endl;
	cout << "# elements    \t\t\t" << firstV.size() << "\t\t" << secondV.size() << endl;



	// print out the minimun number in two different sets
	auto minVector = findMinNum(firstV, secondV);
	v1Min = minVector[0];
	v2Min = minVector[1];
	cout << setprecision(3) << fixed << "minimum \t\t\t" << v1Min << "\t\t" << v2Min << endl;

	// print out the maximum number in two different sets
	auto maxVector = findMaxNum(firstV, secondV);
	v1Max = maxVector[0];
	v2Max = maxVector[1];
	cout << setprecision(3) << fixed << "maximum \t\t\t" << v1Max << "\t\t" << v2Max << endl;

	// calculate mean value
	for (unsigned int i = 0; i < firstV.size(); i++)
	{
		v1Total += firstV[i];
	}

	for (unsigned int i = 0; i < secondV.size(); i++)
	{
		v2Total += secondV[i];
	}
	v1Mean = v1Total / firstV.size();
	v2Mean = v2Total / secondV.size();
	cout << setprecision(3) << fixed << "mean \t\t\t\t" << v1Mean << "\t\t" << v2Mean << endl;

	// find regression 
	// the all data in the first data set will be x
	// the all data in the second data set will be y 
	// We have to calculate the regression line a and b before we sort our data sets
	b = findRegressionLineB(firstV, secondV);
	if (b != 0)
	{
		a = findRegressionLineA(b, v1Mean, v2Mean);
	}
	else
	{
		cout << "Cannot find the regression line b!!!!" << endl;
	}

	// calculate median value
	// we need to sort the vactors before we do this
	//sort array using quickSort()
	qSort(firstV);

	//// Printing the sorted data.
	//cout << "\nSorted Data For First Vector";

	//for (unsigned int i = 0; i < firstV.size(); i++)
	//	cout << "->" << firstV[i];


	////cout << endl;
	qSort(secondV);

	//// Printing the sorted data.
	//cout << "\nSorted Data For Second Vector";

	//for (unsigned int i = 0; i < secondV.size(); i++)
	//	cout << "->" << secondV[i];

	//cout << endl;

	//find middle position of array for calculation median
	int v1MiddleValue = (firstV.size() + 1) / 2;
	int v2MiddleValue = (secondV.size() + 1) / 2;
	// calculate median for the first vector 
	if ((firstV.size() % 2) == 0)
	{
		// if the lenght of array is even
		v1Median = (firstV[v1MiddleValue - 1] + firstV[v1MiddleValue]) / 2;
	}
	else
	{
		// if the lenght of array is odd
		v1Median = firstV[v1MiddleValue - 1];
	}

	// calculate median for the second vector 
	if ((secondV.size() % 2) == 0)
	{
		// if the lenght of array is even
		v2Median = (secondV[v2MiddleValue - 1] + secondV[v2MiddleValue]) / 2;
	}
	else
	{
		// if the lenght of array is odd
		v2Median = secondV[v2MiddleValue - 1];
	}
	cout << setprecision(3) << fixed << "median \t\t\t\t" << v1Median << "\t\t" << v2Median << endl;




	// finding mode in two sets 
	v1Mode = searchVectorForMod(firstV, firstV.size());
	v2Mode = searchVectorForMod(secondV, secondV.size());
	checkMode(v1Mode, v2Mode);


	// calculate variance 
	v1Variance = calcVariance(firstV, v1Mean, firstV.size());
	v2Variance = calcVariance(secondV, v2Mean, secondV.size());
	cout << setprecision(3) << fixed << "variance \t\t\t" << v1Variance << "\t\t" << v2Variance << endl;

	// calculate Standard  Deviation
	v1Deviation = calculateDeviation(v1Variance);
	v2Deviation = calculateDeviation(v2Variance);
	cout << setprecision(3) << fixed << "std. dev. \t\t\t" << v1Deviation << "\t\t" << v2Deviation << endl;

	// calculate mean absoulte deviations about mean
	cout << "mean absoulte deviations" << endl;
	v1MeanAbsoulteDeviationAboutMean = calculateMeanAbsoulteDeviation(firstV, v1Mean, firstV.size());
	v2MeanAbsoulteDeviationAboutMean = calculateMeanAbsoulteDeviation(secondV, v2Mean, secondV.size());
	cout << "->about the mean\t\t" << v1MeanAbsoulteDeviationAboutMean << "\t\t" << v2MeanAbsoulteDeviationAboutMean << endl;

	// calculate mean absoulte deviations about median
	v1MeanAbsoulteDeviationAboutMedian = calculateMeanAbsoulteDeviation(firstV, v1Median, firstV.size());
	v2MeanAbsoulteDeviationAboutMedian = calculateMeanAbsoulteDeviation(secondV, v2Median, secondV.size());
	cout << "->about the median\t\t" << v1MeanAbsoulteDeviationAboutMedian << "\t\t" << v2MeanAbsoulteDeviationAboutMedian << endl;

	// calculate mean absoulte deviations about mode
	checkModeForDeviation(v1Mode, v2Mode, firstV, secondV);


	// print out the regression line a and b
	cout << setprecision(3) << fixed << "regression line\t\t\t" << "a = " << a << "\t" << "b = " << b << endl;

	// calculate outliers 
	cout << "Outliers(2X)\t\t\t";
	find2DOutliers(firstV, v1Mean, v1Deviation);
	find2DOutliers(secondV, v2Mean, v2Deviation);
	cout << endl;

	cout << "Outliers(3X)\t\t\t";
	find3DOutliers(firstV, v1Mean, v1Deviation);
	find3DOutliers(secondV, v2Mean, v2Deviation);
	cout << endl;

}