//
//  NumMethods.cpp
//  NS Code
//
//  Created by Marc Salinas on 1/29/19.
//  Copyright © 2019 Marc Salinas. All rights reserved.
//

#include "NumMethods.hpp"
#include "Conversions.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <omp.h>
#include <vector>

data2 dm;
Convert conv;
using namespace std;

const double pi = 4*atan(1);
const double mE = 0.511;        // Mass of electron (MeV)
const double mU = 2.3;          // Mass of up quark (MeV)
const double mD = 4.8;          // Mass of down quark (MeV)
const double mS = 95.0;         // Mass of strange quark (MeV)
const double mP = 939.0;    // Mass of proton (MeV)
const double mN = 939.0;    // Mass of neutron (MeV)
const double mMU = 105.6583745; // Mass of muon (MeV)

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

// Data Manipulation class
//*********************************************************************************
//*********************************************************************************

// find the number of rows in data
int data2 :: rowcount(string txtfile) {
    int count = 0;                  // Set initial count to 0
    string line;                    // Initialize variable to store line's string temporarily
    ifstream file(txtfile);         // Read the txt file
    if (!file) {                    // Throw error if file isn't found
        cout << "RowCount: Error opening file from path: " << txtfile << endl;
        file.close();
        exit(0);
    }
    
    // Use getline function from ifstream to count non empty lines
    while (getline(file, line)) {
        count++;
    }
    file.close();
    return count; // Return number of lines
}

//----------------------------------------------------------------------------------------

// identify the type of file
string data2 :: filetype(string txtfile) {
    
    ifstream file;
    string line;
    size_t pos;
    file.open(txtfile);
    
    if(!file){                      // Throw error if file isn't found
        cout << "Filetype: Error opening file from path: " << txtfile << endl;
        file.close();
        exit(0);
    }
    
    getline(file,line);             // get line from file
    pos=line.find(",");             // search
    if(pos!=string::npos) {         // string::npos is returned if string is not found
        file.close();
        return "CSV";
    } else {
        file.close();
        return "TXT";
    }
}

//----------------------------------------------------------------------------------------

// find number of columns in data
int data2 :: colcount(string txtfile) {
    string line;
    ifstream file(txtfile);
    if (!file) {                    // Throw error if file isn't found
        cout << "Colcount: Error opening file from path: " << txtfile << endl;
        file.close();
        exit(0);
    } else {
        if (filetype(txtfile) == "TXT") {
            int numcols=0;
            getline(file,line);
            stringstream iss(line);
            do {
                std::string sub;
                iss >> sub;
                if (sub.length()) {
                    ++numcols;
                }
            }
        while(iss);
        return numcols;
        } else {
            int numcols = 0;
            while(getline(file,line)) {
                stringstream linestream(line);
                string value;
                while(getline(linestream,value,',')) {
                    ++numcols;
                }
                return numcols;
            }
        }
    }
    return 0;
}

//----------------------------------------------------------------------------------------

// import data into an array
void data2 :: importdata(string txtfile, double ** &array) {
    int numrows = rowcount(txtfile);            // Get the number of rows
    int numcols = colcount(txtfile);            // Get the number of cols
    string type = filetype(txtfile);            // Get the file type
    ifstream in(txtfile);
    
    // Initialize array size
    array = new double*[numrows];
    for (int i = 0; i<numrows; i++) {
        array[i] = new double[numcols];
    }
    
    if (type == "TXT") {                        // Import data for txt
        for (int j=0; j<numrows; ++j) {
            for (int i=0; i<numcols; ++i) {
                in >> array[j][i];
            }
        }
        in.close();
    }
    if (type == "CSV") {                        // Import data for csv
        ifstream file;
        string line;
        int i=0;
        int j=0;
        file.open(txtfile);
        while(getline(file,line)) {
            stringstream linestream(line);
            string value;
            while(getline(linestream,value,',')) {
                //cout << i << "  " << j << endl;
                array[i][j] = stod(value);
                j=j+1;
            }
            i=i+1;
            j=0;
        }
        in.close();
    }
};

// import data into an array
void data2 :: importdata_string(string txtfile, string ** &array) {
    int numrows = rowcount(txtfile);            // Get the number of rows
    int numcols = colcount(txtfile);            // Get the number of cols
    string type = filetype(txtfile);            // Get the file type
    ifstream in(txtfile);
    
    // Initialize array size
    array = new string*[numrows];
    for (int i = 0; i<numrows; i++) {
        array[i] = new string[numcols];
    }
    
    if (type == "TXT") {                        // Import data for txt
        for (int j=0; j<numrows; ++j) {
            for (int i=0; i<numcols; ++i) {
                in >> array[j][i];
            }
        }
        in.close();
    }
    if (type == "CSV") {                        // Import data for csv
        ifstream file;
        string line;
        int i=0;
        int j=0;
        file.open(txtfile);
        while(getline(file,line)) {
            stringstream linestream(line);
            string value;
            while(getline(linestream,value,',')) {
                //cout << i << "  " << j << endl;
                array[i][j] = stod(value);
                j=j+1;
            }
            i=i+1;
            j=0;
        }
        in.close();
    }
};

//----------------------------------------------------------------------------------------

// Print out an array
void data2 :: print(double** &array, int numrows, int numcols, bool file, string filename) {
    if (file == true) {
        ofstream out(filename);
        for (int i=0; i<numrows; ++i) {
            for (int j=0; j<(numcols-1); ++j) {
                out << scientific << setprecision(15) << array[i][j] << "  ";
            }
        out << scientific << setprecision(15) << array[i][numcols-1] << endl;
        }
    } else {
        for (int i=0; i<numrows; ++i) {
            for (int j=0; j<numcols; ++j) {
                cout << array[i][j] << "  ";
            }
        cout << endl;
        }
    }
};

// Print out an array
void data2 :: convert_array(double** &array, int numrows, int numcols, int col_to_convert, double conv_factor) {
    for (int i=0; i<numrows; ++i) {
        array[i][col_to_convert] = conv_factor*array[i][col_to_convert];
    }
};
//------
//----------------------------------------------------------------------------------------

// Interpolate a value in an array
double data2 :: interpolate(int numrows, int numcols, double** array, double point, int pointcol, int ycol, bool sortflag) {
    
    // Error checking
    if (ycol >= numcols || ycol < 0) {                                              // check to see column exists
        cout << "Interpolation error: column is out of bounds" << endl;
        print(array,numrows,numcols,true,"debug.txt");
        exit(0);
    }

    if (numrows <= 1) {                                                           // check to see if interpolation is possible
        cout << "Interpolation error: not enough points to interpolate" << endl;
        exit(0);
    }

    // Sort data if unsorted
    if (sortflag == false) {
        sortasc(array, pointcol, numrows, numcols);
    }

    double eps = 1e-10*abs(array[numrows-1][pointcol]-array[numrows-2][pointcol]);
    if (point > array[numrows-1][pointcol]) {       // check to see if point is in bounds
        if ( abs(point-array[numrows-1][pointcol]) < eps) {
            return array[numrows-1][pointcol];
        } else {
            cout << "Interpolation error: point is not in bounds (greater): " << array[0][pointcol] << ", " << point << ", " << array[numrows-1][pointcol] << endl;
            exit(0);
        }
    }

    eps = 1e-10*abs(array[0][pointcol]-array[1][pointcol]);
    if (point < array[0][pointcol]) {       // check to see if point is in bounds
        if ( abs(point-array[0][pointcol]) < eps) {
            return array[0][pointcol];
        } else {
            cout << "Interpolation error: point is not in bounds (smaller): " << array[0][pointcol] << ", " << point << ", " << array[numrows-1][pointcol] << endl;
            exit(0);
        }
    }

    // Set up bisection for the upper and lower limits of the point
    int mid = floor(numrows/2);
    double sol, slope, yint;
    int midpoint = 0;
    int upperbound = -1;
    int lowerbound = -1;
    int n;
    
    // base case for an array of 2 or 3 rows
    if (mid == 1) {
        slope = (array[numrows-1][ycol] - array[0][ycol])/(array[numrows-1][pointcol] - array[0][pointcol]);
        yint = array[numrows-1][ycol] - slope*array[numrows-1][pointcol];
        sol = slope*point + yint;
    }
    
    // Bisection to find a bounds
    n = mid;
    while(mid != 0) {
        midpoint = n - 1;   // rownumber to row location in array
        mid = floor(mid/2);
        //cout << mid << " " << array[midpoint][pointcol] << endl; //(For Debug only)
        if (point >= array[midpoint][pointcol]) {
            n = n + mid;
        } else {
            n = n - mid;
        }
    }
    
    // identify the bound as upper or lower
    if (array[midpoint][pointcol] >= point || midpoint == (numrows - 1)) {
        upperbound = midpoint;
        lowerbound = midpoint -1;
        while (array[lowerbound][pointcol] > point) {
            if (lowerbound - 1 == 0) {
                lowerbound = 0;
                upperbound = 1;
                break;
            } else {
                lowerbound = lowerbound - 1;
                upperbound = upperbound - 1;
            }
        }
    } else if (array[midpoint][pointcol] < point || midpoint == 0) {
        lowerbound = midpoint;
        upperbound = midpoint + 1;
        while (array[upperbound][pointcol] < point) {
            if (upperbound + 1 == numrows - 1) {
                lowerbound = numrows - 2;
                upperbound = numrows - 1;
                break;
            } else {
                lowerbound = lowerbound + 1;
                upperbound = upperbound + 1;
            }
        }
    }
    
    
    
    if (upperbound == -1 || lowerbound == -1 || lowerbound > upperbound) {
        cout << "Interpolation error unknown" << endl;
        exit(0);
    }
    
    // Interpolate
    slope = (array[upperbound][ycol] - array[lowerbound][ycol])/(array[upperbound][pointcol] - array[lowerbound][pointcol]);
    yint = array[upperbound][ycol] - slope*array[upperbound][pointcol];
    sol = slope*point + yint;

    return sol;
}


// Interpolate multiple values in an array
int data2 :: interpolate_multi(int numrows, int numcols, double** array, double point, int pointcol, double* &sol_array, bool sortflag) {
    int n_ycols = numcols-1;

    if (numrows <= 1) {                                                           // check to see if interpolation is possible
        cout << "Interpolation error: not enough points to interpolate" << endl;
        exit(0);
    }

    // Sort data if unsorted
    if (sortflag == false) {
        sortasc(array, pointcol, numrows, numcols);
    }

    double eps = 1e-10*abs(array[numrows-1][pointcol]-array[numrows-2][pointcol]);
    if (point > array[numrows-1][pointcol]) {       // check to see if point is in bounds
        if ( abs(point-array[numrows-1][pointcol]) < eps) {
            return array[numrows-1][pointcol];
        } else {
            cout << "Multi Interpolation error: point is not in bounds (greater): " << array[0][pointcol] << ", " << point << ", " << array[numrows-1][pointcol] << endl;
            exit(0);
        }
    }

    eps = 1e-10*abs(array[0][pointcol]-array[1][pointcol]);
    if (point < array[0][pointcol]) {       // check to see if point is in bounds
        if ( abs(point-array[0][pointcol]) < eps) {
            return array[0][pointcol];
        } else {
            cout << "Interpolation error: point is not in bounds (smaller): " << array[0][pointcol] << ", " << point << ", " << array[numrows-1][pointcol] << endl;
            exit(0);
        }
    }

    // Set up bisection for the upper and lower limits of the point
    int mid = floor(numrows/2);
    double* slope; double* yint;
    sol_array = new double[n_ycols];
    slope = new double[n_ycols];
    yint = new double[n_ycols];
    int midpoint = 0;
    int upperbound = -1;
    int lowerbound = -1;
    int n;
    
    // base case for an array of 2 or 3 rows
    if (mid == 1) {
        for (int i=1; i<numcols; ++i) {
            slope[i-1] = (array[numrows-1][i] - array[0][i])/(array[numrows-1][pointcol] - array[0][pointcol]);
            yint[i-1] = array[numrows-1][i] - slope[i-1]*array[numrows-1][pointcol];
            sol_array[i-1] = slope[i-1]*point + yint[i-1];
        }
    }
    
    // Bisection to find a bounds
    n = mid;
    while(mid != 0) {
        midpoint = n - 1;   // rownumber to row location in array
        mid = floor(mid/2);
        //cout << mid << " " << array[midpoint][pointcol] << endl; //(For Debug only)
        if (point >= array[midpoint][pointcol]) {
            n = n + mid;
        } else {
            n = n - mid;
        }
    }
    
    // identify the bound as upper or lower
    if (array[midpoint][pointcol] >= point || midpoint == (numrows - 1)) {
        upperbound = midpoint;
        lowerbound = midpoint -1;
        while (array[lowerbound][pointcol] > point) {
            if (lowerbound - 1 == 0) {
                lowerbound = 0;
                upperbound = 1;
                break;
            } else {
                lowerbound = lowerbound - 1;
                upperbound = upperbound - 1;
            }
        }
    } else if (array[midpoint][pointcol] < point || midpoint == 0) {
        lowerbound = midpoint;
        upperbound = midpoint + 1;
        while (array[upperbound][pointcol] < point) {
            if (upperbound + 1 == numrows - 1) {
                lowerbound = numrows - 2;
                upperbound = numrows - 1;
                break;
            } else {
                lowerbound = lowerbound + 1;
                upperbound = upperbound + 1;
            }
        }
    }
    
    if (upperbound == -1 || lowerbound == -1 || lowerbound > upperbound) {
        cout << "Interpolation error unknown" << endl;
        exit(0);
    }
    
    for (int i=1; i<numcols; ++i) {
        // Interpolate
        slope[i-1] = (array[upperbound][i] - array[lowerbound][i])/(array[upperbound][pointcol] - array[lowerbound][pointcol]);
        yint[i-1] = array[upperbound][i] - slope[i-1]*array[upperbound][pointcol];
        sol_array[i-1] = slope[i-1]*point + yint[i-1];
    }
    delete slope;
    delete yint;
    return 0;
}


void data2 :: cubicspline(double**arr, double**&spline, int nrows, int xcol, int ycol) {
    int n = nrows-2;
    double a[n-1];
    double b[n];
    double c[n-1];
    double d[n];
    double cp[n-1];
    double dp[n];

    spline = new double*[nrows];
    for (int i = 0; i<nrows; i++) {
        spline[i] = new double[3];
    }

    for (int j=0; j<(n-1); ++j) {
        a[j] = 1.0/6.0*(arr[j+2][xcol] - arr[j+1][xcol]);
        b[j] = 1.0/3.0*(arr[j+2][xcol] - arr[j][xcol]);
        c[j] = 1.0/6.0*(arr[j+2][xcol] - arr[j+1][xcol]);
        d[j] = (arr[j+2][ycol] - arr[j+1][ycol])/(arr[j+2][xcol] - arr[j+1][xcol]) - (arr[j+1][ycol] - arr[j][ycol])/(arr[j+1][xcol] - arr[j][xcol]);
    }
    b[n-1] = 1.0/3.0*(arr[n+1][xcol] - arr[n-1][xcol]);
    d[n-1] = (arr[n+1][ycol] - arr[n][ycol])/(arr[n+1][xcol] - arr[n][xcol]) - (arr[n][ycol] - arr[n-1][ycol])/(arr[n][xcol] - arr[n-1][xcol]);

    
    cp[0] = c[0]/b[0];
    dp[0] = d[0]/b[0];
    for (int j=1; j<(n-1); ++j) {
        cp[j] = c[j]/(b[j] - a[j-1]*cp[j-1]);
        dp[j] = (d[j] - a[j-1]*dp[j-1])/(b[j]-a[j-1]*cp[j-1]);
    }
    dp[n-1] = (d[n-1] - a[n-2]*dp[n-2])/(b[n-1]-a[n-2]*cp[n-2]);

    spline[nrows-1][2] = 0.0;
    spline[n][2] = dp[n-1];

    for (int i=n-1; i>0; --i) {
        spline[i][2] = dp[i-1] - cp[i-1]*spline[i+1][2];
    }

    spline[0][2] = 0.0;

    for (int i=0; i<nrows; ++i) {
        spline[i][0] = arr[i][xcol];
        spline[i][1] = arr[i][ycol];
    }
}

double data2 :: splinecalc(double** spline, int x1, int x2, double x) {
    double a,b,c,d,y;
    
    a = (spline[x2][0] - x)/(spline[x2][0] - spline[x1][0]);
    b = 1.0-a;
    c = 1.0/6.0*(pow(a,3.0)-a)*pow(spline[x2][0]-spline[x1][0],2.0);
    d = 1.0/6.0*(pow(b,3.0)-b)*pow(spline[x2][0]-spline[x1][0],2.0);
    y = a*spline[x1][1] + b*spline[x2][1] + c*spline[x1][2] + d*spline[x2][2];
    //cout << spline[x1][0] << "  " << spline[x2][0] << "  " << a << "  " << b << "  " << c << "  " << d << "  " << y << endl;

    return y;
}

double data2 :: transpose_file(string file) {
    double** temp_array;
    double** temp_array_t;
    int nrows = rowcount(file);
    int ncols = colcount(file);
    importdata(file,temp_array);
    create(temp_array_t,ncols,nrows);

    
    for(int i=0; i<nrows; ++i) {
        for (int j=0; j<ncols; ++j) {
            temp_array_t[j][i] = temp_array[i][j];
        }
    }
    cleanup(temp_array,nrows);
    print(temp_array_t,ncols,nrows,true,file);
    cleanup(temp_array_t,ncols);
    return 0;
}

//----------------------------------------------------------------------------------------

// Bubble sort ascending order
bool data2 :: sortasc(double ** array, int col, int numrows, int numcols) {
    
    if (col >= numcols || col < 0) {
        cout << "Sorting error: column is out of bounds" << endl;           // Throw error if col is out of bounds
        cout << "number of cols are: " << numcols << " and index was: " << col << endl;           // Throw error if col is out of bounds
        exit(0);
    }
    
    // set up a temp array to store data
    double temp[numcols];

    int count;
    int numruns = 0;
    bool sorted = false;                                                    // Assume unsorted
    while (sorted == false) {                                               // Run while unsorted
        count = 0;                                                          // Count how many changes were made
        for (int i=0; i<(numrows-1); ++i) {
            for (int j=0; j<numcols; ++j) {                                 // store current row in temp array
                temp[j] = array[i+1][j];
            }
            if (array[i][col] > array[i+1][col]) {                          // compare two columns to see if they are in desc order
                //cout << i << endl;
                for (int j=0; j<numcols; ++j) {                             // Swap them if they are in desc order
                    array[i+1][j] = array[i][j];
                    array[i][j] = temp[j];
                }
                count = count + 1;                                          // add to count
            }
        }
        if (count == 0) {                                                   // set sorted flag to true if no more swaps are necessary
            sorted = true;
        } else {
            numruns = numruns + 1;                                          // number of swaps overall counter
        }
    }
    //cout << "Number of sort runs in col "<< col << ": " << numruns << endl;
    if (numruns == 0) {
        return true;                                                        // return true if array was sorted without any runs
    } else {
        return false;                                                       // return false if sorting was necessary
    }
}

//----------------------------------------------------------------------------------------

// cleanup pointer array
void data2 :: cleanup(double** &array, int nrows) {
    for (int i = 0; i<nrows; i++) {
        delete[] array[i];
    }
    delete[] array;
}

// cleanup pointer array
void data2 :: cleanup_string(string** &array, int nrows) {
    for (int i = 0; i<nrows; i++) {
        delete[] array[i];
    }
    delete[] array;
}

void data2 :: cleanup3d(double*** &array, int nrows, int ncols) {
    for (int i = 0; i<nrows; i++) {
        for (int j=0; j<ncols; j++) {
            delete[] array[i][j];
        }
        delete[] array[i];
    }
    delete[] array;
}

//----------------------------------------------------------------------------------------

// Find maximum of data array
double data2 :: findmax(double** array, int col, int numrows, int numcols) {

    if (col >= numcols || col < 0) {                                        // Throw error if col is out of bounds
        cout << "Find Max error: column is out of bounds" << endl;
        exit(0);
    }

    double max = array[0][col];                                             // find the max of an array by brute force
    for (int i=0; i<numrows; ++i) {
        if (array[i][col] > max) {
            max = array[i][col];
        }
    }
    return max;
}

double data2 :: findmax_vec(vector<vector<double>> array, int col, int numrows, int numcols) {

    if (col >= numcols || col < 0) {                                        // Throw error if col is out of bounds
        cout << "Find Max error: column is out of bounds" << endl;
        exit(0);
    }

    double max = array[0][col];                                             // find the max of an array by brute force
    for (int i=0; i<numrows; ++i) {
        if (array[i][col] > max) {
            max = array[i][col];
        }
    }
    return max;
}

//----------------------------------------------------------------------------------------

// Find minimum of data array
double data2 :: findmin(double** array, int col, int numrows, int numcols) {

    if (col >= numcols || col < 0) {                                    // Throw error if col is out of bounds
        cout << "Find min error: column is out of bounds" << endl;
        exit(0);
    }

    double min = array[0][col];                                         // Find the minimum by brute force
    for (int i=0; i<numrows; ++i) {
        if (array[i][col] < min) {
            min = array[i][col];
        }
    }
    return min;
}

// Find a given point in an array
int data2 :: findvalue(double **array, int nrows, int ncols, double point, int pcol, double tol) {

    if (pcol >= ncols || pcol < 0) {                                  // Throw error if column is out of bounds
        cout << "Find value error: column is out of bounds" << endl;
        exit(0);
    }

    double error = abs(point-array[0][pcol]);                           // set initial error
    double newerror;                                                 // error tolerance of point
    int irow = 0;
    for (int i=0; i<nrows; ++i) {
        newerror = abs(point-array[i][pcol]);                           // calculate how new error
        if (newerror < error) {                                         // smallest error is point within tolerance
            error = newerror;                                           // set the old error as new error
            irow = i;
        }
    }

    if (abs(point - array[irow][pcol])/point > tol) {
        cout << "Value " << point << " not within tolerance: " << endl; // throw error if the value is not close enough within tolerance
        cout << "Closest value is: " << array[irow][pcol] << endl;      // show closest value found regardless
        dm.print(array,nrows,ncols,true,"findvalue.txt");
        exit(0);
    } else {
        return irow;
    }
}

void data2 :: zero(double** array, int nrows, int ncols) {
    for (int i=0; i<nrows; ++i) {
        for (int j=0; j<ncols; ++j) {
            array[i][j] = 0.0;
        }
    }
}

void data2 :: create(double** &array, int nrows, int ncols) {
    array = new double*[nrows];
    for (int i=0; i<nrows; i++) {
        array[i] = new double[ncols];
    }
}

void data2 :: create3d(double*** &array, int nrows, int ncols, int nstacks) {
    array = new double**[nrows];
    for (int i=0; i<nrows; i++) {
        array[i] = new double*[ncols];
        for (int j=0; j<ncols; ++j) {
            array[i][j] = new double[nstacks];
        }
    }
}

void data2 :: invert_order(double** &array, int nrows, int ncols) {
    double** temp_array;
    temp_array = new double*[nrows];
    for (int i=0; i<nrows; i++) {
        temp_array[i] = new double[ncols];
    }
    
    for (int j=0; j<ncols; ++j) {
        for (int i=0; i<nrows; ++i) {
            temp_array[i][j] = array[nrows-i-1][j];
        }
    }

    for (int j=0; j<ncols; ++j) {
        for (int i=0; i<nrows; ++i) {
            array[i][j] = temp_array[i][j];
        }
    }

    cleanup(temp_array,nrows);
}

void data2 :: copy_pointer(double** pointer, double** &array, int nrows, int ncols) {
    for (int i=0; i<nrows; ++i) {
        for (int j=0; j<ncols; ++j) {
            array[i][j] = pointer[i][j];
        }
    }
}

/*
void data2 :: CompOSE_format(string eos) {
    double** eos_array;
    importdata(eos,eos_array);
    int nrows = rowcount(eos_array);

    ofstream out("CompOSE_files/eos.thermo");
    out << mN << "     " << mP << "     " << 1 << endl;
    for (int i=0; i<nrows; ++i) {
        out << 1 << "     " << i+1 << "     " << 
    }
}
*/

//*********************************************************************************
//*********************************************************************************

// Numerical Methods class
//*********************************************************************************
//*********************************************************************************

void nummeth :: twopointD(double **oriarr, double ** &newarr, int ncols, int nrows, int ycol, int xcol) {
    double dydx1, dydx2;
    
    // create new array with an additonal column for the derivative
    newarr = new double*[nrows];
    for (int i = 0; i<nrows; i++) {
        newarr[i] = new double[ncols+1];
    }
    
    // fill the new array with the same data
    for (int j=0; j<ncols; ++j) {
        for (int i=0; i<nrows; ++i) {
            newarr[i][j] = oriarr[i][j];
        }
    }
    
    newarr[0][ncols] = (newarr[1][ycol] - newarr[0][ycol])/(newarr[1][xcol] - newarr[0][xcol]);
    newarr[1][ncols] = (newarr[2][ycol] - newarr[0][ycol])/(newarr[2][xcol] - newarr[0][xcol]);
    
    for (int i=2; i<(nrows-2); ++i) {
        dydx1 = (newarr[i+1][ycol] - newarr[i-1][ycol])/(newarr[i+1][xcol] - newarr[i-1][xcol]);
        dydx2 = (newarr[i+2][ycol] - newarr[i-2][ycol])/(newarr[i+2][xcol] - newarr[i-2][xcol]);
        newarr[i][ncols] = 0.5*(dydx1 + dydx2);
    }
    
    newarr[nrows-2][ncols] = (newarr[nrows-1][ycol] - newarr[nrows-3][ycol])/(newarr[nrows-1][xcol] - newarr[nrows-3][xcol]);
    newarr[nrows-1][ncols] = (newarr[nrows-1][ycol] - newarr[nrows-2][ycol])/(newarr[nrows-1][xcol] - newarr[nrows-2][xcol]);
}

// Find intersection point of 2 arrays
double nummeth :: findintersect(double **arr1, int nrows1, int ncols1, int xcol1, int ycol1, double **arr2, int nrows2, int ncols2, int xcol2, int ycol2) {
    double cand;
    double m1, b1, m2, b2;

    if (xcol1 >= ncols1 || xcol1 < 0) {                                     // Throw error if columns are out of bounds
        cout << "Find Intersect error: column is out of bounds" << endl;
        exit(0);
    }

    if (xcol2 >= ncols2 || xcol2 < 0) {                                     // Throw error if columns are out of bounds
        cout << "Find Intersect error: column is out of bounds" << endl;
        exit(0);
    }

    // Running linear spline interpolations
    for (int i=0; i<(nrows1-1); ++i) {
        for (int j=0; j<(nrows2-1); ++j) {
            m1 = (arr1[i+1][ycol1] - arr1[i][ycol1])/(arr1[i+1][xcol1] - arr1[i][xcol1]);
            m2 = (arr2[j+1][ycol2] - arr2[j][ycol2])/(arr2[j+1][xcol2] - arr2[j][xcol2]);
            b1 = arr1[i+1][ycol1] - m1*arr1[i+1][xcol1];
            b2 = arr2[j+1][ycol2] - m2*arr2[j+1][xcol2];
            cand = (b2 - b1)/(m1 - m2); // candidate intersection point of two lines
            // check if the candidate point is in the spline range
            if (cand >= arr1[i][xcol1] && cand <= arr1[i+1][xcol1] && cand >= arr2[j][xcol2] && cand <= arr2[j+1][xcol2]) {
                return cand;
            }
        }
    }
    cout << "No intersection point found " << endl;
    exit(0);
}

double nummeth :: zbrent(double (*func)(double, double[10], double, double[2]), double x1, double x2, double tol, double kf, double couplings[10], double fields_v[2]) {
    int iter;
    int ITMAX = 100;
    double EPS = 3.0e-8;
    double a=x1, b=x2, c=x2, d,min1,min2;
    double e = 0.0;
    double fa = (*func)(kf,couplings,a,fields_v), fb = (*func)(kf,couplings,b,fields_v), fc,p,q,r,s,tol1,xm;

    if ( (fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
        cout << "Root must be bracketed in zbrent" << endl;
        for (int i=0; i<10; ++i) {
            cout << couplings[i] << "  ";
        }
        cout << endl;
        exit(0);
    }
    fc = fb;
    for (iter=1; iter<=ITMAX; iter++) {
        if ( (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c=a;
            fc=fa;
            e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
        xm = 0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0) return b;
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            s = fb/fa;
            if (a == c) {
                p = 2.0*xm*s;
                q = 1.0-s;
            } else {
                q = fa/fc;
                r = fb/fc;
                p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q = (q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) q = -q;
            p = fabs(p);
            min1 = 3.0*xm*q - fabs(tol1*q);
            min2 = fabs(e*q);
            if (2.0*p < (min1<min2 ? min1:min2)) {
                e=d;
                d=p/q;
            } else {
                d=xm;
                e=d;
            }
        } else {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if (fabs(d) > tol1) {
            b +=d;
        } else {
            b += SIGN(tol1,xm);
        }
        fb=(*func)(kf,couplings,b,fields_v);
    }
    cout << "Maximum number of iterations exceeded in zbrent" << endl;
    return 0.0;
}

void nummeth :: pretovconv(double** &eos, int encol, int prcol, double enconv, double prconv, int nrows) {
    for (int i=0; i<nrows; ++i) {
        eos[i][encol] = eos[i][encol]*enconv;
        eos[i][prcol] = eos[i][prcol]*prconv;
    }
}

// dm/dr TOV equation (unitless)
double nummeth :: dmdr(double r, double p, double m, double en) {
    double eq;
    eq = r*r*en; // Dimensionless
    return eq;
}

// dp/dr TOV equation (unitless)
double nummeth :: dpdr(double r, double p, double m, double en) {
    double eq;
    eq = -0.5*(en+p)*(m+pow(r,3)*p)/(r*(r-m)); //Dimensionless
    return eq;
}

double nummeth :: dvdr(double r, double m, double p) {
    double eq;
    eq = (m+pow(r,3)*p)/(r*(r-m));
    return eq;
}

double nummeth :: dYdr(double r, double Y, double en, double p, double v, double m, double wbar) {
    double eq, jp, j;

    jp = -(en+p)*exp(-v)*r;
    j = exp(-v)*(1.0-m/r);
    eq = -2.0*jp*wbar/(r*j) -(0.5*jp/j + 4.0/r)*Y;
    return eq;
}

double nummeth :: dwdr(double r, double Y) {
    double eq;
    eq = Y;
    return eq;
}

double nummeth :: Rdnldr(double r, double m, double dpde, double en, double p, int l, double nl) {
    double eq, Fr, Qr;
    Fr = (r + 0.5*pow(r,3)*(p-en) )/(r-m);
    Qr = r/(r-m) * ( 0.5*(5.0*en + 9.0*p + (en+p)/dpde) - l*(l+1)/pow(r,2) ) - pow( (m+pow(r,3)*p)/(pow(r,2)-m*r) , 2);
    eq = -pow(nl,2)/r - Fr*nl/r - r*Qr;
    return eq;
}

// Single star icp tov solve
void nummeth :: tovsolve(double icp, double h, double** eos, int nrows, int ncols, int encol, int prcol, int dpdecol, double MR[2], bool print) {
    
    if (prcol >= ncols || prcol < 0 || encol >= ncols || encol < 0) {   // throw error for out of bounds columns
        cout << "TOV error: column is out of bounds" << endl;
        exit(0);
    }

    double r = 1e-20;       // set initial radius close to zero
    double p = icp;         // set central pressure
    double m = 0;           // inital mass to zero
    double en,dpde;
    double k1,k2,k3,k4;
    double l1,l2,l3,l4;
    double pmin = eos[0][prcol];

    if (print == true) {
        ofstream out("SingleICP.txt");
        while (p>pmin) {
            en = dm.interpolate(nrows, ncols, eos, p, prcol, encol, true);  // Interpolate the nuclear/quark eos
            dpde = dm.interpolate(nrows, ncols, eos, p, prcol, dpdecol, true);
            out << scientific << setprecision(15) << en << "  " << p << "  " << dpde << "  " << r << "  " << m << endl;
            // Integration routine (Possible to optimize stepsize)
            k1 = dmdr(r,p,m,en);
            l1 = dpdr(r,p,m,en);
            k2 = dmdr(r+h/2.0,p+l1*h/2.0,m+k1*h/2.0,en);
            l2 = dpdr(r+h/2.0,p+l1*h/2.0,m+k1*h/2.0,en);
            k3 = dmdr(r+h/2.0,p+l2*h/2.0,m+k2*h/2.0,en);
            l3 = dpdr(r+h/2.0,p+l2*h/2.0,m+k2*h/2.0,en);
            k4 = dmdr(r+h,p+l3*h,m+k3*h,en);
            l4 = dpdr(r+h,p+l3*h,m+k3*h,en);
            m = m + 1.0/6.0*h*(k1 + 2.0*k2 + 2.0*k3 + k4);
            p = p + 1.0/6.0*h*(l1 + 2.0*l2 + 2.0*l3 + l4);
            r = r + h;
        }
    } else {
        while (p>pmin) {
            en = dm.interpolate(nrows, ncols, eos, p, prcol, encol, true);  // Interpolate the nuclear/quark eos
            dpde = dm.interpolate(nrows, ncols, eos, p, prcol, dpdecol, true);
            // Integration routine (Possible to optimize stepsize)
            
            k1 = dmdr(r,p,m,en);
            l1 = dpdr(r,p,m,en);
            k2 = dmdr(r+h/2.0,p+l1*h/2.0,m+k1*h/2.0,en);
            l2 = dpdr(r+h/2.0,p+l1*h/2.0,m+k1*h/2.0,en);
            k3 = dmdr(r+h/2.0,p+l2*h/2.0,m+k2*h/2.0,en);
            l3 = dpdr(r+h/2.0,p+l2*h/2.0,m+k2*h/2.0,en);
            k4 = dmdr(r+h,p+l3*h,m+k3*h,en);
            l4 = dpdr(r+h,p+l3*h,m+k3*h,en);

            m = m + 1.0/6.0*h*(k1 + 2.0*k2 + 2.0*k3 + k4);
            p = p + 1.0/6.0*h*(l1 + 2.0*l2 + 2.0*l3 + l4);
            r = r + h;
        }
    }
    
    MR[0] = m;   // save the final star mass
    MR[1] = r;   // save final star radius
    //cout << "icp: " << icp << "  M: " << m << "  R: " << conv.rnonetokm(r) << endl;
}

// Multiple central pressure star configurations
vector<vector<double>> nummeth :: multitov(double h, double** eos, int nrows, int ncols, int encol, int prcol, int dpdecol, int npoints, string filesave, double pr0) {
    double en, icp;
    int space;
    
    if (prcol >= ncols || prcol < 0 || encol >= ncols || encol < 0) { // throw error for out of bounds columns
        cout << "TOV error: column is out of bounds" << endl;
        exit(0);
    }

    if (npoints > nrows) {
        cout << " too many points: calculating each pr point" << endl;
        npoints = nrows;
    }
    
    double MR[2];                                            // initalize the mass radius array
    ofstream out(filesave + "_RM.txt");                         // File name save
    int index = dm.findvalue(eos,nrows,ncols,pr0,prcol,1.0);
    vector<vector<double>> TOV_data;

    space = ceil(1.0*(nrows-(index+1.0))/npoints);
    double dens;
    while (index < nrows) {
        icp = eos[index][prcol];
        tovsolve(icp, h, eos, nrows, ncols, encol, prcol, dpdecol, MR, false);                 // solve tov equations
        en = dm.interpolate(nrows, ncols, eos, icp, prcol, encol, true);
        dens = dm.interpolate(nrows,ncols,eos,icp,prcol,0,true);
        out << dens << "  " << en*conv.energyCONV(1,0) << "  " << icp*conv.energyCONV(1,0) << "  " << conv.rnonetokm(MR[1]) << "  " << MR[0] << "  " << log((en*conv.energyCONV(1,0) + icp*conv.energyCONV(1,0))/(939.0*dens)) << "  " << 32887.0*sqrt(MR[0]/pow(conv.rnonetokm(MR[1]),3.0))*0.9 << endl;
        //cout << to_string(index) + "/" + to_string(nrows) << "  " << "mue: " << dm.interpolate(nrows, ncols, eos, icp, prcol, 0, true) << " icp: " << icp*conv.energyCONV(1,0) << " ice: " << en*conv.energyCONV(1,0) << " RM: " << conv.rnonetokm(MR[0][1]) << " " << MR[0][0] << endl;
        index = index + space;  //spacing of each mass radius point in pressure space
        //cout << "LOVE: " << Nlove("SingleICP.txt", 3, 2, 0, MR[0][1], 1e-3) << endl;
        TOV_data.push_back({dens, en*conv.energyCONV(1,0), icp*conv.energyCONV(1,0), conv.rnonetokm(MR[1]), MR[0]});
    }
    return TOV_data;
}

// relativistic love number, inertia, and quadrupole moment
void nummeth :: RloveMCMC(double** EOS, int dpdecol, int encol, int prcol, double h, double ILQR[4], double icp, int nrows, int ncols) {
    double r0 = 1e-5;   // Set the integration starting point (close to zero)
    double en, cs2, r, v, wbar,Y,vaR;
    double I = 0;
    int l = 2;          // Set the order to 2
    double nl = l*1.0;      // boundary condition for clairaut-radau equation       
    double M, R;
    double r1,r2,r3,r4;
    double p1,p2,p3,p4;
    double m1,m2,m3,m4;
    double b1,b2,b3,b4;
    double a1,a2,a3,a4;
    double l1,l2,l3,l4;
    double pmin = EOS[0][prcol];
    double p = icp;         // set central pressure
    double m = 0;           // inital mass to zero
    
    r = r0;             //integrate from center
    v = 0;
    while (p>pmin) {
        en = dm.interpolate(nrows, ncols, EOS, p, prcol, encol, true);  // Interpolate the nuclear/quark eos
        a1 = dmdr(r,p,m,en);
        l1 = dpdr(r,p,m,en);
        m1 = dvdr(r,m,p);
        a2 = dmdr(r+h/2.0,p+l1*h/2.0,m+a1*h/2.0,en);
        l2 = dpdr(r+h/2.0,p+l1*h/2.0,m+a1*h/2.0,en);
        m2 = dvdr(r+h/2.0,m+a1*h/2.0,p+l1*h/2.0);
        a3 = dmdr(r+h/2.0,p+l2*h/2.0,m+a2*h/2.0,en);
        l3 = dpdr(r+h/2.0,p+l2*h/2.0,m+a2*h/2.0,en);
        m3 = dvdr(r+h/2.0,m+a2*h/2.0,p+l2*h/2.0);
        a4 = dmdr(r+h,p+l3*h,m+a3*h,en);
        l4 = dpdr(r+h,p+l3*h,m+a3*h,en);
        m4 = dvdr(r+h,m+a3*h,p+l3*h);
        m = m + 1.0/6.0*h*(a1 + 2.0*a2 + 2.0*a3 + a4);
        p = p + 1.0/6.0*h*(l1 + 2.0*l2 + 2.0*l3 + l4);
        v = v + 1.0/6.0*h*(m1 + 2.0*m2 + 2.0*m3 + m4);
        r = r+h;
    }
    M = m;
    R = r;
    vaR = v;

    r = r0;
    v = log(1.0-M/R) - vaR;
    wbar = r0;
    Y = r0;
    p = icp;
    m = 0;
    ofstream out("SingleLove.txt");
    while (p>pmin) {
        en = dm.interpolate(nrows, ncols, EOS, p, prcol, encol, true);  // Interpolate the nuclear/quark eos
        cs2  = dm.interpolate(nrows, ncols, EOS, p, prcol, dpdecol, true);   //speed of sound
        
        a1 = dmdr(r,p,m,en);
        b1 = dpdr(r,p,m,en);
        r1 = Rdnldr(r,m,cs2,en,p,l,nl);
        l1 = dwdr(r,Y);
        p1 = dYdr(r,Y,en,p,v,m,wbar);
        m1 = dvdr(r,m,p);

        a2 = dmdr(r+h/2.0, p+b1*h/2.0, m+a1*h/2.0, en);
        b2 = dpdr(r+h/2.0, p+b1*h/2.0, m+a1*h/2.0, en);
        r2 = Rdnldr(r+h/2.0, m+a1*h/2.0, cs2, en, p+b1*h/2.0, l, nl+r1*h/2.0);
        l2 = dwdr(r+h/2.0, Y+p1*h/2.0);
        p2 = dYdr(r+h/2.0, Y+p1*h/2.0,en, p+b1*h/2.0, v+m1*h/2.0, m+a1*h/2.0, wbar+l1*h/2.0);
        m2 = dvdr(r+h/2.0, m+a1*h/2.0, p+b1*h/2.0);

        a3 = dmdr(r+h/2.0, p+b2*h/2.0, m+a2*h/2.0,en);
        b3 = dpdr(r+h/2.0, p+b2*h/2.0, m+a2*h/2.0,en);
        r3 = Rdnldr(r+h/2.0, m+a2*h/2.0, cs2, en, p+b2*h/2.0, l, nl+r2*h/2.0);
        l3 = dwdr(r+h/2.0, Y+p2*h/2.0);
        p3 = dYdr(r+h/2.0, Y+p2*h/2.0,en, p+b2*h/2.0, v+m2*h/2.0, m+a2*h/2.0, wbar+l2*h/2.0);
        m3 = dvdr(r+h/2.0, m+a2*h/2.0, p+b2*h/2.0);

        a4 = dmdr(r+h, p+b3*h, m+a3*h,en);
        b4 = dpdr(r+h, p+b3*h, m+a3*h,en);
        r4 = Rdnldr(r+h, m+a3*h, cs2,en, p+b3*h, l, nl+r3*h);
        l4 = dwdr(r+h, Y+p3*h);
        p4 = dYdr(r+h, Y+p3*h,en, p+b3*h, v+m3*h, m+a3*h, wbar+l3*h);
        m4 = dvdr(r+h, m+a3*h, p+b3*h);

        double C = M/(2.0*R);   // compactness
        double d1 = 2.0*C*(6.0-3.0*nl+15.0*C*nl-24.0*C);
        double d2 = 4.0*pow(C,3)*(13.0-11.0*nl+3.0*C*nl-2.0*C+2.0*pow(C,2)*(nl+1.0));
        double d3 = 3.0*pow(1.0-2.0*C,2)*(2.0-nl+2.0*C*nl-2.0*C)*log(1.0-2.0*C);
        double k2 = 8.0*pow(C,5)/5.0*pow(1.0-2.0*C,2)*(2.0-2.0*C+2.0*nl*C-nl)/(d1+d2+d3);
        
        out << r*2.95324203 << "  " << m << "  " << cs2 << "  " << nl << "  " << k2 << "  " << 2.0/3.0*k2*pow(C,-5) << endl;

        m = m + 1.0/6.0*h*(a1 + 2.0*a2 + 2.0*a3 + a4); 
        p = p + 1.0/6.0*h*(b1 + 2.0*b2 + 2.0*b3 + b4); 
        nl = nl + 1.0/6.0*h*(r1 + 2.0*r2 + 2.0*r3 + r4); // clairaut-radau equation
        wbar = wbar + 1.0/6.0*h*(l1 + 2.0*l2 + 2.0*l3 + l4); 
        Y = Y + 1.0/6.0*h*(p1 + 2.0*p2 + 2.0*p3 + p4);
        v = v + 1.0/6.0*h*(m1 + 2.0*m2 + 2.0*m3 + m4);
        if (cs2 < 0) {
            cout << "bad: " << cs2 << endl;
        }
        r = r+h;
    }
    I = Y*pow(r,4)/(3.0*wbar + Y*r);

    double C = M/(2.0*R);   // compactness
    double d1 = 2.0*C*(6.0-3.0*nl+15.0*C*nl-24.0*C);
    double d2 = 4.0*pow(C,3)*(13.0-11.0*nl+3.0*C*nl-2.0*C+2.0*pow(C,2)*(nl+1.0));
    double d3 = 3.0*pow(1.0-2.0*C,2)*(2.0-nl+2.0*C*nl-2.0*C)*log(1.0-2.0*C);
    double k2 = 8.0*pow(C,5)/5.0*pow(1.0-2.0*C,2)*(2.0-2.0*C+2.0*nl*C-nl)/(d1+d2+d3);
    ILQR[0] = 4.0*I/pow(m,3);         // DIMENSIONLESS I
    ILQR[1] = 2.0/3.0*k2*pow(C,-5);       // Tidal Deformability
    ILQR[2] = 0;
    ILQR[3] = R;
    //cout << "nl: " << nl << endl;
    cout << "mass: " << M << endl;
    //cout << "Compactness: " << C << endl;
    cout << "k2: " << k2 << "  I: " << I << endl;
    cout << "TDeform: " << 2.0/3.0*k2*pow(C,-5) << endl;
    cout << "Radius: " << R*2.95324203 << endl;
}

double nummeth :: Urca_threshold(double** eos, int ncols, int nrows, int Yp_col, int Yp_Urca_col, int dens_col) {
    for (int i=0; i<nrows; ++i) {
        if (eos[i][Yp_col] > eos[i][Yp_Urca_col]) {
            return eos[i][dens_col];
        }
    }
    cout << "Error: No Urca process allowed" << endl;
    return 0;
}

double nummeth :: det(double** A, int n) {
    double det = 1.0;
    for (int i = 0; i < n; i++) {
        int pivot = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(A[j][i]) > abs(A[pivot][i])) {
                pivot = j;
            }
        }
        if (pivot != i) {
            swap(A[i], A[pivot]);
            det *= -1;
        }
        if (A[i][i] == 0) {
            return 0;
        }
        det *= A[i][i];
        for (int j = i + 1; j < n; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i + 1; k < n; k++) {
                A[j][k] -= factor * A[i][k];
            }
        }
    }
    return det;
}
 

bool nummeth :: invertMatrix(double** matrix, double** inverse, int n) {
    // Create an augmented matrix
    double** augmented = new double*[n];
    for (int i = 0; i < n; ++i) {
        augmented[i] = new double[2*n];
        for (int j = 0; j < n; ++j) {
            augmented[i][j] = matrix[i][j];
        }
        for (int j = n; j < 2*n; ++j) {
            augmented[i][j] = (i == j-n) ? 1.0 : 0.0;
        }
    }

    // Perform Gaussian elimination
    for (int i = 0; i < n; ++i) {
        // Find the pivot element
        double maxEl = abs(augmented[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(augmented[k][i]) > maxEl) {
                maxEl = abs(augmented[k][i]);
                maxRow = k;
            }
        }

        // Swap rows
        for (int k = 0; k < 2*n; ++k) {
            double tmp = augmented[maxRow][k];
            augmented[maxRow][k] = augmented[i][k];
            augmented[i][k] = tmp;
        }

        // Make the pivot element equal to 1
        double pivot = augmented[i][i];
        if (pivot == 0) {
            return false; // The matrix is singular and cannot be inverted
        }

        for (int j = 0; j < 2*n; ++j) {
            augmented[i][j] /= pivot;
        }

        // Make all elements in the current column equal to 0 except the pivot
        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = augmented[k][i];
                for (int j = 0; j < 2*n; ++j) {
                    augmented[k][j] -= factor * augmented[i][j];
                }
            }
        }
    }

    // Extract the inverse matrix from the augmented matrix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inverse[i][j] = augmented[i][j+n];
        }
    }

    // Clean up memory
    for (int i = 0; i < n; ++i) {
        delete[] augmented[i];
    }
    delete[] augmented;

    return true;
}