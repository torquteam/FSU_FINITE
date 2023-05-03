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
data2 dm;
Convert conv;
using namespace std;

const double pi = 4*atan(1);
const double mE = 0.511;        // Mass of electron (MeV)
const double mU = 2.3;          // Mass of up quark (MeV)
const double mD = 4.8;          // Mass of down quark (MeV)
const double mS = 95.0;         // Mass of strange quark (MeV)
const double mP = 938.27231;    // Mass of proton (MeV)
const double mMU = 105.6583745; // Mass of muon (MeV)

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