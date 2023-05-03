//
//  NumMethods.hpp
//  NS Code
//
//  Created by Marc Salinas on 1/29/19.
//  Copyright © 2019 Marc Salinas. All rights reserved.
//

#ifndef NumMethods_hpp
#define NumMethods_hpp

#include <stdio.h>
#include <string>
using namespace std;

class data2 {
public:
    int rowcount(string txtfile);
    string filetype(string txtfile);
    int colcount(string txtfile);
    void importdata(string txtfile, double **&array);
    void print(double **&array, int numrows, int numcols, bool file, string filename);
    double interpolate(int numrows, int numcols, double** array, double point, int pointcol, int ycol, bool sortflag);
    void cubicspline(double**arr, double**&spline, int nrows, int xcol, int ycol);
    double splinecalc(double** spline, int x1, int x2, double x);
    bool sortasc(double ** array, int col, int numrows, int numcols);
    void cleanup(double** &array, int nrows);
    void cleanup3d(double*** &array, int nrows, int ncols);
    double findmin(double** array, int col, int numrows, int numcols);
    double findmax(double** array, int col, int numrows, int numcols);
    int findvalue(double **array, int nrows, int ncols, double point, int pcol, double tol);
    void zero(double** array, int nrows, int ncols);
    void create(double** &array, int nrows, int ncols);
    void invert_order(double** &array, int nrows, int ncols);
    void create3d(double*** &array, int nrows, int ncols, int nstacks);
    void copy_pointer(double** pointer, double** &array, int nrows, int ncols);
};

class nummeth {
public:
    void twopointD(double **oriarr, double ** &newarr, int ncols, int nrows, int ycol, int xcol);
    double findintersect(double **arr1, int nrows1, int ncols1, int xcol1, int ycol1, double **arr2, int nrows2, int ncols2, int xcol2, int ycol2);
};

#endif /* NumMethods_h */
