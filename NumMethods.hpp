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
#include <vector>
using namespace std;

class data2 {
public:
    int rowcount(string txtfile);
    string filetype(string txtfile);
    int colcount(string txtfile);
    void importdata(string txtfile, double **&array);
    void importdata_string(string txtfile, string ** &array);
    void print(double **&array, int numrows, int numcols, bool file, string filename);
    double interpolate(int numrows, int numcols, double** array, double point, int pointcol, int ycol, bool sortflag);
    int interpolate_multi(int numrows, int numcols, double** array, double point, int pointcol, double* &y_array, bool sortflag);
    void cubicspline(double**arr, double**&spline, int nrows, int xcol, int ycol);
    double splinecalc(double** spline, int x1, int x2, double x);
    bool sortasc(double ** array, int col, int numrows, int numcols);
    void cleanup(double** &array, int nrows);
    void cleanup3d(double*** &array, int nrows, int ncols);
    void cleanup_string(string** &array, int nrows);
    double findmin(double** array, int col, int numrows, int numcols);
    double findmax(double** array, int col, int numrows, int numcols);
    int findvalue(double **array, int nrows, int ncols, double point, int pcol, double tol);
    void zero(double** array, int nrows, int ncols);
    void create(double** &array, int nrows, int ncols);
    void invert_order(double** &array, int nrows, int ncols);
    void create3d(double*** &array, int nrows, int ncols, int nstacks);
    void copy_pointer(double** pointer, double** &array, int nrows, int ncols);
    void convert_array(double** &array, int numrows, int numcols, int col_to_convert, double conv_factor);
    double transpose_file(string file);
    double findmax_vec(vector<vector<double>> array, int col, int numrows, int numcols);
};

class nummeth {
public:
    void twopointD(double **oriarr, double ** &newarr, int ncols, int nrows, int ycol, int xcol);
    double findintersect(double **arr1, int nrows1, int ncols1, int xcol1, int ycol1, double **arr2, int nrows2, int ncols2, int xcol2, int ycol2);
    double zbrent(double (*func)(double, double[8], double, double[2]), double x1, double x2, double tol, double kf, double couplings[8], double fields_v[2]);
    void pretovconv(double** &eos, int encol, int prcol, double enconv, double prconv, int nrows);
    double dmdr(double r, double p, double m, double en);
    double dpdr(double r, double p, double m, double en);
    double dvdr(double r, double m, double p);
    double dYdr(double r, double Y, double en, double p, double v, double m, double wbar);
    double dwdr(double r, double Y);
    double Rdnldr(double r, double m, double dpde, double en, double p, int l, double nl);
    void tovsolve(double icp, double h, double** eos, int nrows, int ncols, int encol, int prcol, int dpdecol, double MR[2], bool print);
    vector<vector<double>> multitov(double h, double** eos, int nrows, int ncols, int encol, int prcol, int dpdecol, int npoints, string filesave, double pr0);
    void RloveMCMC(double** EOS, int dpdecol, int encol, int prcol, double h, double ILQR[4], double icp, int nrows, int ncols);
    double Urca_threshold(double** eos, int ncols, int nrows, int Yp_col, int Yp_Urca_col, int dens_col);
    bool invertMatrix(double** A, double** inverse, int n);
    double det(double** A, int n);
};

#endif /* NumMethods_h */
