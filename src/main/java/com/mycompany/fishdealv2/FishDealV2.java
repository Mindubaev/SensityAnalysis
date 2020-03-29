/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.fishdealv2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import lpsolve.LogListener;
import lpsolve.LpSolve;
import lpsolve.LpSolveException;
import org.la4j.LinearAlgebra;
import org.la4j.Matrix;
import org.la4j.Vector;
import org.la4j.inversion.MatrixInverter;
import org.la4j.matrix.dense.Basic2DMatrix;

/**
 *
 * @author Artur
 */
public class FishDealV2 {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws LpSolveException {
        solveLinealProg(new File("src/input.txt"));
    }
    
    public static void solveLinealProg(File file){
        try(BufferedReader reader=new BufferedReader(new FileReader(file))){
            Double[] c=getArrayFromString(reader.readLine());//цены
            Double[] a=getArrayFromString(reader.readLine());//ограниче по столбцам
            int k=c.length;
            int m=a.length; 
            Double[][] matrix=getMatrix(reader, m, k);
            LpSolve problem=getLpProblem(matrix, c, a);
            problem.setTrace(false);
            problem.setDebug(false);
            problem.setVerbose(0);
            //problem.setLowbo(2, 1);
            int status=problem.solve();
            System.err.println("status:"+status);
            problem.printLp();
            if (status==0){
                problem.printSolution(1);
                problem.printObjective();
                printSensitivityAnalysis(problem, k);
                printAnalyzeAllAvalableSolution(problem,k);
                printEquals(problem, k, false);
                problem.deleteLp();
            }
        }catch(Exception ex){
            ex.printStackTrace();
        }
    }
    
    
    
    public static LpSolve getLpProblem(Double[][] matrix,Double[] c,Double[] a) throws LpSolveException{
        double[] canonicalC=getCanonicalC(c, a.length*2);
        double[][] canonicalMatrix=getCanonicalMatrix(matrix, a.length,c.length);
        LpSolve problem=LpSolve.makeLp(0, canonicalC.length);
        problem.strSetObjFn(arrToString(canonicalC));
        for (int i=0;i<a.length;i++){
            String row=arrToString(canonicalMatrix[i]);
            problem.strAddConstraint(row, LpSolve.EQ, a[i]);
        }
        return problem;
    }
    
    public static String arrToString(double[] arr){
        String string = "";
        for (double d : arr) {
            string = string + d + " ";
        }
        return string.substring(0,string.length()-1);
    }
    
    public static double[][] getCanonicalMatrix(Double[][] matrix,int constNum,int varNum){
        int newRowSize=varNum+constNum*2;
        double[][] canonicalMatrix=new double [matrix.length][newRowSize];
        for (int i=0;i<matrix.length;i++){
            Double[] row=matrix[i];
            for (int j=0;j<newRowSize;j++){
                if (j<row.length)
                    canonicalMatrix[i][j]=row[j];
                else{
                    if (j==i+row.length)
                        canonicalMatrix[i][j]=-1;
                    else{
                        if (j==i+row.length+constNum)
                            canonicalMatrix[i][j]=1;
                        else
                            canonicalMatrix[i][j]=0;
                    }
                }
            }
        }
        return canonicalMatrix;
    }
    
    public static double[] getCanonicalC(Double[] c,int constNum){
        double[] canonicalC=new double[c.length+constNum];
        for (int i=0;i<canonicalC.length;i++){
            if (i<c.length)
                canonicalC[i]=c[i];
            else{
                if (i>c.length+constNum/2-1)
                    canonicalC[i]=99999999999.0;
                else
                    canonicalC[i]=0;
            }
        }
        return canonicalC;
    }
    
    public static Double[][] getMatrix(BufferedReader reader,int m,int k) throws IOException{
        Double[][] matrix=new Double[m][k];
        for (int j=0;j<m;j++){
            String[] input=reader.readLine().split(" ");
            for (int i=0;i<k;i++){
                matrix[j][i]=Double.parseDouble(input[i]);
            }
        }
         return matrix;
    }
    
    public static Double[] getArrayFromString(String line){
        String[] info=line.split(" ");
        List<Double> c=new ArrayList<>();
        for (int i=0;i<info.length;i++){
            c.add(Double.parseDouble(info[i]));
        }
        return c.toArray(new Double[c.size()]);
    }
    
    public static int[] getBasisIndexes(LpSolve problem) throws LpSolveException{
        int basisSize=problem.getNrows();
        int[] basisIndex=new int[basisSize+1];
        problem.getBasis(basisIndex, false);
        for (int i=1;i<basisIndex.length;i++){
            basisIndex[i]=basisIndex[i]*(-1)-basisSize;
        }
        return basisIndex;
    }
    
    public static Matrix getBasisMatrix(LpSolve problem) throws LpSolveException{
        int basisSize=problem.getNrows();
        int[] basisIndex=getBasisIndexes(problem);
        double[][] basisArr=new double[basisSize][basisSize];
        for (int i=0;i<basisIndex.length-1;i++){
            int index=basisIndex[i+1];
            double[] column=problem.getPtrColumn(index);
            for (int j=1;j<column.length;j++){
                basisArr[i][j-1]=column[j];
            }
        }
        return new Basic2DMatrix(basisArr);
    }
    
    public static Matrix getA(LpSolve problem) throws LpSolveException{
        int rowNum=problem.getNrows();
        int numOfVar=problem.getNcolumns();
        double[][] a=new double[rowNum][numOfVar];
        for (int i=1;i<=rowNum;i++){
            double[] row=problem.getPtrRow(i);
            for (int j=1;j<row.length;j++){
                a[i-1][j-1]=row[j];
            }
        }
        return new Basic2DMatrix(a);
    }
    
    public static Matrix getSimplexTable(LpSolve problem) throws LpSolveException{
        Matrix basisMatrix=getBasisMatrix(problem);
        basisMatrix=basisMatrix.transpose();
        MatrixInverter inverter=basisMatrix.withInverter(LinearAlgebra.InverterFactory.GAUSS_JORDAN);
        Matrix inverseBasis=inverter.inverse();
        Matrix simplex=inverseBasis.multiply(getA(problem));
        return simplex;
    }
    
    public static void printEquals(LpSolve problem,int numOfPrimalVar,boolean printAllVar) throws LpSolveException{
        System.out.println("Система неравенств коэффициентов целевой функции, где dn изменение n коэффициента: ");
        Matrix simplexTable=getSimplexTable(problem);
        double[] c=problem.getPtrRow(0);
        int[] basisIndexes=getBasisIndexes(problem);
        for (int i=0;i<simplexTable.columns();i++){
            Vector column=simplexTable.getColumn(i);
            if (!isBasisVector(column)){
                String equals="";
                double x0=0;
                for (int j=0;j<column.length();j++){
                    int indx=basisIndexes[j+1];
                    x0=x0+column.get(j)*c[indx];
                    if (indx<=numOfPrimalVar || printAllVar)
                        equals=equals+"("+column.get(j)+"*d"+indx+")+";
                }
                x0-=c[i+1];
                if (i+1<=numOfPrimalVar || printAllVar)
                    equals=equals+"("+x0+")-d"+(i+1)+"<=0";
                else
                    equals=equals+"("+x0+")<=0";
                System.out.println(equals);
            }
        }
    }
    
    public static void printSensitivityAnalysis(LpSolve problem,int valNum) throws LpSolveException{
        double[] cUpperBounds=new double[problem.getNcolumns()];
        double[] cLowwerBounds=new double[problem.getNcolumns()];
        problem.getSensitivityObj(cUpperBounds, cLowwerBounds);
        System.out.println("Интервалы стоимости кормов, в рамках которых решение остаётся оптимальным");
        String str="";
        for (int i=0;i<cUpperBounds.length;i++){
            if (i<valNum)
                str=str+"C"+(i+1)+"=["+cUpperBounds[i]+";"+cLowwerBounds[i]+"] ";
        }
        System.out.println(str);
    }
    
    public static boolean isBasisVector(Vector v){
        boolean find=false;
        for (double x:v){
            if ((x!=0 && x!=1) || (x==1 && find))
                return false;
            if (x==1)
                find=true;
        }
        return true;
    }

    public static void printArr(double[] arr) {
        String string = "";
        for (double d : arr) {
            string = string + d + " ";
        }
        System.out.println(string);
    }
    
    public static void printArr(int[] arr) {
        String string = "";
        for (double d : arr) {
            string = string + d + " ";
        }
        System.out.println(string);
    }
    
    public static void printAnalyzeAllAvalableSolution(LpSolve problem,int varNum) throws LpSolveException{
        Interval[] intervals=anylyzeAvalableSolution(problem, 100, 10);
        System.out.println("Интервалы стоимости кормов, в рамках которых существует решение");
            for (int i=0;i<varNum;i++)
                System.out.println("с"+(i+1)+":"+intervals[i]);
    }
    
    public static Interval[] anylyzeAvalableSolution(LpSolve problem,int maxNumOfIterations,double delta) throws LpSolveException{//res-найденное оптимаоьное решение
        double[] defaultK=new double[problem.getNonzeros()];
        problem.getRow(0,defaultK);
        Interval[] kIntervals=new Interval[problem.getNcolumns()];
        for (int i=1;i<=problem.getNcolumns();i++){
            int counter=0;
            double deltarb=delta;
            double deltalb=delta;
            double[] rbk=Arrays.copyOf(defaultK, defaultK.length);
            double[] lbk=Arrays.copyOf(defaultK, defaultK.length);
            boolean findrb=false;
            while (counter<maxNumOfIterations){
                rbk[i]=rbk[i]+deltarb;
                problem.setObjFn(rbk);
                problem.solve();
                double[] newRes=new double[problem.getNonzeros()];
                problem.getVariables(newRes);
                if (newRes==null){
                    findrb=true;
                    rbk[i]=rbk[i]-deltarb;
                    deltarb=deltarb/2;
                }
                lbk[i]=lbk[i]-deltalb;
                problem.setObjFn(lbk);
                problem.solve();
                newRes=new double[problem.getNonzeros()];
                problem.getVariables(newRes);
                if (lbk[i]<0 || newRes==null){
                    lbk[i]=lbk[i]+deltalb;
                    deltalb=deltalb/2;
                }
                counter++;
            }
            if (findrb)
                kIntervals[i-1]=new Interval(lbk[i], rbk[i]);
            else
                kIntervals[i-1]=new Interval(lbk[i]);
        }
        problem.setObjFn(defaultK);
        return kIntervals;
    } 

}
