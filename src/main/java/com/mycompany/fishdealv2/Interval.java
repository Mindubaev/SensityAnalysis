/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.fishdealv2;

/**
 *
 * @author Artur
 */
public class Interval {
    
    private double leftBorder;
    private double rightBorder;
    public boolean rbIsDef=true;

    public Interval() {
    }

    public Interval(double leftBorder, double rightBorder) {
        this.leftBorder = leftBorder;
        this.rightBorder = rightBorder;
    }


    public Interval(double leftBorder) {
        this.leftBorder=leftBorder;
        this.rbIsDef=false;
    }
    

    @Override
    public String toString() {
        String rightBound="";
        if (rbIsDef)
            rightBound=String.valueOf(rightBorder);
        else
            rightBound="+infinity";
        return "["+this.getLeftBorder()+";"+rightBound+"]";
    }    

    public double getLeftBorder() {
        return leftBorder;
    }

    public void setLeftBorder(double leftBorder) {
        this.leftBorder = leftBorder;
    }

    public double getRightBorder() {
        return rightBorder;
    }

    public void setRightBorder(double rightBorder) {
        this.rightBorder = rightBorder;
    }
    
}
