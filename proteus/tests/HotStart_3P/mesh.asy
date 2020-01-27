
unitsize(4.0 inches / 1.000000);
size(11 inches);
real Lx=1.000000;
real Ly=1.000000;
real offset=0.0125Lx;
real x=0.000000;
real y=0.000000;
string strx="$1.00\mbox{m}$";
string stry="$1.00\mbox{m}$";
draw(strx,(x,y-offset)--(x+Lx,y-offset),S,black,Bars,Arrows,PenMargins);
draw(stry,(x-offset,y)--(x-offset,y+Ly),W,black,Bars,Arrows,PenMargins);
import graph;
import palette;
pen[] allPens=Wheel();
pen[] myPens = new pen[3+1];
for(int i=0;i<3+1;++i)
  {
   int iPen = round(i*allPens.length/(3+1));
   myPens[i] = allPens[iPen];
  }
draw((0.000000,0.000000)--(1.000000,0.000000),myPens[2]+linewidth(0.01));
draw((1.000000,0.000000)--(1.000000,1.000000),myPens[1]+linewidth(0.01));
draw((1.000000,1.000000)--(0.000000,1.000000),myPens[3]+linewidth(0.01));
draw((0.000000,1.000000)--(0.000000,0.000000),myPens[0]+linewidth(0.01));
