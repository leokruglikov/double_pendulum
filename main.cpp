#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#define PI 3.14159



double f1(double t,double angle ,double angular_velocity1);
double f2(double t,double angle ,double angular_velocity2);


double g1(double t, double angle_1, double angle_2, double angular_velocity_1, double angular_velocity_2, double m1, double m2, double L1, double L2, double g_acceleration);
double g2(double t, double angle_1, double angle_2, double angular_velocity_1, double angular_velocity_2, double m1, double m2, double L1, double L2, double g_acceleration);




int main(){
std::cout<<"Programm starting"<<std::endl;
    //setting the lenght and the mass
double m1=1.55;
double m2=0.45;
double L1=0.35;
double L2=0.35;
    //the angle that is changing over time
double angle1;
double angle2;
double omega1;
double omega2;
    //the variables necessary for 4th order runge kutta method for grater presition
    //to solve a system of differential equations
double k10;
double k11;
double k12;
double k13;
double l10;
double l11;
double l12;
double l13;
double k20;
double k21;
double k22;
double k23;
double l20;
double l21;
double l22;
double l23;
double t=0;
double tau=20;
double dt=0.01;
double g =9.81;
    double x1,y1,x2,y2;
    //setting the initial conditions
std::vector<double> angular_velocity_1(tau/dt); angular_velocity_1[0]=0; //(y)'' = function1
                                                //(y)'= z=g()
                                                //z'= function of z,y
std::vector<double> angular_velocity_2(tau/dt);      angular_velocity_2[0]=0; //set initial conditions for velocity
std::vector<double> theta_1(angular_velocity_1.size()); theta_1[0]=PI/2;   //set initial conditions for theta
std::vector<double> theta_2(angular_velocity_2.size()); theta_2[0]=PI/4;   //set initial conditions for

    //MAIN LOOP TO CALCULATE THE RUNGE-KUTTA COEFFICIENTS

for(int i = 1; i< angular_velocity_1.size(); i++){
    
    //First order coefficients for R-K. method
    l10=dt*g1(t,theta_1[i-1], theta_2[i-1], angular_velocity_1[i-1],angular_velocity_2[i-1],m1,m2,L1,L2,g);
    k10=dt*f1(t,theta_1[i-1],angular_velocity_1[i-1]);
    l20=dt*g2(t,theta_1[i-1], theta_2[i-1], angular_velocity_1[i-1], angular_velocity_2[i-1], m1, m2, L1, L2, g);
    k20=dt*f2(t,theta_2[i-1], angular_velocity_2[i-1]);
    
    //Secont order coefficients for R-K. method
    l11=dt*g1(t+0.5*dt, theta_1[i-1]+0.5*k10, theta_2[i-1]+0.5*k20, angular_velocity_1[i-1]+0.5*l10, angular_velocity_2[i-1]+0.5*l20,m1,m2,L1,L2,g);
    k11=dt*f1(t+0.5*dt, theta_1[i-1]+0.5*k10, angular_velocity_1[i-1]+0.5*l10);
    k21=dt*f2(t+0.5*dt, theta_2[i-1]+0.5*k20, angular_velocity_2[i-1]+0.5*l20);
    l21=dt*g2(t+0.5*dt, theta_1[i-1]+0.5*k10, theta_2[i-1]+0.5*k20, angular_velocity_1[i-1]+0.5*l10, angular_velocity_2[i-1]+0.5*l20,m1,m2,L1, L2,g );
    k12=dt*f1(t+0.5*dt, theta_1[i-1]+0.5*k11, angular_velocity_1[i-1]+0.5*l11);
    l12=dt*g1(t+0.5*dt, theta_1[i-1]+0.5*k11, theta_2[i-1]+0.5*k21, angular_velocity_1[i-1]+0.5*l11, angular_velocity_2[i-1]+0.5*l21,m1, m2,L1, L2, g );
    
    //Third order coefficients for R-K. method

    k22=dt*f2(t+0.5*dt, theta_2[i-1]+0.5*k21, angular_velocity_2[i-1]+0.5*l21);
    l22=dt*g2(t+0.5*dt, theta_1[i-1]+0.5*k11, theta_1[i-1]+0.5*k21, angular_velocity_1[i-1]+0.5*l11, angular_velocity_2[i-1]+0.5*l21, m1, m2, L1, L2,g);
    k13=dt*f1(t+dt,theta_1[i-1]+k12, angular_velocity_1[i-1]+l12);
    l13=dt*g1(t+dt,theta_1[i-1]+k12, theta_2[i-1]+k22, angular_velocity_1[i-1]+l12, angular_velocity_2[i-1]+l22, m1,m2,L1,L2,g);
    k23=dt*f2(t+dt, theta_2[i-1]+k22, angular_velocity_2[i-1]+l22);
    l23=dt*g2(t+dt, theta_1[i-1]+k12, theta_2[i-1]+k22,angular_velocity_1[i-1]+l12, angular_velocity_2[i-1]+l22, m1,m2,L1,L2,g);
    
    //calculating the new state variables from
    theta_1[i]=theta_1[i-1]+(k10+2*k11+2*k12+k13)/6;
    angular_velocity_1[i]=angular_velocity_1[i-1]+(l10+2*l11+2*l12+l13)/6;
    theta_2[i]=theta_2[i-1]+(k20+2*k21+2*k22+k23)/6;
    angular_velocity_2[i]=angular_velocity_2[i-1]+(l20+2*l21+2*l22+l23)/6;
    
    t=t+dt;
    x1=L1*sin(theta_1[i]);
    x2=x1+L2*sin(theta_2[i]);
    y1= -L1*cos(theta_1[i]);
    y2= y1-L2*cos(theta_2[i]);
        

}
    //veariables for the state evolution depending on time
    std::vector<double> X1(angular_velocity_1.size());
    std::vector<double> Y1(angular_velocity_1.size());
    std::vector<double> X2(angular_velocity_1.size());
    std::vector<double> Y2(angular_velocity_1.size());
    for (int i =0; i<angular_velocity_1.size(); i++) {
        X1[i]=L1*sin(theta_1[i]);
        Y1[i]=-L1*cos(theta_1[i]);
        X2[i]=X1[i]+L2*sin(theta_2[i]);
        Y2[i]=Y1[i]-L2*cos(theta_2[i]);
    }
    
    //putting all the data in file
    
std::ofstream file;
file.open("data.txt");
if(!file.is_open()){
    std::cout<<"The file wasn't opened "<<std::endl;
}
for(int i =0; i<theta_1.size(); i++){
    file<<X1[i]<<"\t"<<Y1[i]<<"\t"<<X2[i]<<"\t"<<Y2[i]<<"\t"<<0<<"\n";
}

return 0;

}



double f1(double t,double angle ,double angular_velocity1){ //several parameters are taken, but only the angular velocity is taken into account
   
    return angular_velocity1;
}
double f2(double t,double angle ,double angular_velocity2){ //several parameters are taken, but only the angular velocity is taken into account
    return angular_velocity2;
}
double g1(double t, double angle_1, double angle_2, double angular_velocity_1, double angular_velocity_2, double m1, double m2, double L1, double L2, double g_acceleration){
    double temporary= (-g_acceleration*(2*m1 + m2)*sin(angle_1)-(m2*g_acceleration*sin(angle_1-2*angle_2))- (2*sin(angle_1-angle_2)*m2*(pow(angular_velocity_2, 2)*L2+
                       pow(angular_velocity_1, 2)*L1*cos(angle_1-angle_2)) ))/(L1*(2*m1+m2-m2*cos(2*angle_1-2*angle_2)));
return temporary;

}
double g2(double t, double angle_1, double angle_2, double angular_velocity_1, double angular_velocity_2, double m1, double m2, double L1, double L2, double g_acceleration){
    double temporary = (2*sin(angle_1-angle_2)*(pow(angular_velocity_1, 2)*L1*(m1+m2)+g_acceleration*(m1+m2)*cos(angle_1)+ pow(angular_velocity_2,2)*L2*m2*cos(angle_1-angle_2)))/(L2*(2*m1+m2-m2*cos(2*angle_1-2*angle_2)));

    return temporary;
}
