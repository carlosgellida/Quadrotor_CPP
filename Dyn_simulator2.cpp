/* This is an efficient program to calculate the dynamic
    of non inertial multibody trees 
    made by: Carlos Gellida   */

#include <istream>
#include <fstream>

#include <iostream>
#include <cmath>
#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor-blas/xlinalg.hpp"
#include "xtensor/xfixed.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xcomplex.hpp"
#include <time.h> 

#include "xtensor/xcsv.hpp"
#include "xtensor/xnpy.hpp"

using namespace std; 
using namespace xt; 
using namespace xt::linalg ; 

/* Declaration of Global variables */ 

const size_t DoF = 13 ; // Input the number of Degree of freedom 
const size_t Nb = 10 ; //Number of bodies (including the ones that are rigidly attached)
const size_t Nframes = 8; //Number of frames (excluding inertial)
const size_t Nj = Nframes - 1 ; // Number of joints
const int pi = 3.14159265 ;

xtensor<float, 1> dot_q({int(DoF)}) ; // Generalized velocities
xtensor<float, 1> q_part({int(Nj) + 3}) ; // Initial joint coordinates
xtensor<float, 2> Chi_base({6, 6}) ; // Pose of base body

xtensor_fixed<double, xshape<Nframes>> lambda = {0, 0, 0, 0, 0, 5, 6} ;  // Parent vector 

xtensor_fixed<float, xshape<Nj, 6>> ExtendedAxisInfo(){
  /* Here the rows of S_0 are the extended axis of 
      action of each joint */
  xtensor_fixed<float, xshape<Nj, 6>> S_0 ; 
  S_0 = zeros<float> ({int(Nj), 6})  ;
  S_0(0, 5) = 1;
  S_0(1, 5) = 1;
  S_0(2, 5) = 1;
  S_0(3, 5) = 1;
  S_0(4, 5) = 1;
  S_0(5, 3) = 1;
  S_0(6, 3) = 1;
  return S_0; 
}

xtensor_fixed<float, xshape<Nj, 6>> InitialPositions(){
  xtensor_fixed<float, xshape<Nj, 6>> P_0; 
  P_0 = {{0.25, -0.25, 0.1},
         {0.25, 0.25, 0.1},
         {-0.25, 0.25, 0.1},
         {-0.25, -0.25, 0.1},
         {0, 0, -0.05},
         //{0, 0, -0.05},
         //{0, float(0.2*sin(0.7854)), -float(0.2*cos(0.7854))}} ; 
         {0, 0, -0.1},
         {0, float(0.2*sin(0.7854)), -float(0.2*cos(0.7854))}} ;
  return P_0 ; 
}

xtensor_fixed<float, xshape<Nb, 3, 3>> InertiaInfoI(){
/* This function provides InfoI, which contains one 3x3 matriz per each

    body wich contains the following: The first row means the inertia over the three main inertia axis, second
    row means the traslation to the body frame (obtained in agreement with the me-
    thodology), third row is the Axis-angle parameters of rotation to the body 
    frame, and  */

  xtensor_fixed<float, xshape<Nb, 3, 3>> Info1 = zeros<float> ({10, 3, 3}) ; 
  view(Info1, 0, 0, all()) = xarray<float> ({0.0008125, 0.0008125, 0.001125}); 

  view(Info1, 1, 0, all()) = xarray<float> ({0.66667e-6, 0.00444417, 0.00444417}) ; 
  view(Info1, 1, 2, all()) = xarray<float> ({0, 0, 0.7854}) ; 

  view(Info1, 2, 0, all()) = xarray<float> ({0.66667e-6, 0.00444417, 0.00444417}) ; 
  view(Info1, 2, 2, all()) = xarray<float> ({0, 0, -0.7854}) ; 

  view(Info1, 3, 0, all()) = xarray<float> ({1.36561e-05, 1.36771e-05, 2.13783e-08}) ; 

  view(Info1, 4, 0, all()) = xarray<float> ({1.36561e-05, 1.36771e-05, 2.13783e-08}) ; 

  view(Info1, 5, 0, all()) = xarray<float> ({1.36561e-05, 1.36771e-05, 2.13783e-08}) ; 

  view(Info1, 6, 0, all()) = xarray<float> ({1.36561e-05, 1.36771e-05, 2.13783e-08}) ; 

  view(Info1, 7, 0, all()) = xarray<float> ({8.41667e-05, 8.41667e-05, 1.66667e-06}) ;  
  view(Info1, 7, 1, all()) = xarray<float> ({0, 0, 0.05}) ; 

  view(Info1, 8, 0, all()) = xarray<float> ({6.68333e-04, 6.68333e-04, 3.33333e-06}) ;  
  //view(Info1, 8, 1, all()) = xarray<float> ({0, 0, 0.1}) ; 
  view(Info1, 8, 1, all()) = xarray<float> ({0, -float(0.1*sin(0.7854)), float(0.1*cos(0.7854))}) ; 
  view(Info1, 8, 2, all()) = xarray<float> ({0.7854, 0, 0}) ;  

  view(Info1, 9, 0, all()) = xarray<float> ({6.68333e-04, 6.68333e-04, 3.33333e-06}) ;  
  //view(Info1, 9, 1, all()) = xarray<float> ({0, 0, 0.1}) ; 
  view(Info1, 9, 1, all()) = xarray<float> ({0, -0.1, 0}) ; 
  view(Info1, 9, 2, all()) = xarray<float> ({1.5708, 0, 0}) ;  
  return Info1 ; 
}

xtensor_fixed<float, xshape<2, Nb>> InertiaInfoII() {
  /* This function provides Info2, a matrix with size 2xNb, which contains the information 
      of which bodies are rigidly joined between them, and their masses */

  xtensor_fixed<float, xshape<2, Nb>> Info2 = {{0, 0, 0, 1, 2, 3, 4, 5, 6, 7},
                   {0.3, 0.1, 0.1, 2.54e-3, 2.54e-3, 2.54e-3, 2.54e-3, 0.1, 0.2, 0.2}}; 

  return Info2 ; 
}

xtensor_fixed<float, xshape<3, 3>> skew3(xtensor_fixed<float, xshape<3, 1>> a) {
  /* This functions calculate the skew-symmetric matrix of a 3x1 array  */
  xtensor_fixed<float, xshape<3, 3>> skew_a = zeros<float> ({3, 3}); 
  skew_a(0, 1) = -a(2, 0) ; 
  skew_a(1, 0) = a(2, 0) ; 
  skew_a(0, 2) = a(1, 0) ; 
  skew_a(2, 0) = -a(1, 0) ; 
  skew_a(1, 2) = -a(0, 0) ; 
  skew_a(2, 1) = a(0, 0) ; 
  return skew_a ; 
}

xarray<float> vector_norm(xarray<float> vector){
  /* It gives the norm two of column arrays */
  vector = vector*vector ; 
  return sqrt(sum(vector)) ; 
}

xtensor_fixed<size_t, xshape<2>> shp(xarray<float> A){
  // It gives the shape of a matrix 
  xtensor_fixed<size_t, xshape<2>> shape ; 
  shape(0) = col(A, 0).size() ; 
  shape(1) = A.size()/shape(0) ; 
  return shape ; 
}
  

xtensor_fixed<float, xshape<Nframes, 6, 6>> EIM(xarray<float> Info1, xarray<float> Info2){
  //  This function generates the Extended Inertia Matrix (EIM) 
  xtensor_fixed<float, xshape<6, 6>> M_int = zeros<float> ({6, 6}) ; 
  xtensor_fixed<float, xshape<Nframes, 6, 6>> M = zeros<float> ({int(Nframes), 6, 6}) ;
  xtensor_fixed<float, xshape<Nb, 3, 3>> I = zeros<float> ({int(Nb), 3, 3}) ; 
  xtensor_fixed<float, xshape<Nb, 3, 3>> R = zeros<float> ({int(Nb), 3, 3}) ;
  xtensor_fixed<float, xshape<Nb, 3>> Pos_cm = zeros<float> ({int(Nb), 3}) ; //Position of center of mass in local frame
  xarray<float> a ; 
  xtensor_fixed<float, xshape<3, 3>> K; // Axis of rotation
  float theta, mass; // Angle of rotation 
  int j = -1 ; 
  
  view(I, all(), 0, 0) = view(Info1, all(), 0, 0) ; 
  view(I, all(), 1, 1) = view(Info1, all(), 0, 1) ; 
  view(I, all(), 2, 2) = view(Info1, all(), 0, 2) ; 

  Pos_cm = view(Info1, all(), 1, all()) ; 

  
  for (int i = 0; i < Nb; i++){
    
    if (Info2(0, i) != Info2(0, i-1)) {
      j = j + 1 ; 
    }


    mass = Info2(1, i) ; 

    theta = vector_norm(view(Info1, i, 2, all() ))(0, 0) ; 
    
    if (allclose(view(Info1, i, 2, all() ), 0)){
      view(R, i, all(), all()) = eye(3) ; 
      K = zeros<float> ({3, 3}) ; 
    }
    else {
      a = (1/theta) * view(Info1, i, 2, all()) ; // Axis of rotation
      K = skew3(a) ; 
      view(R, i, all(), all()) = eye(3) + sin(theta)*K + (1 - cos(theta)) * dot(K, K); //Rodriguez Formula
    }

    view(I, i, all(), all()) = dot( view(R, i, all(), all()), dot(view(I, i, all(), all()), transpose(view(R, i, all(), all()) ) ) ) ; 

    
    view(M_int, range(0, 3), range(0, 3)) = mass *  eye(3) ; 
    view(M_int, range(0, 3), range(3, 6)) = mass *  skew3(view(Pos_cm, i, all())) ; 
    view(M_int, range(3, 6), range(0, 3)) = -mass *  skew3(view(Pos_cm, i, all())) ; 
    view(M_int, range(3, 6), range(3, 6)) = -mass * dot(skew3(view(Pos_cm, i, all())), skew3(view(Pos_cm, i, all()))) + view(I, i, all(), all()); 

    if (Info2(0, i) != Info2(0, i-1)) {
      view(M, j, all(), all()) = M_int; 
    } 
    else {
      view(M, j, all(), all()) = view(M, j, all(), all()) + M_int ;
    }
  } 
  return M ; 
}  

xtensor_fixed<float, xshape<6, 6>> Skew6_I(xtensor_fixed<float, xshape<6>> a){
  xtensor_fixed<float, xshape<6, 6>> Skew_A ; 
  view(Skew_A, range(0, 3), range(0, 3)) = skew3( view(a, range(3, 6))); 
  view(Skew_A, range(0, 3), range(3, 6)) = skew3( view(a, range(0, 3))); 
  view(Skew_A, range(3, 6), range(0, 3)) = zeros<float> ({3, 3}) ; 
  view(Skew_A, range(3, 6), range(3, 6)) = skew3( view(a, range(3, 6))); 
  return Skew_A ; 
}

xtensor_fixed<float, xshape<6, 6>> Skew6_II(xtensor_fixed<float, xshape<6>> a){
  xtensor_fixed<float, xshape<6, 6>> Skew_A ; 
  view(Skew_A, range(0, 3), range(0, 3)) = zeros<float> ({3, 3}) ; 
  view(Skew_A, range(0, 3), range(3, 6)) = skew3( view(a, range(0, 3))); 
  view(Skew_A, range(3, 6), range(0, 3)) = skew3( view(a, range(0, 3))); 
  view(Skew_A, range(3, 6), range(3, 6)) = skew3( view(a, range(3, 6))); 
  return Skew_A ; 
}

void Parent_matrix(xarray<float>& PM){
  // This function construc the parent matrix 
  for(int i = 0; i < Nj; i ++){
    if(lambda(i) == 0){
      view(PM, i, i) = 1 ; 
    } 
    else{
      view(PM, i, all()) = view(PM, (int(lambda(i))-1), all()) ; 
      view(PM, i, i) = 1 ;
    } 
  }
  return ;
}

void Exp_funct(xtensor_fixed<float, xshape<6, 6>>& Exp, xtensor_fixed<float, xshape<6, 6>> Skew_S, xtensor_fixed<float, xshape<6, 1>> S, float q){
  int joint ; 
  if(allclose(view(S, range(0, 3)), 0)){
    joint = 0 ; // Is a rotative joint
  }
  else {
    joint = 1 ; // Is a prismatic joint
    }
  Exp = eye(6) + joint * q * Skew_S + (1 - joint) * sin(q) * Skew_S +  (1 - joint) * (1 - cos(q)) * dot(Skew_S, Skew_S) ; 
}

void Rodriguez_form(xtensor<float, 2>& R, xtensor<float, 1> Axis, float theta){  
  xtensor<float, 2> skew_axis3({3, 3}) ; 
  skew_axis3 = skew3(Axis) ; 
  R = eye(3) + sin(theta) * skew_axis3 +  (1 - cos(theta)) * dot(skew_axis3, skew_axis3) ; 
}

xtensor_fixed<float, xshape<3, 3>> Rodriguez_form2(xtensor_fixed<float, xshape<3>> Axis, xtensor_fixed<float, xshape<3, 3>>& R){   //////////////////////
  xtensor<float, 2> skew_axis3({3, 3}) ; 
  float theta = vector_norm(Axis)(0) ; 
  skew_axis3 = skew3(Axis) ; 
  R = eye(3) + (1/theta)*sin(theta) * skew_axis3 +  (1 - cos(theta))/(theta*theta) * dot(skew_axis3, skew_axis3) ; 
  return R ; 
}

void Chi_make(xtensor<float, 2>& R, xtensor_fixed<float, xshape<3>> P, xtensor_fixed<float, xshape<6, 6>>& Chi){
  view(Chi, range(3, 6), range(0, 3)) = zeros<float> ({3, 3}) ; 
  view(Chi, range(0, 3), range(0, 3)) = R ; 
  view(Chi, range(0, 3), range(3, 6)) = dot(R, skew3(P)) ; 
  view(Chi, range(3, 6), range(3, 6)) = R ; 
}

void Chi_make2(xtensor_fixed<float, xshape<6>>& delta_q_part, xtensor_fixed<float, xshape<6, 6>>& Chi_int){
  xtensor_fixed<float, xshape<3>> Axis_theta = view(delta_q_part, range(3, 6)) ; 
  xtensor_fixed<float, xshape<3>> P = view(delta_q_part, range(0, 3)) ;
  xtensor_fixed<float, xshape<3, 3>> R ; 
  Rodriguez_form2(Axis_theta, R) ; 
  view(Chi_int, range(0, 3), range(0, 3)) = R ; 
  view(Chi_int, range(3, 6), range(3, 6)) = R ; 
  //view(Chi_int, range(0, 3), range(3, 6)) = dot(R, skew3(P)) ;
  view(Chi_int, range(0, 3), range(3, 6)) = dot(R, skew3(P)) ;  
  view(Chi_int, range(3, 6), range(0, 3)) = zeros<float> ({3, 3}) ; 
}

 xtensor_fixed<float, xshape<4>> quat_product(xtensor_fixed<float, xshape<4>>& quat1, xtensor_fixed<float, xshape<4>>& quat2){
  //This function calculates the special product of quaternion
  xtensor_fixed<float, xshape<4>> quat_res = zeros<float> ({4}) ; 
  quat_res(0) = quat1(0)*quat2(0) - dot( view(quat1, range(1, 4)), view(quat2, range(1, 4)))(0) ; 
  view(quat_res, range(1, 4)) = quat1(0)*view(quat2, range(1, 4)) + quat2(0)*view(quat1, range(1, 4)) + dot(skew3(view(quat1, range(1, 4))), view(quat2, range(1, 4)) ) ; 
  return quat_res ; 
} 
  
 xtensor_fixed<float, xshape<3, 3>> quattorot(xtensor_fixed<float, xshape<4>>& quat){
  // This function transform a quaternion into a rotation matrix
  xtensor_fixed<float, xshape<3, 3>> R ; 
  R = eye(3) - 2*quat(0)*skew3(view(quat, range(1, 4))) + 2*dot(skew3(view(quat, range(1, 4))), skew3(view(quat, range(1, 4)))) ; 
  return R ; 
}  

xtensor_fixed<float, xshape<4, 4>> skew4(xtensor_fixed<float, xshape<3>>& omega){
  xtensor_fixed<float, xshape<4, 4>> omega_skew = zeros<float> ({4, 4}) ; 
  
  omega_skew(0, 1) = - omega(0) ;
  omega_skew(0, 2) = - omega(1) ;
  omega_skew(0, 3) = - omega(2) ;
  
  omega_skew(1, 0) = omega(0) ;
  omega_skew(1, 2) = omega(2) ;
  omega_skew(1, 3) = -omega(1) ;
  
  omega_skew(2, 0) = omega(1) ; 
  omega_skew(2, 1) = -omega(2) ;
  omega_skew(2, 3) = omega(0) ;
  
  omega_skew(3, 0) = omega(2) ;
  omega_skew(3, 1) = omega(1) ;
  omega_skew(3, 2) = -omega(0) ; 
  
  return omega_skew ;
  }
  
void Rodriguez_form4(xtensor_fixed<float, xshape<4, 4>>& Exp4, xtensor_fixed<float, xshape<3>> Axis, float theta){  
  xtensor_fixed<float, xshape<4, 4>> skew_axis4 = skew4(Axis) ;  
  xtensor_fixed<float, xshape<4, 4>> R ;
  R = eye(4) + sin(theta) * skew_axis4 +  (1 - cos(theta)) * dot(skew_axis4, skew_axis4) ; 
}  

void  Normalized_quat(xtensor_fixed<float, xshape<4>>& quat){
  float norm = sqrt(dot(quat, quat)(0)) ; 
  quat = quat/norm ; 
} 
  
void inverse_dynamic( xtensor_fixed<float, xshape<6, 6>> Chi_base, xtensor<float, 1> dot_q, 
                      xtensor<float, 1> q_part, xtensor_fixed<float, xshape<Nframes, 6>> P, 
                      xtensor_fixed<float, xshape<Nj, 6>> S_0, xtensor_fixed<float, xshape<Nframes>> lambda, 
                      xtensor<float, 3> M,
                      xtensor<float, 2>& C, xtensor<float, 2>& H)
  {
    //This function perform the inverse dynamic
  
  xtensor<float, 1> V({6}) ; 
  xtensor<float, 1> V_i({6}) ;  //Twist local
  xtensor<float, 1> Mu({6}) ; 
  xtensor<float, 2> J_i({6, int(DoF)}) ;
  xtensor<float, 2> d_J_i({6, int(DoF)}) ;  
  xtensor<float, 2> Chi({6, 6}) ; 
  xtensor<float, 2> H_i({Nframes, Nframes}) ; 
  
  xtensor_fixed<float, xshape<6>> S ; 
  xtensor_fixed<float, xshape<6, 6>> Skew_S ; 
  xtensor_fixed<float, xshape<6, 6>> Exp, Chi_0; 
  xtensor_fixed<float, xshape<Nframes, 6, DoF>> J, d_J; 

  float d_q, q ; 
  int k ; 

  J_i = zeros<float> ({6, int(DoF)}) ; 
  d_J_i = zeros<float> ({6, int(DoF)}) ;
  
  J = zeros<float> ({int(Nframes), 6, int(DoF)}) ; 
  d_J = zeros<float> ({int(Nframes), int(DoF)}) ; 

  view(J_i, range(0, 6), range(0, 6)) = eye(6) ; 
  V = dot(J_i, dot_q) ;  // Twist

  Mu = dot(view(M, 0, all(), all()), V) ; // Momentum 

  view(H, range(0, 6), range(0, 6)) = view(M, 0, all(), all()) ; 
  view(C, range(0, 6), range(0, 6)) = -Skew6_II(Mu) ; 

  view(J, 0, all(), all()) = J_i ; 
  view(d_J, 0, all(), all()) = d_J_i ; 


  for(int i = 0; i < Nj ; i++){
    k = lambda(i) ; 
    d_q = dot_q(i + 6) ; 
    q = q_part(i + 3) ; 
    S = view(S_0, i, all()) ; 
    //S = view(S_0, all(), i) ; 

    Skew_S = Skew6_I(S); 
    Exp_funct(Exp, -Skew_S, view(S_0, 1, all()), q) ; // minus sign added
    Chi_0 = eye(6) - Skew6_I(view(P, i, all())) ; 
    Chi = dot(Exp, Chi_0) ;     

    J_i = view(J, k, all(), all()) ; 
    d_J_i = view(d_J, k, all(), all()) ; 

    J_i = dot(Chi, J_i) ; 
    view(J_i, all(), float(i + 6)) = S ; 

    V = dot(J_i, dot_q) ;

    V_i = d_q * S ; 

    Mu = dot(view(M, (i + 1), all(), all()), V) ; // Momentum 

    d_J_i = - dot(Skew6_I(V_i), dot(Chi, view(J, k, all(), all())) ) +  dot(Chi, view(d_J, k, all(), all()) ); 

    //cout << "J_i = " << J_i << endl << endl ;  
    //cout << "view(M, (i + 1), all(), all()) = " << view(M, (i + 1), all(), all()) << endl << endl ; 
      
    H_i = dot(transpose(J_i), dot(view(M, (i + 1), all(), all()), J_i)) ;

    H = H + H_i ; 
    
    C = C + dot(transpose(J_i), dot(view(M, (i + 1), all(), all()), d_J_i)) - dot(dot(transpose(J_i), Skew6_II(Mu)), J_i);
    /*cout << "Chi con i = " << (i+1) << endl << Chi << endl << endl ; 
    cout << "Exp con i = " << (i+1) << endl << Exp << endl << endl ; 
    cout << "q con i = " << (i+1) << " " << q << endl  ;  
    cout << "Ji con i = " << (i+1) << endl << J_i << endl << endl ;  */

    view(J, (i + 1), all(), all()) = J_i ; 
    view(d_J, (i + 1), all(), all()) = d_J_i ; 
  } 
  return ; 
}

//int main(int argc, char* argv[]){
int main(){

  xarray<float> PM = zeros<float> ({int(Nj), int(Nj)}) ; 

  xarray<float> Noc = zeros<float> ({int(Nj), 1}) ; // Number of columxshape<(Nj + 3)>ns of jacobian of father body

  Parent_matrix(PM) ; //Construction of parent matrix 

  Noc = sum(PM, 1) + 5 ;

  /* Introducing information about system configuration and inertia */

  xtensor_fixed<float, xshape<Nj, 6>>  S_0 = ExtendedAxisInfo() ; 

  xtensor_fixed<float, xshape<Nj, 3>>  P_0 = InitialPositions() ; 

  xtensor_fixed<float, xshape<Nb, 3, 3>> Info1 = InertiaInfoI() ; 
  xtensor_fixed<float, xshape<2, Nb>> Info2 = InertiaInfoII() ;

  xtensor_fixed<float, xshape<Nframes, 6, 6>> M = EIM(Info1, Info2) ; 
  

  xtensor_fixed<float, xshape<Nj, 6>> P = zeros<float> ({6, 1}) ; 

  view(P, all(), range(0, 3)) = P_0 ; 

  /* Introduction of the initial conditions of symulation */

  
  xtensor_fixed<float, xshape<DoF>> dot_q_0 = zeros<float> ({int(DoF)}) ; // Generalized velocities
  //dot_q_0(10) = 0.1 ; 

  xtensor_fixed<float, xshape<(Nj + 3)>> q_part_0 = zeros<float> ({(int(Nj) + 3)}) ; // Initial joint coordinates t = 0

  xtensor_fixed<float, xshape<3, 3>>  R_base_0 = eye<float> (3) ; // Orientation of base body w.r.t. inertial t = 0

  xtensor_fixed<float, xshape<1, 3>>  P_base_0 = ones<float> ({1, 3}) ; //Position of base body w.r.t. inertial t = 0

  xtensor_fixed<float, xshape<6, 6>> Chi_base_0 = zeros<float> ({6, 6}) ; 
  view(Chi_base_0, range(0, 3), range(0, 3)) = R_base_0 ; 
  view(Chi_base_0, range(0, 3), range(3, 6)) = dot(R_base_0, skew3(P_base_0) );
  view(Chi_base_0, range(3, 6), range(3, 6)) = R_base_0 ;
  view(Chi_base_0, range(3, 6), range(0, 3)) = zeros<float> ({3, 3}) ; 

  

  /* Begin the simulation */

  
  //float sim_time = 3e-3 ; // Simulation time (seconds)
  float sim_time = 10 ; // Simulation time (seconds)
  float step = 1e-3 ; // Time step
  int n_step = int(sim_time/step) ; // Number of time steps

  xtensor<float, 1> ddot_q({DoF}) ; // Time derived of dot_q
  xtensor<float, 1> C_vector({DoF}) ; // C * dot_q
  xtensor<float, 1> V({6}) ; // Twist 
  xtensor<float, 1> twist({6}) ; // Angular velocity
  xtensor<float, 1> dot_q_part({Nj + 3}) ; // Generalized velocities excluding base body angular velocity
  xtensor<float, 1> ddot_q_part({Nj + 3}) ; // Generalized aceleration excluding base body angular aceleration
  xtensor<float, 1> Axis({3}) ; // Euler parameters 
  xtensor<float, 1> Axis_unit({3}) ; // Instant axis of action 
  xtensor<float, 2> R({3, 3}) ; // Instant rotation matrix
  xtensor<float, 1> delta_angle({1}) ; // Instant angle of turning
  xtensor_fixed<float, xshape<3>> omega = zeros<float> ({3}) ;
  xtensor_fixed<float, xshape<6>> delta_q_base = zeros<float> ({6}) ;
  xtensor_fixed<float, xshape<6, 6>> Chi_int = zeros<float> ({6, 6}) ; 
  xtensor_fixed<float, xshape<Nj +3>> delta_q_part = zeros<float> ({int(Nj) + 3}) ; 
  xtensor_fixed<float, xshape<3>> delta_p = zeros<float> ({3}) ;

  R = ones<float> ({3, 3}) ; 

  dot_q_part = zeros<float> ({Nj + 3}) ; 

  dot_q = dot_q_0 ; 
  dot_q(10) = 0.1 ; 

  q_part = q_part_0 ; 
  //q_part(8) = pi/4 ;
  //q_part(9) = pi/4 ; 
  Chi_base = Chi_base_0 ; 

  xtensor<float, 2> C ({int(DoF), int(DoF)}) ; // Generalized Coriolis Matrix
  xtensor<float, 2> H ({int(DoF), int(DoF)}) ; // Generalized Inertia Matrix
 

  xarray<float> time_record = zeros<int> ({int(n_step)}) ; // Time vector for plotting
  xtensor<float, 2> q_part_record({size_t(n_step), size_t(Nj + 3) }) ; // Time vector for plotting

  /*
  cout << "Solution with Axis angle and exponential matrix" << endl ; 
  
  double time = double(clock()) ; 

  for (int h = 0; h < n_step; h++){
  
    C = zeros<float> ({int(DoF), int(DoF)}) ; 
    H = zeros<float> ({int(DoF), int(DoF)}) ; 
  
    inverse_dynamic(Chi_base, dot_q, q_part, P, S_0, lambda, M, C, H) ; // Evaluation of inverse dynamic
    
    C_vector = - dot(C, dot_q) ;  //Minus sign added
    ddot_q = solve(H, C_vector) ;   

    dot_q = dot_q + ddot_q*step ; // Euler integration
    
    omega = view(dot_q, range(3, 6)) ; // Angular velocity of base body

    view(dot_q_part, range(0, 3)) = view(dot_q, range(0, 3)) ; 
    view(dot_q_part, range(3, Nj + 3)) = view(dot_q, range(6, DoF)) ; 

    delta_q_part = dot_q_part*step ; 
    q_part = q_part + delta_q_part; // Euler integration
    
    Axis = omega*step ; 
    
    delta_q_base = hstack(xtuple(view(delta_q_part, range(0, 3)), Axis )); 
    
    delta_angle = vector_norm(Axis) ; // Instant angle of rotation
    Axis_unit = Axis/delta_angle(0) ; // Instant unitary axis of rotation import math import math  

    Rodriguez_form(R, Axis, delta_angle(0)) ; // The instant rotation matrix is generated
    
    delta_p = view(delta_q_part, range(0, 3)) ; 

    
    //Chi_make(R, delta_p, Chi_int) ; // The instant Chi matrix of base body
    
    Chi_make2(delta_q_base, Chi_int) ; 
    
    Chi_base = dot(Chi_base, Chi_int) ; // Current Chi matrix of base body    

    time_record(h) = h ; 
    view(q_part_record, h, all()) = q_part ; 
  }

  time = double(clock()) - time ; 
  double time_taken =  time/double(CLOCKS_PER_SEC);
  cout << "The process 1 takes: " << time_taken << "s" << endl ; 
  cout << q_part_record << endl ; 
  
  */ 
  
 
 // Solution using quaternion
  cout << "Solution with quaternion and simple sum" << endl ;
 
  xtensor_fixed<float, xshape<4>> omega_quat = zeros<float> ({4}); 
  xtensor_fixed<float, xshape<4>> dot_quat = zeros<float> ({4}); 
  xtensor_fixed<float, xshape<4>> quat = zeros<float> ({4}); 
  
  double time = double(clock()) ; 

  for (int h = 0; h < n_step; h++){
  
    C = zeros<float> ({int(DoF), int(DoF)}) ; 
    H = zeros<float> ({int(DoF), int(DoF)}) ; 
  
    inverse_dynamic(Chi_base, dot_q, q_part, P, S_0, lambda, M, C, H) ; // Evaluation of inverse dynamic
    
    C_vector = - dot(C, dot_q) ;  //Minus sign added
    ddot_q = solve(H, C_vector) ;   

    dot_q = dot_q + ddot_q*step ; // Euler integration
    
    
    omega = view(dot_q, range(3, 6)) ; // Angular velocity of base body
    

    view(dot_q_part, range(0, 3)) = view(dot_q, range(0, 3)) ; 
    view(dot_q_part, range(3, Nj + 3)) = view(dot_q, range(6, DoF)) ; 

    delta_q_part = dot_q_part*step ; 
    q_part = q_part + delta_q_part; // Euler integration
    
    
    view(omega_quat, range(1, 4)) = omega ; 
    
    dot_quat = 0.5 * quat_product(quat, omega_quat) ; 

    quat = quat + dot_quat*step ; 
    
    Normalized_quat(quat) ; 
    
    //delta_q_base = hstack(xtuple(view(delta_q_part, range(0, 3)), Axis )); 
    
    //delta_angle = vector_norm(Axis) ; // Instant angle of rotation
    //Axis_unit = Axis/delta_angle(0) ; // Instant unitary axis of rotation import math import math  

    //Rodriguez_form(R, Axis, delta_angle(0)) ; // The instant rotation matrix is generated
    
    R = quattorot(quat) ; 
    
    //delta_p = view(delta_q_part, range(0, 3)) ; 

    
    //Chi_make(R, delta_p, Chi_int) ; // The instant Chi matrix of base body
    
    //Chi_make2(delta_q_base, Chi_int) ; 
    
    view(Chi_base, range(0, 3), range(0, 3)) = R ; 
    view(Chi_base, range(3, 6), range(3, 6)) = R ; 
    view(Chi_base, range(0, 3), range(3, 6)) = dot(R, skew3(view(q_part, range(0, 3)) )) ; 
    view(Chi_base, range(3, 6), range(0, 3)) = zeros<float> ({3, 3}) ; 
    
    //Chi_base = dot(Chi_base, Chi_int) ; // Current Chi matrix of base body    

    time_record(h) = h ; 
    view(q_part_record, h, all()) = q_part ; 
  }

  time = double(clock()) - time ; 
  double time_taken =  time/double(CLOCKS_PER_SEC);
  cout << "The process 1 takes: " << time_taken << "s" << endl ; 
  cout << q_part_record << endl ; 
  
  
  
  /* Output to a file, to be plot in Matlab, Python or gnuplot */
  
  fstream my_file;
  my_file.open("q_part_record.txt", ios::out);

  xt::dump_csv(my_file, q_part_record);

  my_file.close();
  
  return 0 ; 
}

