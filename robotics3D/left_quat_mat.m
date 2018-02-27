function [L] = left_quat_mat(q1)

Q1_MATR(1,1) = q1(4,1) ;
Q1_MATR(1,2) = q1(3,1) ;
Q1_MATR(1,3) = -q1(2,1) ;
Q1_MATR(1,4) = q1(1,1) ;

Q1_MATR(2,1) = -q1(3,1) ;
Q1_MATR(2,2) = q1(4,1) ;
Q1_MATR(2,3) = q1(1,1) ;
Q1_MATR(2,4) = q1(2,1) ;

Q1_MATR(3,1) = q1(2,1) ;
Q1_MATR(3,2) = -q1(1,1) ;
Q1_MATR(3,3) = q1(4,1) ;
Q1_MATR(3,4) = q1(3,1) ;

Q1_MATR(4,1) = -q1(1,1) ;
Q1_MATR(4,2) = -q1(2,1) ;
Q1_MATR(4,3) = -q1(3,1) ;
Q1_MATR(4,4) = q1(4,1) ;


L = Q1_MATR;