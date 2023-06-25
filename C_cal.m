function [C_k, summation_ck] = C_cal(ph_k,ph_kprime,dis_factor,k,summation)

first=ph_k* transpose(ph_k-dis_factor*ph_kprime);
summation_ck=first+summation;
C_k=(1/(k+1))*(summation_ck);