clear all
% To plot the speed up
% hold on

% hold on
% s_50_20=[0.045861 0.02349 0.016438 0.012432 0.0105948 0.0085279 0.00853896 0.0066769 0.00589514 0.00580716];
% s_100_100=[0.515961 0.247477 0.213734 0.200553 0.198133 0.161565 0.171224 0.0941899 0.0555649 0.052906];
% s_200_200=[3.11253 1.32781 1.04006 0.943346 0.975984 0.895256 0.898627 0.897025 0.840289 0.84508];
% 
% serial50=0.045861; 
% serial100=0.515961;
% serial200=3.11253;
% speedup200=s_200_200/serial200
% speedup100=s_100_100/serial100;
% speedup50=s_50_20/serial50;
% plot(speedup200,'blue')
% plot(speedup100,'red')
% plot(speedup50,'green')

% To plot the solution:
 m=20;
 n=50;

 load uu.dat
 load vv.dat
 load pp.dat

 % Define grid
         
   for i=1:m
    for j=1:n          
            U(i,j)=uu(i+(j-1)*m,1);
             V(i,j)=vv(i+(j-1)*m,1);
            P(i,j)=pp(i+(j-1)*m,1);
    end
   end

 %mesh(U)
 %mesh(V)
 mesh(P)
