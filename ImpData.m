% This code would import data from lammps log files into matlab 

clc, clear
delimiterIn = ' ';
headerlinesIn = 269 ;
k=0;

% Data column order:
%1Step 2CPU 3PotEng 4KinEng 5Temp 6Lx 7Ly 8Press 9v_xc_x 10v_xc_y 11v_xc_z 
%12c_pe_c60 13c_lennard 14c_ke_c60 15v_vc_x 16v_vc_y 17v_vc_z 
%18v_x1_x 19v_x1_y 20v_x1_z 21v_x2_x 22v_x2_y 23v_x2_z 24c_pe_sub 25c_ke_sub 
%26v_wc_x 27v_wc_y 28v_wc_z 29v_w12_x 30v_w12_y 31v_w12_z 32c_temp_c60 33c_temp_sub

%T=[1,2,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21  ];
T= [1,5,10,20,30,35,50,60,75,100,150,200,250,300,400,500,600,700]; %,800,900

for i=T         
    filename = sprintf('log.G12_1_T%d', i);
    eval(['temp = importdata(filename,delimiterIn,headerlinesIn)']);
    k=k+1;
    imdata(:,:,k)=temp.data;
end
save allT2s.mat imdata



