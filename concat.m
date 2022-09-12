cd 'D:\dwl\MSC\PK\TP\Data months\DWC2020'% Directory of the files
m5=readtable('D1_1_2020.csv');%file name
m5(:,5)=[];
m5=m5{:,:};
for j=1:12
    if j==1 || j==3 || j==5 || j==7 || j==8 || j==10 || j==12
       dmax=31;
    elseif j==2
        dmax=29;
    else
        dmax=30;
    end
     for i=1:dmax
        matrixMay2=readtable(['D' num2str(i) '_' num2str(j) '_2020']);
        matrixMay2(:,5)=[];
        matrixMay2(:,1:3)=[];
        matrixMay2=matrixMay2{:,:};
        matrix3=cat(2,m5,matrixMay2);
        m5=matrix3;
     end
end


            