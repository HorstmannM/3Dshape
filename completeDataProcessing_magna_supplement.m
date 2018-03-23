%Matlab script for shape/form analysis of two to be compared morphotypes A and B
%The conducted analyses are described in the manuscript
%SCAN, EXTRACT, WRAP, COMPUTE – A 3D METHOD TO ANALYSE MORPHOLOGICAL SHAPE DIFFERENCES
%by Horstmann et al. 

%Start of the script
clear                                                   %clear workspace
clc                                                     %clear command window

%download the files necessary to retrace the analysis presented in
%Horstmann et al.
%the files downloaded here can be replaced by any other obj-files
%containing point clouds with equally numbered casts
outfilename_i_01=websave('i_01.obj', 'https://ndownloader.figshare.com/files/9152929?private_link=4981ca6e558c688f5c23');
outfilename_i_02=websave('i_02.obj', 'https://ndownloader.figshare.com/files/9152926?private_link=4981ca6e558c688f5c23');
outfilename_i_03=websave('i_03.obj', 'https://ndownloader.figshare.com/files/9152923?private_link=4981ca6e558c688f5c23');
outfilename_i_04=websave('i_04.obj', 'https://ndownloader.figshare.com/files/9152920?private_link=4981ca6e558c688f5c23');
outfilename_i_05=websave('i_05.obj', 'https://ndownloader.figshare.com/files/9152974?private_link=4981ca6e558c688f5c23');
outfilename_i_06=websave('i_06.obj', 'https://ndownloader.figshare.com/files/9152971?private_link=4981ca6e558c688f5c23');
outfilename_i_07=websave('i_07.obj', 'https://ndownloader.figshare.com/files/9152968?private_link=4981ca6e558c688f5c23');
outfilename_i_08=websave('i_08.obj', 'https://ndownloader.figshare.com/files/9152965?private_link=4981ca6e558c688f5c23');
outfilename_i_09=websave('i_09.obj', 'https://ndownloader.figshare.com/files/9152959?private_link=4981ca6e558c688f5c23');

outfilename_c_1=websave('c_1.obj', 'https://ndownloader.figshare.com/files/9152950?private_link=4981ca6e558c688f5c23');
outfilename_c_2=websave('c_2.obj', 'https://ndownloader.figshare.com/files/9152953?private_link=4981ca6e558c688f5c23');
outfilename_c_3=websave('c_3.obj', 'https://ndownloader.figshare.com/files/9152956?private_link=4981ca6e558c688f5c23');
outfilename_c_4=websave('c_4.obj', 'https://ndownloader.figshare.com/files/9152947?private_link=4981ca6e558c688f5c23');
outfilename_c_5=websave('c_5.obj', 'https://ndownloader.figshare.com/files/9152944?private_link=4981ca6e558c688f5c23');
outfilename_c_6=websave('c_6.obj', 'https://ndownloader.figshare.com/files/9152941?private_link=4981ca6e558c688f5c23');
outfilename_c_7=websave('c_7.obj', 'https://ndownloader.figshare.com/files/9152962?private_link=4981ca6e558c688f5c23');
outfilename_c_8=websave('c_8.obj', 'https://ndownloader.figshare.com/files/9152938?private_link=4981ca6e558c688f5c23');
outfilename_c_9=websave('c_9.obj', 'https://ndownloader.figshare.com/files/9152935?private_link=4981ca6e558c688f5c23');
outfilename_c_10=websave('c_10.obj', 'https://ndownloader.figshare.com/files/9152932?private_link=4981ca6e558c688f5c23');

%count morphotype A (here "c_..."-files)
listA=dir('*c*.obj');
fileListA={listA.name}';
nMorphA=numel(fileListA);

%count morphotype B ("i_0..."-files)
listB=dir('*i*.obj');
fileListB={listB.name}';
nMorphB=numel(fileListB);

%load data into Matlab from downloaded .obj-files in same working directory


for n=1:nMorphA                                             %import data for morphotype A in struct
    changingFileNameA=sprintf('c_%d.obj', n);
    rawDataA(:,n)=importdata(changingFileNameA);
    n=n+1;
end

nPoints=numel(rawDataA(1).data(:,1))
master([1:nPoints],1)=[1:nPoints];                          %create a vector to start with, end is dependent on number of vertices in .obj-file

for l=1:nMorphA                                             %create table from struct
    master=[master, rawDataA(1,l).data];
    l=l+1;
end


for j=1:nMorphB                                             %import data for morphotype B into struct
    changingFileNameB=sprintf('i_0%d.obj', j);
    rawDataB(:,j)=importdata(changingFileNameB);
    j=j+1;
end

for k=1:nMorphB                                             %add morphotype B values to existing morphotype A table
    master=[master, rawDataB(1,k).data];
    k=k+1;
end

master(:,1)=[];                                             %delete first column (numbering)


master1=master;                                             %keep master table and work with copy "master1"
a=1;                

    %exclude points of the grid that have not been projected onto the
    %models surface for all individual point clouds. A threshold of 700 for
    %the z-coordinate was chosen in our case, as we positioned the
    %unprojected 2D grid 700 µm above the animals symmetry plane, which is
    %also the plane in which the origin is located. The threshold needs to
    %be set far of the range, in which the actual specimens are located in
    %3D space.

while (a<=(nMorphA+nMorphB))                                %=<19 due to number of undefended+defended specimens in table = 19; adjust to your respective dataset
     
     z=1;
     nDelete=0;      
     all=0;
     master1=sortrows(master1, 3*a);                        %sort rows to find unprojected points (position deviates hugely in Z-coordinate)
     while (all  <= nPoints)
        
        if ((master1((z),3*a))<=(-700))                     %find number of rows to delete (Z-coordinate > 700)
             nDelete=nDelete+1;                             %count the number of rows to delete
             z=z+1;
        else ((master1((z),3*a))>=(-700))
                        break;
        end
     
     all=all+1;
     end
     
     if (nDelete==0)
          disp(['Nothing to delete!']);
     else 
     master1([1:(nDelete+1)],:)=[];
     end                                                    %delete the unprojected lines
     a=a+1;  
end

    %test plotting to check for further unprojected points
    
q=2;                                                        %which sample shall be exported for check? set q to the animal intended: 0-9 morphotype A, 10-18 morphotype B
                                                            %the .txt-file content can be visualised as a .xyz-file for a check with e.g. Meshlab                                                         %meshlab
dlmwrite('sample.txt',[master1(:, ((1+1*q):(3+1*q)))]);     %if required exchange name "sample" with the actual sample name

dlmwrite('master1.txt',[master1]);                          %save the master table


    %averaging the deformation results for overall displacement analysis

allMorphs=load('master1.txt');                               %create new master table to work without touching the original master table

    %separate the data for morphotypes
A=allMorphs(:, 1:(3*nMorphA));                               %separate, single data for morphotype A
B=allMorphs(:, (3*nMorphA+1):(3*nMorphA+3*nMorphB));         %separate, single data for morphotype B


    %positions of data to be averaged

avgAxNumbers(:,1)=[1:3:3*nMorphA];                           %create vector for x,y,z data position in morphotype A 
avgAyNumbers(:,1)=[2:3:3*nMorphA];
avgAzNumbers(:,1)=[3:3:3*nMorphA];


avgBxNumbers(:,1)=[1:3:3*nMorphB];                           %create vector for x,y,z data position in morphotype B
avgByNumbers(:,1)=[2:3:3*nMorphB];
avgBzNumbers(:,1)=[3:3:3*nMorphB];

    %extract data for averaging (position in table)
i=0;
while (i<numel(master1(:,1)))
    i=i+1;
    
    Ax(i,:)=A(i, avgAxNumbers);                              %extract data for A, B in columns according to morphotype and xyz
    Ay(i,:)=A(i, avgAyNumbers);
    Az(i,:)=A(i, avgAzNumbers);
    
    Bx(i,:)=B(i, avgBxNumbers);
    By(i,:)=B(i, avgByNumbers);
    Bz(i,:)=B(i, avgBzNumbers);
    
    disp(['Already ' num2str(i) ' done!'])
end

    %actual averaging

i=0;
while (i<numel(master1(:,1)))
   i=i+1;
   j=1;
   avgPAx=0;
   avgPAy=0;
   avgPAz=0;
   avgPBx=0;
   avgPBy=0;
   avgPBz=0;
   while j<nMorphA
        avgPAx=avgPAx+Ax(i,j);
        avgPAy=avgPAy+Ay(i,j);
        avgPAz=avgPAz+Az(i,j);
        j=j+1;
   end
   j=1;
   while j<nMorphB
        avgPBx=avgPBx+Bx(i,j);
        avgPBy=avgPBy+By(i,j);
        avgPBz=avgPBz+Bz(i,j);
        j=j+1;
   end
   avgAx(i,:)=avgPAx/nMorphA;
   avgAy(i,:)=avgPAy/nMorphA;
   avgAz(i,:)=avgPAz/nMorphA;
   
   avgBx(i,:)=avgPBx/nMorphB;
   avgBy(i,:)=avgPBy/nMorphB;
   avgBz(i,:)=avgPBz/nMorphB;
   
   disp(['Already ' num2str(i) ' done!'])
end

    %create point clouds of both morphotypes
avgA=[avgAx avgAy avgAz];                                   %combine to resulting averaged vector for morphotype A

avgB=[avgBx avgBy avgBz];                                   %combine to resulting averaged vector for morphotype B


%adjust positioning and orientation for shape comparison
    %call procrustes
[d, Bt, tr]=procrustes(avgA,avgB,'Scaling',false);          %this is for partial Procrustes analysis, for the full Procrustes anylsis, the 'Scaling' argument "false" needs to be replaced by "true"

    %display and write the transformed data into textfile
Bt;                                                         %for display remove semicolon
dlmwrite('Bt.txt',Bt)                                       %save transformed average to text-file

    %also write C average into textfile
dlmwrite('avgA.txt',avgA)

    %calculate displacement for X,Y,Z-axis
displacements=Bt-avgA;

    %calculate overall displacement

displacements(:,4)=sqrt((displacements(:,1)).^2+(displacements(:,2)).^2+(displacements(:,3)).^2);

    %write transformed data of morphotype "I" to textfile
dlmwrite('displacements.txt',displacements)


%plot results for overall displacement and according to coordinate axes
S=5;                                                        %define point size in point cloud
allMorphs=displacements(:,4);                               %enter as argument for displacement: 1 for x, 2 for y, 3 for z, 4 for overall displacement; choosing 3, insert a minus in front of "displacements" to get the right colourcoding (inverted for the z-axis)
XB=Bt(:,1);                                                 %define X, Y and Z based on Procrustes-adjusted shape
YB=Bt(:,2);
ZB=Bt(:,3);
scatter3(XB,YB,ZB,S,allMorphs);                             %actual plotting
axis equal                                                  %equal spacing on axes
xlabel ('X');                                               %attach labels to axes
ylabel ('Y');
zlabel ('Z');
colormap(jet)                                               %choose the colormap jet (blue to red)
colorbar                                                    %insert colorbar

view(0, 90);                                                %define viewing angle (in z-direction here)
caxis ([-0 100])                                            %adjust to deformation range you want to look at, by defining the upper and lower threshold for the colorbar 


    %export figures  
set(gca,'XTick', -1500:1000:300);                           %adjust spacing on axes; necessary if axis are on 
set(gca,'FontSize',20)
axis off                                                    %decide if you like it with or without axes
grid off                                                    %decide if you like it with or without grid

print('dispZ_fig_1', '-djpeg')                              %save to jpeg-file, enter desired filename in the first quotation marks


%Wilcoxon-test between both morphotypes
    %preparation for morphotype A
m=1;
o=0;
for m=1:nMorphA                                             %create individual textfiles for each specimen (A-morphotype)
   changingFileNameA2=sprintf('A_%d.obj', m);
   individualX=A(:, (1+o*3));
   individualY=A(:, (2+o*3));
   individualZ=A(:, (3+o*3));
   individual=[individualX individualY individualZ];
   txtfile=['A_' num2str(m) '.txt'];
   dlmwrite(txtfile ,[individual]);
   m=m+1;
   o=o+1;
    
end

    %preparation for morphotype B
m=1;
o=0;
for m=1:nMorphB                                             %create individual textfiles for each specimen (B-morphotype)
   changingFileNameB2=sprintf('B_%d.obj', m);
   individualX=B(:, (1+o*3));
   individualY=B(:, (2+o*3));
   individualZ=B(:, (3+o*3));
   individual=[individualX individualY individualZ];
   txtfile=['B_' num2str(m) '.txt'];
   dlmwrite(txtfile ,[individual]);
   m=m+1;
   o=o+1;
   end

    %load extracted individual data into workspace
    %for morphotype A
    i=1;
    while i<=nMorphA
        indA=sprintf('A_%d.txt', i);
        if (i==1)
            AArray{1}=load(indA);
        elseif (i==2)
            AArray{2}=load(indA);
        elseif (i==3)
            AArray{3}=load(indA);
        elseif (i==4)
            AArray{4}=load(indA);
        elseif (i==5)
            AArray{5}=load(indA);
        elseif (i==6)
            AArray{6}=load(indA);
        elseif (i==7)
            AArray{7}=load(indA);
        elseif (i==8)
            AArray{8}=load(indA);
        elseif (i==9)
            AArray{9}=load(indA);
        elseif (i==10)
            AArray{10}=load(indA);
        elseif (i==11)
            AArray{11}=load(indA);
        elseif (i==12)
            AArray{12}=load(indA);
        elseif (i==13)
            AArray{13}=load(indA);
         elseif (i==14)
            AArray{14}=load(indA);
         elseif (i==15)
            AArray{15}=load(indA);
        else
            disp('You have more than 15 samples, add futher cases in the script.')
        end
        i=i+1;
    end
               
    %for morphotype B
    i=1;
    while i<=nMorphB
        indB=sprintf('B_%d.txt', i);
        if (i==1)
            BArray{1}=load(indB);
        elseif (i==2)
            BArray{2}=load(indB);
        elseif (i==3)
            BArray{3}=load(indB);
        elseif (i==4)
            BArray{4}=load(indB);
        elseif (i==5)
            BArray{5}=load(indB);
        elseif (i==6)
            BArray{6}=load(indB);
        elseif (i==7)
            BArray{7}=load(indB);
        elseif (i==8)
            BArray{8}=load(indB);
        elseif (i==9)
            BArray{9}=load(indB);
        elseif (i==10)
            BArray{10}=load(indB);
        elseif (i==11)
            BArray{11}=load(indB);
        elseif (i==12)
            BArray{12}=load(indB);
        elseif (i==13)
            BArray{13}=load(indB);
         elseif (i==14)
            BArray{14}=load(indB);
         elseif (i==15)
            BArray{15}=load(indB);
        else
            disp('You have more than 15 samples, add futher cases in the script.')
        end
        i=i+1;
    end 

    %Procrustes-fit all individual data to a selected specimen (A3 here);
    %For your own data select an individual that represents the prototype
    %best
    
    %which model shall represent the prototype? Type model number 1-10 of
    %morphotype A into the curly braces below
    prototype=AArray{3};
       
        %conduct procrustes fits onto prototype; the current setting will
        %produce a partial procrustes analysis, a full procrustes analysis
        %is achieved by replacing the 'scaling' arguments "false" with "true"
    i=1;
    while i<=nMorphA
        [d, procrustesTableA{i}, tr]=procrustes(prototype,AArray{i},'Scaling',false);
        i=i+1;
    end
    
    i=1;
    while i<=nMorphB
        [d, procrustesTableB{i}, tr]=procrustes(prototype,BArray{i},'Scaling',false);
        i=i+1;
    end
       
    %compare single specimens of both morphotypes in overlapping view to
    %check that the Procrustes superimposition worked
figure;
    S=5;
XB=procrustesTableB{6}(:,1);                            %change the treatment and case number to check how the animals are positioned relative to each other 
YB=procrustesTableB{6}(:,2);
ZB=procrustesTableB{6}(:,3);

XA=procrustesTableA{5}(:,1);
YA=procrustesTableA{5}(:,2);
ZA=procrustesTableA{5}(:,3);

hold on                                                 %plot two point clouds in one figure
scatter3(XB,YB,ZB,S,'g');
scatter3(XA,YA,ZA,S,'r');
hold off

axis equal                                              %equal spacing on axes
xlabel ('X');
ylabel ('Y');
zlabel ('Z');
view(0, 90);                                            %viewing angle (in z-direction here)
    %end of comparison

    
    %assemble all procrustes transformed data in one single table
    i=1;
ergMatrix(:,1)=1:numel(procrustesTableA{1}(:,1));
while i<=(nMorphA)
    ergMatrix=[ergMatrix, procrustesTableA{i}];
    i=i+1;
end

i=1;
while (i<=nMorphB)
    ergMatrix=[ergMatrix, procrustesTableB{i}];
    i=i+1;
end
ergMatrix(:,1)=[];

thirdElementAx=[1:3:(3*nMorphA)];                              %this is preparation for Ax; adjust here and below to the number of cases you have per treatment (here 10*3 coordinates for C)
thirdElementAy=[2:3:(3*nMorphA)];                              %this is preparation for Ay
thirdElementAz=[3:3:(3*nMorphA)];                              %this is preparation for Az

thirdElementBx=[(3*nMorphA+1):3:(3*(nMorphA+nMorphB))];        %this is preparation for Bx; since B-data is attached to the A-data, start with 3*10+1 and end with 3*10 + 3*9 (the number of Ind-cases)
thirdElementBy=[(3*nMorphA+2):3:(3*(nMorphA+nMorphB))];        %this is preparation for By
thirdElementBz=[(3*nMorphA+3):3:(3*(nMorphA+nMorphB))];        %this is preparation for Bz


    %extract Procrustes-fitted data for Wilcoxon-test
i=0;
while (i<numel(ergMatrix(:,1)))
    i=i+1;
    
    Ax(i,:)=ergMatrix(i, thirdElementAx);                      %extract data for C and I in columns for these two treatments and separated for xyz
    Ay(i,:)=ergMatrix(i, thirdElementAy);
    Az(i,:)=ergMatrix(i, thirdElementAz);
    
    Bx(i,:)=ergMatrix(i, thirdElementBx);
    By(i,:)=ergMatrix(i, thirdElementBy);
    Bz(i,:)=ergMatrix(i, thirdElementBz);
    
    disp(['Already ' num2str(i) ' done!'])                     %plot progress of calculation
end

    %actual statistical testing with U-tests
    
    
xTestA=zeros(numel(AArray(1,:)),1);                     
xTestA(numel(AArray(1,:)),1)=NaN;                              %create table for morphotype A
yTestA=zeros(numel(AArray(1,:)),1);
yTestA(numel(AArray(1,:)),1)=NaN;
zTestA=zeros(numel(AArray(1,:)),1);
zTestA(numel(AArray(1,:)),1)=NaN;

xTestB=zeros(numel(AArray(1,:)),1);                     
xTestB(numel(AArray(1,:)),1)=NaN;                              %create table for morphotype B
yTestB=zeros(numel(AArray(1,:)),1);
yTestB(numel(AArray(1,:)),1)=NaN;
zTestB=zeros(numel(AArray(1,:)),1);
zTestB(numel(AArray(1,:)),1)=NaN;

i=0;
while (i<numel(ergMatrix(:,1)))
    i=i+1;
    
    xTestA=Ax(i,:)';                                           %transpose the matrices containing the values separated by xyz and treatment and extract (within this loop) values for a single point as preparation for the Wilcoxon-tests 
    xTestB(1:numel(BArray(1,:)),1)=Bx(i,:)';
    
    yTestA=Ay(i,:)';
    yTestB(1:numel(BArray(1,:)),1)=By(i,:)';
    
    zTestA=Az(i,:)';
    zTestB(1:numel(BArray(1,:)),1)=Bz(i,:)';
         
    pux=ranksum(xTestA,xTestB);                                %conduct the Wilcoxon-tests
    puy=ranksum(yTestA,yTestB);
    puz=ranksum(zTestA,zTestB);
   
    pValuesU(i,1)=pux;                                         %concatenate all Wilcoxon-test variables in one variable "pValuesU"
    pValuesU(i,2)=puy;
    pValuesU(i,3)=puz;
    
    disp(['Already ' num2str(i) ' done!'])                     %plot progress of calculation
end

    %how many points were found significant?
pValuesUxs=sort(pValuesU(:,1));                                %sort for number of significantly different x-coordinates
pValuesUys=sort(pValuesU(:,2));                                %sort for number of significantly different y-coordinates
pValuesUzs=sort(pValuesU(:,3));                                %sort for number of significantly different z-coordinates

%for X
i=1;
countX=0;
while i<=numel(pValuesUxs(:))
if pValuesUxs(i)<=0.01                                         %adjust the significance threshold for X by varying 0.01
    countX=countX+1;
    i=i+1;
else
    disp('Not significant!')
    i=i+1;
end
end

%for Y
i=1;
countY=0;
while i<=numel(pValuesUys(:))
if pValuesUys(i)<=0.01                                  %adjust the significance threshold for Y by varying 0.01
    countY=countY+1;
    i=i+1;
else
    disp('Not significant!')
    i=i+1;
end
end

%for Z
i=1;
countZ=0;
while i<=numel(pValuesUzs(:))
if pValuesUzs(i)<=0.01                                  %adjust the significance threshold for Z by varying 0.01
    countZ=countZ+1;
    i=i+1;
else
    disp('Not significant!')
    i=i+1;
end
end


    %decide, for which direction results shall be plotted
i=0;
a=2;                                    %enter 1, 2 or 3 as value for a to get x-z in the following plotting
while (i<numel(ergMatrix(:,1)))
    i=i+1;
coloursPs_u(i,1)=pValuesU(i,a);
end

    %plot Wilcoxon-test results 
S=5;
Colours=coloursPs_u;
XB=Bt(:,1);
YB=Bt(:,2);
ZB=Bt(:,3);
scatter3(XB,YB,ZB,S,Colours);

axis equal                                              %equal spacing on axes
xlabel ('X');                                           %plot axis labels
ylabel ('Y');
zlabel ('Z');

view(0, 90);                                            %define viewing angle (in z-direction here)

cm=flipud(jet);                                         %inverted colourmap (intuitive red = high significance)
colormap(cm);
caxis ([0 0.02])                                        %(adjust to significance range you want to look at)
caxis([0 0.1])

    %export images for posters, presentation etc.
set(gca,'XTick', -1500:500:300);                        %adjust the axes intervals
set(gca,'FontSize',20)
colorbar                                                %plot colorbar
axis off                                                %decide if you want to disable axes
grid off                                                %decide if you want to disable grid

    %create a figure in the current folder
print('morphotype_utest_0.01_x', '-djpeg')


    % plot results of the Wilcoxon-test in blue, red, yellow, green as indicators of
    % significance level (additional plotting with uncontinuous colour regime)
pTest=pValuesU(:,3);
Colours=[0 0 0];

for index=1:length(pTest)

  if(pTest(index)<0.001)
      Colours(index,:)=[0 255 0];                       %coloursPz is the new vector with only colour tags 
  elseif(pTest(index)<0.01)
      Colours(index,:)=[255 255 0]; 
  elseif(pTest(index)<0.05)
      Colours(index,:)=[255 0 0];
  else
      Colours(index,:)=[0 0 255];
  end
end

    %actual plotting
S=5;
Colours;
XB=Bt(:,1);
YB=Bt(:,2);
ZB=Bt(:,3);
scatter3(XB,YB,ZB,S,Colours);
axis equal                                              %equal spacing on axes
xlabel ('X');
ylabel ('Y');
zlabel ('Z');
view(0, 90);                                            %define viewing angle (in z-direction here)

%-------------------------------------------------------------------------------------------
    %check test results for multiple testing (Storey and Tibshirani (2003))
    %export the p-values as textfile 
dlmwrite('pValues.txt',pValuesU);

    %calculate the FDR-based approach (after Storey and Tibshirani (2003))
[FDR_x, Q_x, Pi0_x]= mafdr(pValuesU(:,1), 'showplot', true);
[FDR_y, Q_y, Pi0_y]= mafdr(pValuesU(:,2), 'showplot', true);
[FDR_z, Q_z, Pi0_z]= mafdr(pValuesU(:,3), 'showplot', true);

    %plot in command, in how many of the conducted tests there is in fact a
    %difference between morphotype A and B
PiA_x=(1-Pi0_x)*100;
PiA_y=(1-Pi0_y)*100;
PiA_z=(1-Pi0_z)*100;
disp(['For the x-direction, at least ' num2str(PiA_x) '% of the 127087 conducted tests are expected to follow the alternative hypothesis, i.e. alteration is going on.']);
disp(['For the y-direction, at least ' num2str(PiA_y) '% of the 127087 conducted tests are expected to follow the alternative hypothesis, i.e. alteration is going on.']);
disp(['For the z-direction, at least ' num2str(PiA_z) '% of the 127087 conducted tests are expected to follow the alternative hypothesis, i.e. alteration is going on.']);

    %plot the results of the correction
S=5;
    %select direction to plot with q-values
Colours=Q_x;                                            %execute this line to plot the q-values for x
Colours=Q_y;                                            %execute this line to plot the q-values for y
Colours=Q_z;                                            %execute this line to plot the q-values for z

XB=Bt(:,1);
YB=Bt(:,2);
ZB=Bt(:,3);
scatter3(XB,YB,ZB,S,Colours);

axis equal                                              %equal spacing on axes
xlabel ('X');                                           %plot axis labels
ylabel ('Y');
zlabel ('Z');

view(0, 90);                                            %define viewing angle (in z-direction here)

cm=flipud(jet);                                         %inverted colourmap (intuitive red = high significance)
colormap(cm);
colorbar
caxis([0 0.01])

    %print results of q-values to figures
set(gca,'XTick', -1500:500:300);                        %adjust the axes intervals
set(gca,'FontSize',20)
colorbar                                                %plot colorbar
axis off                                                %decide if you want to disable axes
grid off                                                %decide if you want to disable grid
print('q-value_001_x', '-djpeg')
print('q-value_001_y', '-djpeg')
print('q-value_001_z', '-djpeg')


    %plot p-value against q-value
    %for x
plot(sort(pValuesU(:,1)), sort(Q_x), 'k');
xlabel ('p-Values');                                     %plot axis labels
ylabel ('q-Values');
ylim([0 0.35]);
set(gca,'FontSize',20)
print('p-value-q-value-plot_x', '-djpeg')                %export image of the p-q-plot

    %for y
plot(sort(pValuesU(:,2)), sort(Q_y), 'k');
xlabel ('p-Values');                                     %plot axis labels
ylabel ('q-Values');
ylim([0 0.35]);
set(gca,'FontSize',20)
print('p-value-q-value-plot_y', '-djpeg')                %export image of the p-q-plot

    %for z
plot(sort(pValuesU(:,3)), sort(Q_z), 'k');
xlabel ('p-Values');                                     %plot axis labels
ylabel ('q-Values');
ylim([0 0.35]);
set(gca,'FontSize',20)
print('p-value-q-value-plot_z', '-djpeg')                %export image of the p-q-plot

%----------------------------------------------------------------------------------------------------

    %confidence blops-comparison

    %calculate the averaged point position, variance and confidence intervals for both morphotypes

%t-table values for 95 % confidence intervals
tTable95=[12.71 4.303 3.182	2.776 2.571	2.447 2.365	2.306 2.262	2.228 2.201	2.179 2.16 2.145 2.131 2.12	2.11 2.101 2.093 2.086 2.08	2.074 2.069	2.064 2.06 2.056 2.052 2.048 2.045 2.042];

    %for morphotype A
i=0;
while (i<numel(ergMatrix(:,1)))
    i=i+1;
    
    nXA=length(Ax(i,:));                                    %length of point vector
    MittelXA(i,1)=mean(Ax(i,:));                            %average x-positions of all points
    VarianzXA(i,1)=var(Ax(i,:));                            %calculate the respective variance
    XCIA(i,1)=tTable95(nMorphA) * sqrt(VarianzXA(i,1)/nXA); %calculate single confidence intervals for morphotype undefended (c) in x-direction, t-values for CI's taken from table, according to number of samples  
    
    
    
    nYA=length(Ay(i,:));                                    %length of point vector
    MittelYA(i,1)=mean(Ay(i,:));                            %average y-positions of all points
    VarianzYA(i,1)=var(Ay(i,:));                            %calculate the respective variance
    YCIA(i,1)=tTable95(nMorphA) * sqrt(VarianzYA(i,1)/nYA); %calculate single confidence intervals for morphotype undefended (c) in y-direction
    
    nZA=length(Az(i,:));                                    %length of point vector
    MittelZA(i,1)=mean(Az(i,:));                            %average z-positions of all points
    VarianzZA(i,1)=var(Az(i,:));                            %calculate the respective variance
    ZCIA(i,1)=tTable95(nMorphA) * sqrt(VarianzZA(i,1)/nZA); %calculate single confidence intervals for morphotype undefended (c) in z-direction
    
    disp(['Already ' num2str(i) ' done!'])                  %plot progress of calculation
    
end

    %for morphotype B
i=0;
while (i<numel(ergMatrix(:,1)))
    i=i+1;
    
    nXB=length(Bx(i,:));                                    %length of point vector
    MittelXB(i,1)=mean(Bx(i,:));                            %average x-positions of all points
    VarianzXB(i,1)=var(Bx(i,:));                            %calculate the respective variance
    XCIB(i,1)=tTable95(nMorphB) * sqrt(VarianzXB(i,1)/nXB); %calculate single confidence intervals for morphotype defended (i) in x-direction, t-values for CI's taken from table, according to number of samples  
    
    nYB=length(By(i,:));                                    %length of point vector
    MittelYB(i,1)=mean(By(i,:));                            %average y-positions of all points
    VarianzYB(i,1)=var(By(i,:));                            %calculate the respective variance
    YCIB(i,1)=tTable95(nMorphB) * sqrt(VarianzYB(i,1)/nYB); %calculate single confidence intervals for morphotype defended (i) in y-direction
    
    nZB=length(Bz(i,:));                                    %length of point vector
    MittelZB(i,1)=mean(Bz(i,:));                            %average z-positions of all points
    VarianzZB(i,1)=var(Bz(i,:));                            %calculate the respective variance
    ZCIB(i,1)=tTable95(nMorphB) * sqrt(VarianzZB(i,1)/nZB); %calculate single confidence intervals for morphotype defended (i) in z-direction
    
    disp(['Already ' num2str(i) ' done!'])                  %plot progress of calculation
end

    % test whether the confidence intervals overlap (this step takes a lot of computing time!)
    %to compute the overlap with the function ellipsoid(), the
    %ELLIPSOIDAL TOOLBOX by Peter Gagarinov, Alex Kurzhanskiy
    %needs to be installed/referenced in path list
    %<<http://systemanalysisdpt-cmc-msu.github.io/ellipsoids/doc/chap_install.html?highlight=download>>

i=0;
while (i<numel(ergMatrix(:,1)))
    i=i+1;
       
    QA=[(XCIA(i,1)^(-2)) 0 0; 0 (YCIA(i,1)^(-2)) 0; 0 0 (ZCIA(i,1)^(-2))];         %create ellipsoid matrix for C
    qA=[MittelXA(i,1); MittelYA(i,1); MittelZA(i,1)];                              %calculate center for confidence ellipsoid of respective point (c)
    
    QB=[(XCIB(i,1)^(-2)) 0 0; 0 (YCIB(i,1)^(-2)) 0; 0 0 (ZCIB(i,1)^(-2))];      %create ellipsoid matrix for I
    qB=[MittelXB(i,1); MittelYB(i,1); MittelZB(i,1)];                           %calculate center for confidence ellipsoid of respective point (ind)
    
    ellipsoidA=ellipsoid(qA,inv(QA));                                     %create both ellipsoids step by step within the loop
    ellipsoidB=ellipsoid(qB,inv(QB));
    
    distABBlobs(i,1)=ellipsoidA.distance(ellipsoidB);             %calculate distance between ellipsoids
    
    if(distABBlobs(i,1)>0)                              %associate a distance larger than 0 (= no overlap) with 0
        check(i,1)=0;
    elseif(distABBlobs(i,1)<=0)                         %associate a distance smaller/equal 0 (= overlap) with 1
        check(i,1)=1;
    else
        check(i,1)=2;                                   %indicates error, since distCIBlob can usually be either greater, smaller or equal 0
    end
    disp(['Already ' num2str(i) ' done!'])              %plot progress of calculation
end

    %create a colorscheme, in which red indicates no overlapping, blue overlapping

Colours=[0 0 0];
for index=1:length(check)
                                                        %Colours is the new vector with only colour tags
    if(check(index)==0)
Colours(index,:)=[255 0 0];                             %red
    elseif(check(index)==1)
         Colours(index,:)=[0 0 255];                    %blue
    else
        Colours(index,:)=[155 0 155];                   %violet (indicates errors in the calculation, if present)
    end
end


    %plot results
scatter3(XB,YB,ZB,S,Colours, 'filled');

axis equal                                              %equal spacing on axes
xlabel ('X');                                           %add axis labels
ylabel ('Y');
zlabel ('Z');
view(0, 90);                                            %define viewing angle (in z-direction here)

    %export figures
set(gca,'XTick', -1500:300:300);                        %adjust grid
set(gca,'FontSize',20)                                  %adjust font size

print('confidenceBlops_fig', '-djpeg')                  %save as jpeg-file


%final saving of workspace
save('workspaceCompleteDataProcessing.mat');            %save the workspace with the variable values
load('workspaceCompleteDataProcessing.mat');            %load existing workspace (for later changes on the script)